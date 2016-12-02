# This script has routines to calculate various statistics from the survivor output file.
# Will be combined into a single file with other contributors' routines when done.

use File::Basename;
use List::Util qw(min max);


# extract the executable directory so can find the data structures code

my $data_struc_code = 'data_structures.pl';
my $exe_dirname = dirname(__FILE__);
require "$exe_dirname/$data_struc_code";

usage() unless(@ARGV > 1);
my ($vcf_file, $tech_file, $outfile_prefix) = ($ARGV[0], $ARGV[1], $ARGV[2]);
my $min_callers=4; #to do: make this an argument

my $techniques = map_technologies($tech_file);

# Load the data
print STDERR "loading data by event...\n";
my $data_by_event_href = generate_data_by_event($vcf_file,$tech_file);
print STDERR "loading data by caller...\n";
my $data_by_caller_href = generate_data_by_caller($vcf_file,$tech_file);
print STDERR "processing...\n";

# total number of events
$total_events = scalar ( keys %$data_by_event_href ) ;

# --------------------------------------------------------------
# set up the file names for the output files before running each analysis.
# filenames are the prefix input on the command line, with an analysis-specific suffix.
# the file extensions are stored in the suffix hash declared here.
# --------------------------------------------------------------

my %suffix;  

# --------------------------------------------------------------
# set up and summarize event filters
# --------------------------------------------------------------

# Set up values to use for filtering events
($filters_href,$total_failed_events) = &compute_filter_values($data_by_event_href, $min_callers);

# Print statistics on the number of events failing each filter
$suffix{'filters'} = '.filters.tsv';
&print_filters($filters_href,$total_events,$total_failed_events,$outfile_prefix.$suffix{'filters'});


# --------------------------------------------------------------
# Calculate the total number of events called in each callset (ie by a particular calling algorithm).
# If a call set has more than one variant within the range of the event, count as one variant.
# Events are unfiltered.
# --------------------------------------------------------------

$suffix{'callset_event_count'} = '.callset_event_count.tsv';
&callset_event_count( $data_by_caller_href, $outfile_prefix.$suffix{'callset_event_count'});


# --------------------------------------------------------------
# Count the number of variants of each type in each callset.
# If there are multiple variants per event, count them all.
# Print 0's for types not present in a call set.
# Events are unfiltered.
# --------------------------------------------------------------

$suffix{'callset_total_count_by_type'} = '.callset_total_count_by_type.tsv';
&callset_total_count_by_type( $data_by_caller_href, $outfile_prefix.$suffix{'callset_total_count_by_type'});


# --------------------------------------------------------------
# Create histogram of number of callsets per event (ie number of callers that called this SV)
# ignoring the variant type.
# --------------------------------------------------------------

$suffix{'callsets_per_event'} = '.callsets_per_event.tsv';
&callsets_per_event($data_by_event_href,$outfile_prefix.$suffix{'callsets_per_event'});


# --------------------------------------------------------------
# Calculate the consistency of the breakpoints.
# Limit the analysis to those events that are called by
# at least min_caller callers for which all called variants are of the same type.
# Calculate distance of the start and stop for each caller to 
# the median for the event.
# --------------------------------------------------------------
$suffix{'breakpoint_variation'} = '.breakpoint_variation.tsv';
&breakpoint_variation( $data_by_event_href, $filters_href, $outfile_prefix.$suffix{'breakpoint_variation'} );


# --------------------------------------------------------------
# Calculate basic statistics on the gene annotations in the demo vcf file
# --------------------------------------------------------------

$suffix{'annotation_stats'} = '.annotation_counts.tsv';
@strings_to_match=('\bENSG','\bAllRepeats_');
&summarize_annotations($data_by_event_href,$outfile_prefix.$suffix{'annotation_stats'},\@strings_to_match);

exit 0;

#==================================================
# SUBROUTINES
#==================================================


#------------------------------------------
# histogram of the number of callers "supporting" a variant.
# Counts a caller if there is any variant called at the same location  -
# ignores the variant type called by each caller.
#------------------------------------------
sub callsets_per_event {
  my ($data_href,$outfile) = @_;
  my %hist; # number of variants called by the indicated number of callers;

  while(my ($event_id,$event_info) = each(%$data_href)) {
    my $caller_count = $event_info->{'caller_count'};
    $hist{$caller_count}++;
  }

  # print the result to the output file. Header first, then data.
  open(OUT,">$outfile");
  # create a header line for the output file
  my @header = qw(NUM_CALLERS NUM_VARIANTS);
  print OUT join "\t",@header; print OUT "\n";

  # now print the data
  foreach my $k (keys %hist) {
    print OUT "$k\t$hist{$k}\n";
  }
  close(OUT);
  return;
}


#------------------------------------------
# total number of events with at least one variant called in each callset
# ignore variant type
#------------------------------------------
sub callset_event_count{
  my( $data_href, $outfile) = @_;

  # prepare the output file and header
  open(OUT,">$outfile");
  my @header = qw(CALLER N_EVENTS);
  print OUT join "\t",@header; print OUT "\n";

  while(my ($caller,$records) = each(%$data_href)) {
    print OUT $caller,"\t",$records->{"caller_count"},"\n";
  }

  close(OUT);
  return;
}


#------------------------------------------
# total number of variants in each callset by variant type.
# if there is more than one call for a callset for an event, count each sub-event separately.
#------------------------------------------
sub callset_total_count_by_type{
  my( $data_href, $outfile) = @_;

  # prepare the output file and header
  open(OUT,">$outfile");
  my @header = qw(CALLER VAR_TYPE N_VARIANTS);
  print OUT join "\t",@header; print OUT "\n";

  my %counts; #count of each type for each caller
  my %all_types; #track all the types in the file so print complete table with 0's

  while(my ($caller,$records) = each(%$data_href)) {
    while(my ($record_id,$regions) = each(%$records)) {
      my $events = $regions->{"events"};
      my $types = $events->{"types"};
      foreach my $type (@$types) {
	$counts{$caller}{$type}++;
	$all_types{$type} = 1;
      }
    }
  }

  # print the result for this caller
  foreach my $caller (keys %counts) {
    foreach my $type (sort keys %all_types) {
      if (!defined $counts{$caller}{$type}) { $counts{$caller}{$type} = 0; }
      print OUT join ("\t", $caller,$type,$counts{$caller}{$type}); 
      print OUT "\n";
    }
  }

  close(OUT);
  return;
}


#------------------------------------------
# For each event, determine if :
# - it has support from at least num_callers callers (LOW_CALL_SUPPORT)
# - all calls made by each caller are the same type (for those callers with multiple calls by one caller for an event)	 (MULTIPLE_TYPES_PER_CALLER)
# - all calls across callers are the same type (MULTIPLE TYPES ACROSS CALLERS)
# Save data both by event and by filter.
#------------------------------------------
sub compute_filter_values {
  my( $data_href, $min_callers)=@_;
  my %event_type;
  my %event_filters;
  my %filters;

  # loop over each row=event
  while(my ($event_id,$event_info) = each(%$data_href)) {

    my $event_type = $event_info->{"alt"};
    # strip off the brackets
    $event_type =~ s/[\<\>]//g;
    my $num_supporting_callers = $event_info->{"supp"};


    # is this event supported by at least min_caller callers
    if ($num_supporting_callers < $min_callers) {
       $filters{'LOW_CALL_SUPPORT'}{$event_id}=1;
    }

    # are all the svs for this event of the same type
    my $callers_href = $event_info->{"callers"};
    my $events_href;
    my $types_aref;

    # loop over all callers for this event
    while(my($caller_id, $regions_href) = each(%$callers_href)) {
       $events_href = $regions_href->{"events"};
       $types_aref = $events_href->{"types"};
       @uniq_types = &uniq(@{$types_aref});
       # are all svs called by this one caller for this event the same type
       if ( scalar @uniq_types >1 ) {
         $filters{'MULTIPLE_TYPES_ACROSS_CALLERS'}{$event_id} = 1;
         $filters{'MULTIPLE_TYPES_PER_CALLER'}{$event_id} = 1;
       } elsif ($uniq_types[0] ne $event_type)  {
         # count number of callers with a type not matching the event type
         # if there are multiple types for a caller then by definition one 
	 # doesn't match and this is assigned above
         $filters{'MULTIPLE_TYPES_ACROSS_CALLERS'}{$event_id}++;
       }
    }


    # save number of events failing any filter
    # events passing filters are undefined
    if ( defined $filters{'LOW_CALL_SUPPORT'}{$event_id}  ||
         defined $filters{'MULTIPLE_TYPES_PER_CALLER'}{$event_id}   ||
         defined $filters{'MULTIPLE_TYPES_ACROSS_CALLERS'}{$event_id}   ) {
         $total_failed_events++;
    }

  } #end event loop

  return (\%filters,$total_failed_events);

}


#------------------------------------------
# print statistics - counts of number of events failing each filter
#------------------------------------------

sub print_filters {
  my($filters_href,$total_events,$total_failed_events,$outfile)=@_;

  # set up output file
  open(OUT,">$outfile");

  # print total number of events
  print OUT join "\t", ("EVENTS_IN",$total_events);
  print OUT "\n";

  foreach my $f (keys %$filters_href) {
    my $n = scalar (keys %{$filters_href->{$f}});
    print OUT join "\t", ($f,$n);
    print OUT "\n";
  }

  # print TOTAL_FAILED_EVENTS here
  print OUT join "\t", ('TOTAL_FAILED_EVENTS', $total_failed_events);
  print OUT "\n";

  # final count of all events less those failing any filter
  print OUT join "\t", ("EVENTS_OUT",$total_events - $total_failed_events);
  print OUT "\n";

  close(OUT);

  return;
}

#------------------------------------------
# For events (SVs/ rows) 
#  - supported by at least num_callers callers
#  - for which all calls in any caller are the same type
# For each event, compute the distance of the start and stop to the median of the calls
#------------------------------------------
sub breakpoint_variation {
  my( $data_href, $filters_href, $outfile)=@_;
  my %event_type;
  my %caller_count;

  # set up output files
  open(OUT,">$outfile");
  my @header = qw(EVENT_ID TYPE LENGTH NUM_SUPPORTING_CALLERS CALLER START_DIST_TO_MEDIAN STOP_DIST_TO_MEDIAN);
  print OUT join "\t",@header; print OUT "\n";

  # loop over each row=event
  while(my ($event_id,$event_info) = each(%$data_href)) {

    # skip the events that fail the filters 
    
    if ( defined $filters_href->{'LOW_CALL_SUPPORT'}{$event_id}  || 
         defined $filters_href->{'MULTIPLE_TYPES_PER_CALLER'}{$event_id} || 
         defined $filters_href->{'MULTIPLE_TYPES_ACROSS_CALLERS'}{$event_id}  ) {
	 next;
     }

     # process the events that pass the filters

     # some info about the event to save for the output

     $event_length{$event_id} = $event_info->{'stop'} - $event_info->{'start'} + 1;
     $caller_count{$event_id} = $event_info->{'caller_count'};
     $num_supporting_callers{$event_id} = $event_info->{"supp"};
     $event_type{$event_id} = $event_info->{"alt"};
     # strip off the brackets
     $event_type{$event_id} =~ s/[\<\>]//g;

     # Calculate the median start and stop across the callers for this event.
     # For each caller, calculate the abs distance for each start and stop to the median.
     # If there are multiple events for one caller, use the min start and max stop.
     # Assume the coordinates are all on the same chromosome.


     # collect the starts and stops to use for the median calculation
     my $callers_href = $event_info->{"callers"};
     my @estarts;  my @estops; my @ecallers;
     while(my($caller_id, $regions_href) = each(%$callers_href)) {
 
        $events_href = $regions_href->{"events"};
        my $coords_href = $events_href->{"coords"};
        while(my ($type, $regions_aref) = each %$coords_href) {
 	 # find the min start and max stop for this event for this caller
 	 my @cstarts; my @cstops;
         foreach my $coord (@$regions_aref) {
 	   push @cstarts, $coord->{"start"};
 	   push @cstops, $coord->{"stop"};
         }
         $min_start = min @cstarts;
 	 $max_stop = max @cstops;
 	 # save the starts and stops for this event per caller
 	 push @estarts, $min_start;
 	 push @estops, $max_stop;
 	 # save the caller names in the same order
 	 push @ecallers, $caller_id;
        }

      }
      # compute the median event start and stop
      $median_start = &median(@estarts);
      $median_stop = &median(@estops);


      # for each caller's starts and stops, compute the distance to the median
      for(my $i=0; $i<scalar @estarts; $i++) {
        $start_dist = $estarts[$i] - $median_start;
        $stop_dist  = $estops[$i]  - $median_stop; 
        print OUT join "\t", ($event_id,$event_type{$event_id},$event_length{$event_id},$caller_count{$event_id},$ecallers[$i],$start_dist,$stop_dist);
        print OUT "\n";
      }
 
 } #end while event
  close(OUT);
  return;
}


#---------------------------------------------
# Count numbers of SVs/events with predefined annotation types
# Serves as an example that can be modified for different annotations.
# Here, count genes (ENSG), AllRepeats -> found by regular expression match using
# values in @strings_to_match
# For the test file, the overlapped vcf entries represent overlap with 1000 genomes svs.
#---------------------------------------------
sub summarize_annotations {
  my( $data_href, $outfile, $strings_to_match_aref)=@_;
  my %annot_count;

  #initialize so 0 values are printed
  $annot_count{'1000_genomes_SVs'} = 0;
  foreach my $s (@$strings_to_match_aref) {
    $annot_count{$s} = 0;
  }

  while(my ($event_id,$event_info) = each(%$data_href)) {
    if ($event_id =~ /DEL00137688SUR/) { print "processing DEL00137688SUR\n"; }
    my $annotations_href = $event_info->{"annotations"};
    while(my($annot_data_type, $annot_data_values) = each(%$annotations_href)) {
       # fetch the set of annotations for this event
       if ($annot_data_type eq 'overlapped_VCF') {
	 if ($annot_data_values >0) {
	   $annot_count{'1000_genomes_SVs'}++;
	 }
       } elsif ($annot_data_type eq 'overlapped_Annotations') {
	 # does this SV have the annotations represented by the string match
	 foreach my $s (@$strings_to_match_aref) {
	   $match = &find_string_match( (join ' ',@$annot_data_values) ,$s);
	   if ($match == 1) {
	     $annot_count{$s}++;
	   }
	 }
       }
    }
  }

  # print the results

  open(OUT,">$outfile");
  print OUT join "\t",('ANNOT_TYPE', 'SV_COUNT'); print OUT "\n";
  foreach my $s (keys %annot_count) {
    print OUT join "\t",($s, $annot_count{$s}); print OUT "\n";
  }
  close(OUT);

  return;
}

#---------------------------------------------
# Test whether the annotation string passed has the reg ex representing an annotation type
#---------------------------------------------
sub find_string_match {
  my($annot_data_values,$s)=@_;
  my $match=0;

  if ( $annot_data_values =~ /$s/) { return 1; }
  return $match;
}

#---------------------------------------------
# compute median. 
# from http://www.perlmonks.org/index.pl?node_id=474564
# reimplement here to reduce dependencies
#---------------------------------------------
sub median {
  my @vals = sort {$a <=> $b} @_;
  my $len = @vals;
  if($len%2) #odd?
  {
     return $vals[int($len/2)];
  }
  else #even
  {
     return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
  }
}


#---------------------------------------------
# replacement for uniq from List:MoreUtils
# used to reduce dependencies.
#---------------------------------------------
sub uniq {
  my(@ar) = @_;
  my %h;
  foreach $a (@ar) {
    $h{$a} = 1;
   }
   return (keys %h);
}

1;

