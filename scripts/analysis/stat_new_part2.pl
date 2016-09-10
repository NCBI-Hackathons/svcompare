# This script has routines to calculate various statistics from the survivor output file.
# Will be combined into a single file with other contributors' routines when done.

use File::Basename;
use List::Util qw(min max);
use List::MoreUtils qw(apply uniq);


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


# --------------------------------------------------------------
# set up the file names for the output files before running each analysis.
# filenames are the prefix input on the command line, with an analysis-specific suffix.
# the file extensions are stored in the suffix hash declared here.
# --------------------------------------------------------------

my %suffix;  


# --------------------------------------------------------------
# Calculate the total number of events called in each callset (ie by a particular calling algorithm).
# If a call set has more than one variant within the range of the event, count as one variant.
# --------------------------------------------------------------

$suffix{'callset_event_count'} = '.callset_event_count.tsv';
&callset_event_count( $data_by_caller_href, $outfile_prefix.$suffix{'callset_event_count'});


# --------------------------------------------------------------
# Count the number of variants of each type in each callset.
# If there are multiple variants per event, count them all.
# Print 0's for types not present in a call set.
# --------------------------------------------------------------

$suffix{'callset_total_count_by_type'} = '.callset_total_count_by_type.tsv';
callset_total_count_by_type( $data_by_caller_href, $outfile_prefix.$suffix{'callset_total_count_by_type'});


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
$suffix{'breakpoint_variation_id_counts'} = '.breakpoint_variation.id_counts.tsv';
&breakpoint_variation( $data_by_event_href, $min_callers, $outfile_prefix.$suffix{'breakpoint_variation'},  $outfile_prefix.$suffix{'breakpoint_variation_id_counts'} );



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
# For events (SVs/ rows) 
#  - supported by at least num_callers callers
#  - for which all calls in any caller are the same type
# For each event, compute the distance of the start and stop to the median of the calls
# Two output files: 
# - $outfile: details per event per caller
# - $count_outfile: summary statistics from event filtering
#------------------------------------------
sub breakpoint_variation {
  my( $data_href, $min_callers, $outfile,  $count_outfile )=@_;
  my %event_type;
  my %caller_count;

  # counters
  my $num_events_low_call_support=0;
  my $num_events_multiple_types_per_caller=0;
  my $num_events_multiple_types_across_callers=0;
  my $num_events_in=0;
  my $num_events_out=0;

  # set up output files
  open(OUT,">$outfile");
  open(COUT,">$count_outfile");
  my @header = qw(EVENT_ID TYPE LENGTH NUM_SUPPORTING_CALLERS CALLER START_DIST_TO_MEDIAN STOP_DIST_TO_MEDIAN);
  print OUT join "\t",@header; print OUT "\n";

  # loop over each row=event
  while(my ($event_id,$event_info) = each(%$data_href)) {
    $caller_count{$event_id} = $event_info->{'caller_count'};
    $event_length{$event_id} = $event_info->{'stop'} - $event_info->{'start'} + 1;
    $num_events_in++;

    # is this event supported by at least min_caller callers
    if ($caller_count{$event_id} < $min_callers) {
      $num_events_low_call_support++;
      $failed{$event_id} = 1;
      next;
    }

    # are all the svs for this event of the same type
    my $callers_href = $event_info->{"callers"};
    my $events_href;
    my $types_aref;
    while(my($caller_id, $regions_href) = each(%$callers_href)) {
       if (defined $failed{$event_id} ) { next; }
       $events_href = $regions_href->{"events"};
       $types_aref = $events_href->{"types"};
       @uniq_types = uniq(@{$types_aref});
       # check if multiple svs called by this one caller for this event are all the same type
       if ( scalar @uniq_types >1 ) {
         $num_events_multiple_types_per_caller++;
         $failed{$event_id} = 1;
         next;
       }

       # the first time see an event, assign the type
       if (!defined $event_type{$event_id}) { $event_type{$event_id} = $types_aref->[0]; }

       # skip events with a different type in any caller
       # already know that this event in this caller only has a single type
       if ($types_aref->[0] ne $event_type{$event_id}) {
         $num_events_multiple_types_across_callers++;
         $failed{$event_id} = 1;
	 next;
       }
    }

    # This event passes the filters:
    # minimum number of supporting callers, and only a single type both within and across callers


    # Calculate the median start and stop across the callers for this event.
    # For each caller, calculate the abs distance for each start and stop to the median.
    # If there are multiple events for one caller, use the min start and max stop.
    # Assume the coordinates are all on the same chromosome.

    # collect the starts and stops to use for the median calculation
    # make sure this is not an event to skip
    if (defined $failed{$event_id}) { next; }
    $num_events_out++;
    my @estarts;  my @estops; my @ecallers;
    while(my($caller_id, $regions_href) = each(%$callers_href)) {

       $events_href = $regions_href->{"events"};
       my $coords_href = $events_href->{"coords"};
       while(my ($type, $regions_aref) = each %$coords_href) {
	 # find the min start and max stop for the this event for this caller
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
     $median_start = &median(\@estarts);
     $median_stop = &median(\@estops);

     # for each caller's starts and stops, compute the distance to the median
     for(my $i=0; $i<scalar @estarts; $i++) {
       $start_dist = $estarts[$i] - $median_start;
       $stop_dist = $estops[$i] - $median_stop;
       print OUT join "\t", ($event_id,$event_type{$event_id},$event_length{$event_id},$caller_count{$event_id},$ecallers[$i],$start_dist,$stop_dist);
       print OUT "\n";
     }
   }

  # print run statistics
  print COUT join "\t",('INPUT_EVENTS', $num_events_in); print COUT "\n";
  print COUT join "\t",('LOW_CALL_SUPPORT', $num_events_low_call_support); print COUT "\n";
  print COUT join "\t",('MULTIPLE_TYPES_PER_CALLER', $num_events_multiple_types_per_caller); print COUT "\n";
  print COUT join "\t",('MULTIPLE_TYPES_ACROSS_CALLERS', $num_events_multiple_types_across_callers); print COUT "\n";
  print COUT join "\t",('TOTAL_FAILED_EVENTS', scalar keys %failed); print COUT "\n";
  print COUT join "\t",('OUTPUT_EVENTS', $num_events_out); print COUT "\n";

  # clean up
  close(OUT);
  close(COUT);
  return;
}


#---------------------------------------------
# compute median. pass array by reference
# from http://www.perlmonks.org/index.pl?node_id=474564
#---------------------------------------------
sub median { $_[0]->[ @{$_[0]} / 2 ] }
