# This script has routines to calculate various statistics from the survivor output file.
# Will be combined into a single file with other contributors' routines when done.

use File::Basename;

# extract the executable directory so can find the data structures code

my $data_struc_code = 'data_structures.pl';
my $exe_dirname = dirname(__FILE__);
require "$exe_dirname/$data_struc_code";

usage() unless(@ARGV > 1);
my ($vcf_file, $tech_file, $outfile_prefix) = ($ARGV[0], $ARGV[1], $ARGV[2]);

my $techniques = map_technologies($tech_file);

# Load the data
my $data_by_event_href = generate_data_by_event($vcf_file,$tech_file);
my $data_by_caller_href = generate_data_by_caller($vcf_file,$tech_file);


# --------------------------------------------------------------
# set up the file names for the output files before running each analysis.
# filenames are the prefix input on the command line, with an analysis-specific suffix.
# the file extensions are stored in the suffix hash declared here.
# --------------------------------------------------------------

my %suffix;  


# --------------------------------------------------------------
# Calculate the total number of variants called in each callset (ie by a particular calling algorithm).
# If a call set has more than one variant within the range of the event, count as one variant.
# --------------------------------------------------------------

$suffix{'callset_total_count'} = '.callset_total_count.tsv';
&callset_total_count( $data_by_caller_href, $outfile_prefix.$suffix{'callset_total_count'});


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
# total number of variants in each callset
# ignore variant type
#------------------------------------------
sub callset_total_count{
  my( $data_href, $outfile) = @_;

  # prepare the output file and header
  open(OUT,">$outfile");
  my @header = qw(CALLER N_VARIANTS);
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
