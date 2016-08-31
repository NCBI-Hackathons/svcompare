require 'data_structures.pl'

usage() unless(@ARGV > 1);
my ($vcf_file, $tech_file, $outfile_prefix) = ($ARGV[0], $ARGV[1], $ARGV[2]);


my $techniques = map_technologies($tech_file);

my @headers = ();

my $data = load_data($vcf_file);
my $statistics = create_stats();

# set up the file names for the output files before running each analysis.
# filenames are the prefix input on the command line, with an analysis-specific suffix.
my %suffix;  

$suffix{'breakpoint_variation'} = '.breakpoint_variation.tsv';
$suffix{'breakpoint_variation_id_counts'} = '.breakpoint_variation.id_counts.tsv';
breakpoint_variation( $outfile_prefix.$suffix{'breakpoint_variation'},  $outfile_prefix.$suffix{'breakpoint_variation_id_counts'} );

$suffix{'total_callset_count'} = '.total_callset_count.tsv';
total_callset_count( $outfile_prefix.$suffix{'total_callset_count'});

$suffix{'total_callset_count_by_type'} = '.total_callset_count_by_type.tsv';
total_callset_count_by_type( $outfile_prefix.$suffix{'total_callset_count_by_type'});
# callsets_per_event();


exit 0;
