#!/bin/bash

# Command line usage for generating statistics for the test data
# perl generate_statistics.pl <annotated_vcf> <table_mapping_callers_to_technology> <output_filename_prefix> 

# To generate statistics on the test data, unzip the vcf and then run:
perl generate_statistics.pl ../../test_data/GIAB_H002_081716_anno.vcf ../../test_data/caller_tech.map.txt giab
