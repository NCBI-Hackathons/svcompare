# SVCompare
Comparison and estimation of breakpoints for SVs across different technologies and resources.
Summaries of SV annotations.

**************************************
## Introduction:

This github represents a collection of the methods that our group developed to analyze, annotate and more accurately predict breakpoints for SVs on the hackathon (08/15-08/17 @ NIH).

The methods developed consists of:

- SURVIVOR_ant: a comparison and annotate method for a given vcf file.
- generate_statistics.pl: provides statics of callsets in which an event is supported.
- run_all: wrapper for generate_statistics.pl

************************************
## Installation
svcompare tool can be in installed at the root system level or at the user local.

1. clone the repository
2. gunzip svcompare.zip
3. place svcompare folder in a desired location

************************************
## Repo content

### test_data: 
- GIAB_H002_081716.vcf.zip: includes the merged calls of GIAB HG002 based on 16 call sets.
- GIAB_H002_081716.vcf_anno.zip: The resulting file from SURVIVOR_ant based on the merged file of GIAB HG002.
- caller_tech.map.txt: table linking callers and technologies to use in statistics generation.
    

************************************
## Citation:
please cite our article if you use parts of our methods or data sets. (will soon be announced) 



