# Usage Example:
# python dbvar_vcf_to_sites.py < estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.noheader > dbvar_estd219_sites.vcf

#This scripts takes a dbVar VCF file consisting of sample-level variants  callsfrom many samples  and
#produces a unique set of variants found in any sample in the input  as an output VCF


import sys

result_dict = {}

for line in sys.stdin:
    row = line.rstrip().split("\t")
    info_list = row[7].split(";")
    info_dict = {}
    for x in info_list:
       if "=" in x:
        key, val = x.split("=")
        info_dict[key] = val

    info_field = "SVTYPE=" + info_dict["SVTYPE"] + ";" + "END=" + info_dict["END"] + ";" \
         "REGION=" + info_dict["REGION"] 

    result = row[0] + "\t" + row[1] + "\t" + info_dict["REGION"] + "\t" + row[3] + "\t" \
                   + row[4] + "\t.\t.\t" +  info_field

    result_dict[result] = 1              

for item in result_dict.keys():
	print item








