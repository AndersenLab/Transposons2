#!/usr/bin/env python 
# this script takes a WB element->family conversion file, a gff file (or any file) with element names, and the name of an output file
# it replaces the element names in the gff file with the specified family name in the conversion file  
# USE:correct_names_2.py <element->family conversion file> <gff_file>(element names) <name_of_outpus_file>
import sys
import re
import os
file1 = sys.argv [1]
round2_names_for_gff={}
INPUT_WB_FAMILIES = open(file1,"r")
for line in INPUT_WB_FAMILIES:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	wb_name = items[0]
	wb_family = items[1]
	round2_names_for_gff[wb_name] = wb_family

INPUT_WB_FAMILIES.close()


gff_file = sys.argv[2]
GFF= open(gff_file, "r")

out_file = sys.argv[3]
R2_GFF = open(out_file, "w")
for line in GFF:
	line = line.rstrip('\n')
	for key in round2_names_for_gff.keys(): #replace with more efficient method later
		value = round2_names_for_gff[key]
		if re.search(key, line):
			line = line.replace(key, value)
	R2_GFF.write("{line}\n".format(**locals()))
R2_GFF.close()
