#!/usr/bin/env python 
# this script takes the gff of the WB elements and repalces the element names with the appropriate family names
# USE:correct_names.py in the appropriate directory
import sys
import re
import os

round2_names_for_gff={}
INPUT_WB_FAMILIES = open("round2_WB_familes.txt","r")
for line in INPUT_WB_FAMILIES:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	wb_name = items[0]
	wb_family = items[1]
	round2_names_for_gff[wb_name] = wb_family

INPUT_WB_FAMILIES.close()

gff_file="/lscr2/andersenlab/kml436/sv_files/TEMP_TELOCATE_RETROSEQ_sim_files/WB_all_dups_renamed.gff"
GFF= open(gff_file, "r")
R2_GFF = open("round2_WB_tes.gff", "w")
for line in GFF:
	line = line.rstrip('\n')
	for key in round2_names_for_gff.keys(): #replace with more efficient method later
		value = round2_names_for_gff[key]
		if re.search(key, line):
			line = line.replace(key, value)
	R2_GFF.write("{line}\n".format(**locals()))
R2_GFF.close()