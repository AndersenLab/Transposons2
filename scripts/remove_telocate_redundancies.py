#!/usr/bin/env python 
# this script takes a telocate nonredundant output file and removes redundant calls
# USE: remove_telocate_redundancies.py <closest.txt>
import sys
import re
import os
from subprocess import Popen, PIPE

in_file = sys.argv[1]
IN_FILE= open(in_file, "r")
OUT_FILE= open("closest_prepped.txt", "w")

#in_file2 = sys.argv[2]
#IN_FILE2= open(in_file2, "r")

#count_T=0
#count_F=0

TEs={}
full_TE_info={}
for line in IN_FILE:
	items= re.split("[\t]",line)
	transposon_ID=items[9]
	distance=items[10]
	if transposon_ID in TEs.keys():
		prev_distance=TEs[transposon_ID]
		if int(distance) > int(prev_distance): #take farther away one to test worst case
			TEs[transposon_ID]=distance
			full_TE_info[transposon_ID]=line

	else:
		TEs[transposon_ID]=distance
		full_TE_info[transposon_ID]=line
		
IN_FILE.close()
for key in full_TE_info.keys():
	value=full_TE_info[key]
	OUT_FILE.write(value)


IN_FILE.close()
OUT_FILE.close