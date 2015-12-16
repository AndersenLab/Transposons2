#!/usr/bin/env python 
# this script checks to see if there are any positions at which telocate called a reference transposon whereas temp made an absence call
# in such a case, the 2 callers are contradictory
# outputs:
# 1) files containing all contradictory calls for the strains (bed info)
# 2) summary files in order of: strain, total number of contradictory calls, total reference calls, and total absence calls
# USE: temp_telocate_compatibility.py <data directory> <sample_list file>

import sys
import re
from subprocess import Popen, PIPE

data_dir= sys.argv[1]#"/lscr2/andersenlab/kml436/git_repos2/Transposons2/data"
sample_list= sys.argv[2]#"/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/full_sample_list.txt"

CONTRADICTORY_CALLS=open("contradictory_calls.txt", "w")
CONTRADICTORY_CALLS_SUMMARY=open("contradictory_calls_summary.txt", "w")

SAMPLE_LIST=open(sample_list, "r")
for line in SAMPLE_LIST:
	strain=line.rstrip('\n')
	telocate="{data_dir}/{strain}/final_results/org_{strain}_telocate_nonredundant.bed_org".format(**locals())
	temp="{data_dir}/{strain}/final_results/org_{strain}_temp_absence_nonredundant.bed_org".format(**locals())

	TELOCATE=open(telocate, "r")
	TEMP=open(temp, "r")
	reference_count=0
	absence_count=0
	contradictory_count=0
	reference_calls={}

	# put all positions of the reference counts into a dictionary
	for line in TELOCATE:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		chromosome=items[0]
		start=items[1]
		transposon=items[3]
		reference_calls["{chromosome}_{start}".format(**locals())]=line
		#count total number of reference calls with telocate
		reference_count+=1
	TELOCATE.close()

	# check if any positions of the absence calls are in the dictionary of the reference calls
	for line in TEMP:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		chromosome=items[0]
		start=items[1]
		transposon=items[3]
		temp_support=items[4]
		orientation=items[5]
		position="{chromosome}_{start}".format(**locals())
		#count total number of absence calls with temp
		absence_count+=1

		if position in reference_calls.keys():
			value = reference_calls[position]
			contradictory_count+=1

			CONTRADICTORY_CALLS.write("{strain}\t{value}\t{transposon}\t{temp_support}\t{orientation}\n".format(**locals()))
	CONTRADICTORY_CALLS_SUMMARY.write("{strain}\t{contradictory_count}\t{reference_count}\t{absence_count}\n".format(**locals()))

	TEMP.close()


CONTRADICTORY_CALLS.close()
CONTRADICTORY_CALLS_SUMMARY.close()







