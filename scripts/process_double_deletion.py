#!/usr/bin/env python 
# this script takes the output from a transposon caller and splits one call involving 2 or more TEs into individual calls
# applicable with temp absence caller in which a deletion spanning 2 TEs and the region between them is counted as 1 call rather than 2
# reports any reference TE found within the span of the multiple TE deletion range determined by TEMP
# USE: process_double_deletion.py <nonredundant.bed_output_file> <position_file_of_all_tes>
# called on by run_RSVSim.sh

import sys
import re
import os

te_calls=sys.argv[1]
TE_CALLS = open(te_calls, "r")

te_pos=sys.argv[2]

out = "tmp_double_deletion.txt"
OUT = open(out, "w")

single=0
single_line = "tmp_single_line.txt"
SINGLE_LINE = open(single_line, "w")

for line in TE_CALLS:
	items = re.split("[\t]",line)
	transposons = items[3]
	if re.search(",", transposons):
		SINGLE_LINE.write(line)
		single=1
	else:
		OUT.write(line)

SINGLE_LINE.close()
OUT.close()
if single==1:
	os.system("intersectBed -wb -a tmp_single_line.txt -b {te_pos} > only_doubles.txt".format(**locals()))
	os.system("""cat only_doubles.txt | awk '{match($4,"_[a-z]+_temp_[a-z]+_[0-9]+",a)}{print $1"\t"$8"\t"$9"\t"$10a[0]"\t"$5"\t"$6}' > int""")
	os.system("cat int tmp_double_deletion.txt > tmp && mv tmp tmp_double_deletion.txt")