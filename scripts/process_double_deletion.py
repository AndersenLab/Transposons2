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
		single=0

SINGLE_LINE.close()
OUT.close()
if single==1:
	os.system("intersectBed -wb -a tmp_single_line.txt -b {te_pos} > only_doubles.txt".format(**locals()))
	os.system("""cat only_doubles.txt | awk '{match($4,"_[a-z]+_temp_[a-z]+_[0-9]+",a)}{print $1"\t"$8"\t"$9"\t"$10a[0]"\t"$5"\t"$6}' > int""")
	os.system("cat int tmp_double_deletion.txt > tmp && mv tmp tmp_double_deletion.txt")







		#match = re.search("(.*)_([a-zA-Z]+_[\d]+)", transposons)
		#identifier = match.group(2) # match the individual transposon identifier (usually rp_#)
		#double_transposon = re.split(",", transposons)
		#for i in range(0,len(double_transposon)):
		#	TE=double_transposon[i]
	#		TE = re.sub(r'_[a-zA-Z]+_[\d]+',"",TE) # get rid of identifier if present
	#		OUT.write('\t'.join(items[0:3]))
#			OUT.write("\t{TE}_{identifier}.{i}\t".format(**locals()))
#			OUT.write('\t'.join(items[4:6]))
#	else:
#		OUT.write(line)

#
#TE_CALLS.close()
#OUT.close()


##alternative version:
#te_calls=sys.argv[1]
#TE_CALLS = open(te_calls, "r")

#out = "tmp_double_deletion.txt"
#OUT = open(out, "w")
#for line in TE_CALLS:
#	items = re.split("[\t]",line)
#	transposons = items[3]
#	if re.search(",", transposons):
#		match = re.search("(.*)_([a-zA-Z]+_[\d]+)", transposons)
#		identifier = match.group(2) # match the individual transposon identifier (usually rp_#)
#		double_transposon = re.split(",", transposons)
#		for i in range(0,len(double_transposon)):
#			TE=double_transposon[i]
#			TE = re.sub(r'_[a-zA-Z]+_[\d]+',"",TE) # get rid of identifier if present
#			OUT.write('\t'.join(items[0:3]))
#			OUT.write("\t{TE}_{identifier}.{i}\t".format(**locals()))
#			OUT.write('\t'.join(items[4:6]))
#	else:
#		OUT.write(line)
#
#
#TE_CALLS.close()
#OUT.close()