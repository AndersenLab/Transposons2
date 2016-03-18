#!/usr/bin/env python
# this script calculates the coverage of a sample at each interval given in the interval file
# USE: calculate_cov_at_intervals.py <sample name>

import re
import sys
from subprocess import Popen, PIPE

OUT=open("sample_coverages_and_positions.txt", 'aw')
bam_dir="/lscr2/andersenlab/dec211/WI/bam"
interval_file="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/TE_matrix/coverage_intervals.txt"
sample=sys.argv[1]

cov_dict={}

with open(interval_file,'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		key,av_pos,chromosome,start,end=items[0:5]
		ID="_".join(items[0:5])
		result, err = Popen(["""samtools depth {bam_dir}/{sample}.bam -r {chromosome}:{start}-{end}|datamash mean 3""".format(**locals())],stdout=PIPE, stderr=PIPE, shell=True).communicate()
		if result:
			coverage=round(float(result),2)
		else:
			coverage=0
		cov_dict[ID]=coverage

OUT.write(sample)
for key in sorted(cov_dict.keys()):
	value=cov_dict[key]
	OUT.write("\t{key}:{value}".format(**locals()))
OUT.write("\n")
OUT.close()