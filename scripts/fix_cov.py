#!/usr/bin/env python
# this script output coerage per strain in a file that can be merged with T_kin_C_matrix_full.txt
# USE: cov_fix.py

import re

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/T_kin_C_matrix_full.txt"
cov_file="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/coverage_and_te_counts.txt"
outfile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/T_coverage.txt"


coverages={}
with open(cov_file, 'r') as IN:
	for line in IN:
		line = line.rstrip('\n')
		items = re.split("[\t]",line)
		sample,absence,insertion,reference,coverage=items[0:len(items)]
		coverages[sample]=coverage


strains=list()
with open(infile, 'r') as IN:
	header = next(IN)
	header = header .rstrip('\n')
	items = re.split("[\t]", header )
	for i in items[1:len(items)]:
		strains.append(i)


OUT=open(outfile, 'w')
OUT.write("trait")

for i in strains:
	OUT.write("\t" + i)
OUT.write('\ncoverage')

for i in strains:
	coverage=coverages[i]
	OUT.write("\t" + coverage)

OUT.write("\n")
OUT.close()