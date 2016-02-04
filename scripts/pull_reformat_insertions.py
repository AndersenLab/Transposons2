#!/usr/bin/env python

import re

infile ="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/kin_matrix_full.txt"
outfile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/insertion_positions.txt"

OUT=open(outfile, 'w')

with open(infile, 'r') as IN:
	next(IN)
	for line in IN:	
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		TE=items[0]
		if re.search("_non-reference",TE):
			match=re.search("(\w+)_(\d+)_(.*)_non-reference",TE)
			chromosome=match.group(1)
			position=match.group(2)
			family=match.group(3)
			OUT.write(chromosome + '\t' + position + '\t' + family + '\n')


OUT.close()