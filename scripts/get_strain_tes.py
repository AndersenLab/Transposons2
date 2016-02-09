#!/usr/bin/env python
import re

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/data_for_figures/CtCp_all_nonredundant.txt"

with open(infile, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromosome=items[0]
		pos=items[1]
		strain=items[6]



		if chromosome == "IV":

			#3794443-3794543
			if int(pos)>= 3794443-1000 and int(pos) <= 3794543+1000:
				if strain=="N2" or strain== "CB4856":
					print chromosome
					print pos
					print strain