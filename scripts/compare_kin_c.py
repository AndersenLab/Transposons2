#!/usr/bin/env python

import re
import sys
from collections import defaultdict

old="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results_feb/final_results/data_for_figures/T_kin_C_matrix_full.txt"
new="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/data_for_figures/T_kin_C_matrix_full.txt"


comp=defaultdict(list)
with open(old, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)

		trait=items[0]
		for index,value in enumerate(items[1:153]):
			tt=(str(index),str(value))
			tog=":".join(tt)
			comp[trait].append(tog)

OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/compare_differences.txt", 'w')
OUT1=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/compare_values.txt", 'w')
with open(new, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		trait=items[0]
		old_dat=comp[trait]
		OUT.write(trait)
		OUT1.write(trait)
		for index,value in enumerate(items[1:153]):
			old_info=old_dat[index]
			old_items=re.split(":",old_info)
			old_index=old_items[0]
			old_value=old_items[1]
			if value != "NA" and old_value !="NA":
				value_diff=float(value)-float(old_value)
				OUT.write('\t' + str(value_diff))
				OUT1.write('\t' +str(old_value) + "/" + str(value))
			elif value == "NA" and old_value =="NA":
				value_diff=0
				OUT.write('\t' + str(value_diff))
				OUT1.write('\t' +str(old_value) + "/" + str(value))

			else:
				value_diff="X"
				OUT.write('\t' + str(value_diff))
				OUT1.write('\t' +str(old_value) + "/" + str(value))
			

		OUT.write('\n')
		OUT1.write('\n')

OUT.close()
OUT1.close()
