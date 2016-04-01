#!/usr/bin/env python

import re

in_one="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full_reduced.txt"
OUT_ONE=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full_id_reduced.txt", 'w')
MASTER_ONE=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/key_T_kin_C_matrix_full_id_reduced.txt", 'w')


count=1


with open(in_one, 'r') as IN:
	header=next(IN)
	OUT_ONE.write(header)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		TE=items[0]
		MASTER_ONE.write(TE + '\t' + "T" + str(count) + '\n')
		OUT_ONE.write("T"+str(count))
		for i in items[1:len(items)]:
			OUT_ONE.write('\t' + i)
		OUT_ONE.write('\n')
		count +=1
	
OUT_ONE.close()
MASTER_ONE.close()