#!/usr/bin/env python

import re

in_one="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full.txt"
in_two="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/kin_matrix_full.txt"
in_three="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_Full_Results_Activity.txt"

OUT_ONE=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full_id.txt", 'w')
OUT_TWO=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/kin_matrix_full_id.txt", 'w')
OUT_THREE=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_Full_Results_Activity_id.txt", 'w')

MASTER_ONE=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/key_T_kin_C_matrix_full_id.txt", 'w')
MASTER_TWO=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/key_kin_matrix_full_id.txt", 'w')
MASTER_THREE=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/key_T_Full_Results_Activity_id.txt", 'w')

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

with open(in_two, 'r') as IN:
	header=next(IN)
	OUT_TWO.write(header)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		TE=items[0]
		MASTER_TWO.write(TE + '\t' + "T" + str(count) + '\n')
		OUT_TWO.write("T"+str(count))
		for i in items[1:len(items)]:
			OUT_TWO.write('\t' + i)
		OUT_TWO.write('\n')
		count +=1
OUT_TWO.close()
MASTER_TWO.close()


with open(in_three, 'r') as IN:
	header=next(IN)
	OUT_THREE.write(header)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		TE=items[0]
		MASTER_THREE.write(TE + '\t' + "T" + str(count) + '\n')
		OUT_THREE.write("T"+str(count))
		for i in items[1:len(items)]:
			OUT_THREE.write('\t' + i)
		OUT_THREE.write('\n')
		count +=1
OUT_THREE.close()
MASTER_THREE.close()
