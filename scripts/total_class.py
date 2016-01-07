#!/usr/bin/env python
# this script counts the total number of dnatransposon,retrotransposon, and transposons of unknown classification found for each strain
# USE: total_class.py
import re
from collections import defaultdict
from collections import Counter
totals=defaultdict(list)
infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/CtCp_all_nonredundant.txt"
strain_matrix="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full.txt"
strain_list=list()

with open(infile) as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]',line)
		method=items[5]
		strain=items[6]
		classTE=items[7]
		ID=method+"_TRANS_total_"+ classTE +"_C"
		totals[ID].append(strain)


# read in strains
with open(strain_matrix) as IN:
	header=next(IN)
	header=header.rstrip('\n')
	strains=re.split("[\t]",header)
	strains=strains[1:len(strains)+1]
	for i in strains:
		strain_list.append(i)



OUT=open("total_classes.txt","w")
for i in totals.keys():
	countedTEs=Counter(totals[i])
	OUT.write(i)
	for strain in strain_list:
		if strain in countedTEs.keys():
			te_count=countedTEs[strain]
			OUT.write("\t" + str(te_count))
		else:
			OUT.write("\t0")

	OUT.write('\n')

OUT.close()


