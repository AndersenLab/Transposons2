#!/usr/bin/env python

import re

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full.txt"
AA={}
RR={}

with open(infile, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		ID=items[0]
		if ID != "coverage":
			items2=re.split("_TRANS_",ID)
			method=items2[0]
			family=items2[1]
			counts=items[1:len(items)]
			counts = [int(x) for x in counts if x != "NA"]
			max_count = max(counts)
			if method == "absent":
				AA[family]=max_count
			elif method == "reference":
				RR[family]=max_count
print len(AA)
print len(RR)
#test= set(AA.keys()) & set (RR.keys())
#print len(test)

AA={key:value for key,value in AA.items() if value <=1}
RR={key:value for key,value in RR.items() if value <=1}


PA= set(AA.keys()) & set(RR.keys())
