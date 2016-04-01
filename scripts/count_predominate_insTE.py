#!/usr/bin/env python
# this script calculates the prodominate TE family isnertions across each strain and cumulatively
# USE: count_predominate.py
import re
from collections import defaultdict
from collections import Counter
infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full_reduced.txt"

predominate=defaultdict(list)
strains={}
with open(infile, 'r') as IN:
	header=next(IN)
	header=header.rstrip('\n')
	header=re.split('\t',header)
	for x,y in enumerate(header[1:],start=1):
		x=int(x)
		strains[x]=y
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		trait=items[0]
		match=re.search("ONE_new_TRANS_(.*)",trait)
		if match and not re.search("total",trait):
			family=match.group(1)
			for i,value in enumerate(items[1:],start=1):
				if value != "NA":
					value=int(value)
					strain=strains[int(i)]
					for _ in range(value): #append the family names the same number of times as the value
						predominate[strain].append(family)
full_predominate=Counter()
for strain in sorted(predominate.keys()):
	predominate_count=Counter(predominate[strain])
	first=predominate_count.most_common(1)
	print strain
	print first
	print '\n'
	full_predominate=full_predominate + predominate_count


print full_predominate
