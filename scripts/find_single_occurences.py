#!/usr/bin/env python

import re
from collections import Counter

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/kin_matrix_full.txt"
outfile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/singles.txt"
OUT=open(outfile, 'w')

occurences=list()

#nons=list()

with open(infile, 'r') as IN:
	next(IN)
	for line in IN:
		line = line.rstrip('\n')
		items = re.split("[\t]", line)
		ID = items[0]
		if not re.search("_new_",ID):
			match=re.search('\w+_\d+_(.*)',ID)
			family=match.group(1)
			occurences.append(family)
		#else:
		#	match = re.search("([A-Z]+)_(\d+)_(([A-Za-z\d+_-]+))_((\w+-)?reference)", ID) #([A-Za-z\d+_])_((\w+-)?reference)\w+_\d+_\d+
		#	te = match.group(4)
		#	nons.append(te)



occurence_count=Counter(occurences)
#nons_count=Counter(nons)
print len(occurence_count)
#print len(nons_count)
#test= list(set(occurence_count.keys() + nons_count.keys()))
#test= set(test)
#print len(test)

singles={family:count for family,count in occurence_count.items() if count <=1}

OUT.write('family\n')
for i in singles.keys():
	OUT.write('absent_TRANS_' + i + '_C' + '\n')
	OUT.write('reference_TRANS_' + i + '_C' + '\n')

OUT.close()