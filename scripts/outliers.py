#!/usr/bin/env python
import re
import statistics


infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full.txt"
OUT_TOTAL = open ("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/outliers_totals.txt", 'w')
OUT_FAMILY = open ("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/outliers_families.txt", 'w')

OUT_TOTAL.write("trait\tstrain\tcount\taverage\tSD\tdifference\n")
OUT_FAMILY.write("trait\tstrain\tcount\taverage\tSD\tdifference\n")

from collections import defaultdict
from collections import Counter
totals=defaultdict(list)
strain_totals=defaultdict(list)

with open(infile) as IN:
	header=next(IN)
	header=header.rstrip('\n')
	strains=re.split("[\t]",header)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		trait=items[0]
		#if re.search("total",trait): 
		#totals[trait].extend(items[1:(len(items)+1)])
		for index,count in enumerate(items[1:len(items)],start=1):
			strain=strains[index]
			strain_totals[trait].append(strain +":" + count)
			if count != "NA":
				totals[trait].append(count)

for i,values in totals.items():
	values=[float(x) for x in values]
	av=statistics.mean(values)
	sd=statistics.pstdev(values)
	#test_difference=abs(av-sd)
	#print test_difference
	for strain_info in strain_totals[i]:
		info=re.split(":", strain_info)
		strain,count=info
		if count != "NA":
			strain_difference=abs(float(av)-float(count))
			if strain_difference>sd:
				if re.search("total",i):
					OUT_TOTAL.write("{i}\t{strain}\t{count}\t{av}\t{sd}\t{strain_difference}\n".format(**locals()))
				else:
					OUT_FAMILY.write("{i}\t{strain}\t{count}\t{av}\t{sd}\t{strain_difference}\n".format(**locals()))

OUT_TOTAL.close()
OUT_FAMILY.close()