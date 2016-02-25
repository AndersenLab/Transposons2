#!/usr/bin/env python
import re
import statistics


infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full.txt"
OUT_TOTAL = open ("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/outliers_totals.txt", 'w')
OUT_FAMILY = open ("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/outliers_families.txt", 'w')

OUT_TOTAL.write("trait\tstrain\tcount\taverage\tSD\tdifference\n")
OUT_FAMILY.write("trait\tstrain\tcount\taverage\tSD\tdifference\n")

OUT_HEADS=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/heads.txt", 'w')
OUT_HEADS.write("Trait\tStrain\tStrain_Count\tStrain_SD\tMEAN\tSD\n".format(**locals()))
OUT_TAILS=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/tails.txt", 'w')
OUT_TAILS.write("Trait\tStrain\tStrain_Count\tStrain_SD\tMEAN\tSD\n".format(**locals()))

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
	unique_values=sorted(set(values))

	no_unique_values= len(unique_values) # get number of unique values
	if no_unique_values>10:
		process_tails="TRUE"
		slice_head=unique_values[0:10] # the first ten values
		slice_tail=unique_values[-10:] # the last ten values
	else:
		process_tails="FALSE"

	



	for strain_info in strain_totals[i]:
		info=re.split(":", strain_info)
		strainN,count=info
		if count != "NA":
			strain_difference=abs(float(av)-float(count))
			if strain_difference> (2*sd):
				if re.search("total",i):
					OUT_TOTAL.write("{i}\t{strain}\t{count}\t{av}\t{sd}\t{strain_difference}\n".format(**locals()))
				else:
					OUT_FAMILY.write("{i}\t{strain}\t{count}\t{av}\t{sd}\t{strain_difference}\n".format(**locals()))

			if process_tails=="TRUE":
				if float(count) in slice_head:
					OUT_HEADS.write("{i}\t{strainN}\t{count}\t{strain_difference}\t{av}\t{sd}\n".format(**locals()))
				if float(count) in slice_tail:
					OUT_TAILS.write("{i}\t{strainN}\t{count}\t{strain_difference}\t{av}\t{sd}\n".format(**locals()))

		#if

OUT_TOTAL.close()
OUT_FAMILY.close()
OUT_TAILS.close()
OUT_HEADS.close()