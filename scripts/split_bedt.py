#!/usr/bin/env python
# this script reformats the TE absence, insertion, reference files that are output by kin_mean.py
# USE: split_bedt.py
import re

#X_9470448_9470448_Tc1A_reference_temp_rp_103_152_+_absent_JU782
#III_11497122_11497122_PALTTAA3_CE_non-reference__temp_rp_4521_9_-_new_EG4347


methods=('absences','insertions','references')
for method in methods:
	OUT=open("bedfile_{method}.bed".format(**locals()),'w')
	with open("T_Samples_{method}_bedt.txt".format(**locals()), 'r') as IN:
		next(IN) #skip first line
		for line in IN:
			items = re.split("[\t]", line)
			info = items[0]
			match = re.search("([A-Z]+)_(\d+)_(\d+)_(([A-Za-z\d+_-]+)_((\w+-)?reference)_+([a-z]+)_([a-z]+)_(\d+)_([\d\.]+)_([\+-])_(\w+)_(\w+))", info) #([A-Za-z\d+_])_((\w+-)?reference)\w+_\d+_\d+
			chromosome = match.group(1)
			start = match.group(2)
			start2 = match.group(3)
			TE = match.group(4)
			RS = match.group(11)
			orient = match.group(12)
			sample = match.group(14)

			OUT.write("{chromosome}\t{start}\t{start2}\t{chromosome}_{start}_{start2}_{TE}\t{RS}\t{orient}\n".format(**locals()))





	OUT.close()