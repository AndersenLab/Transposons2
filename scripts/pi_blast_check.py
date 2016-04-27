#!/usr/bin/env python
# this script checks if the piRNA within the CI of the QTL for a given trait aligned to the TE elements/family for that trait
# also checks if knwon Tc3 seqs from literature align to the Tc3 sequences and  to the genome
#USE: pi_align_check.py

import re
import pickle
from collections import defaultdict
from Bio import SeqIO
from subprocess import Popen, PIPE


def search_blasted(mismatch):
	OUT=open("intervalPIs_{mismatch}.txt".format(**locals()),'w') 
	qtl_overlap="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/QLT_piRNA_overlap.txt"

	in_interval=defaultdict(list)
	blasted_in_interval=defaultdict(list)

	with open("strict_blasts_{mismatch}.txt".format(**locals()), "rb") as fp: #load information on which piRNAs blasted to which TEs
		blasts = pickle.load(fp)

	with open(qtl_overlap, 'r') as IN:
		next(IN)
		for line in IN:
			line = line.rstrip('\n')
			items=re.split('\t', line)
			
			TE=items[0]
			transcript=items[8]
			
			match=re.search('(.*)\(.*\)',TE)
			family=match.group(1)

			if family in blasts.keys():
				check=True
				pis_blasted=blasts[family]
				if transcript in pis_blasted:
					blasted_in_interval[TE].append(transcript)
			else:
				check=False
				print family


			in_interval[TE].append(transcript)

	for i in in_interval.keys():
		no_in_interval=len(in_interval[i])

		if check==True:
			if i in blasted_in_interval.keys():
				no_blasted_in_interval=len(blasted_in_interval[i])
			else:
				no_blasted_in_interval=0
		else:
			no_blasted_in_interval=0
		OUT.write("{i}\t{mismatch}\t{no_blasted_in_interval}\t{no_in_interval}\n".format(**locals()))
	OUT.close()

search_blasted('zero')
search_blasted('one')
search_blasted('two')
search_blasted('three')


########################


