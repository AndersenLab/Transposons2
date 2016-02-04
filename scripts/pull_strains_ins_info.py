#!/usr/bin/env python
# this script takes 3-4 letter gene names as input and pulls out the line of info from essentiality_nonredundant_strain_info.txt if
# the gene names matches those that were supplied
# USE: pull_strains_info_ins.py

import re
import sys
from collections import defaultdict

queries=sys.argv[1:]
in_file="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/essentiality_nonredundant_strain_info.txt"
exon_overlap="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/exon_overlap.txt"
OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/genes_of_interest.txt", 'w')


exons=defaultdict(list)

with open(in_file, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		chromosome=items[0]
		start=items[1]
		te=items[3]
		region=items[4]
		transcript_name=items[5]
		print transcript_name
		gene_name=items[6]
		for query in queries:
			if gene_name == query:
				OUT.write(line)
				if region == "exon":
					if gene_name in exons.keys():
						print "check for duplicates"
						sys.exit()
					else:
						exons[gene_name].extend([chromosome,start,te,transcript_name])

OUT.close()


with open(exon_overlap, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromosome_e,wb,region,start_e,end_e,orient,parent,te_location,te_info=items[0:9]
		print(te_info)


		#I	WormBase	exon	5075431	5079303	+	Parent=Transcript:C30F8.7	5075776	I_5075776_Tc1_non-reference__temp_sr_961_63_+_new_CB4851	new	dnatransposon



