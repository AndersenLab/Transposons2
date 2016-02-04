#!/usr/bin/env python
# this script adds the numbers and names of strains associated with each event in essentiality_nonredundant_GO_only_new.txt
# USE: insertion_distribution.py

import re
from collections import defaultdict

all_nonredundant = "/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/CtCp_all_nonredundant.txt"
ess = "/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/essentiality_nonredundant_GO_only_new.txt"

OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/essentiality_nonredundant_strain_info.txt",'w')
OUT.write("Chromosome\tTE_start\tMethod\tTE\tRegion\tTranscript_Name\tGene_Name\tBiotype\tPhenotype\tGO_Annotation\tNo_Strains\tStrains\n")


strains=defaultdict(list)
with open(all_nonredundant, 'r') as IN:
	for line in IN:
		line = line.rstrip('\n')
		items = re.split("[\t]",line)
		chromosome,start,start2,te_info,orient,method,strain,te_class=items[0:len(items)]
		match = re.search(".*_\d+_(.*)_(\w+-)?reference",te_info)
		family = match.group(1)
		ID="{chromosome}_{start}_{method}_{family}".format(**locals())
		strains[ID].append(strain)


with open(ess, 'r') as IN:
	next(IN) #skip header
	for line in IN:
		line = line.rstrip('\n')
		items = re.split("[\t]",line)
		chromosome,start,method,family,region,transcript,gene,biotype,phenp,GO=items[0:len(items)]
		ins_ID="{chromosome}_{start}_{method}_{family}".format(**locals())
		strain_list=strains[ins_ID]
		no_strains=len(strain_list)
		strain_list=','.join(strain_list)	
		OUT.write("{line}\t{no_strains}\t{strain_list}\n".format(**locals()))

OUT.close()


