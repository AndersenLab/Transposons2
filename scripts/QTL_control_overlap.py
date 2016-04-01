#!/usr/bin/env python
# this script checks if any of the QTL confidence intervals overlap with locations of the TE control genes as cited in literature
# USE: QTL_control_overlap.py
# NOTE: must first transfer "Peak_Table.txt"
import re

lit_genes={}
lit="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/lit_TE_control_genes.txt"
with open(lit, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		lit_genes[line]=0


OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/lit_TE_control_positions.bed", 'w')
infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_gene_positions.gff"
element_info={}
with open(infile, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		chromosome,wb,type_of,start,end,na,orient,na2,info=items[0:9] #actaully bed file, 0-based
		match=re.search("ID=Gene:(\w+\d+);",info)
		element=match.group(1)

		if element in lit_genes.keys():
			start=int(start)-1
			OUT.write("{chromosome}\t{start}\t{end}\t{element}\t.\t{orient}\n".format(**locals()))
OUT.close()

OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/peak_table_pos.bed", 'w')
peak_table="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/Peak_Table.txt"
with open(peak_table, 'r') as IN:
	next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		trait=items[0]
		trait=re.sub('\s+','',trait)
		peak_id=items[1]
		chromosome=items[2]
		left_CI=items[5]
		right_CI=items[6]
		ID=trait+"_"+peak_id
		print trait
		print chromosome
		print peak_id
		print left_CI
		print right_CI
		print ID
		OUT.write("{chromosome}\t{left_CI}\t{right_CI}\t{ID}\t.\t.\n".format(**locals()))
OUT.close()


#check for overlap 

#check for overlap lit genes
cmd="bedtools intersect -wo -a /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/lit_TE_control_positions.bed -b /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/peak_table_pos.bed > /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/QLT_lit_overlap.txt"
from subprocess import Popen, PIPE
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()
#check for overlap piRNA
cmd="bedtools intersect -wo -a /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_piRNA_positions.gff -b /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/peak_table_pos.bed > /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/QLT_piRNA_overlap.txt"
from subprocess import Popen, PIPE
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()