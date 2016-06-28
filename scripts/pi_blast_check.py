#!/usr/bin/env python
# this script checks if the piRNA within the CI of the QTL for a given trait blasted to the TE elements/family for that trait
#USE: pi_blast_check.py

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
# TC3 control
# make blast database of TC3 sequences 
cmd="/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/makeblastdb -in TC3_TE_seqs.fasta  -dbtype nucl -out TC3_database".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

cmd="/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/blastn -db TC3_database -query /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/known_Tc3_piRNAs.fasta -evalue .5 -word_size 5 -outfmt '6 qseqid sseqid pident qlen length mismatch gapopen evalue bitscore qstart qend btop' -max_target_seqs 100 -out TC3_blast.txt -num_threads 10".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

cmd="cat intervalPIs_* > intervalPIs_all.txt"
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()

########################
# check traits of interest
trait_list=[]
with open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/TOI.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		trait_list.append(line)

OUT=open("TOIs_matching_intervalPIs.txt",'w')  # file for TOIs with piRNAs found in the CI of the QTL blasting to their repective TEs  
with open("intervalPIs_all.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		trait,mismatch,found,possible=items[0:]
		found=int(found)
		if trait in trait_list:
			if found !=0:
				OUT.write(line +'\n')
OUT.close()


