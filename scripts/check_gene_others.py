#!/usr/bin/env python

import re
import sys

original_gene_giff="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/gene.gff"
gene_transcripts={}
with open(original_gene_giff, 'r') as IN:
	for line in IN:

		line=line.rstrip('\n')
		items=re.split('[\t]', line)
		WB=items[1]
		if WB=="WormBase":
			chromosome,WB,Gene_type,WB_start,WB_end,NA1,WB_orient,NA2,transcript=items[0:len(items)]
			match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", transcript) 
			transcript_name=match.group(1)	
			match2=re.search("biotype=(.*);", transcript) 

			if match2:
				biotype=match2.group(1)
				gene_transcripts[transcript_name]=biotype
	
				#print transcript_name
				#print biotype



#print gene_transcripts["Y74C9A.6"]
infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/essentiality_nonredundant_GO_only_new.txt"
with open(infile, 'r') as IN:
	for line in IN:

		line=line.rstrip('\n')
		items=re.split('[\t]', line)
		Chr,start,method,TE,region,transcript_name,gene_name,biotype,phenotype,GO=items[0:len(items)]
		if region == "gene":
			
			if gene_transcripts[transcript_name] == "protein_coding":
				print gene_transcripts[transcript_name]
				print line

