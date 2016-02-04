#!/usr/bin/env python

import re
from Bio import SeqIO

all_WB_seqs="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/Wb-TC8.fasta"
HL_gff="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names_alias.bed"
family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"
TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"

elements_to_keep={}
families_to_keep={}


HLOUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round4/WB_pos_element_names_alias_round4.bed", 'w')
TEOUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round4/round2_consensus_set2_round4.fasta",'w')

# Regenerate the HL GFF file
with open(all_WB_seqs, 'r') as IN:
	for line in IN:	
		line=line.rstrip("\n")
		if re.search("^>", line):
			element=line.replace(">","")
			elements_to_keep[element]=0
			families_to_keep[element]=0

with open(HL_gff, 'r') as IN:
	for line in IN:
		items=re.split("[\t]",line)
		WB_ID=items[2]
		if WB_ID in elements_to_keep.keys():
			HLOUT.write(line)

HLOUT.close()

# Regenerate the consensus file
with open(family_renames, 'r') as IN:
	for line in IN:

		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		element,family=items[0:2]
		families_to_keep[element]=0
		families_to_keep[family]=0


fasta_sequences = SeqIO.parse(open(TE_consensus),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq) 
	
	if name in families_to_keep.keys():
		print name
		TEOUT.write(">{name}\n{sequence}\n".format(**locals()))

