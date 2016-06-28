#!/usr/bin/env python
# this script creates the supplementary table of the fasta sequences and family renames
# USE: fasta_supp_table.py <TE_consensus> <family_renames>

import re
import sys
from Bio import SeqIO

TE_consensus=sys.argv[1]
family_renames=sys.argv[2]

OUT=open("TE_seqs.txt", 'w')
OUT.write("Element\tFamily\tSequence\n")
renames={}
with open(family_renames, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		element,family=items[0:2]
		renames[element]=family

consensus_families={}
fasta_sequences = SeqIO.parse(open(TE_consensus),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq) 

	#determine family name
	if name in renames.keys():
		family=renames[name]
	else:
		family=name

	family_short=re.sub("_CE$","",family)
	family_short=re.sub("WBTransposon","WBT",family_short)
	name_short=re.sub("_CE$","",name)
	name_short=re.sub("WBTransposon","WBT",name_short)

	OUT.write("{name_short}\t{family_short}\t{sequence}\n".format(**locals()))

OUT.close()