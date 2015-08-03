#!/usr/bin/env python 
# this script takes a fasta file and outputs non-WB transposons
# USE: pull_REPBASE_fasta.py <fasta_file>
# example: python ../scripts/pull_REPBASE_fasta.py all_transposons.fasta

from Bio import SeqIO
import sys
import re

all_seqs = sys.argv[1] 

repb_seqs ="REPB_all_seqs.fasta"
REPB_SEQs = open(repb_seqs, "w")
fasta_sequences = SeqIO.parse(open(all_seqs),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, fasta.seq.tostring()
	print name
	if re.search("WBTransposon",name):
		print "WB Transposon"
	else:
		print "REPB Transposon"
		REPB_SEQs.write(">" + name + "\n" +  sequence + "\n")