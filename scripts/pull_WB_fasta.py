#!/usr/bin/env python 
# this script takes a fasta file as input and extracts all fasta seqs 
# originating from Wormbase to a new file
# USE: pull_WB_fasta.py <fasta_file>
# example: python ../scripts/pull_WB_fasta.py all_transposons.fasta

from Bio import SeqIO
import sys
import re

all_seqs = sys.argv[1] 

wb_seqs ="WB_all_seqs.fasta"
WB_SEQs = open(wb_seqs, "w")
fasta_sequences = SeqIO.parse(open(all_seqs),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, fasta.seq.tostring()
	print name
	if re.search("WBTransposon",name):
		print "YES"
		WB_SEQs.write(">" + name + "\n" +  sequence + "\n")
