#!/usr/bin/env python 
# this script takes a fasta file as input and extracts all fasta seqs 
# originating from Wormbase to a new file
# USE: pull_WB_fasta.py <fasta_file>
# example: python ../scripts/pull_WB_pos_from_consensus.py /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta

from Bio import SeqIO
import sys
import re

all_seqs = sys.argv[1] 

wb_seqs ="WB_all_seqs_jul30.fasta"
WB_SEQs = open(wb_seqs, "w")
fasta_sequences = SeqIO.parse(open(all_seqs),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, fasta.seq.tostring()
	print name
	if re.search("WBTransposon",name):
		print "YES"
		WB_SEQs.write(">" + name + "\n" +  sequence + "\n")