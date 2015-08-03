#!/usr/bin/env python 
# this script takes a fasta file as input and outputs the reverse complement of the sequence
# USE: reverse_comple.py <fasta_file>
import sys
import re
import os
from subprocess import Popen, PIPE

input_file = sys.argv[1]
OUT = open("WB_all_tes_plus_strand.fasta", "w")

from Bio import SeqIO
from Bio.Seq import Seq
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	print sequence
	seq = Seq(sequence)
	print "\n"
	rev = seq.reverse_complement()
	OUT.write(">" + name + "\n" +  str(rev) + "\n")
OUT.close() 

