#!/usr/bin/env python 
# old check script...ignore
import sys
import re
import os

in_file = sys.argv[1]
out_file = "RSV_te_postions.txt"
OUT_FILE = open(out_file, "w")
from Bio import SeqIO
fasta_sequences = SeqIO.parse(open(in_file),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	length_seq = len(sequence)
	#don't print out the actual chromosomees
	if not re.match(r'(^I$)|(^II$)|(^III$)|(^IV$)|(^V$)|(^X$)',name):
		OUT_FILE.write("{name}\t1\t{length_seq}\t{name}\n".format(**locals()))
OUT_FILE.close()