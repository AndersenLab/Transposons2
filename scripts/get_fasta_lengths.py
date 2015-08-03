#!/usr/bin/env python 
# this script takes a fasta files as input and outputs a "length.txt" specifying the length of each sequence in the fastsa file
# USE: get_fasta_lengths.py <fasta_file>
import sys
import re
import os
from os.path import basename
from subprocess import Popen, PIPE
from  collections import defaultdict

in_fasta = sys.argv[1]
OUT_LENGTH = open("lengths.txt", "w")

from Bio import SeqIO
#from Bio.Seq import Seq
fasta_sequences = SeqIO.parse(open(in_fasta),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	length_seq = len(sequence)
	OUT_LENGTH.write("{name}\t{length_seq}\n".format(**locals()))