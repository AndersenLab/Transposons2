#!/usr/bin/env python 
# this script generates the fastsa consens and corresponding length file for round3 of the TE caller simulations (just TEMP)
# USE: generate_round3_consensus_fasta_file.py <family
# example:  python ../../scripts/generate_round3_consensus_fasta_file.py ../new_WB_famlies.txt  round2_consensus_fasta.fasta
import sys
import re
import os
from os.path import basename
from subprocess import Popen, PIPE
from  collections import defaultdict

fam_file = sys.argv[1]
in_fasta = sys.argv[2]
NEW_FAM_FILE = open(fam_file, "r")
OUT_LENGTH = open("round3_lengths.txt", "w")
OUT_FASTA = open("round3_consensus_fasta.fasta", "w")

fam_to_keep={}

for line in NEW_FAM_FILE:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	WB_ID = items[0]
	WB_family = items[1]
	fam_to_keep[WB_family] =0
	print WB_family
NEW_FAM_FILE.close()

from Bio import SeqIO
#from Bio.Seq import Seq
fasta_sequences = SeqIO.parse(open(in_fasta),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	length_seq = len(sequence)
	if name in fam_to_keep.keys():
		OUT_FASTA.write(">{name}\n{sequence}\n".format(**locals()))
		OUT_LENGTH.write("{name}\t{length_seq}\n".format(**locals()))
OUT_FASTA.close()
OUT_LENGTH.close()
