#!/usr/bin/env python 
# this script removes the headers of the TE families that were inserted in the "genome_rearranged" file from RSVSIM (these headers no longer have sequences associated
# with them because they were cut and pasted elsewhere in the genome but the header naems still remain (identified by the presence of "copy" in the header name))
# USE: edit_rearranged_genome.py <genome_rearranged.fasta> <run_ID>
# called on by run_RSVSim.sh
import sys
import re
import os
import random
from Bio import SeqIO


in_file = sys.argv[1]
run_no = sys.argv[2]
fasta_out = "new_genome_{run_no}.fasta".format(**locals())
FASTA_OUT = open(fasta_out, "w")

fasta_sequences = SeqIO.parse(open(in_file),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq) 
	length = len(sequence)
	if not re.search ("copy", name):
		print name
		FASTA_OUT.write(">{name}\n{sequence}\n".format(**locals()))


