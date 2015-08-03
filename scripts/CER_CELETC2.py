#!/usr/bin/env python 
# this script pulls the fasta seqs of WB elements in the CER6 or CELETC2 families
# USE:CER_CELETC2.py in appropriate folder
import sys
import re
import os
from subprocess import Popen, PIPE

OUT_CER6 = open("WB_CER6_fastas.fasta", "w")
OUT_CELETC2 = open("WB_CELETC2_fastas.fasta", "w")
cer6 = {}
celetc2 = {}

input_file = "WB_all_tes_plus_strand.fasta"
FAM_FILE = open("WB_familes.txt", "r")

for line in FAM_FILE:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	WB_ID = items[0]
	WB_family = items[1]
	if WB_family == "CER6":
		cer6[WB_ID] = 0
	if WB_family == "CELETC2":
		celetc2[WB_ID] = 0
FAM_FILE.close()

from Bio import SeqIO
from Bio.Seq import Seq
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	if name in cer6.keys():
		OUT_CER6.write(">" + name + "\n" +  sequence + "\n")
	if name in celetc2.keys():
		OUT_CELETC2.write(">" + name + "\n" +  sequence + "\n")

OUT_CER6.close() 
OUT_CELETC2.close() 