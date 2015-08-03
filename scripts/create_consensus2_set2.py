#!/usr/bin/env python 
#this script combines the REPbase consensus seqs with the individual WB element seqs while omitting the TC8 family
#USE: create_consensus2_set2.py <file specifying te renames> <WB seqs (plus strand)>
##python /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/create_consensus2_set2.py /lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_familes_set2.txt /lscr2/andersenlab/kml436/git_repos2/Transposons/files/WB_all_tes_plus_strand.fasta
from Bio import SeqIO
from subprocess import Popen, PIPE
import sys
import re
import os
wb_fam_file = sys.argv[1] 
WB_FAM_FILE= open(wb_fam_file, "r")

wb_seqs = sys.argv[2] 
WB_SEQS= open(wb_seqs, "r")

out_file = "Wb-TC8.fasta"
OUT_FILE = open(out_file, "w")

TC8_tes = {}
for line in WB_FAM_FILE:
	line = line.rstrip('\n')
	items = re.split("[\t]",line)
	WB_te = items[0]
	WB_family = items[1]
	if WB_family == "TC8":
		TC8_tes[WB_te] = 0


fasta_sequences = SeqIO.parse(open(wb_seqs),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	if name not in TC8_tes.keys(): #do not include TC8 transposons
		print name
		OUT_FILE.write(">{name}\n{sequence}\n".format(**locals()))
OUT_FILE.close()


cmd = """cat {out_file} /lscr2/andersenlab/kml436/git_repos2/Transposons/files/REPB_all_seqs.fasta > round2_consensus_set2.fasta""".format(**locals())
print cmd
result, err = Popen([cmd], stdout=PIPE, stderr=PIPE, shell=True).communicate()
WB_FAM_FILE.close()