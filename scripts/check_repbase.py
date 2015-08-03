#!/usr/bin/env python 
# this script checks which repbase names in repbase_header.txt are not in REPB_all_seqs.fasta (renames)
# USE: check_repbase.py
import re

REPB_all="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/REPB_all_seqs.fasta"
repbase="/lscr2/andersenlab/kml436/repbase_header.txt"


from Bio import SeqIO
from Bio.Seq import Seq
IDs={}
fasta_sequences = SeqIO.parse(open(REPB_all),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	IDs[name]=0
	#don't want it formatted as a string 
	#name, sequence = fasta.id, fasta.seq.tostring()


REPBASE=open(repbase, "r")
for line in REPBASE:
	line=line.rstrip('\n')
	items=re.split("[\t]", line)
	family=items[0]
	family=family.replace(">","")
	if family not in IDs.keys():
		print family