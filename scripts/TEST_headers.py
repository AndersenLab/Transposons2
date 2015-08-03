#!/usr/bin/env python 
# this script takes the fasta file and outputs the transposon names and as the header of a results file
# USE: TEST_headers.py
from Bio import SeqIO
results_file="/lscr2/andersenlab/kml436/git_repos2/Transposons/FINAL_RESULTS.txt"
consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/SET2/AB-PR/consensus_wTC8.fasta" 

RESULTS_FILE=open(results_file, "aw")
consensus_families={}
fasta_sequences = SeqIO.parse(open(consensus),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq) 
	consensus_families[name]=0


for key in sorted(consensus_families.keys()):
	count = consensus_families[key]
	RESULTS_FILE.write("{key}\t".format(**locals()))
RESULTS_FILE.write("\n")
