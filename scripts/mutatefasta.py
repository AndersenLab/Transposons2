#!/usr/bin/env python 
# this script introduces a specified percent of mutations(SNPs) into a fasta sequence
# USE: mutatefasta.py <fasta_file> <percent_of_seq_to_mutate> <element->family_conversion_file>
# example: python /lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/mutatefasta.py /lscr2/andersenlab/kml436/git_repos2/Transposons/files/WB_all_seqs.fasta .05 /lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_familes_set2.txt

#NOTE: make sure to prepare emboss before running this script
#NOTE: allows for the possibility of performing 2 mutations on the same nucleotide position

from Bio import SeqIO
from subprocess import Popen, PIPE
import sys
import re
import os
seq = sys.argv[1] 
per = sys.argv[2]
wb_fam_file = sys.argv[3] 
WB_FAM_FILE= open(wb_fam_file, "r")



TC8_tes = {}
for line in WB_FAM_FILE:
	line = line.rstrip('\n')
	items = re.split("[\t]",line)
	WB_te = items[0]
	WB_family = items[1]
	if WB_family == "TC8":
		TC8_tes[WB_te] = 0






fasta_sequences = SeqIO.parse(open(seq),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	length = len(sequence)
	print name
	print length
	print per
	#mutate percent * length of seq and round to nearest integer
	num_mutations = int(round(float(per)*int(length)))
	print num_mutations
	if name not in TC8_tes.keys():
		out = "ind_{name}.fasta".format(**locals())
		OUT = open(out, "w")
		OUT.write(">{name}\n{sequence}\n".format(**locals()))
		OUT.close()

		#run emboss to introduce SNPs into the fasta seqs
		cmd = """msbar -sequence {out} -snucleotide -count {num_mutations} -point 4 -block 0 -codon 0 -outseq mut_{name}.fasta""".format(**locals())
		print cmd
		result, err = Popen([cmd], stdout=PIPE, stderr=PIPE, shell=True).communicate()

#concatenate the mutated fastas
cmd = """cat mut_* > per{per}_mut_WB_seqs""".format(**locals())
print cmd
result, err = Popen([cmd], stdout=PIPE, stderr=PIPE, shell=True).communicate()