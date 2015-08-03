#!/usr/bin/env python 
# this script processes transposon caller output to calculate how many transposons of each family  and totals were found in a strain
# USE: te_totals.py <bam_name_temp_nonredundant.bed> <consensus_renamed> <bam_name> <name_results_file> <method> <no_of_end_tes>
# called on by ProcessTransposonCallerOutput.py

import sys
import re
from Bio import SeqIO
from datetime import date
import string
import os


te_output = sys.argv[1] ## the new_CT..<file> for temp insertion
consensus = sys.argv[2] ##want parred down consensus file here!!!!!!
bam_name = sys.argv[3]
results_file = sys.argv[4]
method = sys.argv[5]
end_tes = sys.argv[6]
##COUNT TOTAL
file_name=os.path.splitext(results_file)[0]
print file_name
RESULTS_FILE=open(results_file, "aw")
RESULTS_FILE_LF=open("{file_name}_LF.txt".format(**locals()), "aw")
date = date.today()


#add names of all TE families to a dictionary
consensus_families={}
fasta_sequences = SeqIO.parse(open(consensus),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq) 
	consensus_families[name]=0

total_tes = 0
TE_OUTPUT = open(te_output, "r")
for line in TE_OUTPUT:
	line = line.rstrip('\n')
	items = re.split("[\t]", line)
	transposon = items[3]
	match = re.search("(.*)_(\w+-)?reference", transposon)
	family = match.group(1)
	if family in consensus_families.keys():
		consensus_families[family] +=1
		total_tes +=1
	else:
		print "ERROR: The detected family '{family}' does not exist in the consensus file, exiting...".format(**locals()) #produce error if family does not exist in the consensus file
		RESULTS_FILE.write("{date}\t{bam_name}\t{method}\tERROR\n".format(**locals()))
		sys.exit()

TE_OUTPUT.close()
RESULTS_FILE.write("{date}\t{bam_name}\t{method}\t{end_tes}\t{total_tes}\t".format(**locals()))
RESULTS_FILE_LF.write("{date}\t{bam_name}\t{method}_TRANS_total\t{total_tes}\n".format(**locals()))
RESULTS_FILE_LF.write("{date}\t{bam_name}\t{method}_TRANS_end_tes\t{end_tes}\n".format(**locals()))


for key in sorted(consensus_families.keys()):
	count = consensus_families[key]
	RESULTS_FILE.write("{count}\t".format(**locals()))
	RESULTS_FILE_LF.write("{date}\t{bam_name}\t{method}_TRANS_{key}\t{count}\n".format(**locals()))
RESULTS_FILE.write("\n")



