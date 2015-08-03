#this script ........
import sys
import re
import os
from os.path import basename
from subprocess import Popen, PIPE
from Bio import SeqIO

my_directory = sys.argv[1]
in_fasta = sys.argv[2]
NO_MATCH = open("round4_WB_fastas.fasta", "w")
print my_directory
fam_to_keep={}
#my_directory = sys.argv[4]
for results_file in os.listdir(my_directory):
	if re.search('(.*)_fastas.fasta$', results_file):
		match = re.search('(.*)_fastas.fasta$', results_file)
		file_name = match.group(1)
		fam_to_keep[file_name] = 0
		if not re.search('(^CONSENSUS)|(^check)', results_file):
			print results_file
			fasta_sequences = SeqIO.parse(open("{my_directory}/{results_file}".format(**locals())),'fasta')
			for fasta in fasta_sequences:
				name, sequence = fasta.id, str(fasta.seq)
				NO_MATCH.write(">{name}\n{sequence}\n".format(**locals()))


OUT_LENGTH = open("round4_lengths.txt", "w")
OUT_FASTA = open("round4_consensus_fasta.fasta", "w")
fasta_sequences = SeqIO.parse(open(in_fasta),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	length_seq = len(sequence)
	if name in fam_to_keep.keys():
		OUT_FASTA.write(">{name}\n{sequence}\n".format(**locals()))
		OUT_LENGTH.write("{name}\t{length_seq}\n".format(**locals()))
OUT_FASTA.close()
OUT_LENGTH.close()



