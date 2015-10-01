#!/usr/bin/env python
# this script checks that all into/exon/UTR transcript names are also found in the gene.gff file and if not exits with an error
# USE: gene_check.py
import re


gene_file="all_gene_insertions_and_excisions.txt"
file_list=("all_introns_insertions_and_excisions.txt",\
"all_exons_insertions_and_excisions.txt",\
"all_fiveUTR_insertions_and_excisions.txt",\
"all_threeUTR_insertions_and_excisions.txt")

other_names={}
gene_names={}

for i in file_list:
	with open(i, 'r') as IN:
		for line in IN:
			items=re.split("[\t]|;", line) #split by tab or semicolon
			name=items[7]
			#match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w.]+)", name)
			match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", name)
			ID = match.group(1)
			other_names[ID]=0

with open(gene_file,'r') as IN:
	for line in IN:
		items=re.split("[\t]", line)
		name=items[7]
		match = re.search("sequence_name=([\w|\d]+.\d+)", name)
		ID = match.group(1)
		gene_names[ID]=0

for i in other_names.keys():
	if i not in gene_names.keys():
		print "NOT FOUND"
		print "ERROR: Transcript '{i}' does not exist in the gene file, exiting...".format(**locals())
		sys.exit()



