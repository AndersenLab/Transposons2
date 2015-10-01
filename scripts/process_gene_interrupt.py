#!/usr/bin/env python
# this script extracts gene/transcript ID and TE family info from an insertions/excisions results file and outputs a rearranged file
# USE: process_gene_interrupt.py <all_insertions_and_excisions_file>
import sys
import re

in_file = sys.argv[1]
IN_FILE = open(in_file, "r")

out_file = "final_{in_file}".format(**locals())
OUT_FILE = open(out_file, "w")
for line in IN_FILE:
	line = line.rstrip('\n')
	items = re.split("[\t]", line)
	chromosome = items[2]
	category = items[3]
	start_pos = items[5]
	te_start = items[8]
	te_info = items[9]
	te_orient = items[11]
	match = re.search("(.*)_(\w+-)?reference", te_info)
	family = match.group(1)

	ID_info = items[7]
	match2 = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", ID_info) #jsut pull gene name, remove splice info
	ID = match2.group(1)

	other_fields = '\t'.join(items[0:7])

	OUT_FILE.write("{other_fields}\t{ID}\t{te_start}\t{family}\t{te_orient}\n".format(**locals()))


IN_FILE.close()