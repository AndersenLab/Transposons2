#!/usr/bin/env python 
# this script takes the positions of the reference TEs and simulated TEs and adjusts the bed file
# to refelect the new posititons in the rearranged genome
# USE:adjust_RSVSIM_positions.py <bed_files of te positions>
# called on by run_RSVSim.sh

import sys
import re
import os
import random

chr_add_length = {}
chr_add_length["I"] = 0
chr_add_length["II"] = 0
chr_add_length["III"] = 0
chr_add_length["IV"] = 0
chr_add_length["V"] = 0
chr_add_length["X"] = 0

in_file = sys.argv[1]
IN_FILE = open(in_file, "r")
out_file = "adjusted_pos_{in_file}".format(**locals())
OUT_FILE = open(out_file, "w")
for line in IN_FILE:
	line = line.strip('\n')
	items = re.split("\t",line)
	chromosome = items[0]
	start_pos = int(items[1])
	end_pos = int(items[2])
	te_family = items[3]
	ref = items[4]
	orientation = items[5]
	length = end_pos - start_pos + 1  ##need to add the one
	print start_pos
	print end_pos
	print length
	if ref == "R":
		new_start = start_pos + chr_add_length[chromosome]
		new_end = end_pos + chr_add_length[chromosome] 
		print chr_add_length[chromosome]
	else:
		new_start = start_pos + chr_add_length[chromosome]
		new_end = end_pos + chr_add_length[chromosome]
		chr_add_length[chromosome] += length #add length if a simulated (not original reference) te
		print chr_add_length[chromosome]
	OUT_FILE.write("{chromosome}\t{new_start}\t{new_end}\t{te_family}\t{ref}\t{orientation}\n".format(**locals()))

	####TAKE INTO ACCOUNT ORIENTATION????


IN_FILE.close