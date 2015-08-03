#!/usr/bin/env python 
# this script checks to make sure each strain has values for the 3 transposon caller methods and prints out the strain name and method if not found
# USE: find_error.py <FINAL_RESULTS.txt>


import sys
import re
absent={}
new={}
reference={}
strain_dict={}
te_input = sys.argv[1]
TE_INPUT = open(te_input, "r")
for line in TE_INPUT:
	line = line.rstrip('\n')
	items = re.split("[\t]",line)
	strain_name = items[1]
	method = items[2]
	strain_dict[strain_name]=1
	if method=="absent":
		absent[strain_name]=1
	if method=="new":
		new[strain_name]=1
	if method=="reference":
		reference[strain_name]=1


for strain in strain_dict.keys():
	if strain not in absent.keys():
		print "absent"
		print strain
	if strain not in new.keys():
		print "new"
		print strain
	if strain not in reference.keys():
		print "reference"
		print strain
TE_INPUT.close()