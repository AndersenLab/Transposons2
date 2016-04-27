#!/usr/bin/env python
# this script reformats grantham score information
# USE: one_to_three.py

import sys
import re
import os

three_letter ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', 'Y':'TYR',    \
'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA',    \
'G':'GLY', 'P':'PRO', 'C':'CYS'}

infile="grantham_info.txt"
OUT=open("grantham_scores.txt",'w')

with open(infile) as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=line.split()
		before=items[0]
		after=items[1]
		items[0]=three_letter[before].capitalize()
		items[1]=three_letter[after].capitalize()

		OUT.write(items[0] + "2" + items[1] + '\t' + items[2] + '\n')
OUT.close()