#1/usr/bin/env python
import re


infile="NA_counts_at_positions.txt"

keeps={}
with open(infile, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		position,count,frac=items
		if float(frac) <.10:
			keeps[position]=frac

OUT=open("kin_matrix_ins_reduced.txt", 'w')
org="kin_matrix_ins.txt"

with open(org, 'r') as IN:
	header=next(IN)
	OUT.write(header)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		pos=items[0]
		if pos in keeps.keys():
			OUT.write(line +'\n')

OUT.close()
