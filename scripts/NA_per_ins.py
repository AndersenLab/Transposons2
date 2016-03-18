#!/usr/bin/env python
import re
infile="kin_matrix_ins.txt"
OUT=open("NA_counts_at_positions.txt", 'w')
with open(infile, 'r') as IN:
	next(IN)
	for line in IN:
		NA_count=0
		line=line.rstrip('\n')
		items=re.split('\t',line)
		trait=items[0]
		len_it=len(items)+1
		strain_no=len_it-2

		for i in items[1:len_it]:
			if i=="NA":
				NA_count +=1
		NA_fraction=float(NA_count)/strain_no
		NA_fraction= round(NA_fraction,2)
		OUT.write("{trait}\t{NA_count}\t{NA_fraction}\n".format(**locals()))

OUT.close()