#!/usr/bin/env python
# this script clipa the WB gene positions by 10 base pairs on each side and output a "clipped" gff
# USE: modify_gff.py
import re

OUT=open("clipped_WB_gene_positions.gff",'w')
with open("WB_gene_positions.gff", 'r') as IN:
	for line in IN:
		items=re.split("[\t]", line)
		items[3]=str(int(items[3])+10)
		items[4]=str(int(items[4])-10)

		new_line ='\t'.join(items[0:9])
		OUT.write(new_line)

OUT.close()