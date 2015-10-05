#!/usr/bin/env python
# this script classifies TE position as genic or intergenic
# it also outputs the sequeunce name if a TE was found in a gene and whether or not that TE was in the the "border" region of a gene(within 10bp from the end)
# or if it is in an "internal" region
# USE: separate_gene_assignments.py

import re

internal_TEs={}
full_ins_TEs={}
all_TEs={}


'''
chr X
start_TE X
end_TE x
start_gene/NA X
TE X
orient X
RS X
part X
gene_name/NA X
gene class protein_coding/pseduogene/NA X
''' 



# put TEs overlapped with genes into dict
with open("all_window.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		TE=items[12]
		
		#gene info:
		gene_start=items[3]
		gene_end=items[4]
		gene_info=items[8]
		match=re.search("sequence_name=([A-za-z\d\.]+);", gene_info)
		gene_name=match.group(1)
		match=re.search("biotype=([A-za-z]+);", gene_info)
		gene_class=match.group(1)
		#TE info:
		match = re.search("([A-Z]+)_(\d+)_(\d+)_(([A-Za-z\d+_-]+))_((\w+-)?reference)_+([a-z]+)_([a-z]+)_(\d+)_([\d\.]+)_([\+-])_([A-Za-z]+)_(\w+)", TE) #([A-Za-z\d+_])_((\w+-)?reference)\w+_\d+_\d+
		chromosome = match.group(1)
		start = match.group(2)
		start2 = match.group(3)
		te = match.group(4)
		RS = match.group(11)
		orient = match.group(12)
		method=match.group(13)
		sample = match.group(14)

		new_info="{chromosome}\t{start}\t{start2}\t{method}\t{gene_start}\t{gene_end}\t{te}\t{orient}\t{RS}\t{gene_name}\t{gene_class}".format(**locals())
		internal_TEs[TE]=new_info


		
# put insertion TEs overlapped with  full genes giff into dict
with open("insertions_full_window.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		TE=items[12]
		
		#gene info:
		gene_start=items[3]
		gene_end=items[4]
		gene_info=items[8]
		match=re.search("sequence_name=([A-za-z\d\.]+);", gene_info)
		gene_name=match.group(1)
		match=re.search("biotype=([A-za-z]+);", gene_info)
		gene_class=match.group(1)
		#TE info:
		match = re.search("([A-Z]+)_(\d+)_(\d+)_(([A-Za-z\d+_-]+))_((\w+-)?reference)_+([a-z]+)_([a-z]+)_(\d+)_([\d\.]+)_([\+-])_([A-Za-z]+)_(\w+)", TE) #([A-Za-z\d+_])_((\w+-)?reference)\w+_\d+_\d+
		chromosome = match.group(1)
		start = match.group(2)
		start2 = match.group(3)
		te = match.group(4)
		RS = match.group(11)
		orient = match.group(12)
		method=match.group(13)
		sample = match.group(14)

		new_info="{chromosome}\t{start}\t{start2}\t{method}\t{gene_start}\t{gene_end}\t{te}\t{orient}\t{RS}\t{gene_name}\t{gene_class}".format(**locals())
		full_ins_TEs[TE]=new_info


# put all TEs into dict
with open("all.bed", 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		TE=items[3]

		#TE info:
		match = re.search("([A-Z]+)_(\d+)_(\d+)_(([A-Za-z\d+_-]+))_((\w+-)?reference)_+([a-z]+)_([a-z]+)_(\d+)_([\d\.]+)_([\+-])_([A-Za-z]+)_(\w+)", TE) #([A-Za-z\d+_])_((\w+-)?reference)\w+_\d+_\d+
		chromosome = match.group(1)
		start = match.group(2)
		start2 = match.group(3)
		te = match.group(4)
		RS = match.group(11)
		orient = match.group(12)
		method=match.group(13)
		sample = match.group(14)
		gene_start="NA"
		gene_end="NA"
		gene_name="NA"
		gene_class="NA"

		new_info="{chromosome}\t{start}\t{start2}\t{method}\t{gene_start}\t{gene_end}\t{TE}\t{orient}\t{RS}\t{gene_name}\t{gene_class}".format(**locals())
		all_TEs[TE]=new_info


OUT=open("TE_gene_interrupt_output.txt", 'w')
for key, value in all_TEs.items():
	if key in full_ins_TEs.keys() and key not in internal_TEs.keys():
		part="border"
		overall="Genic"
		value=full_ins_TEs[key]
	elif key in internal_TEs.keys():
		part="internal"
		overall="Genic"
		value=internal_TEs[key]
	else:
		part="intergenic"
		overall="Intergenic"
	value="{value}\t{part}\t{overall}".format(**locals())
	OUT.write(value)
	OUT.write('\n')
OUT.close()




