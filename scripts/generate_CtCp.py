#!/usr/bin/env python
# this script generates the CtCp_all_nonredundant_file
# USE: generate_CtCp.py <kin_matrix_full.txt>

import sys
import re

OUT=open("all_nonredundant.txt",'w')
infile=sys.argv[1]
with open(infile, 'r') as IN:
	header=next(IN)
	header=header.rstrip('\n')
	strains=re.split('\t', header)
	for line in IN:
		line =line.rstrip('\n')
		items=re.split('[\t]', line)
		

		for index,i in enumerate(items[1:len(items)],start=1):
			ID=index
			strain=strains[ID]
			if i!="NA":
				TE_info=items[0]
				items2=re.split("_", TE_info)

				chromosome,start=items2[0:2]
				#for insertion calls
				if re.search("non-reference", TE_info):
					if i=="1": #count only positive confirmation of  TE insertion
						method="new"
					else:
						break
				elif i=="1":
					TE_info=TE_info+"_reference" # a reference TE position
					method="reference"
				elif i=="0":
					TE_info=TE_info+"_reference" # a reference TE position
					method="absent"
				else:
					print("ERROR, cannot determine call, exiting...")
					sys.exit()
#

					#output "start" twice for placeholder
				OUT.write("{chromosome}\t{start}\t{start}\t{TE_info}\tNA\t{method}\t{strain}\n".format(**locals()))  ####unidents this by one!!!!
OUT.close()

