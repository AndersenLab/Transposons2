#!/usr/bin/env python
# this script generates the CtCp_all_nonredundant_file
# USE: generate_CtCp.py <kin_matrix_full.txt>

import sys
import re


# get orientations of reference transposons
refs="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names.gff"
ref_tes={}
with open(refs, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('[\t]', line)
		chromosome,start,stop,ref,R,orient=items[0:7]
		ID=chromosome+"_"+start+"_"+ref
		ref_tes[ID]=orient

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

				#I_220561_TURMOIL1_non-reference__temp_rp_7_25_-_new_ED3011
				#for insertion calls
				if re.search("non-reference", TE_info):
					if i=="1": #count only positive confirmation of  TE insertion
						method="new"
						match=re.search('non-reference__temp_[A-Za-z]+_[\d\.]+_[\d\.]+_(\+|-)_', TE_info)
						orient=match.group(1)
						OUT.write("{chromosome}\t{start}\t{start}\t{TE_info}\t{orient}\t{method}\t{strain}\n".format(**locals()))  
				elif i=="1":
					if TE_info in ref_tes.keys():
						orient=ref_tes[TE_info]
					TE_info=TE_info+"_reference" # a reference TE position
					method="reference"
					OUT.write("{chromosome}\t{start}\t{start}\t{TE_info}\t{orient}\t{method}\t{strain}\n".format(**locals()))  
				elif i=="0":
					if TE_info in ref_tes.keys():
						orient=ref_tes[TE_info]
					TE_info=TE_info+"_reference" # a reference TE position
					method="absent"
					OUT.write("{chromosome}\t{start}\t{start}\t{TE_info}\t{orient}\t{method}\t{strain}\n".format(**locals()))  
				else:
					print("ERROR, cannot determine call, exiting...")
					sys.exit()


OUT.close()

