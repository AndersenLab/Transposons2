#!/usr/bin/env python
# this script scores the reference/absence calls and checks that NA counts match
# USE: kin_hash_AF.py <sample_list>



import re
from collections import defaultdict
import sys
from subprocess import Popen, PIPE
import os
import sys

#file_dir="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files"
WB_fam_pos="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_family_names.gff"
sample_list=sys.argv[1]
data_dir="/lscr2/andersenlab/kml436/git_repos2/Transposons2/data"
NA_file="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/non_Nas.txt"

SAMPLES={}
SAMPLES=defaultdict(list)

total_TEs=0 # the number fo reference TEs
#put all reference TE positions into a dictionary as the keys:
with open(WB_fam_pos, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		chromosome, start, end, TE, R, orient =items[0:(len(items))]
		info=(chromosome, start, TE)
		if TE != "TC8": #ignore TC8s
			ID="_".join(info)
			SAMPLES[ID] #add key to dict
			total_TEs +=1



#put counts of calls for each strain into a dictionary
non_NA={}
with open(NA_file, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\s]", line)
		strain=items[0]
		calls=items[1]
		non_NA[strain]=int(calls)

with open(sample_list, 'r') as SAMPLE_LIST:
	for sample in SAMPLE_LIST:
		pos_list=list()
		NA_count=0
		sample=sample.rstrip('\n')
		#if sample has an absence call at a reference position, score it as a zero
		with open("{data_dir}/{sample}/final_results/{sample}_temp_absence_nonredundant.bed".format(**locals())) as IN:
			for line in IN:
				line=line.rstrip('\n')
				items=re.split("[\t]", line)
				chromosome, start, end, TE, R, orient=items[0:(len(items))]
				match=re.search("(.*)_(\w+-)?reference", TE)
				family=match.group(1)
				info=(chromosome, start, family)
				ID="_".join(info)
				SAMPLES[ID].append("{sample}:0".format(**locals()))
				pos_list.append(ID)

		#if sample has an reference call at a reference position, score it as a 1	
		with open("{data_dir}/{sample}/final_results/{sample}_telocate_nonredundant.bed".format(**locals())) as IN:
			for line in IN:
				line=line.rstrip('\n')
				items=re.split("[\t]", line)
				chromosome, start, end, TE, R, orient=items[0:(len(items))]
				match=re.search("(.*)_(\w+-)?reference", TE)
				family=match.group(1)
				info=(chromosome, start, family)
				ID="_".join(info)
				SAMPLES[ID].append("{sample}:1".format(**locals()))
				pos_list.append(ID)

		#if no absence or reference call, score the sample with "NA" and check that "NA" counts here match NA counts calculated with a previous script 	
		for key in SAMPLES.keys():
			if key not in pos_list:
				SAMPLES[key].append("{sample}:NA".format(**locals()))
				NA_count +=1
				#calculate non_NA count:
		non_nas=int(total_TEs)-NA_count
		if int(non_nas)!=int(non_NA[sample]):
		#if int(non_nas)!=int(non_NA[sample] +1):
			print "Inconsistency in NA counts for sample {sample}: {non_nas}, exiting...".format(**locals())
			sys.exit()



#print kinship matrix
KIN_MATRIX=open("kin_matrix_AF.txt", 'w')
KIN_MATRIX.write("TE")

#print headers
value=SAMPLES.values()[1] #pull first value from dictionary
value=sorted(value)
for i in value:
	strain=(re.split(":", i))[0]
	KIN_MATRIX.write("\t{strain}".format(**locals()))
KIN_MATRIX.write('\n')

#output the final matrix
for key,value in SAMPLES.items():
	value=sorted(value)
	KIN_MATRIX.write(str(key)) #or output key ID here?
	for i in value:
		score=(re.split(":", i))[1]
		KIN_MATRIX.write("\t{score}".format(**locals()))
	KIN_MATRIX.write('\n')
KIN_MATRIX.close()

result, err = Popen(["""head -n1 kin_matrix_AF.txt > tmp"""],stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""cat kin_matrix_AF.txt |sed 1d| sort -t"_" -k1,1 -k2,2n >>tmp """],stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""mv tmp kin_matrix_AF.txt """],stdout=PIPE, stderr=PIPE, shell=True).communicate()






