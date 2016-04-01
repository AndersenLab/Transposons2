#!/usr/bin/env python
# this scripts compare the absence/reference calls of TEMP/TELOCATE for the CER1 transposon and compares them to a "truth set"
# from "Chromosome-scale selective sweeps shape Caenorhabditis elegans genomic diversity" and output a file of the strain, TE call, truth call, 
# and whether the TE call was correct,incorrect, or "NA" if a comparison could not be made due to missing data
# USE: find_cers.py in directory containing the "Full_Results.txt" and "cer_PA.txt" files

import re

calls="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/kin_C_matrix_full.txt"
CALLS=open(calls,"r")

absences={}
references={}
call={}

first_line=True
for line in CALLS:
	line=line.rstrip('\n')
	items=re.split("[\t]",line)
	strain=items[0]
	if first_line:
		ab=items.index('absent_TRANS_CER1_C') #get index of CER1 absence column
		ref=items.index('reference_TRANS_CER1_C') #get index of CER1 reference column
		first_line=False
	else:
		absences[strain]=items[ab] #key:strain, value: absent_TRANS_CER1 column
		references[strain]=items[ref] #key:strain, value: reference_TRANS_CER1 column

for i in absences.keys():
	if absences[i] == references[i]:
		print "ERROR: TE cannot be both present and absent" # not necessarily error, could both be NAs, adjust if needed
	if absences[i]=="1":
		call[i]="absent"
	elif references[i]=="1":
		call[i]="present"
	else:
		call[i]="NA"
CALLS.close()

#########################################################################################################
ng="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/cer_PA.txt" 
NG=open(ng, "r")

truth_call={}
for line in NG:
	line=line.rstrip('\n')
	items=re.split("[\t]",line)
	name=items[0]
	pa=items[1]
	if name in call.keys():
		if pa=="1":
			truth_call[name]="absent"
		elif pa=="0":
			truth_call[name]="present"
		else:
			truth_call[name]="NA"
#########################################################################################################
outfile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/cer_comparison.txt"
OUTFILE=open(outfile, "w")
for i in sorted(call.keys()):
	found=call[i]
	if i in truth_call.keys():
		truth=truth_call[i]
		if call[i] =="NA" or truth_call[i]=="NA":
			comparison="NA" #if either the TE callers or truth caller output an "NA" result, no comparison can be made
		elif call[i]==truth_call[i]:
			comparison="CORRECT" #if the TEMP/TELOCATE transposon caller result equals the truth call, it was correct
		else:
			comparison="INCORRECT" #if the TEMP/TELOCATE transposon caller result equals the truth call, it was incorrect
	else:
		truth="NA"
		comparison="NA"
	OUTFILE.write("{i}\t{found}\t{truth}\t{comparison}\n".format(**locals()))
