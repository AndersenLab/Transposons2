#!/usr/bin/env python
# this script splits the files contain all positions of called TEs among the strains into separate files depending on the TE family and method caller
# USE: kin.py <all_nonredundant.txt>
# ex: kin.py /lscr2/andersenlab/kml436/git_repos2/Transposons/data/all_nonredundant.txt

import sys
import re
import os
from  collections import defaultdict
from subprocess import Popen, PIPE


kin_step2="/lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/kin_step2.sh"
##NEED TO EDIT SAMPLE FILE IN BELOWSCRIPT TOO
kin_step3="/lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/kin_temp.py"
transpose_matrix="/lscr2/andersenlab/kml436/git_repos2/Transposons/scripts/transpose_matrix.sh"
#
#
#
#
#change smaple list to full one later
sample_list="/lscr2/andersenlab/kml436/git_repos2/Transposons/data/full_sample_list.txt"
####
#
#
#
#
os.system("mkdir TE_matrix")
dir=os.getcwd() # get current directory
os.chdir("{dir}/TE_matrix".format(**locals()))

def rreplace(s, old, new, occurrence):
	li = s.rsplit(old, occurrence)
	return new.join(li)

##########################################################
# PULL POSITIONS BASED ON FAMILY AND METHOD
##########################################################
all_nonredundant=sys.argv[1]
ALL_NONREDUNDANT=open(all_nonredundant, "r")

output_files={}
output_files=defaultdict(list)

# create dictonary(key=file name of TE family and method, value=list of detected transposon info)
for line in ALL_NONREDUNDANT:
	line=line.rstrip('\n')
	items=re.split("[\t]",line)
	te_info=items[3]
	match=re.search("(.*)_(\w+-)?reference", te_info)
	family=match.group(1)
	method=items[6]
	file_name="{method}_{family}".format(**locals())
	output_files[file_name].append(line)

ALL_NONREDUNDANT.close()

##########################################################
# COLLAPSE INTO UNIQUE POSITIONS
##########################################################

for i in output_files.keys():
	value=output_files[i]
	OUT_FILE=open(i,"w")
	for te in value:
		OUT_FILE.write(te + '\n')
	OUT_FILE.close()
	#sort file by position
	#os.system("""sort -k1,1 -k2,2n {i} > tmp && mv tmp {i}""".format(**locals()))
	result, err = Popen(["""sort -k1,1 -k2,2n {i} > tmp && mv tmp {i}""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
	first_line = True
	collapsed_transposons={}
	te_ID=0
	
	IND_FILE=open(i, "r")

	for line in IND_FILE:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		chromosome=items[0]
		start_pos=items[1]
		end_pos=items[2]
		method=items[6]

		# ADD IN STRANDEDNESS

		if first_line == False:
			if chromosome == prev_chromosome:
				if method == "new": # allow for bp differences in insertion calls
					distance = 50
				else:
					distance = 0 # do not allow for bp differences in absence or reference call--should already be correct/exact
				if (int(start_pos)-int(prev_end_pos)) <= int(distance) :
					line=prevLine
					prevLine = rreplace(prevLine, "\t{prev_end_pos}".format(**locals()), "\t{end_pos}".format(**locals()), 1)  # replace last occurence...need to avoid number after
					#
					#
					#
					#
					#
					#
					#
					#REDO
					#collapsed_transposons[te_ID] = prevLine
					#reset prev end position
					prev_end_pos=end_pos
					#don't increase te_ID
				else:
					te_ID+=1
					prev_chromosome=chromosome
					prev_start_pos=start_pos
					prev_end_pos=end_pos
					collapsed_transposons[te_ID]=line
					prevLine=line
			else:
				te_ID+=1
				prev_chromosome=chromosome
				prev_start_pos=start_pos
				prev_end_pos=end_pos
				collapsed_transposons[te_ID]=line
				prevLine=line
		else:
			prev_chromosome=chromosome
			prev_start_pos=start_pos
			prev_end_pos=end_pos
			collapsed_transposons[te_ID]=line
			prevLine=line
			first_line=False
	#print collapsed transposons to a new file
	IND_FILE.close()
	final_out = "final" + "_" + i
	FINAL_OUT=open(final_out, "w")
	for ID in collapsed_transposons.keys():
		###TAKE  MIDPOINT HERE?!?!?!??!?!?!
		##
		#
		#
		#
		#
		TE = collapsed_transposons[ID]
		FINAL_OUT.write(TE + '\n')
	FINAL_OUT.close()
	result, err = Popen(["""sort -k1,1 -k2,2n {final_out} > tmp && mv tmp {final_out}""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()

##########################################################
# SORT POSITION FILES
##########################################################
result, err = Popen(["""cat final_* > cleaned_positions.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""sort -k1,1 -k2,2n cleaned_positions.gff > tmp && mv tmp cleaned_positions.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()

result, err = Popen(["""cat final_new* > cleaned_positions_new.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""sort -k1,1 -k2,2n cleaned_positions_new.gff > tmp && mv tmp cleaned_positions_new.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()

result, err = Popen(["""cat final_reference* > cleaned_positions_reference.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""sort -k1,1 -k2,2n cleaned_positions_reference.gff > tmp && mv tmp cleaned_positions_reference.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()

result, err = Popen(["""cat final_absent* > cleaned_positions_absent.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""sort -k1,1 -k2,2n cleaned_positions_absent.gff > tmp && mv tmp cleaned_positions_absent.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()

##########################################################
# RUN BEDTOOLS WINDOW ON ALL SAMPLES
##########################################################
result, err = Popen(["""bash {kin_step2} {sample_list}""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#output files from above step are named "insertions_bedt.txt" "absences_bedt.txt" "references_bedt.txt"

##########################################################
# GENERATE KINSHIP MATRIX
##########################################################
#ensure that the family found in bedtools matches that of the unique postion in the gff and output a matrix:
result, err = Popen(["""python {kin_step3} insertions_bedt.txt cleaned_positions_new.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""python {kin_step3} references_bedt.txt cleaned_positions_reference.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""python {kin_step3} absences_bedt.txt cleaned_positions_absent.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#transpose matrices:
result, err = Popen(["""bash {transpose_matrix} Samples_insertions_bedt.txt T_Samples_insertions_bedt.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""bash {transpose_matrix} Samples_references_bedt.txt T_Samples_references_bedt.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""bash {transpose_matrix} Samples_absences_bedt.txt T_Samples_absences_bedt.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()









