#!/usr/bin/env python
# this script:
# 1) outputs each TE insertion call into new files based on family
# 2) collapses TEs of same famly within 50 base pairs of one another
# 3) outputs all the unqiue TE positions to a new file
# 4) calculates the coverage for each sample at each insertion postion +/- 25 base pairs
# 5) for each unique position scores strains where a TE was not found as a 0 if coverage was above 8, or "NA" if coverage below 8
# 6) outputs a kinship matrix
# USE: kin_hash.py <CtCp_ll_nonredundant.txt>
# ex: kin_hash.py /lscr2/andersenlab/kml436/git_repos2/Transposons2/kin_hash/CtCp_all_nonredundant.txt

import sys
import re
import os
import statistics
from  collections import defaultdict
from subprocess import Popen, PIPE


#python ../scripts/kin_mean.py /lscr2/andersenlab/kml436/git_repos2/Transposons2/kintest/CtCp_all_nonredundant.txt 


kin_step2="/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/kin_step2.sh"
##NEED TO EDIT SAMPLE FILE IN BELOWSCRIPT TOO
kin_step3="/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/kin_temp.py"
transpose_matrix="/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/transpose_matrix.sh"

bam_dir="/lscr2/andersenlab/dec211/RUN/v2_snpset/bam"
#
#
#
#
#change smaple list to full one later
sample_list="/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/full_sample_list.txt"
#sample_list="/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/test_list.txt"
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
	if method=="new": #only need this method for the insertions
		output_files[file_name].append(line)

ALL_NONREDUNDANT.close()

##########################################################
# COLLAPSE INTO UNIQUE POSITIONS
##########################################################
FINAL_SAMPLES={}
FINAL_SAMPLES=defaultdict(list)
FINAL_POSITIONS={}
FINAL_POSITIONS=defaultdict(list)
MATRIX={}
MATRIX=defaultdict(list)
POSITION_COUNTER=0
for i in output_files.keys():

	TE_positions={}
	TE_positions=defaultdict(list)
	TE_samples={}
	TE_samples=defaultdict(list)
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
		sample=items[7]

		# ADD IN STRANDEDNESS

		if first_line == False:
			if chromosome == prev_chromosome:
				if method == "new": # allow for bp differences in insertion calls
					distance = 50
				else:
					distance = 0 # do not allow for bp differences in absence or reference call--should already be correct/exact
				if (int(start_pos)-int(prev_end_pos)) <= int(distance) :
					line=prevLine
					TE_positions[te_ID].append(start_pos)
					TE_samples[te_ID].append(sample)
					###add another dictiionaty
					#prevLine = rreplace(prevLine, "\t{prev_end_pos}".format(**locals()), "\t{end_pos}".format(**locals()), 1)  # replace last occurence...need to avoid number after
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
					TE_positions[te_ID].append(start_pos)
					TE_samples[te_ID].append(sample)
					prevLine=line
			else:
				te_ID+=1
				prev_chromosome=chromosome
				prev_start_pos=start_pos
				prev_end_pos=end_pos
				collapsed_transposons[te_ID]=line
				TE_positions[te_ID].append(start_pos)
				TE_samples[te_ID].append(sample)
				prevLine=line
		else:
			prev_chromosome=chromosome
			prev_start_pos=start_pos
			prev_end_pos=end_pos
			collapsed_transposons[te_ID]=line
			TE_positions[te_ID].append(start_pos)
			TE_samples[te_ID].append(sample)
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

		TE = collapsed_transposons[ID]
		items=re.split("[\t]", TE)
		chromosome=items[0]
		info='\t'.join(items[3:8]) # don't include end position in here because it should be the same as the start position

		
		TE_positions[ID] = map(float, TE_positions[ID]) #convert strings in list to floating point number
		# take the mean of the start positions,round it, and use the integer value

		average_start=int(round(statistics.mean(TE_positions[ID])))
		FINAL_OUT.write("{chromosome}\t{average_start}\t{average_start}\t{info}\n".format(**locals()))

		#add to the parent hashes because above block is only for individual TE families
		TE_samples[ID]=[i+":1" for i in TE_samples[ID]] # mark "1" for samples where that TE was found
		FINAL_SAMPLES[POSITION_COUNTER]=TE_samples[ID]
		FINAL_POSITIONS[POSITION_COUNTER].extend((chromosome, average_start, info))
		POSITION_COUNTER+=1

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
#result, err = Popen(["""bash {kin_step2} {sample_list}""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#output files from above step are named "insertions_bedt.txt" "absences_bedt.txt" "references_bedt.txt"

##########################################################
# GENERATE KINSHIP MATRIX
##########################################################
#ensure that the family found in bedtools matches that of the unique postion in the gff and output a matrix:
#result, err = Popen(["""python {kin_step3} insertions_bedt.txt cleaned_positions_new.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#result, err = Popen(["""python {kin_step3} references_bedt.txt cleaned_positions_reference.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#result, err = Popen(["""python {kin_step3} absences_bedt.txt cleaned_positions_absent.gff""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#transpose matrices:
#result, err = Popen(["""bash {transpose_matrix} Samples_insertions_bedt.txt T_Samples_insertions_bedt.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#result, err = Popen(["""bash {transpose_matrix} Samples_references_bedt.txt T_Samples_references_bedt.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#result, err = Popen(["""bash {transpose_matrix} Samples_absences_bedt.txt T_Samples_absences_bedt.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()

#for each unique transposon postion, in the coverage interval fil record its chromomse and interval 25 base pair up and down stream
FINAL_INTERVALS={}
FINAL_INTERVALS=defaultdict(list)
COVERAGE_INTERVALS=open("coverage_intervals.txt", "w")
for key,value in FINAL_POSITIONS.items():
	chromosome,av_pos=value[0:2]
	#av_pos=items[1]
	upstream=int(int(av_pos)-25)
	downstream=int(int(av_pos)+25)
	FINAL_INTERVALS[key].extend((av_pos,chromosome,upstream,downstream))
	COVERAGE_INTERVALS.write("{key}\t{av_pos}\t{chromosome}\t{upstream}\t{downstream}\n".format(**locals()))
COVERAGE_INTERVALS.close()

#sort the coverage interval file
result, err = Popen(["cat coverage_intervals.txt | sort -k3,3 -k4,4n > tmp && mv tmp coverage_intervals.txt"], stdout=PIPE, stderr=PIPE, shell=True).communicate()

#for each sample, calculate the coverage at each of those intervals
SAMPLE_COV=open("sample_coverages_and_positions.txt", 'w')
with open(sample_list,'r') as IN:
	for line in IN:
		sample=line.rstrip('\n')
		for key,value in FINAL_INTERVALS.items():
			av_pos,chromosome,start,end=value[0:4]
			result, err = Popen(["""samtools depth {bam_dir}/{sample}.bam -r {chromosome}:{start}-{end}|datamash mean 3""".format(**locals())],stdout=PIPE, stderr=PIPE, shell=True).communicate()
			if result:
				coverage=round(float(result),2)
			else:
				coverage=0

			#find the strains that have already been scored for that TE within a 1, these need to further scoring so analyzing the coerage at that position can be skipped
			strain_list=[(re.split(":", x))[0] for x in FINAL_SAMPLES[key]]

			# if the sample has not already been scored and the coverage is greater than 8, that sample is marked 0, if less than 8 it's marked with "NA"
			if sample not in strain_list:
				if coverage >= 8:
					FINAL_SAMPLES[key].append("{sample}:0".format(**locals()))
				elif coverage <8:
					FINAL_SAMPLES[key].append("{sample}:NA".format(**locals()))
				else:
					print("coverage not found...exiting....")
					sys.exit

			SAMPLE_COV.write("{key}\t{chromosome}\t{av_pos}\t{sample}\t{coverage}\n".format(**locals()))
SAMPLE_COV.close()


KIN_MATRIX=open("kin_matrix_ins.txt", 'w')
KIN_MATRIX.write("TE")

#print headers
value=FINAL_SAMPLES[1]
value=sorted(value)
for i in value:
	strain=(re.split(":", i))[0]
	KIN_MATRIX.write("\t{strain}".format(**locals()))
KIN_MATRIX.write('\n')

#output the final matrix
for key,value in FINAL_SAMPLES.items():
	value=sorted(value)

	full_info=FINAL_POSITIONS[key]
	full_info=map(str, full_info) #convert the integers to strings
	te_info='_'.join(full_info)
	te_info = te_info.replace("\t", "_")

	KIN_MATRIX.write(str(te_info)) #or output key ID here?
	for i in value:
		score=(re.split(":", i))[1]
		KIN_MATRIX.write("\t{score}".format(**locals()))
	KIN_MATRIX.write('\n')

KIN_MATRIX.close()

result, err = Popen(["""head -n1 kin_matrix_ins.txt  > tmp"""],stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""cat kin_matrix_ins.txt  |sed 1d| sort -t"_" -k1,1 -k2,2n >>tmp """],stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""mv tmp kin_matrix_ins.txt  """],stdout=PIPE, stderr=PIPE, shell=True).communicate()


#error check for smae TE at same position
#why getting only zeros
