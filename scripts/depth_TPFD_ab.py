#!/usr/bin/env python 
# this script goes through a simulation set and for every detected position made by the transposon caller outputs the interval 20 bp upstream and downstream of that position
# with the average coverage acorss that interval and wherther the call was a true or false positive
# USE: depth_TPFD.py **change "my_dir" to specify the simulation set
# NOTE: change "my_dir" to the directory of interest

import os
import re
from subprocess import Popen, PIPE
import statistics
my_dir="/lscr2/andersenlab/kml436/RSV_simulations_Aug31_half/RSV_simulations_Aug31"

means=[]

###CHANGE TO 1 through 8
run_list=[1,2,3,4,5,6,7,8]
for i in run_list:
	print "Processing Run:"
	print i

	run=i
	bam_file="{my_dir}/run_{run}/merged_bam_{run}.sorted.bam".format(**locals()) 
	call_file="{my_dir}/run_{run}/run_{run}_filter_results_temp/intermediates_RECALCULATED/new_CT_run_{run}_temp_nonredundant.bed".format(**locals())
	final_positions="{my_dir}/run_{run}/run_{run}_filter_results_temp/adjusted_pos_SIM_tes_{run}.bed".format(**locals())###
	sc_file="{my_dir}/run_{run}/run_{run}_filter_results_temp/intermediates_RECALCULATED/SC_new_CT_run_{run}_temp_nonredundant.bed_to_adjusted_pos_SIM_tes_{run}.bed_compare.bed".format(**locals())
	SC=open(sc_file, "r")
	true_positives={}
	true_pos_position={}

	# add true positive calls (within 20 base and with a matching family to a dictionary )
	for line in SC:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		distance=items[12]
		true_fam=items[3]
		te=items[9]
		match = re.search("(.*)_(\w+-)?reference", te)
		found_fam=match.group(1)
		te_info="_".join(items[6:10])
		pos_info="_".join(items[0:4]) ###
		if int(distance)<20: #make sure distance is within 20 bp
			if found_fam==true_fam:	# make sure families match
				true_positives[te_info]=0
				true_pos_position[pos_info]=0###
	SC.close()


	CALL_FILE=open(call_file, "r")
	OUT=open("{run}_depth_TPDF.bed".format(**locals()), "w")

	#add coverage and true of false positive info to the new_CT_temp_nonredundant file
	for line in CALL_FILE:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		chromosome=items[0]
		start=items[1]
		#get position 20 bp up and downstream of the start position
		# dependent on chromosome NOT being at the very ends of the chromosome..change later if this becomes an issue
		left_interval= int(start)-20
		right_interval = int(start)+20
		end=items[2]
		te=items[3]
		N=items[4]
		orient=items[5]
		te_info2="_".join(items[0:4])
		result, err = Popen(["""samtools depth {bam_file} -r {chromosome}:{left_interval}-{right_interval}| datamash mean 3""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
		 # result is the mean coverage over the interval
		if result=="": # if samtools depth returns no coverage, set result equal to zero
			result='0\n'

		means.append(float(result))


		if te_info2 in true_positives.keys():
			TPFD="TP"
		else:
			TPFD="FD"
		OUT.write("{run}\t{chromosome}\t{left_interval}\t{right_interval}\t{te}\t{N}\t{orient}\t{TPFD}\t{result}".format(**locals()))

	CALL_FILE.close()
	OUT.close()
	FINAL_POSITIONS=open(final_positions, "r")###
	OUT2=open("{run}_depth_TPDF.bed_FN".format(**locals()), "w")
	print "Calculating interval coverages..."
	for line in FINAL_POSITIONS:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		chromosome=items[0]
		start=items[1]
		#get position 20 bp up and downstream of the start position
		# dependent on chromosome NOT being at the very ends of the chromosome..change later if this becomes an issue
		left_interval= int(start)-20
		right_interval = int(start)+20
		end=items[2]
		te=items[3]
		N=items[4]
		orient=items[5]
		pos_info2="_".join(items[0:4])
		result, err = Popen(["""samtools depth {bam_file} -r {chromosome}:{left_interval}-{right_interval}| datamash mean 3""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
		 # result is the mean coverage over the interval
		if result=="": # if samtools depth returns no coverage, set result equal to zero
			result='0\n'
		
		
		
		if pos_info2 not in true_pos_position.keys():
			TPFD="FN"
			OUT2.write("{run}\t{chromosome}\t{left_interval}\t{right_interval}\t{te}\t{N}\t{orient}\t{TPFD}\t{result}".format(**locals()))


	FINAL_POSITIONS.close()
	OUT2.close()

final_mean=statistics.mean(means)
final_SD=statistics.stdev(means)
one_SDs=final_mean+(1*final_SD)
two_SDs=final_mean+(2*final_SD)
three_SDs=final_mean+(3*final_SD)
four_SDs=final_mean+(4*final_SD)


MEAN_COVERAGE=open("mean_coverage_and_sd.txt","w")
MEAN_COVERAGE.write("Mean\tOneSD\tTwoSD\tThreeSD\tFourSD\n")
MEAN_COVERAGE.write("{final_mean}\t{one_SDs}\t{two_SDs}\t{three_SDs}\t{four_SDs}\n".format(**locals()))


header="chromosome\tleft_interval\tright_interval\tTE\tN\torient\tcall\tcoverage"
result, err = Popen(["""cat *.bed > summary_depth_TPFD_abs.txt"""], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""cat *.bed_FN > summary_depth_TPFD_FN.txt"""], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(["""echo '{header}' | cat - summary_depth_TPFD_abs.txt > temp && mv temp summary_depth_TPFD_abs.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen([""" cat summary_depth_TPFD_abs.txt summary_depth_TPFD_FN.txt > temp && mv temp summary_depth_TPFD_abs.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
