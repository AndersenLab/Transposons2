#!/usr/bin/env python 
# this script adds TP and FD annotation to the contradictory call file for a simulation
# USE: add_TP_FD_to_contradictory_sim_calls.py <summary_depth_TPFD_abs.txt> <summary_depth_TPFD_abs.txt> <contradictory_calls.txt> 

###python /lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts/add_TP_FD_to_contradictory_sim_calls.py ../RSV_simulations_Aug31_tenth/RSV_simulations_Aug31/depth_TPFD_abs/summary_depth_TPFD_abs.txt ../RSV_simulations_Aug31_tenth/RSV_simulations_Aug31/depth_TPFD_ref/summary_depth_TPFD_ref.txt contradictory_calls.txt 

import sys
import re 
tpfd_a=sys.argv[1] 
tpfd_r=sys.argv[2] 
cc=sys.argv[3]

abs_dict={}
ref_dict={}

#go through absence depth TPFD file first
with open(tpfd_a, 'r') as TPFD_A:
	for line in TPFD_A:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		num=items[0] #the run ID number
		te=items[4] #the transposon ID
		annotation=items[7] #whether the call was TP or FD
		if num != "chromosome": #skip first line
			ID= "run_" + num+ "_" + te
			abs_dict[ID]=annotation


#go through reference depth TPFD file first
with open(tpfd_r, 'r') as TPFD_R:
	for line in TPFD_R:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		num=items[0] #the run ID number
		te=items[4] #the transposon ID
		annotation=items[7] #whether the call was TP or FD
		if num != "chromosome": #skip first line
			ID= "run_" + num+ "_" + te
			ref_dict[ID]=annotation


CONTRA_CALLS=open("contra_call_TPFD.txt",'w')


with open(cc, "r") as CC:
	for line in CC:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		run=items[0]
		temp=items[4]
		telocate=items[7]

		temp_ID = run + "_" + temp
		telocate_ID = run + "_" + telocate
		print temp_ID
		print telocate_ID

		if temp_ID in abs_dict.keys() and telocate_ID in ref_dict.keys():
			temp_value=abs_dict[temp_ID]
			telocate_value=ref_dict[telocate_ID]
			CONTRA_CALLS.write("{line}\t{temp_value}\t{telocate_value}\n".format(**locals()))




#with open(cc) as CC: