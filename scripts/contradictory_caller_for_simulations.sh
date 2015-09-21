#!/bin/bash
# this script goes through the RSV simulations for the various coverage subsets and outputs contradictory call file
# with information on whether the TEMP or the TELCOATE call was a TP or FD
# USE: contradictory_caller_for_simulations.sh in directory of choice

scripts_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts
samples=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/run_list.txt
data_dir_tenth=/lscr2/andersenlab/kml436/RSV_simulations_Aug31_tenth/RSV_simulations_Aug31
data_dir_quarter=/lscr2/andersenlab/kml436/RSV_simulations_Aug31_quarter/RSV_simulations_Aug31
data_dir_half=/lscr2/andersenlab/kml436/RSV_simulations_Aug31_half/RSV_simulations_Aug31
data_dir_full=/lscr2/andersenlab/kml436/RSV_simulations_Aug31

directories=($data_dir_tenth $data_dir_quarter $data_dir_half $data_dir_full)
for data_dir in "${directories[@]}";

do
	python ${scripts_dir}/temp_telocate_compatibility_sim_version.py $data_dir $samples
	python ${scripts_dir}/add_TP_FD_to_contradictory_sim_calls.py  ${data_dir}/depth_TPFD_abs/summary_depth_TPFD_abs.txt ${data_dir}/depth_TPFD_ref/summary_depth_TPFD_ref.txt contradictory_calls.txt 
# write R script

	name=`echo $data_dir| awk -F/ '{print $(NF-1)}'`
	echo $name #full will get name "kml436"

	mv contra_call_TPFD.txt ${name}_contra_call_TPFD.txt
	mv contradictory_calls_summary.txt ${name}_contradictory_calls_summary.txt
	mv contradictory_calls.txt ${name}_contradictory_calls.txt
done


