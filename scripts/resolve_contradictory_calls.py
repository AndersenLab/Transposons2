#!/usr/bin/env python 

import re

####CHANGE BELOW LATeR
data_dir="/lscr2/andersenlab/kml436/git_repos2/Transposons2/data"
sample_file="/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/full_sample_list.txt"
SAMPLE_FILE=open(sample_file, "r")
for line in SAMPLE_FILE:
	strain =line.rstrip('\n')


	telocate_removals={}
	temp_removals={}
	contra_calls="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/contradictory_calls.txt"
	CONTRA_CALLS=open(contra_calls, "r")
	for line in CONTRA_CALLS:
		line =line.rstrip('\n')
		items=re.split("[\t]", line)
		strain_name=items[0]
		if strain_name==strain:
			telocate_support=int(items[5])
			temp_support=int(items[8])
			tel_id='_'.join(items[1:5])
			temp_id='_'.join(items[1:4])
			temp_id = temp_id + "_" + items[7]

			if temp_support>telocate_support:
				telocate_removals[tel_id]=0

			elif temp_support==telocate_support:
				telocate_removals[tel_id]=0
			else:
				temp_removals[temp_id]=0

	CONTRA_CALLS.close()

	#for i in telocate_removals.keys():
	#	print i
	fileIN="{data_dir}/{strain}/final_results/org_{strain}_telocate_nonredundant.bed_org".format(**locals())
	fileOUT="{data_dir}/{strain}/final_results/{strain}_telocate_nonredundant.bed".format(**locals())

	ORG_TEL=open(fileIN, "r")
	NEW_TEL=open(fileOUT, "w")

	for line in ORG_TEL:
		line =line.rstrip('\n')
		items=re.split("[\t]", line)
		tel_call='_'.join(items[0:4])
		if tel_call not in telocate_removals.keys():
			NEW_TEL.write(line)
			NEW_TEL.write("\n")
		#else:
			#print tel_call

	ORG_TEL.close()
	NEW_TEL.close()

	fileIN2="{data_dir}/{strain}/final_results/org_{strain}_temp_absence_nonredundant.bed_org".format(**locals())
	fileOUT2="{data_dir}/{strain}/final_results/{strain}_temp_absence_nonredundant.bed".format(**locals())

	ORG_TEMP=open(fileIN2, "r")
	NEW_TEMP=open(fileOUT2, "w")

	for line in ORG_TEMP:
		line =line.rstrip('\n')
		items=re.split("[\t]", line)
		temp_call='_'.join(items[0:4])
		if temp_call not in temp_removals.keys():
			NEW_TEMP.write(line)
			NEW_TEMP.write("\n")
		#else:
			#print temp_call


	ORG_TEMP.close()
	NEW_TEMP.close()

