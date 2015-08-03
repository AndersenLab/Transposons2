#!/usr/bin/env python 
# this script takes the base transposon traits for each strain and outputput additional measures of transposon counts and activity
# column names in correct order will be output to "column_names_activity.txt"
# USE: reformat_te_traits.py <Full_Results.txt>

import sys
import re


trait_list={}
family={}
te_input = sys.argv[1]

first_line=True
TE_INPUT = open(te_input, "r")
OUT_FILE = open("Full_Results_Activity.txt", "w")
COLUMN_NAMES = open("column_names_activity.txt", "w")

for line in TE_INPUT:
	strain_traits={}
	line = line.rstrip('\n')
	items = re.split("[\t]",line)
	strain_name = items[0]
	activity_strain_traits = {}
	if first_line:
		#nomber_traits=len(items)
		for i in range(1,(len(items))): ##CHECK!!!..include last and avoid first "trait"
			#put te family names into a dictionary
			TE=re.split("_TRANS_",items[i])
			TE_family=TE[1]
			family[TE_family]=0

			trait_list[i]=items[i]  #remove 1,triat value later
	else:
		for i in range(1,(len(items))): 
			#print i
			trait_name = trait_list[i]
			#print trait_name
			strain_traits[trait_name] = items[i]

		#calcualte  activity measurements as above per family
		for i in family.keys(): #this will already include total
			ref = strain_traits["reference_TRANS_{i}".format(**locals())]
			new = strain_traits["new_TRANS_{i}".format(**locals())]
			absent = strain_traits["absent_TRANS_{i}".format(**locals())]

			# assign activity values to a new dictionary
			activity_strain_traits["no_ref_plus_new_TRANS_{i}".format(**locals())] = int(ref) + int(new)
			activity_strain_traits["no_abs_plus_new_TRANS_{i}".format(**locals())] = int(absent) + int(new)
			activity_strain_traits["no_abs_minus_new_TRANS_{i}".format(**locals())] = int(absent) - int(new)
			activity_strain_traits["no_ref_plus_new_D2_TRANS_{i}".format(**locals())] = (float(ref) + float(new))/2 # divide by 2
			activity_strain_traits["no_abs_plus_new_D2_TRANS_{i}".format(**locals())] = (float(absent) + float(new))/2 # divide by 2
			activity_strain_traits["no_abs_minus_new_AV_TRANS_{i}".format(**locals())] = abs(int(absent) - int(new)) # absolute value

			#family count/total tes called for that particular method
			activity_strain_traits["no_ref_fam_over_total_TRANS_{i}".format(**locals())] = (float(ref)/float(strain_traits["reference_TRANS_total"]))*100
			activity_strain_traits["no_new_fam_over_total_TRANS_{i}".format(**locals())] = (float(new)/float(strain_traits["new_TRANS_total"]))*100
			activity_strain_traits["no_abs_fam_over_total_TRANS_{i}".format(**locals())] = (float(absent)/float(strain_traits["absent_TRANS_total"]))*100

			#family presence/absence data 
			if int(ref) > 0:
				activity_strain_traits["no_ref_PA_TRANS_{i}".format(**locals())] = 1 # presence/absence
			else:
				activity_strain_traits["no_ref_PA_TRANS_{i}".format(**locals())] = 0 # presence/absence

			if int(new) > 0:
				activity_strain_traits["no_new_PA_TRANS_{i}".format(**locals())] = 1 # presence/absence
			else:
				activity_strain_traits["no_new_PA_TRANS_{i}".format(**locals())] = 0 # presence/absence

			if int(absent) > 0:
				activity_strain_traits["no_abs_PA_TRANS_{i}".format(**locals())] = 1 # presence/absence
			else:
				activity_strain_traits["no_abs_PA_TRANS_{i}".format(**locals())] = 0 # presence/absence

		OUT_FILE.write(strain_name)
		for i in sorted(activity_strain_traits.keys()):
			value = activity_strain_traits[i]
			OUT_FILE.write("\t{value}".format(**locals()))
		OUT_FILE.write('\n')

	
		#print trait_name
		#print items[i]

	first_line = False
COLUMN_NAMES.write("trait")
for i in sorted(activity_strain_traits.keys()):
	COLUMN_NAMES.write("\t{i}".format(**locals()))
COLUMN_NAMES.write('\n')
print "DONE"
TE_INPUT.close()




		#no_ref_plus_new = int(ref) + int(new)
		#no_abs_plus_new = int(absent) + int(new)
		#no_abs_minus_new = int(absent) - int(new)
		#no_ref_plus_new_D2 = (int(ref) + int(new))/2 # divide by 2
		#no_abs_plus_new_D2 = (int(absent) + int(new))/2 # divide by 2
		#no_abs_minus_new_AV = abs(int(absent) - int(new)) # absolute value
		#activity_strain_traits[no_ref_plus_new] = no_ref_plus_new
		#activity_strain_traits[no_abs_plus_new] = no_abs_plus_new
		#activity_strain_traits[no_abs_minus_new] = no_abs_minus_new
		#activity_strain_traits[no_ref_plus_new_D2] = no_ref_plus_new_D2
		#activity_strain_traits[no_abs_plus_new_D2] =no_abs_plus_new_D2 
		#activity_strain_traits[no_abs_minus_new_AV] = no_abs_minus_new_AV