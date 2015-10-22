#!/usr/bin/env python 
# this script takes the base transposon traits for each strain and outputput additional measures of transposon counts and activity
# column names in correct order will be output to "column_names_activity.txt"
# USE: activity_calculator.py <kin_C_matrix_full.txt> prviously:<Full_Results.txt>

import sys
import re
import os


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
	print strain_name
	if first_line:
		#nomber_traits=len(items)
		for i in range(1,(len(items))): ##CHECK!!!..include last and avoid first "trait"
			#put te family names into a dictionary
			TE=re.split("_TRANS_",items[i])
			TE_family=TE[1]
			family[TE_family]=0

			trait_list[i]=items[i]  #remove 1,trait value later
	else:
		for i in range(1,(len(items))): 
			#print i
			trait_name = trait_list[i]
			#print trait_name
			strain_traits[trait_name] = items[i]

		#calcualte  activity measurements as above per family
		for i in family.keys(): #this will already include total

			if "reference_TRANS_{i}".format(**locals()) in strain_traits.keys(): # reference TE wont have all repbase fam names
				ref = strain_traits["reference_TRANS_{i}".format(**locals())]
			else:
				ref="NA"
			if "ONE_new_TRANS_{i}".format(**locals()) in strain_traits.keys():
				new = strain_traits["ONE_new_TRANS_{i}".format(**locals())]
			else:
				new="NA"
			if "absent_TRANS_{i}".format(**locals()) in strain_traits.keys():
				absent = strain_traits["absent_TRANS_{i}".format(**locals())]
			else:
				absent="NA"

			# assign activity values to a new dictionary
			if ref =="NA" and new =="NA":
				activity_strain_traits["no_ref_plus_ONE_new_TRANS_{i}".format(**locals())] = "NA"
				activity_strain_traits["no_ref_plus_ONE_new_D2_TRANS_{i}".format(**locals())] = "NA" # divide by 2
			elif ref=="NA": #as long as one value is not NA, can set the NA value to zero in order to do the calculations
				activity_strain_traits["no_ref_plus_ONE_new_TRANS_{i}".format(**locals())] = int(0) + int(new)
				activity_strain_traits["no_ref_plus_ONE_new_D2_TRANS_{i}".format(**locals())] = (float(0) + float(new))/2 # divide by 2
			elif new=="NA":
				activity_strain_traits["no_ref_plus_ONE_new_TRANS_{i}".format(**locals())] = int(ref) + int(0)
				activity_strain_traits["no_ref_plus_ONE_new_D2_TRANS_{i}".format(**locals())] = (float(ref) + float(0))/2 # divide by 2
			else: # calcualte normally if neither value is NA
				activity_strain_traits["no_ref_plus_ONE_new_TRANS_{i}".format(**locals())] = int(ref) + int(new)
				activity_strain_traits["no_ref_plus_ONE_new_D2_TRANS_{i}".format(**locals())] = (float(ref) + float(new))/2 # divide by 2


			if absent =="NA" and new =="NA":
				activity_strain_traits["no_abs_plus_ONE_new_TRANS_{i}".format(**locals())] = "NA"
				activity_strain_traits["no_abs_minus_ONE_new_TRANS_{i}".format(**locals())] = "NA"
				activity_strain_traits["no_abs_plus_ONE_new_D2_TRANS_{i}".format(**locals())] = "NA" # divide by 2
				activity_strain_traits["no_abs_minus_ONE_new_AV_TRANS_{i}".format(**locals())] = "NA"
			elif absent=="NA":
				activity_strain_traits["no_abs_plus_ONE_new_TRANS_{i}".format(**locals())] = int(0) + int(new)
				activity_strain_traits["no_abs_minus_ONE_new_TRANS_{i}".format(**locals())] = int(0) - int(new)
				activity_strain_traits["no_abs_plus_ONE_new_D2_TRANS_{i}".format(**locals())] = (float(0) + float(new))/2 # divide by 2
				activity_strain_traits["no_abs_minus_ONE_new_AV_TRANS_{i}".format(**locals())] = abs(int(0) - int(new)) # absolute value
			elif new=="NA":
				activity_strain_traits["no_abs_plus_ONE_new_TRANS_{i}".format(**locals())] = int(absent) + int(0)
				activity_strain_traits["no_abs_minus_ONE_new_TRANS_{i}".format(**locals())] = int(absent) - int(0)
				activity_strain_traits["no_abs_plus_ONE_new_D2_TRANS_{i}".format(**locals())] = (float(absent) + float(0))/2 # divide by 2
				activity_strain_traits["no_abs_minus_ONE_new_AV_TRANS_{i}".format(**locals())] = abs(int(absent) - int(0)) # absolute value
			else: # calcualte normally if neither value is NA
				activity_strain_traits["no_abs_plus_ONE_new_TRANS_{i}".format(**locals())] = int(absent) + int(new)
				activity_strain_traits["no_abs_minus_ONE_new_TRANS_{i}".format(**locals())] = int(absent) - int(new)
				activity_strain_traits["no_abs_plus_ONE_new_D2_TRANS_{i}".format(**locals())] = (float(absent) + float(new))/2 # divide by 2
				activity_strain_traits["no_abs_minus_ONE_new_AV_TRANS_{i}".format(**locals())] = abs(int(absent) - int(new)) # absolute value

			if ref != "NA":
				if int(ref)==0:
					activity_strain_traits["no_ref_fam_over_total_TRANS_{i}".format(**locals())] = 0
					activity_strain_traits["no_ref_PA_TRANS_{i}".format(**locals())] = 0
				else:
					activity_strain_traits["no_ref_fam_over_total_TRANS_{i}".format(**locals())] = round((float(ref)/float(strain_traits["reference_TRANS_total_C"]))*100,2)
					activity_strain_traits["no_ref_PA_TRANS_{i}".format(**locals())] = 1
			else:
				activity_strain_traits["no_ref_fam_over_total_TRANS_{i}".format(**locals())] = "NA"
				activity_strain_traits["no_ref_PA_TRANS_{i}".format(**locals())] ="NA"

			if new != "NA":
				if int(new)==0:
					activity_strain_traits["no_ONE_new_fam_over_total_TRANS_{i}".format(**locals())] = 0
					activity_strain_traits["no_ONE_new_PA_TRANS_{i}".format(**locals())] = 0
				else:
					activity_strain_traits["no_ONE_new_fam_over_total_TRANS_{i}".format(**locals())] = round((float(new)/float(strain_traits["ONE_new_TRANS_total_C"]))*100,2)
					activity_strain_traits["no_ONE_new_PA_TRANS_{i}".format(**locals())] = 1
			else: 
				activity_strain_traits["no_ONE_new_fam_over_total_TRANS_{i}".format(**locals())]= "NA"
				activity_strain_traits["no_ONE_new_PA_TRANS_{i}".format(**locals())] = "NA"

			if absent !="NA":
				if int(absent)==0:
					activity_strain_traits["no_abs_fam_over_total_TRANS_{i}".format(**locals())] = 0
					activity_strain_traits["no_abs_PA_TRANS_{i}".format(**locals())] = 0
				else:
					activity_strain_traits["no_abs_fam_over_total_TRANS_{i}".format(**locals())] = round((float(absent)/float(strain_traits["absent_TRANS_total_C"]))*100,2)
					activity_strain_traits["no_abs_PA_TRANS_{i}".format(**locals())] = 1
			else:
				activity_strain_traits["no_abs_fam_over_total_TRANS_{i}".format(**locals())] ="NA"
				activity_strain_traits["no_abs_PA_TRANS_{i}".format(**locals())] = "NA"


			#activity_strain_traits["no_ref_plus_new_TRANS_{i}".format(**locals())] = int(ref) + int(new)
			#activity_strain_traits["no_ref_plus_new_D2_TRANS_{i}".format(**locals())] = (float(ref) + float(new))/2 # divide by 2
			#activity_strain_traits["no_abs_plus_new_TRANS_{i}".format(**locals())] = int(absent) + int(new)
			#activity_strain_traits["no_abs_minus_new_TRANS_{i}".format(**locals())] = int(absent) - int(new)
			#activity_strain_traits["no_abs_plus_new_D2_TRANS_{i}".format(**locals())] = (float(absent) + float(new))/2 # divide by 2
			#activity_strain_traits["no_abs_minus_new_AV_TRANS_{i}".format(**locals())] = abs(int(absent) - int(new)) # absolute value



			#family count/total tes called for that particular method
			#activity_strain_traits["no_ref_fam_over_total_TRANS_{i}".format(**locals())] = (float(ref)/float(strain_traits["reference_TRANS_total"]))*100
			#activity_strain_traits["no_new_fam_over_total_TRANS_{i}".format(**locals())] = (float(new)/float(strain_traits["new_TRANS_total"]))*100
			#activity_strain_traits["no_abs_fam_over_total_TRANS_{i}".format(**locals())] = (float(absent)/float(strain_traits["absent_TRANS_total"]))*100

			#family presence/absence data 
			#if int(ref) > 0:
			#	activity_strain_traits["no_ref_PA_TRANS_{i}".format(**locals())] = 1 # presence/absence
			#else:
			#	activity_strain_traits["no_ref_PA_TRANS_{i}".format(**locals())] = 0 # presence/absence

			#if int(new) > 0:
			#	activity_strain_traits["no_new_PA_TRANS_{i}".format(**locals())] = 1 # presence/absence
			#else:
			#	activity_strain_traits["no_new_PA_TRANS_{i}".format(**locals())] = 0 # presence/absence

			#if int(absent) > 0:
			#	activity_strain_traits["no_abs_PA_TRANS_{i}".format(**locals())] = 1 # presence/absence
			#else:
			#	activity_strain_traits["no_abs_PA_TRANS_{i}".format(**locals())] = 0 # presence/absence

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
COLUMN_NAMES.close()
OUT_FILE.close()

os.system(" cat column_names_activity.txt Full_Results_Activity.txt > tmp && mv tmp Full_Results_Activity.txt")






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