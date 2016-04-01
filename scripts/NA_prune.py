#!/usr/bin/env python
import re
import sys
from collections import defaultdict

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/kinship/T_kin_C_matrix_NAs_reduced.txt"
results_file="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full_reduced.txt"
OUT=open("pruned_data.txt",'w')

counts=defaultdict(list)
with open(infile, 'r') as IN:
	head_infile=next(IN)

	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		trait=items[0]
		trait= trait=re.sub('_NA$','',trait)
		for strain in items[1:]:
			counts[trait].append(strain) #this will be zero through 151 strains

with open(results_file, 'r') as IN:
	head_result_file=next(IN)
	OUT.write(head_result_file)
	if head_result_file != head_infile: # make sure that headers between the two input files matched
		sys.exit("HEADERS DO NOT MATCH")
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t', line)
		trait_name=items[0]
		trait_name=re.sub('_C$','',trait_name)
		previous_sites=False
		if trait_name in counts.keys(): #insertion traits
			OUT.write(trait_name)
			for a,strain_value in enumerate(items[1:],start=0):

				NA_perc=counts[trait_name]
				info=(counts[trait_name])[a]
				
				match=re.search("(.*)\(\d+/(\d+)\)",info)
				percentage=float(match.group(1))
				total_sites=int(match.group(2))

				if previous_sites==True:
					if total_sites != previous_sites:
						sys.exit("Total potential sites do not match between strains...exiting...")
					previous_sites=total_sites
				else:
					previous_sites=total_sites
				#check that all toals are the same per line

				if total_sites >10:
					if percentage<.10:
						strain_value=strain_value
					else:
						strain_value="NA"

				else: # when total_strain is <=10
					if percentage==0.0:
						strain_value=strain_value
					else:
						strain_value="NA" #strains with low coverage across trait sites and all strain values that were already NAs will be "NA"
					print strain_value
				OUT.write('\t' + strain_value)
			OUT.write('\n')
		elif re.search("total", trait_name) or \
		re.search("absent",trait_name) or \
		re.search("reference", trait_name) or \
		re.search("coverage",trait_name):
			pass
			OUT.write(line + '\n')
		else:
			print "Check trait consistency...."


OUT.close()


# Rules:
#Traits with a count  <= 10 should remove all strains with any NAs.
#Traits with a count >10 should remove all strains with NA counts greater than or equal to 10% rounded up. 

				






