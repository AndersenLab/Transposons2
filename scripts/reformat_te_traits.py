#!/usr/bin/env python 
# this script takes the base transposon traits for each strain and outputput additional measures of transposon counts and activity
# USE: reformat_te_traits.py <FINAL_RESULTS_LF>  # NEED TO EDIT DESCRIPTION OF THIS SCRIPT

import sys
import re
import os
from subprocess import Popen, PIPE

trait_list=()
te_input = sys.argv[1]
file_name=os.path.splitext(te_input)[0]
RESULTS_FILE_REFORMATTED=open("{file_name}_ReF.txt".format(**locals()), "w")
result, err = Popen(["""cat %s | cut -f2|uniq""" %(te_input)], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result=result.split()
result=list(set(result))
te_dict={}
# initiate a dictionary for all traits
for i in result:
	TE_INPUT = open(te_input, "r")
	for line in TE_INPUT:
		line = line.rstrip('\n')
		items = re.split("[\t]",line)
		strain_name = items[1]
		trait = items[2]
		result = items[3]
		if strain_name==i:
			te_dict[trait]=result

		#exec """{trait}[{strain}]={result}""".format(**locals())
	TE_INPUT.close()
	RESULTS_FILE_REFORMATTED.write(i)
	for TE in sorted(te_dict.keys()):
		value=te_dict[TE]
		print i
		print TE
		print value
		RESULTS_FILE_REFORMATTED.write("\t{value}".format(**locals()))
	RESULTS_FILE_REFORMATTED.write('\n')



	#trait[strain] = result
	#trait_list.append(trait)



#for  i in trait_list:
#	print i

#TE_INPUT.close()





#trait_list=()
#te_input = sys.argv[1]
#result, err = Popen(["""cat %s | cut -f3|uniq""" %(te_input)], stdout=PIPE, stderr=PIPE, shell=True).communicate()
#result=result.split()
# initiate a dictionary for all traits
#for i in result:
#	i={}

#new_TRANS_total={}

#TE_INPUT = open(te_input, "r")
#for line in TE_INPUT:
#	line = line.rstrip('\n')
#	items = re.split("[\t]",line)
#	strain = items[1]
#	trait = items[2]
#	result = items[3]

#	exec """{trait}[{strain}]={result}""".format(**locals())

#TE_INPUT.close()

	#trait[strain] = result
	#trait_list.append(trait)



#for  i in trait_list:
#	print i

#TE_INPUT.close()