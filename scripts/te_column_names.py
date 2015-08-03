#!/usr/bin/env python 
# this script takes the base transposon traits for each strain and outputput additional measures of transposon counts and activity
# USE: reformat_te_traits.py <FINAL_RESULTS_LF>

import sys
import re
import os
from subprocess import Popen, PIPE

trait_list=()
te_input = sys.argv[1]
RESULTS_FILE=open("column_names.txt", "w")
result, err = Popen(["""cat %s | cut -f2|uniq""" %(te_input)], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result=result.split()
te_dict={}
# initiate a dictionary for all traits

TE_INPUT = open(te_input, "r")
for line in TE_INPUT:
	line = line.rstrip('\n')
	items = re.split("[\t]",line)
	strain = items[1]
	trait = items[2]
	result = items[3]
	te_dict[trait]=result

	#exec """{trait}[{strain}]={result}""".format(**locals())
TE_INPUT.close()

for TE in sorted(te_dict.keys()):
	RESULTS_FILE.write("{TE}\t".format(**locals()))