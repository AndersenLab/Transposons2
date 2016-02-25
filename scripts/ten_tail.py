#!/usr/bin/env python

import re

infile="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full.txt"

with open(infile, 'r') as IN:
	header=next(IN)
	header=header.rstrip('\n')
	strains=re.split('\t',header)
	print strains
	for i,j in enumerate(strains):
		strain="ID"+str(i)



	for line in IN:
		line=line.rstrip('\n')
		line=re.split('\t',line)
		trait=line[0]
		print trait
		trait={}
		for i,value in enumerate(line[1:len(line)],start=1):
			ID="ID"+str(i)
			trait[ID]=value

		# NEED THIS TO BE A NUENRICAL SORT

		#trait = dict((k,int(v)) for k,v in trait.iteritems())
		for key,value in trait.items():
			trait[key]=int(value)
		t=sorted(zip(trait.values(),trait.keys()))
		for i in t[0:10]:
			print i

		#just find tenth largest and tenth smallesy




#dict of families