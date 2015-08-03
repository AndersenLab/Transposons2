#!/usr/bin/env python

import sys
import re

reference_in=sys.argv[1]
absence_in=sys.argv[2]

REFERENCE_IN=open(reference_in, "r")
ABSENCE_IN=open(absence_in, "r")

reference_counts={}
absence_counts={}
#count the reference instances
first_line=True
for line in REFERENCE_IN:
	count_ones=0
	count_zeros=0
	if first_line==False: #don't count first line
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		number=len(items)
		te_id=items[0]
		te_items=re.split("_", te_id)
		te_id='_'.join(te_items[0:4])
		#print te_id
		for i in items[1:126]:
			if i=="1":
				count_ones+=1
			if i=="0":
				count_zeros+=1
		total=int(count_ones) + int(count_zeros)
		if int(count_ones) + int(count_zeros) != 124: # make sure all 124 samples are accounted for each reference transposition event
			print "Error..samples information is missing"
		reference_counts[te_id]=int(count_ones)
	first_line=False
REFERENCE_IN.close()

#coutn the absence instances
first_line=True
for line in ABSENCE_IN:
	count_ones=0
	count_zeros=0
	if first_line==False: #don't count first line
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		number=len(items)
		te_id=items[0]
		te_items=re.split("_", te_id)
		te_id='_'.join(te_items[0:4])

		for i in items[1:126]:
			if i=="1":
				count_ones+=1
			if i=="0":
				count_zeros+=1
		total=int(count_ones) + int(count_zeros)
		if int(count_ones) + int(count_zeros) != 124: # make sure all 124 samples are accounted for each reference transposition event
			print "Error..samples information is missing"
		absence_counts[te_id]=int(count_ones)
	first_line=False

ABSENCE_IN.close()

found_already={}

for i in reference_counts.keys():
	#ref=reference_counts[i]
	if i not in absence_counts.keys():
		absence_counts[i]=0
		#print "YES"
		#ab=absence_counts[i]
	#else:
		#absence_counts[i]=0
	found_already[i]=0

for i in absence_counts.keys():
	if i not in found_already.keys():
		reference_counts[i]=0

OUT=open("total_ref_ab_counts.txt", "w")
for i in found_already.keys():
	ref=reference_counts[i]
	ab=absence_counts[i]
	OUT.write("{i}\t{ref}\t{ab}\n".format(**locals()))


