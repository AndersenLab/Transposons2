#!/usr/bin/env python
# this script add GO term annotation to a file that has the transcript ID in the 6th column
#USE: add_GO_terms.py <WB_GO.txt> <infile>

import sys
import re
import os
from collections import defaultdict

annotations=defaultdict(list)

infile=sys.argv[1]
infile2=sys.argv[2]

with open(infile, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		GO=items[4]
		#transcript name sometime in columnn 11, bu if not, it's in column 3
		if items[10] !='':
			transcript=items[10]
		else:
			transcript=items[2]
		annotations[transcript].append(GO)

#for i,value in annotations.items():
#	print str(value)

pref=os.path.splitext(infile2)[0]
outfile=pref+"_GO.txt"
OUT=open(outfile, 'w')
with open(infile2, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		transcriptID=items[5]
		if transcriptID	in annotations.keys():
			GO_annotations=annotations[transcriptID]
			GO_annotations=','.join(GO_annotations)	
		else:
			GO_annotations="NA"

		OUT.write("{line}\t{GO_annotations}\n".format(**locals()))
OUT.close()