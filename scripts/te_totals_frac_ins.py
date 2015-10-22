#!/usr/bin/env python

import sys
import re
from collections import defaultdict


infile=sys.argv[1]
indexed_fam={}
all_families={}

OUT=open("kin_frac_matrix_ins.txt", 'w') #fractions
OUTC=open("kin_C_matrix_ins.txt", 'w') #counts

with open(infile, 'r') as IN:
	headers=next(IN)
	headers=headers.rstrip('\n')
	headers=re.split("[\t]",headers)

	for ID,i in enumerate(headers[1:len(headers)],start=1):
		match=re.search("\w+_\d+_(.*)_(\w+-)?reference.*", i)
		family=match.group(1)
		indexed_fam[ID]=family
		all_families[family]=0

	OUT.write("Strain")
	OUTC.write("Strain")
	#out insertion names
	for i in sorted(all_families.keys()):
		nameN="new_TRANS_"+i
		OUT.write("\tONE_"+nameN)
		OUT.write("\tZERO_"+nameN)
		OUTC.write("\tONE_"+nameN+"_C")
		OUTC.write("\tZERO_"+nameN+"_C")
	OUTC.write("\tONE_new_TRANS_total_C")
	OUT.write('\n')
	OUTC.write('\n')

	for line in IN:
		line=line.rstrip('\n')
		items=re.split("[\t]",line)
		strain=items[0]
		OUT.write(strain)
		OUTC.write(strain)
		fam_coverage={}
		fam_coverage=defaultdict(int)
		ONE_total={}
		ONE_total=defaultdict(int)
		ZERO_total={}
		ZERO_total=defaultdict(int)
		NA_total={}
		NA_total=defaultdict(int)
		ONE_C=0 # count the total number of insertions

		for ID,value in enumerate(items[1:len(items)],start=1):
			TE=indexed_fam[ID]
			if value=="1":
				fam_coverage[TE] +=1
				ONE_total[TE] +=1
				ONE_C +=1
			elif value=="0":
				fam_coverage[TE] +=1
				ZERO_total[TE] +=1
			elif value=="NA":
				NA_total[TE] +=1
			else:
				print "ERROR: No coverage found, exiting..."
				sys.exit()

		#out calls
		for i in sorted(all_families.keys()):
			if i not in fam_coverage.keys():
				ONE_frac="NA"
				ZERO_frac="NA"

				ONE="NA"
				ZERO="NA"
			else:
				ONE=ONE_total[i]
				ZERO=ZERO_total[i]
				NA=NA_total[i]

				frac=float(ONE)/float(ONE+ZERO+NA)	
				ONE_frac=round(frac,2)

				frac2=float(ZERO)/float(ONE+ZERO+NA)	
				ZERO_frac=round(frac2,2)

				frac3=float(NA)/float(ONE+ZERO+NA)	
				NA_frac=round(frac3,2)	

			#print(str(ONE)+'\t'+ str(ZERO)+'\t'+str(NA))
			#print(str(frac)+'\t'+ str(frac2)+'\t'+str(frac3))
			#print(str(ONE_frac)+'\t'+ str(ZERO_frac)+'\t'+str(NA_frac))

			OUT.write("\t{ONE_frac}".format(**locals()))
			OUT.write("\t{ZERO_frac}".format(**locals()))
			OUTC.write("\t{ONE}".format(**locals()))
			OUTC.write("\t{ZERO}".format(**locals()))
		OUTC.write("\t{ONE_C}".format(**locals()))
		OUT.write('\n')
		OUTC.write('\n')
OUT.close()
OUTC.close()



#ALTERNATIVE:
		#out calls
		#for i in sorted(all_families.keys()):
		#	if i in fam_coverage.keys():
		#		coverage=fam_coverage[i]
		#		ONE=ONE_total[i]
		#		ZERO=ZERO_total[i]
		#		frac=float(ONE)/float(coverage)	
		#		ONE_frac=round(frac,2)	
		#		frac2=float(ZERO)/float(coverage)	
		#		ZERO_frac=round(frac2,2)				
		#	else:
		#		NA=NA_total[i] #not using this yet, placeholder
		#		ONE_frac="NA"
		#		ZERO_frac="NA"
