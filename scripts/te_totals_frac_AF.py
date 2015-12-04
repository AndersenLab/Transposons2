#!/usr/bin/env python

import sys
import re
from collections import defaultdict


infile=sys.argv[1]
indexed_fam={}
all_families={}

OUT=open("kin_frac_matrix_AF.txt", 'w')
OUTC=open("kin_C_matrix_AF.txt", 'w')

with open(infile, 'r') as IN:
	headers=next(IN)
	headers=headers.rstrip('\n')
	headers=re.split("[\t]",headers)

	for ID,i in enumerate(headers[1:len(headers)],start=1):
		match=re.search("\w+_\d+_(.*)", i)
		family=match.group(1)
		indexed_fam[ID]=family
		all_families[family]=0

	OUT.write("trait")
	OUTC.write("trait")
	#out reference and absence names
	for i in sorted(all_families.keys()):
		nameR="reference_TRANS_"+i
		OUT.write("\t"+nameR)
		OUTC.write("\t"+nameR+"_C")
		nameA="absent_TRANS_"+i
		OUT.write("\t"+nameA)
		OUTC.write("\t"+nameA+"_C")
	OUTC.write("\treference_TRANS_total_C")
	OUTC.write("\tabsent_TRANS_total_C")
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
		REF_total={}
		REF_total=defaultdict(int)
		ABS_total={}
		ABS_total=defaultdict(int)
		NA_total={}
		NA_total=defaultdict(int)
		REF_C=0 # count the total number of references
		ABS_C=0 # count the total number of absences

		for ID,value in enumerate(items[1:len(items)],start=1):
			TE=indexed_fam[ID]
			if value=="1":
				fam_coverage[TE] +=1
				REF_total[TE] +=1
				REF_C +=1
			elif value=="0":
				fam_coverage[TE] +=1
				ABS_total[TE] +=1
				ABS_C +=1
			elif value=="NA":
				NA_total[TE] +=1
			else:
				print "ERROR: No coverage found, exiting..."
				sys.exit()

		#reference calls
		for i in sorted(all_families.keys()):
			if i not in fam_coverage.keys():
				REF_frac="NA"
				ABS_frac="NA"

				REF="NA"
				ABS="NA"
			else:
				REF=REF_total[i]
				ABS=ABS_total[i]
				NA=NA_total[i]

				frac=float(REF)/float(REF+ABS+NA)	
				REF_frac=round(frac,2)

				frac2=float(ABS)/float(REF+ABS+NA)	
				ABS_frac=round(frac2,2)

				frac3=float(NA)/float(REF+ABS+NA)	
				NA_frac=round(frac3,2)	



			#print(str(REF)+'\t'+ str(ABS)+'\t'+str(NA))
			#print(str(frac)+'\t'+ str(frac2)+'\t'+str(frac3))
			#print(str(REF_frac)+'\t'+ str(ABS_frac)+'\t'+str(NA_frac))

			OUT.write("\t{REF_frac}".format(**locals()))
			OUT.write("\t{ABS_frac}".format(**locals()))
			OUTC.write("\t{REF}".format(**locals()))
			OUTC.write("\t{ABS}".format(**locals()))
		OUTC.write("\t{REF_C}".format(**locals()))
		OUTC.write("\t{ABS_C}".format(**locals()))
		OUT.write('\n')
		OUTC.write('\n')
OUT.close()
OUTC.close()


#ALTERNATIVE:
# 		#reference calls
# 		for i in sorted(all_families.keys()):
# 			if i in fam_coverage.keys():
# 				coverage=fam_coverage[i]
# 				total=fam_total[i]
# 				frac=float(total)/float(coverage)	
# 				rounded_frac=round(frac,2)				
# 			else:
# 				rounded_frac="NA"
# 			OUT.write("\t{rounded_frac}".format(**locals()))

# 		#absence calls
# 		for i in sorted(all_families.keys()):
# 			if i in fam_coverage.keys():
# 				coverage=fam_coverage[i]
# 				total=fam_total[i]
# 				frac=1-(float(total)/float(coverage))	
# 				rounded_frac=round(frac,2)		
# 			else:
# 				rounded_frac="NA"
# 			OUT.write("\t{rounded_frac}".format(**locals()))
