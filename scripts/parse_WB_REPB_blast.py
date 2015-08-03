#this script ........
import sys
import re
import os
from subprocess import Popen, PIPE

top_blast_hits = sys.argv[1]
TOP_BLAST_HITS = open(top_blast_hits, "r")

out = "WB_REPB_blast_family_comparison.txt"
OUT = open(out, "w")
OUT.write("WB_TE\tWB_Family\tBlast_Hit\tPercent_Identity\tAlng_Length\tMismatches\tGap_Open\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_Value\tBit_Score\n")

IDs_in_blast={}
for line in TOP_BLAST_HITS:
	items= re.split("[\t]",line)
	transposon_ID= items[0]
	blast_fields= '\t'.join(items[1:12])
	print blast_fields
	IDs_in_blast[transposon_ID] = blast_fields


wb_tes = sys.argv[2]
WB_TES = open(wb_tes, "r")
for line in WB_TES:
	line= line.rstrip('\n')
	items = re.split("[\t]",line)
	WB_ID = items[0]
	family = items[1]

	if WB_ID in IDs_in_blast.keys():
		value = IDs_in_blast[WB_ID]
		OUT.write("{WB_ID}\t{family}\t{value}".format(**locals()))
	else:
		OUT.write("{WB_ID}\t{family}\tNo_Hit\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t9999\tNA\n".format(**locals()))