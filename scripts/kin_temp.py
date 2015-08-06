#!/usr/bin/env python
# this script determines if takes bedtools window output, checks that the family matches between the given and found TE
# and output presence/absence info for each sample and unique TE (chromosome_startPos_endPos_TE)
# USE: kin_temp.py <bedt_file> <cleaned_positions_file>
# ex: kin_temp.py insertions_bedt.txt cleaned_positions_new.gff
import re
import sys
from subprocess import Popen, PIPE

file_name=sys.argv[1] #insertions_bedt.txt
file_pos=sys.argv[2] #cleaned_positions_new.gff
sample_file="/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/full_sample_list.txt"
result, err = Popen(["""cat {file_pos} | awk '{{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8}}'""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate() ###too out$7 to $6

insertion_tes={}
result=result.split()
for i in result:
	print i
	insertion_tes[i]=0


result, err = Popen(["""cat {sample_file} | cut -f1""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
samples={}
result=result.split()
for i in result:
	i=i.rstrip('\n')
	samples[i]=0

OUT_SAMPLE=open("Samples_" + file_name, "w")
OUT_SAMPLE.write("sample\t")
for i in sorted(insertion_tes.keys()):
	OUT_SAMPLE.write("{i}\t".format(**locals()))
OUT_SAMPLE.write('\n')
###FOR EACH STRAIN
for sample in samples.keys():
	print sample
	OUT_SAMPLE.write(sample)
	BEDT_FILE=open(file_name, "r")
	sampleD={}
	for line in BEDT_FILE:
		line=line.rstrip('\n')
		items=re.split("[\t]", line)
		strain=items[0]
		given_fam_info=items[5]
		found_fam_info=items[13] #changed here
		if sample==strain:
			#print "YES"
			#print line
			match=re.search("(.*)_(?:rp|sr)_(\d)+",given_fam_info)
			#print given_fam_info
			given_fam=match.group(1)
			#print given_fam
			#print "NEXT"
			#print found_fam_info
			match2=re.search("(.*)_(?:rp|sr)_(\d)+",found_fam_info)
			found_fam=match2.group(1)

			#print found_fam
			if given_fam==found_fam:


				full_info="_".join(items[2:10])

				#full_info=full_info + "_" + items[8]
				#print full_info
				sampleD[full_info]=0
		#for k in sampleD.keys():
			#print k
	BEDT_FILE.close()
	for i in sorted(insertion_tes.keys()):
		if i in sampleD.keys():
			print "YESSSSSSSSSS"
			OUT_SAMPLE.write("\t1")
		else:
			OUT_SAMPLE.write("\t0")
	OUT_SAMPLE.write("\n")

