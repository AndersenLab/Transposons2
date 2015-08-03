 #!/usr/bin/env python
 # this script calculates the TPR and FDR for each family 
 # USE: family_TPR_FDR.py <consensus file used for simulations> <position files of known spike ins> <bed file of correct found families> <bed file of incorrectly found families>

import sys
import re
import os
from subprocess import Popen, PIPE
from Bio import SeqIO

consensus = sys.argv[1]
known_positions = sys.argv[2]
correct_file = sys.argv[3]
incorrect_file = sys.argv[4]
out_name = sys.argv[5]

CORRECT_FILE = open(correct_file,"r")
INCORRECT_FILE = open(incorrect_file,"r")
OUT_NAME = open("FAMILY_TFPN_{out_name}".format(**locals()), "w")


###NEED TO DO FOR EACH METHOD

#set up dictionaries
all_families={}
known_positions_count={} ####NEED TO SET ALL TO ZERO
correct_families={}
incorrect_families={}

#add all families in the consensus file to a dictionary
fasta_sequences = SeqIO.parse(open(consensus),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	all_families[name] = 0

for key in all_families.keys():
	known_positions_count[key] = 0
	correct_families[key] = 0
	incorrect_families[key] = 0

#count how many times each familiy was simulated and should have been found
result, err = Popen(["cat {known_positions} | cut -f4".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result=result.split()
for i in result:
	if i in all_families.keys():
		#print i
		known_positions_count[i] +=1

#count how many of each family was found correctly(True Positives)
for line in CORRECT_FILE:
	items = re.split("[\t]",line)
	transposon = items[3]
	match = re.search("(.*)_(\w+-)?reference", transposon) #match = re.search("(.*)_non-reference", transposon)
	family = match.group(1)
	#print family
	correct_families[family] +=1
CORRECT_FILE.close

#count how many of each family was found incorrectly(False Positives)
for line in INCORRECT_FILE:
	items = re.split("[\t]",line)
	transposon = items[4] # output is shifted one from comm in this file
	match = re.search("(.*)_(\w+-)?reference", transposon) #match = re.search("(.*)_non-reference", transposon)
	family = match.group(1)
	#print family
	incorrect_families[family] +=1
INCORRECT_FILE.close()

#output TPR and FDR for each family
for family in all_families.keys():
	#calculate TPR
	TP = int(correct_families[family])
	FP = int(incorrect_families[family])
	simulated = int(known_positions_count[family])

	if simulated==0:
		TPR = "NA" # no true positives if nothing simulated 
	else:
		TPR = float(TP)/(simulated) * 100

	if FP is not 0:
		if TP is 0:
			positives = int(FP)
			FDR = float(FP)/(positives) * 100
		else:
			positives= int(TP+FP)
			FDR = float(FP)/(positives) * 100
	else:
		FDR = 0

	print "FAMILY:{family}\tSIMUALTED:{simulated}\tTP:{TP}\tFP:{FP}\tTPR:{TPR}\tFDR:{FDR}\n".format(**locals())

	OUT_NAME.write("{out_name}\t{family}\t{simulated}\t{TP}\t{FP}\t{TPR}\t{FDR}\n".format(**locals()))



	

