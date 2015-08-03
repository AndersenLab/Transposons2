#!/usr/bin/env python 
# this script randomly samples the TE families and outputs 100 of them to be used in the RSVSIm simulations
# USE: create_RSVSim_genome.py <consensus_family_file> <run_number>
# called on by run_RSVSim.sh

import sys
import re
import os
import random
from Bio import SeqIO
from subprocess import Popen, PIPE
##loop this 10x over later
in_file = sys.argv[1]
run_no = sys.argv[2]
#read in familie names
result, err = Popen(["""cat %s| awk '$0 ~/^>/ {print $0}'""" %(in_file)], stdout=PIPE, stderr=PIPE, shell=True).communicate()
family_count={}

fasta_out = "RSVSIM_ins_seqs_{run_no}.txt".format(**locals())
FASTA_OUT = open(fasta_out, "w")
te_pos= "RSVSIM_te_pos_{run_no}.txt".format(**locals())
TE_POS= open(te_pos, "w")

result = result.split('\n')
#get rid of ">" in family names
result = [i.replace('>', '') for i in result]
for i in result: ######GET RID OF THE BALNK ONE!!!!!!
	print i
	family_count[i] = 0
    #if i != '': #get rid of balnk key
       # unique_distances[i]=0
for i in range(0,100): 
	random_family = (random.choice(result))

	fasta_sequences = SeqIO.parse(open(in_file),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq) 
		length = len(sequence)
		if name == random_family:
			family_count[name] += 1
			copy_number = family_count[name] #needs to have unique names, so if the same transposon is simulated more than once, add a copy number to the name
			FASTA_OUT.write(">{name}_copy{copy_number}\n{sequence}\n".format(**locals()))
			TE_POS.write("{name}_copy{copy_number}\t1\t{length}\t{name}_copy{copy_number}\n".format(**locals()))
FASTA_OUT.close()
TE_POS.close()

reference = "/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa"
os.system("cat {fasta_out} {reference} > RSVSIM_genome_{run_no}.fasta".format(**locals()))