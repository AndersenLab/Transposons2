#!/usr/bin/env python 
# this script calculates the mean and standard deviation of the TPR and FDR rate per transposon family
# USE:family_TFPN_average4.py <directory> <consensus_fasta_file>
# NOTE: comment/uncomment line 90-92 depending on the method that needs to analzed

import sys
import re
import os
from subprocess import Popen, PIPE
from  collections import defaultdict
from Bio import SeqIO
import statistics

out ="BEDCOMPARE_FAMILY_MEANS.txt"
OUT = open(out, "w")
OUT.write ("M1\tfam\tTPR\tFDR\tTPR_error\tFDR_error\n")

all_families={}
consensus = sys.argv[2]
fasta_sequences = SeqIO.parse(open(consensus),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    all_families[name] = 0

###START FUNCTION HERE
def calculate_mean(TE_program):
    my_directory = sys.argv[1] 
    files=''
    items = []
    for results_file in os.listdir(my_directory):
        match = re.findall("(FAMILY_TFPN_ALL*)",results_file)
        if len(match) >0:
            print "yes"
            print results_file
            files+=str(" {results_file}".format(**locals()))
    TPR_fam={}
    FDR_fam={}
    TPR_fam = defaultdict(list)
    FDR_fam = defaultdict(list)

    fam_found = {}

    #get rid of leading space...next tiem append space after file:
    files= files[1:]
    files_to_test = files.split(' ')
    for sim_file in files_to_test:
        OPEN_SIM_FILE = open(sim_file, "r")
        for line in OPEN_SIM_FILE:
            if re.search(TE_program,line):
                line = line.rstrip('\n')
                items= re.split("[\t]",line)   # WILL NEED TO CHANGE THESE
                M1 = items[0]
                fam = items[1]
                TPR = items[5]
                FDR = items[6]

                if TPR !="NA":
                    TPR_fam[fam].append(TPR)
                    fam_found[fam] = 0
                FDR_fam[fam].append(FDR)
    for key in FDR_fam.keys():
        print key
        print FDR_fam[key]

    for key in sorted(all_families.keys()):
        if key in fam_found.keys():
            TPR_fam[key] = map(float, TPR_fam[key]) #convert strings in list to integers
            mean_TPR = statistics.mean(TPR_fam[key])
            standard_deviation_TPR = statistics.pstdev(TPR_fam[key])
        else:
            mean_TPR = "NA"
            standard_deviation_TPR = "NA"

        FDR_fam[key] = map(float, FDR_fam[key]) 
        print key
        print FDR_fam[key]
        mean_FDR = statistics.mean(FDR_fam[key])
        standard_deviation_FDR = statistics.pstdev(FDR_fam[key])
        print "The mean_TPR is {mean_TPR}".format(**locals())
        print "The standard deviation TPR is {standard_deviation_TPR}".format(**locals())
        OUT.write ("{M1}\t{key}\t{mean_TPR}\t{mean_FDR}\t{standard_deviation_TPR}\t{standard_deviation_FDR}\n".format(**locals()))
        ##OUT.write ("{M1}_error\t{key}\t{standard_deviation_TPR}\t{standard_deviation_FDR}\n".format(**locals()))


#calculate_mean # ADD IN ALL OPTION!!!!!!!!!!!!!!!!!!
#!!!!!!!!
#!!!!!!!!
#XXXXXXXX
#AND UNCOMMENT BELOW
calculate_mean("temp") #last setting
calculate_mean("retroseq")
calculate_mean("telocate") #last setting


OUT.close()
