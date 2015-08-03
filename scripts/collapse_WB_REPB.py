#!/usr/bin/env python 
# this script writes out WB to Repbase blast matches and takes those with no hit to repbase and blasts them agianst one another
# USE:collapse_WB_REPB.py <blast_comparison_file>
# example:/lscr2/andersenlab/kml436/git_repos2/Transposons/files$ python ../scripts/collapse_WB_REPB.py WB_REPB_blast_family_comparison.txt
import sys
import re
import os
from subprocess import Popen, PIPE

blast_compare_file = sys.argv[1]
BLAST_COMPARE_FILE = open(blast_compare_file, "r")

INT1 = open("intermediate_family_matches_e0.txt", "w")
INT2 = open("intermediate_e0.txt", "w")
INT1.write("WB_ID\tWB_family\tREPB_family\tevalue\n")
INT2.write("WB_ID\tWB_family\tREPB_family\tevalue\n")

INT3 = open("intermediate_non_family_match.txt", "w")
INT3.write("WB_TE\tWB_Family\tBlast_Hit\tPercent_Identity\tAlng_Length\tMismatches\tGap_Open\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_Value\tBit_Score\n")

repb_te_families ={}
wb_te_families ={}
wb_no_match={}
wb_no_match_familes={}
firstline = True
IDs_in_blast={}

for line in BLAST_COMPARE_FILE:
	if firstline:    #skip first line
		firstline = False
		continue
	items= re.split("[\t]",line)
	WB_ID = items[0]
	WB_family = items[1]
	REPB_family = items[2]
	evalue = float(items[11])

	#write out blast hits where best hit Repbase family matched the Wormbase family and the evale was zero
	if WB_family == REPB_family and evalue == 0.0: ##AND IF EVALUSE IS ZERO
		print "YES"
		print line
		repb_te_families[REPB_family] = 0
		wb_te_families[WB_family] = 0
		INT1.write("{WB_ID}\t{WB_family}\t{REPB_family}\t{evalue}\n".format(**locals()))

	#write out blast hits where the families did not match but the evalue was zero (cases where Wormabse family is represented by a Repabse family by a differnet name)
	elif evalue == 0.0:
		print "evalue is {evalue}".format(**locals())
		print line
		repb_te_families[REPB_family] = 0
		wb_te_families[WB_family] = 0
		INT2.write("{WB_ID}\t{WB_family}\t{REPB_family}\t{evalue}\n".format(**locals()))

	#add to dictionaries all hits that did not the above 2 conditions
	else:
		wb_no_match_familes[WB_ID]= WB_family
		wb_no_match[WB_ID]=line

#for each key in dictionary, ignore those that have a family already represented in repbase or whose family matches a wormbase family already represented 
final_no_matches={}
for key in wb_no_match_familes.keys():
	family = wb_no_match_familes[key]
	if family in repb_te_families.keys() or family in wb_te_families.keys():
		print "Family Already Represented"

	#write all other blast hits not matching the above condition to a new file and a new dictionary
	else:
		value = wb_no_match[key]
		INT3.write(value)
		final_no_matches[key] = 0

#write the WB transposon fastas with no mathces to a new file
NO_MATCH = open("no_match.fasta", "w")
from Bio import SeqIO
fasta_sequences = SeqIO.parse(open("WB_all_seqs.fasta"),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	#name, sequence = fasta.id, fasta.seq.tostring()
	if name in final_no_matches.keys():
		NO_MATCH.write(">" + name + "\n" +  sequence + "\n")
NO_MATCH.close()

#make blast databse of these seqeucnes and blast these against one another
#os.system("/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/makeblastdb -in no_match.fasta -dbtype nucl -out no_match_database")
result, err = Popen(['/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/makeblastdb -in no_match.fasta -dbtype nucl -out no_match_database'], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(['/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/blastn -db no_match_database -query no_match.fasta -evalue 1 -outfmt 6 -max_target_seqs 10 -out no_match_blast.txt -num_threads 22'], stdout=PIPE, stderr=PIPE, shell=True).communicate()


	#blast_fields= '\t'.join(items[1:12])
	#print blast_fields
	#IDs_in_blast[transposon_ID] = blast_fields