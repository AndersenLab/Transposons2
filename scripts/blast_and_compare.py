
#this script ........
import sys
import re
import os
from os.path import basename
from subprocess import Popen, PIPE
from  collections import defaultdict

top_blast_hits = sys.argv[1]
TOP_BLAST_HITS = open(top_blast_hits, "r")

###############################################################
###############################################################
# COMPARE WORMBASE FAMILIES TO THE FAMILY SUGGESTED BY BLAST
###############################################################
###############################################################
NEW_WB_FAMILIES = open("new_WB_famlies.txt", "w")

out = "WB_REPB_blast_family_comparison.txt"
OUT = open(out, "w")
OUT.write("WB_TE\tWB_Family\tBlast_Hit\tPercent_Identity\tAlng_Length\tMismatches\tGap_Open\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_Value\tBit_Score\n")

IDs_in_blast={}
WB_family_dct={}
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
	WB_family_dct[WB_ID]=family

	if WB_ID in IDs_in_blast.keys():
		value = IDs_in_blast[WB_ID]
		OUT.write("{WB_ID}\t{family}\t{value}".format(**locals()))
	else:
		OUT.write("{WB_ID}\t{family}\tNo_Hit\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t9999\tNA\n".format(**locals()))
OUT.close()

###############################################################
# CHECK FOR WHICH TRANSPOSONS SHOUDL BE FLIPPED TO REVERSE COMPLEMET LATER
###############################################################
###############################################################
orientation={}
wb_tes_bed = sys.argv[3]
WB_TES_BED= open(wb_tes_bed, "r")
for line in WB_TES_BED:
	items= re.split("[\t]",line)
	WB_TE = items[3]
	WB_Orientation = items[5]
	if WB_Orientation == "-":
		orientation[WB_TE]=WB_Orientation
		print WB_TE
		print "ORIENTATION"
		print WB_Orientation
WB_TES_BED.close()


###############################################################
###############################################################
# SEPARATE TRANSPOSONS BASE ON CORRECT FAM, RENMAE FAM, or NO MATCH
###############################################################
###############################################################

BLAST_COMPARE_FILE = open("WB_REPB_blast_family_comparison.txt", "r")

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
		NEW_WB_FAMILIES.write("{WB_family}\t{REPB_family}\n".format(**locals()))
	#write out blast hits where the families did not match but the evalue was zero (cases where Wormabse family is represented by a Repabse family by a differnet name)
	elif evalue == 0.0:
		print "evalue is {evalue}".format(**locals())
		print line
		repb_te_families[REPB_family] = 0
		wb_te_families[WB_family] = 0
		INT2.write("{WB_ID}\t{WB_family}\t{REPB_family}\t{evalue}\n".format(**locals()))
		NEW_WB_FAMILIES.write("{WB_family}\t{REPB_family}\n".format(**locals()))

	#add to dictionaries all hits that did not match the above 2 conditions
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
		final_no_matches[key] = WB_family_dct[key]


#write the WB transposon fastas with no mathces to a new file
NO_MATCH = open("no_match.fasta", "w")
from Bio import SeqIO
from Bio.Seq import Seq
fasta_sequences = SeqIO.parse(open("WB_all_seqs.fasta"),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	#don't want it formatted as a string 
	seq = Seq(sequence)
	###IF STATEMENT
	rev=seq.reverse_complement()
	#name, sequence = fasta.id, fasta.seq.tostring()
	if name in final_no_matches.keys():
		if name in orientation.keys():
			print name
			print "OUTPUTTTTTTTTTTTTTTTTTT"
			print orientation[name]
			NO_MATCH.write(">" + name + "\n" +  str(rev) + "\n")
		else:
			NO_MATCH.write(">" + name + "\n" +  sequence + "\n")
NO_MATCH.close()
####HERE SWITCH
###############################################################
###############################################################
# BLAST NO MATCHES AGAINST EACH OTHER
###############################################################
###############################################################

#make blast database of these seqeucnes and blast these against one another
result, err = Popen(['/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/makeblastdb -in no_match.fasta -dbtype nucl -out no_match_database'], stdout=PIPE, stderr=PIPE, shell=True).communicate()
result, err = Popen(['/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/blastn -db no_match_database -query no_match.fasta -evalue 1 -outfmt 6 -max_target_seqs 10 -out no_match_blast.txt -num_threads 22'], stdout=PIPE, stderr=PIPE, shell=True).communicate()

IN_NO_MATCH_BLAST = open("no_match_blast.txt", "r")
OUT_NO_MATCH_BLAST = open("no_match_blast_comparison.txt", "w")
WB_IDs_found={}
for line in IN_NO_MATCH_BLAST:
	items= re.split("[\t]",line)
	WB_ID1 = items[0]
	family_ID1 = WB_family_dct[WB_ID1]
	WB_ID2 = items[1]
	family_ID2 = WB_family_dct[WB_ID2]
	blast_info='\t'.join(items[2:12])
	WB_IDs_found[WB_ID1] = 0
	OUT_NO_MATCH_BLAST.write("{WB_ID1}\t{family_ID1}\t{WB_ID2}\t{family_ID2}\t{blast_info}".format(**locals()))
for no_hit in final_no_matches.keys():
	no_hit_fam = final_no_matches[no_hit]
	if no_hit not in WB_IDs_found.keys():
		OUT_NO_MATCH_BLAST.write("{no_hit}\t{no_hit_fam}\tNo_Hit\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t9999\tNA\n".format(**locals()))
OUT_NO_MATCH_BLAST.close()

result, err = Popen(['cat no_match_blast_comparison.txt |sort -k1,1 > tmp && mv tmp no_match_blast_comparison.txt'], stdout=PIPE, stderr=PIPE, shell=True).communicate()

seen={}
consensus_families={}
#consensus_families= defaultdict(list)
final={}
BLAST = open("no_match_blast_comparison.txt", "r")
for line in BLAST:
	items= re.split("[\t]",line)
	WB_ID1 = items[0]
	WB_ID1_family = items[1]
	WB_ID2 = items[2]
	WB_ID2_family = items[3]
	evalue = float(items[12])
	if evalue==0.0 and re.search("WBTransposon",WB_ID1_family):
		if not re.search("WBTransposon", WB_ID2_family):
			#instance of renaming
			print 'RENAME'
			final[WB_ID1] = WB_ID2_family
			seen[WB_ID1_family]=0
			NEW_WB_FAMILIES.write("{WB_ID1_family}\t{WB_ID2_family}\n".format(**locals()))
	if evalue==0.0 and re.search("WBTransposon",WB_ID1_family):
		if re.search("WBTransposon", WB_ID2_family):
			if WB_ID1_family not in seen.keys():
				seen[WB_ID1_family]=0
				seen[WB_ID2_family]=0
				print "FOLLOWING NOT SEEN:"
				print WB_ID1_family
				final[WB_ID2_family]= WB_ID1_family
				#consensus_families[WB_ID1_family].append(WB_ID2_family)
			elif WB_ID1_family in seen.keys() and WB_ID2_family not in seen.keys():
				seen[WB_ID2_family]=0
				final[WB_ID2_family]= WB_ID1_family
				#consensus_families[WB_ID1_family].append(WB_ID2_family)
			elif WB_ID1_family in seen.keys() and WB_ID2_family  in seen.keys(): ##what about the fourth condition?
				print "already seen"
	elif evalue==9999 and re.search("NA", WB_ID2_family):
		print "no hit to self"

	else:
		if not re.search("WBTransposon", WB_ID1_family):
			final[WB_ID1] = WB_ID1_family

for key in final.keys():
	value=final[key]
	VALUE = open("{value}_fastas.fasta".format(**locals()), "w")
	VALUE.close()

###############################################################
###############################################################
# OUTPUT FASTA SEQUENCES TO THEIR CORRECT FILES
###############################################################
###############################################################

fasta_sequences = SeqIO.parse(open("no_match.fasta"),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	#name, sequence = fasta.id, fasta.seq.tostring()
	if name in final.keys():
		value=final[name]
		NEW_WB_FAMILIES.write("{name}\t{value}\n".format(**locals()))
		VALUE = open("{value}_fastas.fasta".format(**locals()), "a")
		VALUE.write(">" + name + "\n" +  sequence + "\n")
		VALUE.close()
NEW_WB_FAMILIES.close()

result, err = Popen(['cat new_WB_famlies.txt |sort -k1,1| sort -u -k1,1> tmp && mv tmp new_WB_famlies.txt'], stdout=PIPE, stderr=PIPE, shell=True).communicate()
IN_NEW = open("new_WB_famlies.txt", "r")
family_renames={}
for line in IN_NEW:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	old_name = items[0]
	new_name = items[1]
	family_renames[old_name] = new_name
IN_NEW.close()

round2_names_for_gff={}
ORIGINAL_WB_FAMILIES = open("WB_familes.txt", "r")
INPUT_WB_FAMILIES = open("round2_WB_familes.txt", "w")
for line in ORIGINAL_WB_FAMILIES:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	wb_name = items[0]
	wb_family = items[1]
	if wb_family in family_renames.keys():
		value = family_renames[wb_family]
		INPUT_WB_FAMILIES.write("{wb_name}\t{value}\n".format(**locals()))
		round2_names_for_gff[wb_name] = value
	else:
		INPUT_WB_FAMILIES.write("{wb_name}\t{wb_family}\n".format(**locals()))
		round2_names_for_gff[wb_name] = wb_family

INPUT_WB_FAMILIES.close()

#ALGIN AND CALL CONSENSUS

###############################################################
###############################################################
# CREATE NEW GFF FILE
###############################################################
###############################################################

gff_file="/lscr2/andersenlab/kml436/sv_files/TEMP_TELOCATE_RETROSEQ_sim_files/WB_all_dups_renamed.gff"
GFF= open(gff_file, "r")
R2_GFF = open("round2_WB_tes.gff", "w")
for line in GFF:
	line = line.rstrip('\n')
	for key in round2_names_for_gff.keys(): #replace with more efficient method later
		value = round2_names_for_gff[key]
		if re.search(key, line):
			line = line.replace(key, value)
	R2_GFF.write("{line}\n".format(**locals()))

###############################################################
###############################################################
# ALIGN FASTA FILES
###############################################################
###############################################################
#CLEAN THIS SECTION UP LATER
my_directory=os.getcwd()
print my_directory
os.system("mkdir alignments")
os.system("mkdir final_consensus_fastas")
#my_directory = sys.argv[4]
for results_file in os.listdir(my_directory):
	if re.search("fastas.fasta", results_file):
		print results_file
		result, err = Popen(["/opt/mafft-6.240/bin/mafft {my_directory}/{results_file}> {results_file}.aln".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
		grep_result, err = Popen(['grep -c ">" %s' %results_file], stdout=PIPE, stderr=PIPE, shell=True).communicate()
		print grep_result
		#grep_result = grep_result.rstrip('\n')
		#move files with only one fasta sequence to the final consensus folder
		if int(grep_result)==1:
			print "ONE"
			os.system("mv {results_file} final_consensus_fastas".format(**locals()))
		#otherwise move the alignment and fasta to the alignment folder for further processing
		else:
			os.system("mv {results_file}.aln alignments".format(**locals()))
			os.system("mv {results_file} alignments".format(**locals()))
		#match = re.search("(\d+)\s",result)

#"/opt/mafft-6.240/bin/mafft {dir}/{results_file}> {results_file}.aln".format(**locals())
