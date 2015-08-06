#!/usr/bin/env python
# this scripts adds  "dnatransposon, retrotransposon, or unknown" information to the transposon calls
# USE: assign_cut_or_copy_paste.py <repbase_fasta> <consensus_fasta> <all_nonredundant file>

import sys
import re
import os

repbase_fasta=sys.argv[1] #"/lscr2/andersenlab/kml436/repbase.fasta"
consensus_fasta=sys.argv[2] #"/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta" 
all_nonredundant=sys.argv[3] #"/lscr2/andersenlab/kml436/git_repos2/Transposons2/data/all_nonredundant.txt"

repbase_family={}

REPBASE=open(repbase_fasta, "r")
for line in REPBASE:
	line=line.rstrip('\n')
	if re.search('^>',line):
		items=re.split("[\t]",line)
		TE=items[0]
		TE=TE.replace(">","")
		repbase_family[TE]=0

WB={}
other={}
CONSENSUS=open(consensus_fasta, "r")
for line in CONSENSUS:
	line=line.rstrip('\n')
	if re.search('^>',line):
		items=re.split("[\t]",line)
		consensus_TE=items[0]
		consensus_TE=consensus_TE.replace(">","")
		if consensus_TE not in repbase_family.keys():
			if re.search("WBTransposon",consensus_TE):
				WB[consensus_TE]=0
			else:
				other[consensus_TE]=0
#for transposon in other.keys():
#	print transposon
REPBASE.close()
CONSENSUS.close()


#from repbase
rep_CC={'Chapaev':'dnatransposon',
'CR1':'retrotransposon',
'DNA transposon':'dnatransposon',
'EnSpm/CACTA':'dnatransposon',
'Gypsy':'retrotransposon',
'Harbinger':'dnatransposon',
'hAT':'dnatransposon',
'Helitron':'dnatransposon',
'LTR Retrotransposon':'retrotransposon',
'Mariner/Tc1':'dnatransposon',
'Merlin':'dnatransposon',
'MuDR':'dnatransposon',
'NeSL':'retrotransposon',
'RTE':'retrotransposon',
'SINE':'retrotransposon',
'SINE2/tRNA':'retrotransposon',
'Transposable Element':'dnatransposon', #only one classified as this is TC6, which is a dna transposon
'Vingi':'retrotransposon'}



#classification of wormbase families (also represent repbase renames to WB)
classes={'CER6':'retrotransposon',
'MARINER3':'dnatransposon',
'RTE1':'retrotransposon',
'Tc9':'dnatransposon',
'Tc6':'dnatransposon',
'Tc5':'dnatransposon',
'Tc2':'dnatransposon',
'Tc1':'dnatransposon',
'Tc5A':'dnatransposon',
'Tc5B':'dnatransposon',
'Tc1A':'dnatransposon',
'TC8':'dnatransposon',
'Tc4v':'dnatransposon',
'MARINER4':'dnatransposon',
'MARINER5':'dnatransposon',
'MARINER2':'dnatransposon',
'NTc2A':'dnatransposon',
'LINE2F':'retrotransposon',
'LINE2D':'retrotransposon',
'LINE2E':'retrotransposon',
'LINE2B':'retrotransposon',
'LINE2C':'retrotransposon',
'LINE2A':'retrotransposon'
}



#classes={}
#partly repetitive of above....below added in after manual checks
REPBASE=open(repbase_fasta, "r")
for line in REPBASE:
	line=line.rstrip('\n')
	if re.search('^>',line):
		items=re.split("[\t]",line)
		TE=items[0]
		class_TE=items[1]
		TE=TE.replace(">","")
		repbase_family[TE]=0
		classes[TE]=rep_CC[class_TE] # assign transposon name a transposon class


# add WB familt cut/copy info from WB elements that were originally assigned a WB family but were renamed to an element name
classes.update({
		'WBTransposon00000046':'dnatransposon', #Tc10 dnatransposon
		'WBTransposon00000047':'dnatransposon', #Tc10 dnatransposon
		'WBTransposon00000074':'dnatransposon', #Tc3-CeIIa dnatransposon
		'WBTransposon00000084':'dnatransposon', #Tc3-CeIIb dnatransposon
		'WBTransposon00000416':'dnatransposon', #CELETC2 dnatransposon
		'WBTransposon00000425':'retrotransposon', #CER6 retrotransposon
		'WBTransposon00000428':'retrotransposon', #CER6 retrotransposon
		'WBTransposon00000586':'retrotransposon', #CER1 retrotransposon
		'WBTransposon00000588':'retrotransposon', #LINE2E retrotransposon
		'WBTransposon00000637':'dnatransposon' #Polinton-1_CB dnatransposon
		})

#ADD IN COPT/CUT INFO:
new_out=os.path.basename(all_nonredundant)
ALL_NONREDUNDANT=open(all_nonredundant, "r")
CTCP_ALL_NONREDUNDANT=open("CtCp_{new_out}".format(**locals()), "w")
unknowns={}
for line in ALL_NONREDUNDANT:
	line=line.rstrip('\n')
	items=re.split("[\t]", line)
	te_info=items[3]
	match=re.search("(.*)_(\w+-)?reference", te_info)
	family=match.group(1)
	if family in classes.keys():
		CtCp=classes[family]
	else:
		CtCp="unknown"
		unknowns[family]=0
	CTCP_ALL_NONREDUNDANT.write("{line}\t{CtCp}\n".format(**locals()))

# check that only WB elements renamed to antoher WB element are included in the "unkown" category
# this doesnt mean it doesnt have a family...but those that do were manually added above
CHECK_FILE=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB3.txt", "r")
for line in CHECK_FILE:
	line=line.rstrip('\n')
	items=re.split("[\t]", line)
	WBelement=items[0]
	#if WBelement in unknowns.keys():
	#	print line

REPBASE.close()
ALL_NONREDUNDANT.close()
CHECK_FILE.close()


##checks have passed
