#!/usr/bin/env python 
# this script takes the fasta seqs and checks that each seq name is represetned in the element->family conversion file, position gff file, length files, fasta list file, and bed list file
# if not, an error occured while generating one of the above listed files and the missing seq name will be printed
# USE:CORRECTIONS_WBT_tel.py
import re
from Bio import SeqIO
names={}
consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/round2_consensus_fasta.fasta"
fasta_sequences = SeqIO.parse(open(consensus),'fasta')
for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
	names[name] = 0

wbts="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_familes.txt"
print wbts
WBTS= open(wbts, "r")
for line in WBTS:
	line = line.rstrip('\n')
	items= re.split("[\t]",line)
	old_name = items[0]
	new_name = items[1]
	if new_name not in names.keys():
		print new_name
WBTS.close()

wbts="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/round2_WB_tes.gff"
print wbts
#wbts="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/new_WB_famlies.txt"
WBTS= open(wbts, "r")
for line in WBTS:
        line = line.rstrip('\n')
        items= re.split("[\t]",line)
        new_name = items[2]
        if new_name not in names.keys():
                print new_name
WBTS.close()

wbts="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/lengths.txt"
print wbts
WBTS= open(wbts, "r")
for line in WBTS:
        line = line.rstrip('\n')
        items= re.split("[\t]",line)
        new_name = items[0]
        if new_name not in names.keys():
                print new_name
WBTS.close()

wbts="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/fasta_list.txt"
print wbts
WBTS= open(wbts, "r")
for line in WBTS:
        line = line.rstrip('\n')
        items= re.split("[\t]",line)
        new_name = items[0]
        if new_name not in names.keys():
                print new_name
WBTS.close()



wbts="/lscr2/andersenlab/kml436/git_repos2/Transposons/files/CORRECTIONS/bed_list.txt"
print wbts
WBTS= open(wbts, "r")
for line in WBTS:
        line = line.rstrip('\n')
        items= re.split("[\t]",line)
        new_name = items[0]
        if new_name not in names.keys():
                print new_name
WBTS.close()

