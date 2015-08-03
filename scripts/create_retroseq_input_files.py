#!/usr/bin/env python 
# this script takes a fastsa and gff file as input and generates the files that are necessary as inputs for retroseg
# such as an individual file for each fasta seq and a fasta list containing the paths to such files
# USE:create_retroseq_input_files.py <fasta_files> <gff_file> 
# example: python /exports/people/andersenlab/kml436/scripts/create_retroseq_input_files.py /lscr2/andersenlab/kml436/sv_files/TEMP_TELOCATE_RETROSEQ_sim_files/WB_all_seqs_no_dups.fasta /lscr2/andersenlab/kml436/sv_files/TEMP_TELOCATE_RETROSEQ_sim_files/WB_all_dups_renamed.gff

from Bio import SeqIO
import sys
import re
import os
all_seqs = sys.argv[1] 
pos = sys.argv[2]

os.system("mkdir fastas")
os.system("mkdir beds")
fasta_list = "fasta_list.txt"
FASTA_LIST = open(fasta_list, "w")
bed_list = "bed_list.txt"
BED_LIST = open(bed_list, "w")

dirname = os.path.dirname(os.path.realpath("fastas"))

#PRINT THE FASTAS TO THEIR OWN FILES
fasta_sequences = SeqIO.parse(open(all_seqs),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, fasta.seq.tostring()
	print name
	#uncomment below line and tab others over if only want to pull out WormBase transposons
	#if re.search("WBTransposon",name):  #
	print "YES"
	out = "c_elegans_reference_{name}.fasta".format(**locals())
	OUT = open(out, "w")
	OUT.write(">" + name + "\n" +  sequence + "\n")
	OUT.close()
	os.system("mv {out} fastas".format(**locals()))
	FASTA_LIST.write("{name}\t{dirname}/fastas/{out}\n".format(**locals()))

#PRINT THE POSITIONS TO THEIR OWN FILES
POS = open(pos, "r")
for line in POS:
	line = line.rstrip('\n')
	items = re.split("[\t]",line)
	chromosome = items[0]
	WB_name = items[2]
	start = items[3]
	end = items[4]
	orientation = items[6]
	out_pos = "c_elegans_reference_{WB_name}.bed".format(**locals())
	OUT_POS = open(out_pos, "w")
	OUT_POS.write("{chromosome}\t{start}\t{end}\tinstance{WB_name}\t.\t{orientation}\n".format(**locals()))
	OUT_POS.close()
	os.system("mv {out_pos} beds".format(**locals()))
	BED_LIST.write("{WB_name}\t{dirname}/beds/{out_pos}\n".format(**locals()))

print dirname

