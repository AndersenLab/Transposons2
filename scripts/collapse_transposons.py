#!/usr/bin/env python 
# this script collapses adjacent transposons if they belong to the same family and lie within the length of the corresponding TE family
# USE:collapse_transposons.py <bed_files of te positions>

import sys
import re
import os

in_file= sys.argv[2]

def rreplace(s, old, new, occurrence):
	li = s.rsplit(old, occurrence)
	return new.join(li)


TE_length_file = sys.argv[1]
TE_length_file = open(TE_length_file, "r")

transposon_lengths={}
for line in TE_length_file:
	items= re.split("[\t]",line)
	transposon_ID= items[0]
	length= items[1]
	if transposon_ID in transposon_lengths.keys():
		if length > transposon_lengths[transposon_ID]:
			transposon_lengths[transposon_ID] = length
	else:
		transposon_lengths[transposon_ID]=length


def collapse_transposons(bed_file):

	inp_bed_file = open(bed_file, "r")
	out_file_CT = "CT_{bed_file}".format(**locals())
	my_out_file_CT = open(out_file_CT, "w")


	collapsed_transposons={}
	first_line = True
	for line in inp_bed_file:
		items= re.split("[\s+]",line)
		chromosome = items[0]
		start_pos = items[1]
		end_pos = items[2]
		full_trans_name = items[3]
		read_support = items [4]

		#match = re.search("(.*)_(reference)|(.*)_(non-reference)",full_trans_name)
		match = re.search("(.*)_(\w+-)?reference",full_trans_name) #####fix this part
		if match is not None:
			#print line
			family_found = match.group(1)
			#print family_found
		if first_line == False:
			if chromosome == prev_chromosome and family_found == prev_family_found:
				difference = (int(end_pos)-int(prev_start_pos))
				difference2 = (int(difference) - int(transposon_lengths[family_found]))

				if int(difference2) < 0:
					print "COLLAPSE"
					print difference
					print difference2
					print family_found
					print transposon_lengths[family_found]

					print prevLine
					print line
					print prevLine
					print line
					print transposon_lengths[family_found]

					if float(read_support)> float(prev_read_support):
						prevLine = line


					#prevLine = rreplace(prevLine, "\t{prev_end_pos}".format(**locals()), "\t{end_pos}".format(**locals()), 1)  # replace last occurence...need to avois number after
						collapsed_transposons[prev_full_trans_name] = prevLine # DO NOT RESET prev_full_trans_name.....this is used to keep tract of the family cluster
						prev_end_pos=end_pos
						prev_start_pos=start_pos
						prev_read_support=read_support
						print prevLine
					else:
						print "previous read support was higher..."

				else:
					collapsed_transposons[full_trans_name] = line
					collapsed_transposons[full_trans_name] = line
					prevLine = line
					prev_chromosome = chromosome
					prev_start_pos = start_pos
					prev_end_pos = end_pos
					prev_full_trans_name = full_trans_name
					prev_family_found = family_found
					prev_read_support = read_support
			else:
				collapsed_transposons[full_trans_name] = line
				prevLine = line
				prev_chromosome = chromosome
				prev_start_pos = start_pos
				prev_end_pos = end_pos
				prev_full_trans_name = full_trans_name
				prev_family_found = family_found
				prev_read_support = read_support

		else:
			collapsed_transposons[full_trans_name] = line
			prevLine = line
			prev_chromosome = chromosome
			prev_start_pos = start_pos
			prev_end_pos = end_pos
			prev_full_trans_name = full_trans_name
			prev_family_found = family_found
			prev_read_support = read_support
		first_line = False
	for key in collapsed_transposons.keys():
		my_out_file_CT.write(collapsed_transposons[key])

	#print results_file
	my_out_file_CT.close()


	os.system("cat {out_file_CT}| sort -k1,1 -k2,2n > new_{out_file_CT}_exact_positions".format(**locals()))
	os.system("""cat new_{out_file_CT}_exact_positions| awk '{{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6}}'> new_{out_file_CT}""".format(**locals()))

collapse_transposons(in_file)