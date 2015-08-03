#backup copy of script original version of script in bamsurgeon file
import sys
import re
import os
from subprocess import Popen, PIPE

#set up files
my_distance = sys.argv[1]
my_directory = sys.argv[2]
my_name = sys.argv[5]
my_distance =int(my_distance)
my_out_file = "BEDCOMPARE_{my_name}.txt".format(**locals())
my_out_file_family = "BEDCOMPARE_FAMILY_{my_name}.txt".format(**locals())
my_out_file = open(my_out_file, "w")
my_out_file_family = open(my_out_file_family, "w")
os.system("mkdir intermediates")

##TAKE OUT LOG TEST LATER
log_test = "LOG_test.txt"
log_test= open(log_test, "w")

#write headers to output files
my_out_file.write( "M1\tM2\tConcordance\tTotal_M1\tFraction_in_M1_Found\tFraction_CorrectFam_in_M1_Found\tFraction_in_M1_NOT_Found\tTotal_M2\tFraction_in_M2_Found\tFraction_in_M2_NOT_Found\n".format(**locals()) )
my_out_file_family.write( "M1\tM2\tFamily\tTotal_Known\tNumber_Known_Found\tPercent_Known_Found\tNumber_Found_Correct\tPercent_Found_Correct\n".format(**locals()) )

my_out_file_summary = "BEDCOMPARE_SUMMARY_{my_name}.txt".format(**locals())
my_out_file_summary = open(my_out_file_summary, "w")
#write headers to output files
my_out_file_summary.write("M1\tM2\tDistanceCutoff\tfam\tTPR\tFDR\n")
#####################################################################################
def rreplace(s, old, new, occurrence):
	li = s.rsplit(old, occurrence)
	return new.join(li)


TE_length_file = sys.argv[3]
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

####################################################################################

def compare_beds(file1,file2):
	#find the closest festure in file2 to each feature in file 1 with bedtools
	os.system("closestBed -a {file1} -b {file2} -d -t first > {file2}_to_{file1}_compare.bed".format(**locals()))
	#set dictionaries and variables
	family_total={}
	family_count={}
	family_matched_count={}
	count = 0
	total_transposons =0
	total_family_matched_count=0
	#get total counts of the features found in file1
	result, err = Popen(['wc -l %s' %file1], stdout=PIPE, stderr=PIPE, shell=True).communicate()
	match = re.search("(\d+)\s",result)
	total_transposons = int(match.group(1))
	print total_transposons
	result, err = Popen(['wc -l %s' %file2], stdout=PIPE, stderr=PIPE, shell=True).communicate()
	match = re.search("(\d+)\s",result)
	total_transposons2 = int(match.group(1))
	unique_transposons = {}
	final_transposons = {}
	#####account for ties
	#open the ouput file from bed closest tool
	my_in_file = open("{file2}_to_{file1}_compare.bed".format(**locals()), "r")
	for line in my_in_file:
		items= re.split("[\t]",line)
		ID1 = items[3]
		ID2 = items[9]
		distance = items[12]
		distance = distance.rstrip('\n')
		distance =int(distance)

		# count a feature in file1 as matched if the distance of the closeset feature in file2 is greater than zero and less than the max distance specified in input
		if distance >=0:
			if ID2 in unique_transposons.keys():
				difference = (int(distance) - int(unique_transposons[ID2]))
				if int(difference) < 0:
					unique_transposons[ID2] = distance
					final_transposons[ID2] = line

			else:
				unique_transposons[ID2] = distance
				final_transposons[ID2] = line


	for key in unique_transposons.keys():
		log_test.write(str(key) + "\t" + str(unique_transposons[key]) + "\n")

	SC_file = "SC_{file2}_to_{file1}_compare.bed".format(**locals())
	SC_file = open(SC_file, "w")
	for key in final_transposons.keys():
		line = final_transposons[key]
		SC_file.write(line)
		items= re.split("[\t]",line)
		ID1 = items[3]
		ID2 = items[9]
		distance = items[12]
		distance = distance.rstrip('\n')
		distance =int(distance)
		#create a key for each transposons family if it does not already exists
		if ID1 in family_total.keys():
			family_total[ID1] += 1
		else:
			family_total[ID1] = 1
			family_count[ID1] = 0
			family_matched_count[ID1] = 0
		# count a feature in file1 as matched if the distance of the closeset feature in file2 is greater than zero and less than the max distance specified in input
		if distance <= my_distance and distance >=0:
			count += 1

			#print "{ID1}\t{ID2}".format(**locals())
			#increase the matched count if the closest feature is within X base pairs distance
			family_count[ID1] += 1
			if re.search("(temp)|(telocate)|(retroseq)",file2):
				print "asfjkhdskjfhjkdskfllllllllllllllllllllllllllllllllllllllllllllllllhk"
				match = re.search("(.*)_(\w+-)?reference",ID2)
				family_found = match.group(1)
			else:
				family_found= ID2
			#print family_found
			#if the family of the closest feature found matches what it should be, increase the matched correctly count
			if family_found == ID1:
				print "FOUND!"
				family_matched_count[ID1] += 1
				total_family_matched_count += 1
				#print family_matched_count[ID1]
			else:
				print "WRONG ONE!"
			print "{ID1}\t{ID2}\t{family_found}".format(**locals())








	not_found = total_transposons-count
	fraction_not_found = float (not_found)/(total_transposons) *100
	fraction = float (count)/(total_transposons) *100
	fraction_fam_correct = float (total_family_matched_count)/(total_transposons) *100

	not_found2 = total_transposons2-count
	fraction_not_found2= float (not_found2)/(total_transposons2) *100
	fraction2 = float (count)/(total_transposons2) *100

	my_out_file.write( "{file1}\t{file2}\t{count}\t{total_transposons}\t{fraction}%\t{fraction_fam_correct}%\t{fraction_not_found}%\t{total_transposons2}\t{fraction2}%\t{fraction_not_found2}%\n".format(**locals()) )
	for key in family_count.keys():
		#total feature for each family in file1
		total = family_total[key]
		#total number of above feature found in file2
		num_found = family_count[key]
		#total number of correctly matched features out of those that were matched
		num_correct_from_found = family_matched_count[key]
		#percent matched(regardess if correct) of feature in file1 families
		percent_family_found = float(family_count[key])/(family_total[key]) *100
		#print the results, if a family wasn't found, print "none" in the relevant column
		if family_count[key] >0:
			#percent correctly matched of matched feature in file1 families
			percent_family_found_from_matched = float(num_correct_from_found)/(num_found) *100
			my_out_file_family.write("{file1}\t{file2}\t{key}\t{total}\t{num_found}\t{percent_family_found}\t{num_correct_from_found}\t{percent_family_found_from_matched}\n".format(**locals())) # add overall percent family found
		else:
			my_out_file_family.write( "{file1}\t{file2}\t{key}\t{total}\t{num_found}\t{percent_family_found}\t{num_correct_from_found}\tnone\n".format(**locals()))
	#my_in_file.close()


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
					prevLine = rreplace(prevLine, "\t{prev_end_pos}".format(**locals()), "\t{end_pos}".format(**locals()), 1)  # replace last occurence...need to avois number after
					collapsed_transposons[prev_full_trans_name] = prevLine
					prev_end_pos=end_pos
					print prevLine
				else:
					collapsed_transposons[full_trans_name] = line
					collapsed_transposons[full_trans_name] = line
					prevLine = line
					prev_chromosome = chromosome
					prev_start_pos = start_pos
					prev_end_pos = end_pos
					prev_full_trans_name = full_trans_name
					prev_family_found = family_found
			else:
				collapsed_transposons[full_trans_name] = line
				prevLine = line
				prev_chromosome = chromosome
				prev_start_pos = start_pos
				prev_end_pos = end_pos
				prev_full_trans_name = full_trans_name
				prev_family_found = family_found

		else:
			collapsed_transposons[full_trans_name] = line
			prevLine = line
			prev_chromosome = chromosome
			prev_start_pos = start_pos
			prev_end_pos = end_pos
			prev_full_trans_name = full_trans_name
			prev_family_found = family_found
		first_line = False
	for key in collapsed_transposons.keys():
		my_out_file_CT.write(collapsed_transposons[key])

	#print results_file
	my_out_file_CT.close()


	os.system("cat {out_file_CT}| sort -k1,1 -k2,2n > new_{out_file_CT}".format(**locals()))

def TFPN(file1,file2):
		#find the closest festure in file2 to each feature in file 1 with bedtools
	#get rid of double counts
	#######
	######
	######
	#set dictionaries and variables
	os.system("cat SC_{file2}_to_{file1}_compare.bed| cut -f13 | sort  -n|uniq > intermediate_file.txt".format(**locals()))
	intermediate = "intermediate_file.txt"
	intermediate = open(intermediate, "r")
	cutoff_distances={}
	for line in intermediate:
		line = line.rstrip('\n')
		cutoff_distances[line]=0


	for num in sorted(cutoff_distances.keys()):
		my_in_file = open("SC_{file2}_to_{file1}_compare.bed".format(**locals()), "r")
		result, err = Popen(['wc -l %s' %file1], stdout=PIPE, stderr=PIPE, shell=True).communicate()
		match = re.search("(\d+)\s",result)
		total_transposons = int(match.group(1))
		result, err = Popen(['wc -l %s' %file2], stdout=PIPE, stderr=PIPE, shell=True).communicate()
		match = re.search("(\d+)\s",result)
		total_transposons2 = int(match.group(1))
		family_total={}
		family_count={}
		family_matched_count={}
		count = 0
		total_family_matched_count=0
		#get total counts of the features found in file1

	#####account for ties


		print num
		my_distance = int(num)
		print num
		print my_distance
		aa= my_distance
		family_total={}
		family_count={}
		family_matched_count={}
		count = 0
		total_family_matched_count=0
		for line in my_in_file:
			items= re.split("[\t]",line)
			ID1 = items[3]
			ID2 = items[9]
			distance = items[12]
			distance = distance.rstrip('\n')
			distance =int(distance)
			#create a key for each transposons family if it does not already exists
			if ID1 in family_total.keys():
				family_total[ID1] += 1
			else:
				family_total[ID1] = 1
				family_count[ID1] = 0
				family_matched_count[ID1] = 0
			# count a feature in file1 as matched if the distance of the closeset feature in file2 is greater than zero and less than the max distance specified in input
			if distance <= my_distance and distance >=0:
				print "FOUND"
				print "{distance}\t{my_distance}".format(**locals())

				count += 1
				###################################
				###################################
				###################################


				#increase the matched count if the closest feature is within X base pairs distance
				family_count[ID1] += 1
				if re.search("(temp)|(telocate)|(retroseq)",file2):
					print "THE FILE IS {file2}".format(**locals())
					match = re.search("(.*)_(\w+-)?reference",ID2)
					family_found = match.group(1)
				else:
					family_found= ID2
				#if the family of the closest feature found matches what it should be, increase the matched correctly count
				if family_found == ID1:
					print "THE FAMILY WAS FOUND"
					family_matched_count[ID1] += 1
					total_family_matched_count += 1
				else:
					print "THE FAMILY{family_found} did not match {ID1}".format(**locals())

		not_found = total_transposons-count
		fraction_not_found = float (not_found)/(total_transposons) *100
		fraction = float (count)/(total_transposons) *100
		fraction_fam_correct = float (total_family_matched_count)/(total_transposons) *100

		not_found2 = total_transposons2-count
		fraction_not_found2= float (not_found2)/(total_transposons2) *100
		fraction2 = float (count)/(total_transposons2) *100

		not_found2_fam = total_transposons2-total_family_matched_count
		fraction_not_found2_fam= float (not_found2_fam)/(total_transposons2) *100
		fraction2_fam = float (total_family_matched_count)/(total_transposons2) *100






		my_out_file_summary.write("{file1}\t{file2}\t{my_distance}\tfamily_aware\t{fraction_fam_correct}\t{fraction_not_found2_fam}\n".format(**locals()))
		my_out_file_summary.write("{file1}\t{file2}\t{my_distance}\toverall\t{fraction}\t{fraction_not_found2}\n".format(**locals()))

		my_in_file.close()

#################################################################################################
#################################################################################################




#run the function for each nonredudant output bed file in the directory
for results_file in os.listdir(my_directory):
	#match = re.findall("(temp|telocate|retroseq)_nonredundant.bed$",results_file)
	match = re.findall("(temp|telocate|retroseq)_nonredundant.bed$",results_file)
	if len(match) >0:
		my_file = sys.argv[4]
		print my_file
		os.system("cat {results_file} | sed s'/chrX/X/g'| sed '/^track/d' > {results_file}_temp && mv {results_file}_temp {results_file}".format(**locals()))
		os.system("cat {my_file} | sed s'/chrX/X/g'|sed '/^track/d' > {my_file}_temp && mv {my_file}_temp {my_file}".format(**locals()))
		file_name = match[0]
		collapse_transposons(results_file)
		results_file = "new_CT_{results_file}".format(**locals())
		compare_beds(my_file,results_file)
		TFPN(my_file,results_file)
		print match
os.system("mv *CT* intermediates")
#os.system("cat {results_file} | sed s'/chrX/X/g'| sed '/^track/d' | awk '$4 ~ /_non-reference/ {print $0}' > {results_file}_temp && mv {results_file}_temp {results_file}".format(**locals()))
