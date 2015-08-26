#!/usr/bin/env python 
import os
import re
from subprocess import Popen, PIPE
import statistics
my_dir="/lscr2/andersenlab/kml436/round21_Aug19/round19_Aug13"

means=[]
SDs=[]
###CHANGE TO 1 through 8
run_list=[1,2,3,4,5,6,7,8]
for i in run_list:
	print "Processing Run:"
	print i

	run="run_{i}".format(**locals()) 
	bam_file="{my_dir}/{run}_N2/{run}_N2.sorted.bam".format(**locals()) 
	result, err = Popen(["""samtools depth {bam_file}| datamash mean 3 sstdev 3""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()
	result=result.split('\t')
	mean=float(result[0])
	SD=float(result[1])
	means.append(mean)
	SDs.append(SD)
	print mean
	print SD

mean_of_means=statistics.mean(means)
mean_of_SDs=statistics.mean(SDs)
two_SDS=2*mean_of_SDs
three_SDS=3*mean_of_SDs
four_SDS=4*mean_of_SDs

print mean_of_means
print mean_of_SDs	
print two_SDS
print three_SDS
print four_SDS

