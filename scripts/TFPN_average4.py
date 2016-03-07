import sys
import re
import os
from subprocess import Popen, PIPE
from  collections import defaultdict
import statistics

out ="BEDCOMPARE_MEANS.txt"
OUT = open(out, "w")
OUT.write ("M1\tM2\tDistanceCutoff\tfam\tTPR\tFDR\n")



###START FUNCTION HERE
def calculate_mean(TE_program,raw_or_filter,fam):
    my_directory = sys.argv[1] #argv[0] would be script name
    files=''
    items = []
    for results_file in os.listdir(my_directory):
        match = re.findall("(BEDCOMPARE_SUMMARY*)",results_file)
        if len(match) >0:
            print "yes"
            print results_file
            files+=str(" {results_file}".format(**locals()))
    result, err = Popen(["""cat %s | awk '{if($4 ~/%s/  && $2 ~/%s/ && $2 ~/%s/) {print $3}}'|uniq""" %(files,fam,TE_program,raw_or_filter)], stdout=PIPE, stderr=PIPE, shell=True).communicate()
    print result
    average_TPR={}
    average_FDR={}
    average_TPR = defaultdict(list)
    average_FDR = defaultdict(list)

    unique_distances = {}
    result=result.split('\n')
    for i in result:
        if i != '': #get rid of balnk key
            unique_distances[i]=0


    #get rid of leading space...next tiem append space after file:
    files= files[1:]
    files_to_test = files.split(' ')
    for sim_file in files_to_test:
        sim_TPR = {}
        sim_FDR = {}
        OPEN_SIM_FILE = open(sim_file, "r")
        for line in OPEN_SIM_FILE:
            if re.search(TE_program,line):
                if re.search(raw_or_filter,line):
                    #condense this next line
                    #CHANGE TO FAMILY AWARE LATER
                    #CHANGEEEEEEEEEEEE
                    if re.search (fam, line):
                        items= re.split("[\t]",line)
                        M1 = items[0]
                        M2 = items[1]
                        dist = items[2]
                        TPR = items[4]
                        FDR = items[5]
                        sim_TPR[dist]=TPR
                        sim_FDR[dist]=FDR


        if "0" in sim_TPR.keys():
            print "YES 0 value is present"
        else:
            TPR_value=0
            FDR_value=0
        for uk in sorted(unique_distances.keys(),key=int): #CHANGE TO KEY=INT
            if uk in sim_TPR.keys():
                print "THE KEY IS {uk}".format(**locals())
                TPR_value=sim_TPR[uk]
                FDR_value=sim_FDR[uk]
                average_TPR[uk].append(TPR_value)
                average_FDR[uk].append(FDR_value)

            else:
                average_TPR[uk].append(TPR_value)
                average_FDR[uk].append(FDR_value)
                #print "{uk} is not in the keys".format(**locals())

                ##get rid of the zero later
    for key in sorted(average_TPR.keys(),key=int):
        print key
        print average_TPR[key]
        average_TPR[key] = map(float, average_TPR[key]) #convert strings in list to integers
        average_FDR[key] = map(float, average_FDR[key]) 
        mean_TPR = statistics.mean(average_TPR[key])
        standard_deviation_TPR = statistics.pstdev(average_TPR[key])
        mean_FDR = statistics.mean(average_FDR[key])
        standard_deviation_FDR = statistics.pstdev(average_FDR[key])
        print "The mean_TPR is {mean_TPR}".format(**locals())
        print "The standard deviation TPR is {standard_deviation_TPR}".format(**locals())
        OUT.write ("{M1}\t{M2}\t{key}\t{fam}\t{mean_TPR}\t{mean_FDR}\n".format(**locals()))
        OUT.write ("{M1}\t{M2}_error\t{key}\t{fam}\t{standard_deviation_TPR}\t{standard_deviation_FDR}\n".format(**locals()))

#calculate_mean("temp")
#calculate_mean("temp","_NF","overall")
#calculate_mean("temp","_F","overall") #last setting
#calculate_mean("temp","_NF","family_aware")
calculate_mean("temp","_F","family_aware") #last setting
#calculate_mean("telocate")
#calculate_mean("telocate","_NF","overall")
#calculate_mean("telocate","_F","overall") #last setting
#calculate_mean("telocate","_NF","family_aware")
calculate_mean("telocate","_F","family_aware") #last setting
#calculate_mean("retroseq")
#calculate_mean("retroseq","_NF","overall")
#calculate_mean("retroseq","_F","overall")
#calculate_mean("retroseq","_NF","family_aware")
calculate_mean("retroseq","_F","family_aware")
OUT.close()
