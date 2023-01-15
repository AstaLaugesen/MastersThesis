#!/usr/bin/env python3

import sys,re
from statistics import mean

#Add to list of bacteria that we are interested in (top 50 read amounts)
with open(sys.argv[1]) as file_in:
	BactList = []
	for line in file_in:
		BactList.append(line.strip())
        
#get total length in dict
total_length_dict = {}
length_file = open(sys.argv[2])
for line in length_file:
	key, value = line.split()
	total_length_dict[key] = int(value)

#get total reads in dict
total_reads_dict = {}
reads_file = open(sys.argv[3])
for line in reads_file:
	key, value = line.split()
	total_reads_dict[key] = int(value)

#Read the infile (defense finder files) once
DefSys = []
for line in sys.stdin:
	DefSys.append(line.strip())


defsys_reads_dict={}
defsys_length_dict={}
#find data for the defense system in question
for bacteria in BactList:
	defsys_reads_dict[bacteria] = [0]
	defsys_length_dict[bacteria] = [0]
	for line in DefSys:
		##                                                 CHANGE DEFENSE SYSTEM NAME HERE !!!
		if ((re.search('(\s)'+bacteria+'(\s)',line)) and (re.search(r'\sRM\s',line))):
			if defsys_reads_dict[bacteria]==[0]:
				defsys_reads_dict[bacteria]= [int( re.search(r'\s(\d+)$', line).group(1) )]
				defsys_length_dict[bacteria]= [int( re.search(r'\s(\d+)\s\d+$', line).group(1) )]
			else:
				defsys_reads_dict[bacteria].append( int( re.search(r'\s(\d+)$', line).group(1) ) )
				defsys_length_dict[bacteria].append( int( re.search(r'\s(\d+)\s\d+$', line).group(1) ) )


#find results
result_dict={}
for bacteria in BactList:
	if ( (defsys_length_dict[bacteria]==[0]) or (defsys_reads_dict[bacteria]==[0]) ):
		result_dict[bacteria]=0
	else:
		result_dict[bacteria]= ( mean(defsys_reads_dict[bacteria])/mean(defsys_length_dict[bacteria]) ) / ( total_reads_dict[bacteria]/total_length_dict[bacteria])


                                           
#Setting name of sample
SampleName=0
#print(sys.argv[4])
SampleName=re.search(r'results/(\S+)',sys.argv[4]).group(1)

#Printing results:
print(SampleName, end="\t"), print(*result_dict.values(), sep='\t')

