# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

import argparse
import os
import numpy as np
#from statsmodels.robust.scale import mad
		
def GetContigLen():
	DetailFile = open('contigLength.txt','r')
	Contiglen_list = {}
	for line in DetailFile:
		detail = line.strip().split("\t")
		ctg_n = int(detail[0][6:])
		Contiglen_list[ctg_n]=int(detail[1])
	DetailFile.close()
	return Contiglen_list

def CalMPGap(record,InsertSize):# 	<—— ingap ——>
	global contigLength_dict
	c1=int(record[0][6:])
	c2=int(record[3][6:])
	p1 = int(record[1])
	m1 = int(record[2])
	p2 = int(record[4])
	m2 = int(record[5])
	d1 = d2 =0
	if p1>m1:
		d1 = contigLength_dict[c1]-m1	
	else:
		d1 = m1
	if p2<m2:
		d2 = m2
	else:
		d2 = contigLength_dict[c2]-m2
	return InsertSize-d1-d2

def CalPEGap(record,InsertSize): # ——> ingap <—— 
	global contigLength_dict
	c1=int(record[0][6:])
	c2=int(record[3][6:])
	p1 = int(record[1])
	m1 = int(record[2])
	p2 = int(record[4])
	m2 = int(record[5])
	d1 = d2 =0
	if p1<m1:
		d1 = contigLength_dict[c1]-p1	
	else:
		d1 = p1
	if p2>m2:
		d2 = p2
	else:
		d2 = contigLength_dict[c2]-p2
	return (InsertSize-d1-d2)

def FindGap_In_Dict(gap_name):
	for libname in RawGap:
		if gap_name in RawGap[libname]:
			gap_list = RawGap[libname][gap_name]
			sd = (library_info[libname]['std'])/np.sqrt(len(gap_list))
			return(libname,(np.median(gap_list),int(sd)))
	return(False,0)

def InferFromLargerGap(FinalGap_dict):
	global contigLength_dict
	for name,item in FinalGap_dict.items():
		if item[1]==0:
			ctg1,ctg2 = list(map(lambda x: int(x), name.split('-')))
			ctg3 = ctg2+1
			ctg0 = ctg1-1
			name_01 = str(ctg0)+'-'+str(ctg1)
			name_02 = str(ctg0)+'-'+str(ctg2)
			name_13 = str(ctg1)+'-'+str(ctg3)
			name_23 = str(ctg2)+'-'+str(ctg3)
			gap_1 = gap_2 = 0
			if1, gap_01 = FindGap_In_Dict(name_01)
			if2, gap_02 = FindGap_In_Dict(name_02)
			if if1 and if2:
				gap_1 = (gap_02[0]-gap_01[0]-contigLength_dict[ctg1], gap_02[1]+gap_01[1])

			if3, gap_13 = FindGap_In_Dict(name_13)
			if4, gap_23 = FindGap_In_Dict(name_23)
			if if3 and if4:
				#print(gap_13,gap_23)
				gap_2 = (gap_13[0]-gap_23[0]-contigLength_dict[ctg2], gap_13[1]+gap_23[1])
			
			if gap_1!=0 and gap_2!=0:
				FinalGap_dict[name]=int(((gap_1[0]+gap_2[0])/2)),int(np.sqrt(((gap_1[1])**2+(gap_2[1])**2)/2))
			elif gap_2!=0:
				FinalGap_dict[name]=gap_2
			elif gap_1!=0:
				FinalGap_dict[name]=gap_1

descript="This program re-estimate gap sizes based on tabfiles.\n"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-p', type=str, default='../prerequisite', help='A formatted file specifying the code path, reads directory, and library information. [Prerequisite]')
parser.add_argument('-sd', type=int, default=100, help="Default standard deviation of gap size. [100]")
args = parser.parse_args()

#libtype = args.b
#libmean = args.m
contigLength_dict = GetContigLen()
library_info={}
with open(args.pre,'r') as fin: #'prerequisite'
	longest_ins = 100
	longest_lib = ' '
	for line in fin:
		record = line.strip().split()
		if len(record) == 7:
			library_info[record[0]]={'filename1':record[1],'filename2':record[2], 'mean':int(record[3]), 'std':int(record[4]),'direction':record[5],'extension':record[6]}
			if int(record[3])>longest_ins:
				longest_ins = int(record[3])
				longest_lib = record[0]
		else:
			reads_directory=line.strip().split(':')[1].strip()
			if reads_directory.endswith('/'):
				reads_directory=reads_directory[:-1]

RawGap={}
#ShortGap_dict={}
#LongGap_dict={}
for libname,library in library_info.items():
	tab = open('tab%s.txt'%libname,'r')
	print ("open tab%s.txt"%libname)
	Gap_dict = {}
	for tabline in tab:
		record=tabline.strip().split('\t')
		contig1=int(record[0][6:])
		contig2=int(record[3][6:])
		if abs(int(contig1)-int(contig2))==1 or abs(int(contig1)-int(contig2))==2:
			gap_name = str(min(contig1,contig2))+'-'+str(max(contig1,contig2))
			if library['direction']=='FR':
				gap = CalPEGap(record,library['mean'])
			elif library['direction']=='RF':
				gap = CalMPGap(record,library['mean'])
			if gap_name in Gap_dict:
				Gap_dict[gap_name]+=[gap]
			else:
				Gap_dict[gap_name]=[gap]
	tab.close()
	RawGap[libname] = Gap_dict

FinalGap_dict = {'0':0}
gap_name = '0'
evidence = open("evidence.txt",'r')
for line in evidence:
	evidence_record=line.strip().split()
	if len(evidence_record)<4:
		if len(evidence_record)==1:
			FinalGap_dict.pop(gap_name)
		continue
	gap_name = str(evidence_record[0])+'-'+str(int(evidence_record[0])+1)
	FinalGap_dict[gap_name] = (int(evidence_record[4]),args.sd) ##default: gap = initial gap, sd =100
FinalGap_dict.pop(gap_name)

for gap in FinalGap_dict:
	gap_list = []
	sdsum = 0
	for libname in RawGap:
		Gap_dict = RawGap[libname]
		if gap in Gap_dict:
			gap_list +=Gap_dict[gap]
			sdsum +=len(Gap_dict[gap])*(library_info[libname]['std'])**2
	if len(gap_list)>0:
		final_median = np.median(gap_list)
		sd = np.sqrt(sdsum)/len(gap_list)
		FinalGap_dict[gap] = (int(final_median),int(sd))

# DictFile = open("Gap_dict",'w')
# DictFile.write(str(FinalGap_dict))
# DictFile.close()
# os.system(r"sed -i 's/),/)\n/g' Gap_dict")

InferFromLargerGap(FinalGap_dict)

# DictFile = open("Gap_dict_iter1",'w')
# DictFile.write(str(FinalGap_dict))
# DictFile.close()
# os.system(r"sed -i 's/),/)\n/g' Gap_dict_iter1")

InferFromLargerGap(FinalGap_dict)

with open("evidence_estimate_sd.txt",'w') as fgap:
	evidence = open("evidence.txt",'r')
	scaf_end = False
	for line in evidence:
		evidence_record=line.strip().split()
		if len(evidence_record)<4:
			print(line, file=fgap, end='')
		else:
			gap_name = str(evidence_record[0])+'-'+str(int(evidence_record[0])+1)
			if gap_name in FinalGap_dict:
				gap,sd=FinalGap_dict[gap_name]
				print(evidence_record[0],evidence_record[1],evidence_record[2],evidence_record[3], int(gap), int(sd), sep='\t', file=fgap)
			else:
				print(evidence_record[0],evidence_record[1],evidence_record[2],evidence_record[3], 0, 0, sep='\t', file=fgap) #line, file=fgap,end='')
	evidence.close()
fgap.close()
	
# DictFile = open("FinalGap_dict",'w')
# DictFile.write(str(FinalGap_dict))
# DictFile.close()
# os.system(r"sed -i 's/),/)\n/g' FinalGap_dict")
