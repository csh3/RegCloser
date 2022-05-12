# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

import sys
import re
from Bio.Seq import Seq
import os
import argparse

def CalcSeqLen(cigar):
    sections=re.split('([IMDSH])', cigar)
    length=0
    for k in range(len(sections)):
        if sections[k]=='M' or sections[k]=='D':
            length+=int(sections[k-1])
    return(length)

def ParseRecord(record):
    sections=re.split('([IMDSH])', record[5])
    leftClip=rightClip=0
    mapped=False
    if sections[1] in ['S','H']:
        leftClip+=int(sections[0])
    if sections[-2] in ['S','H']:
        rightClip+=int(sections[-3])
    if leftClip<=maxClip and rightClip<=maxClip:
        mapped=True
    if record[2][-1] in ['A','L','R']:
        contigName=record[2][:-2]
    else:
        contigName=record[2]
    start = int(record[3])-leftClip
    end = int(record[3])+CalcSeqLen(record[5])-1+rightClip
    if record[2].endswith('R'):
        start = Length[contigName]-trimLength+start
        end = Length[contigName]-trimLength+end
    return(mapped, contigName, start, end)


def HighCoverage(contigName, start, end):
    highCoverage=False
    if contigName in highCoverageRegions:
        for region in highCoverageRegions[contigName]:
            if not (end < region[0] or start > region[1]): # start>region[0] and end<region[1]:   #
                highCoverage=True
    return(highCoverage)


def CollectGapReads(contigName1, flag1, start1, end1, record2):
    if not HighCoverage(contigName1, start1, end1):
        if direction=='FR' and flag1 in ['0','256','2048'] and Length[contigName1]-start1+1<=(mean+ml*std):
            if record2[1] in ['0','4']:
                seq=str(Seq(record2[9]).reverse_complement())
                phred=record2[10][::-1]
            else:
                seq=record2[9]
                phred=record2[10]
            with open('FILEs/%s_R.fastq'%contigName1,'a') as fout:
                fout.write('@%s:%d:%d\n'%(record2[0], mean-(Length[contigName1]-start1+1), std))
                fout.write(seq+'\n')
                fout.write('+\n')
                fout.write(phred+'\n')
        elif direction=='FR' and flag1 in ['16','272','2064'] and end1<=(mean+ml*std):
            if record2[1] in ['0','4']:
                seq=record2[9]
                phred=record2[10]
            else:
                seq=str(Seq(record2[9]).reverse_complement())
                phred=record2[10][::-1]
            with open('FILEs/%s_L.fastq'%contigName1,'a') as fout:
                fout.write('@%s:%d:%d\n'%(record2[0], mean-end1, std))
                fout.write(seq+'\n')
                fout.write('+\n')
                fout.write(phred+'\n')
        elif direction=='RF' and flag1 in ['0','256','2048'] and end1<=(mean+ml*std):
            if record2[1] in ['0','4']:
                seq=str(Seq(record2[9]).reverse_complement())
                phred=record2[10][::-1]
            else:
                seq=record2[9]
                phred=record2[10]
            with open('FILEs/%s_L.fastq'%contigName1,'a') as fout:
                fout.write('@%s:%d:%d\n'%(record2[0], mean-end1, std))
                fout.write(seq+'\n')
                fout.write('+\n')
                fout.write(phred+'\n')
        elif direction=='RF' and flag1 in ['16','272','2064'] and Length[contigName1]-start1+1<=(mean+ml*std):
            if record2[1] in ['0','4']:
                seq=record2[9]
                phred=record2[10]
            else:
                seq=str(Seq(record2[9]).reverse_complement())
                phred=record2[10][::-1]
            with open('FILEs/%s_R.fastq'%contigName1,'a') as fout:
                fout.write('@%s:%d:%d\n'%(record2[0], mean-(Length[contigName1]-start1+1), std))
                fout.write(seq+'\n')
                fout.write('+\n')
                fout.write(phred+'\n')


descript="This program scans sam files to Collect reads in the gap regions for local assembly and make tab files of linking information between contigs.\n"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-l', required=True, help='Reads library, required.')
parser.add_argument('-d', required=True, help='Orientation of library, FR (paired-end) or (RF) mate-pair, required.')
parser.add_argument('-m', type=int, required=True, help='Average insert size, required.')
parser.add_argument('-s', type=int, required=True, help='Standard deviation of the insert size, required.')
parser.add_argument('-ml', type=int, default=3,help="Contig end with insert size + ml * std bases were cut out for mapping reads. [3]")
parser.add_argument('-tl', type=int, required=True, help='Cut out length of the contig end, required.')
parser.add_argument('-e', required=True, help='Whether the library is used for local assembly, Y or N, required.')
parser.add_argument('-nr', action="store_true", help='Not re-map anchored reads to the whole draft genome. [null]')
parser.add_argument('-c', type=int, default=5, help='Maximum number of soft-clipped bases on either end of a mapped read. [5]')
parser.add_argument('-mq', type=int, default=60, help='Mapped Reads with mapping quality greater than this value will be identified as anchored reads. [60]')

args = parser.parse_args()

maxClip = args.c
library = args.l
direction = args.d
mean = args.m
std = args.s
extension = args.e
trimLength = args.tl
ml = args.ml
mq = args.mq

Length={}
with open('contigLength.txt') as fin:
    for line in fin:
        record = line.strip().split()
        contigName = record[0]
        contigLen = int(record[1])
        Length[contigName]=contigLen

highCoverageRegions={}
with open('highCoverage.bed') as fin:
    for line in fin:
        record = line.strip().split()
        contigName = record[0]
        start = int(record[1])
        end = int(record[2])
        if contigName in highCoverageRegions:
            highCoverageRegions[contigName].append([start, end])
        else:
            highCoverageRegions[contigName]=[[start, end]]

with open('%s_remap_R1.sam'%library) as fin_1:
    with open('%s_remap_R2.sam'%library) as fin_2:
        with open('tab%s.txt'%library,'w') as fout1:
            read1_backup=''; read2_backup=''
            record1_backup_list=[]; record2_backup_list=[]
            while(True):
                while(True):
                    record1=fin_1.readline().strip().split()
                    if len(record1)==0:
                        break
                    elif not record1[0].startswith('@') and record1[0] != read1_backup:
                        break
                    elif not record1[0].startswith('@'):
                        record1_backup_list.append(record1)

                while(True):
                    record2=fin_2.readline().strip().split()
                    if len(record2)==0:
                        break
                    elif not record2[0].startswith('@') and record2[0] != read2_backup:
                        break
                    elif not record2[0].startswith('@'):
                        record2_backup_list.append(record2)

                if record1_backup_list == [] or record2_backup_list == []:
                    read1_backup=record1[0]; read2_backup=record2[0]
                    record1_backup_list=[record1]; record2_backup_list=[record2]
                    continue

                mappedNum1=0
                for record1_backup in record1_backup_list:
                    if record1_backup[1]!='4' and int(record1_backup[4])>=mq:
                        mapped, contigName, start, end = ParseRecord(record1_backup)
                        if mapped:
                            mappedNum1 += 1
                            flag1 = record1_backup[1]
                            contigName1, start1, end1=contigName, start, end

                mappedNum2=0
                for record2_backup in record2_backup_list:
                    if record2_backup[1]!='4' and int(record2_backup[4])>=mq:
                        mapped, contigName, start, end = ParseRecord(record2_backup)
                        if mapped:
                            mappedNum2 += 1
                            flag2 = record2_backup[1]
                            contigName2, start2, end2=contigName, start, end

                if mappedNum1==1 and mappedNum2==1:
                    if (not HighCoverage(contigName1, start1, end1)) and (not HighCoverage(contigName2, start2, end2)):
                        if flag1 in ['0','2048']:
                            fout1.write(contigName1+'\t'+str(start1)+'\t'+str(end1)+'\t')
                        else:
                            fout1.write(contigName1+'\t'+str(end1)+'\t'+str(start1)+'\t')
                        if flag2 in ['0','2048']:
                            fout1.write(contigName2+'\t'+str(start2)+'\t'+str(end2)+'\n')
                        else:
                            fout1.write(contigName2+'\t'+str(end2)+'\t'+str(start2)+'\n')

                if extension == 'Y' and not args.nr:
                    if mappedNum1==1:
                        CollectGapReads(contigName1, flag1, start1, end1, record2_backup_list[0])
                    if mappedNum2==1:
                        CollectGapReads(contigName2, flag2, start2, end2, record1_backup_list[0])

                if len(record1)==0 and len(record2)==0:
                    break
                else:
                    read1_backup=record1[0]; read2_backup=record2[0]
                    record1_backup_list=[record1]; record2_backup_list=[record2]

if extension == 'Y' and args.nr:
    with open('%s_R1.sam'%library) as fin_1:
        with open('%s_R2.sam'%library) as fin_2:
            read1_backup=''; read2_backup=''
            record1_backup_list=[]; record2_backup_list=[]

            while(True):

                while(True):
                    record1=fin_1.readline().strip().split()
                    if len(record1)==0:
                        break
                    elif not record1[0].startswith('@') and record1[0] != read1_backup:
                        break
                    elif not record1[0].startswith('@'):
                        record1_backup_list.append(record1)

                while(True):
                    record2=fin_2.readline().strip().split()
                    if len(record2)==0:
                        break
                    elif not record2[0].startswith('@') and record2[0] != read2_backup:
                        break
                    elif not record2[0].startswith('@'):
                        record2_backup_list.append(record2)

                if record1_backup_list == [] or record2_backup_list == []:
                    read1_backup=record1[0]; read2_backup=record2[0]
                    record1_backup_list=[record1]; record2_backup_list=[record2]
                    continue

                if record1_backup_list[0][1]!='4':
                    for record1_backup in record1_backup_list:
                        if (int(record1_backup_list[0][4])>=mq and int(record1_backup[4])>=mq) or (int(record1_backup_list[0][4])==0 and int(record1_backup[4])==0):
                            mapped1, contigName1, start1, end1 = ParseRecord(record1_backup)
                            if mapped1:
                                CollectGapReads(contigName1, record1_backup[1], start1, end1, record2_backup_list[0])

                if record2_backup_list[0][1]!='4':
                    for record2_backup in record2_backup_list:
                        if (int(record2_backup_list[0][4])>=mq and int(record2_backup[4])>=mq) or (int(record2_backup_list[0][4])==0 and int(record2_backup[4])==0):
                            mapped2, contigName2, start2, end2 = ParseRecord(record2_backup)
                            if mapped2:
                                CollectGapReads(contigName2, record2_backup[1], start2, end2, record1_backup_list[0])
 
                if len(record1)==0 and len(record2)==0:
                    break
                else:
                    read1_backup=record1[0]; read2_backup=record2[0]
                    record1_backup_list=[record1]; record2_backup_list=[record2]