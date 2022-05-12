# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

import sys
import re
import argparse
from Bio.Seq import Seq


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
    start = int(record[3])
    end = int(record[3])+CalcSeqLen(record[5])-1
    if leftClip<=maxClip and rightClip<=maxClip:
        mapped=True
    return(mapped)


descript="This program screens out reads mapped to the contig ends.\n"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-f1', required=True, help='Fastq file to output one end of paired reads, required.')
parser.add_argument('-f2', required=True, help='Fastq file to output the other end of paired reads, required.')
parser.add_argument('-l', required=True, help='Reads library, required.')
parser.add_argument('-d', required=True, help='Orientation of library, FR (paired-end) or (RF) mate-pair, required.')
parser.add_argument('-c', type=int, default=5, help='Maximum number of soft-clipped bases on either end of a mapped read. [5]')
parser.add_argument('-mq', type=int, default=60, help='Mapped Reads with mapping quality greater than this value will be identified as anchored reads. [60]')

args = parser.parse_args()

maxClip = args.c
library = args.l
direction=args.d
f1=args.f1
f2=args.f2
mq=args.mq

readNum = 0

refLen={}
with open('%s_R1.sam'%library) as fin_1:
    with open('%s_R2.sam'%library) as fin_2:
        with open(f1,'w') as fout_1:
            with open(f2,'w') as fout_2:

                read1_backup=''; read2_backup=''
                record1_backup_list=[]; record2_backup_list=[]

                while(True):

                    while(True):
                        record1=fin_1.readline().strip().split()
                        if len(record1)==0:
                            break
                        elif record1[0].startswith('@') and len(record1)==3:
                            refLen[record1[1][3:]]=int(record1[2][3:])
                        elif not record1[0].startswith('@') and record1[0] != read1_backup:
                            break
                        elif not record1[0].startswith('@'):
                            record1_backup_list.append(record1)

                    while(True):
                        record2=fin_2.readline().strip().split()
                        if len(record2)==0:
                            break
                        elif record2[0].startswith('@') and len(record2)==3:
                            refLen[record2[1][3:]]=int(record2[2][3:])
                        elif not record2[0].startswith('@') and record2[0] != read2_backup:
                            break
                        elif not record2[0].startswith('@'):
                            record2_backup_list.append(record2)

                    if record1_backup_list == [] or record2_backup_list == []:
                        read1_backup=record1[0]; read2_backup=record2[0]
                        record1_backup_list=[record1]; record2_backup_list=[record2]
                        continue

                    record1_backup = record1_backup_list[0]; record1_mapped = False
                    if (direction=='FR' and record1_backup[1]=='0' and not record1_backup[2].endswith('L')) or (direction=='FR' and record1_backup[1]=='16' and not record1_backup[2].endswith('R')) or (direction=='RF' and record1_backup[1]=='0' and not record1_backup[2].endswith('R')) or (direction=='RF' and record1_backup[1]=='16' and not record1_backup[2].endswith('L')):
                        mapped = ParseRecord(record1_backup)
                        if mapped and int(record1_backup[4])>=mq:
                            record1_mapped=True
                    
                    record2_backup = record2_backup_list[0]; record2_mapped = False
                    if (direction=='FR' and record2_backup[1]=='0' and not record2_backup[2].endswith('L')) or (direction=='FR' and record2_backup[1]=='16' and not record2_backup[2].endswith('R')) or (direction=='RF' and record2_backup[1]=='0' and not record2_backup[2].endswith('R')) or (direction=='RF' and record2_backup[1]=='16' and not record2_backup[2].endswith('L')):
                        mapped = ParseRecord(record2_backup)
                        if mapped and int(record2_backup[4])>=mq:
                            record2_mapped=True

                    if record1_mapped or record2_mapped:
                        readNum+=1
                        if record1_backup[1] in ['0','4']:
                            seq=record1_backup[9]
                            quality=record1_backup[10]
                        else:
                            seq=str(Seq(record1_backup[9]).reverse_complement())
                            quality=record1_backup[10][::-1]
                        fout_1.write('@%d\n'%readNum)
                        fout_1.write(seq+'\n')
                        fout_1.write('+\n')
                        fout_1.write(quality+'\n')

                        if record2_backup[1] in ['0','4']:
                            seq=record2_backup[9]
                            quality=record2_backup[10]
                        else:
                            seq=str(Seq(record2_backup[9]).reverse_complement())
                            quality=record2_backup[10][::-1]
                        fout_2.write('@%d\n'%readNum)
                        fout_2.write(seq+'\n')
                        fout_2.write('+\n')
                        fout_2.write(quality+'\n')

                    if len(record1)==0 and len(record2)==0:
                        break
                    else:
                        read1_backup=record1[0]; read2_backup=record2[0]
                        record1_backup_list=[record1]; record2_backup_list=[record2]
