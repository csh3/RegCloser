# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

from Bio import SeqIO
import re
import argparse
# from multiprocessing import Pool,Manager
from scipy.stats import mode

descript="This program identifies high coverage regions in the contig ends.\n"

parser = argparse.ArgumentParser(description=descript)

# parser.add_argument('-t', type=int, default=1, help='number of threads for parallelism')
parser.add_argument('-hs', nargs='+', help='Input sam files used for identifying high coverage regions, required.')
parser.add_argument('-tl', type=int, required=True, help='Cut out length of the contig end, required.')
parser.add_argument('-c', type=int, default=5, help='Maximum number of soft-clipped bases on either end of a mapped read. [5]')
parser.add_argument('-mq', type=int, default=60, help='Mapped Reads with mapping quality greater than this value will be identified as anchored reads. [60]')
parser.add_argument('-ra', type=float, default=1.8, help='Consecutive bases with coverage higher then ra * mode coverage will be marked as high coverage regions. [1.8]')

args = parser.parse_args()

def CalcSeqLen(cigar):
    sections=re.split('([IMDSH])', cigar)
    length=0
    for k in range(len(sections)):
        if sections[k]=='M' or sections[k]=='D':
            length+=int(sections[k-1])
    return(length)

def ParseRecord(record):
    contigName=record[2][:-2]
    start = int(record[3])
    end = int(record[3])+CalcSeqLen(record[5])-1
    sections=re.split('([IMDSH])', record[5])
    leftClip=rightClip=0
    mapped=False
    if sections[1] in ['S','H']:
        leftClip+=int(sections[0])
    if sections[-2] in ['S','H']:
        rightClip+=int(sections[-3])
    if leftClip<=maxClip and rightClip<=maxClip:
        mapped=True
    return(mapped, contigName, start, end)

def countCoverage(filename, Coverage):
    refLen={}
    with open(filename) as fin:
        read_backup=''; record_backup_list=[]
        while(True):

            while(True):
                record=fin.readline().strip().split()
                if len(record)==0:
                    break
                elif record[0].startswith('@') and len(record)==3:
                    refLen[record[1][3:]]=int(record[2][3:])
                elif not record[0].startswith('@') and record[0] != read_backup:
                    break
                elif not record[0].startswith('@'):
                    record_backup_list.append(record)

            if record_backup_list == []:
                read_backup=record[0]
                record_backup_list=[record]
                continue

            if record_backup_list[0][1]!='4':
                for record_backup in record_backup_list:
                    if (int(record_backup_list[0][4]) >= mq and int(record_backup[4]) >= mq) or (int(record_backup_list[0][4]) == 0 and int(record_backup[4]) == 0):
                        mapped, contigName, start, end = ParseRecord(record_backup)
                        if mapped:
                            if record_backup[2].endswith('A') or record_backup[2].endswith('L'):
                                for i in range(start-1,end):
                                    Coverage[contigName][i]+=1
                            elif record_backup[2].endswith('R'):
                                for i in range(start-trimLength-1,end-trimLength):
                                    Coverage[contigName][i]+=1

            if len(record)==0:
                break
            else:
                read_backup=record[0]
                record_backup_list=[record]


# thread = args.t
maxClip = args.c
trimLength = args.tl
samFiles= args.hs
ratio = args.ra
mq = args.mq

Coverage={}
Length={}
with open('contigLength.txt') as fin:
    for line in fin:
        record = line.strip().split()
        contigName = record[0]
        contigLen = int(record[1])
        Length[contigName]=contigLen
        if contigLen<=2*trimLength:
            Coverage[contigName]=[0]*contigLen
        else:
            Coverage[contigName]=[0]*2*trimLength

# m=Manager()
# Coverage=m.dict(Coverage)
# pool=Pool(thread)
# for filename in samFiles:
#     pool.apply_async(func=countCoverage, args=(filename, Coverage))
# pool.close()
# pool.join()

for filename in samFiles:
    print("open", filename)
    countCoverage(filename, Coverage)

coverageList=[]
for contigName in Coverage:
    coverageList+=Coverage[contigName]
threshold=mode(coverageList)[0][0]*ratio
if threshold<5:
    threshold=100
print('mode',mode(coverageList)[0][0],"ratio",ratio,'threshold:',threshold)

with open('highCoverage.bed','w') as fout:
    for contigName in Coverage:
        flag=False
        coverage=Coverage[contigName]
        contigLen=Length[contigName]
        if contigLen<=2*trimLength:
            l=contigLen
        else:
            l=trimLength
        for i in range(l):
            if coverage[i]>threshold and not flag:
                start = i+1
                flag=True
            elif coverage[i]<=threshold and flag:
                end = i
                flag=False
                fout.write('%s\t%d\t%d\n'%(contigName,start,end))
        if flag:
            end = l
            fout.write('%s\t%d\t%d\n'%(contigName,start,end))
            flag=False
        if contigLen>2*trimLength:
            for i in range(trimLength):
                if coverage[i+trimLength]>threshold and not flag:
                    start = contigLen-trimLength+i+1
                    flag = True
                elif coverage[i+trimLength]<=threshold and flag:
                    end = contigLen-trimLength+i
                    flag=False
                    fout.write('%s\t%d\t%d\n'%(contigName,start,end))
            if flag:
                end = contigLen
                fout.write('%s\t%d\t%d\n'%(contigName,start,end))
                flag=False

# with open('PrimCoverage.txt','w') as fout:
#     for contig in Coverage:
#         fout.write('>%s\n'%contig)
#         line=int(len(Coverage[contig])/100)
#         for k in range(line):
#             for d in Coverage[contig][100*k:100*(k+1)]:
#                 fout.write('%d\t'%d)
#             fout.write('\n')
#         if len(Coverage[contig])>line*100:
#             for d in Coverage[contig][100*line:]:
#                 fout.write('%d\t'%d)
#             fout.write('\n')