# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

import os
import argparse
import numpy as np
from multiprocessing import Pool, Manager
import networkx as nx
import MultiAlignment_func_python
import scipy.sparse as sp
import collections
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align
import math

def importContigs():
    contigSeq={}
    contigInd=[]
    for seq_record in SeqIO.parse('initial_contig.fa', 'fasta'):
        contigSeq[seq_record.id] = str(seq_record.seq).upper()
        contigInd.append(seq_record.id)
    return(contigSeq, contigInd)

def importEvidence():
    scaffold_dic={}
    with open('evidence_estimate_sd.txt') as fin:
        for line in fin:
            record = line.strip().split()
            if len(record) == 3:
                scaffold_num=int(record[0])
                scaffold_dic[scaffold_num]=[{'contig':'start', 'gap':0, 'std':0, 'reads':[]}]
            elif len(record) == 6:
                # scaffold_dic[scaffold_num].append({'name': 'contig'+str(record[0]), 'orientation':record[1], 'gap':int(record[4]), 'lower':int(record[5])-100, 'upper':int(record[6])+100})
                contigName=contigInd[int(record[0])-1]
                orientation=record[1]
                if orientation == 'r':
                    contigSeq[contigName]=str(Seq(contigSeq[contigName]).reverse_complement())
                endSeqL=contigSeq[contigName][-100:]
                reads=[{'name':'endSeqL', 'seq':endSeqL, 'phred':[38]*len(endSeqL), 'location':0, 'std':0, 'position':'L'}]
                scaffold_dic[scaffold_num].append({'contig':contigName, 'orientation':orientation, 'gap':int(record[4]), 'std':int(record[5]), 'reads':reads})
                endSeqR=contigSeq[contigName][:100]
                scaffold_dic[scaffold_num][-2]['reads'].append({'name':'endSeqR', 'seq':endSeqR, 'phred':[38]*len(endSeqR), 'location':scaffold_dic[scaffold_num][-2]['gap']+len(endSeqR), 'std':0, 'position':'R'})
    return(scaffold_dic)

def importReads(scaffold_dic, minPhred=3):
    for filename in os.listdir('FILEs'):
        side=filename[-7]
        contigName=filename[:-8]
        for scaffold_num in scaffold_dic:
            for i in range(len(scaffold_dic[scaffold_num])):
                if scaffold_dic[scaffold_num][i]['contig']==contigName:
                    orientation = scaffold_dic[scaffold_num][i]['orientation']
                    for seq_record in SeqIO.parse('FILEs/%s'%filename, 'fastq'):
                        seq_id=seq_record.id.split(':')
                        name=seq_id[-3]
                        location=int(seq_id[-2])
                        std=int(seq_id[-1])
                        seq=str(seq_record.seq).upper()
                        phred = seq_record.letter_annotations["phred_quality"]
                        if orientation == 'r':
                            side = ('R' if side == 'L' else 'L')
                            seq = str(Seq(seq).reverse_complement())
                            phred = phred[::-1]
                        while len(phred)>0:
                            if phred[0] <= minPhred:
                                phred = phred[1:]
                                seq = seq[1:]
                                if side == 'L':
                                    location-=1
                            else:
                                break
                        while len(phred)>0:
                            if phred[-1] <= minPhred:
                                phred = phred[:-1]
                                seq = seq[:-1]
                                if side == 'R':
                                    location-=1
                            else:
                                break
                        if len(seq)>=30:
                            expected_errors = 0
                            for q in phred:
                                expected_errors+=10**(-q/10)
                            if expected_errors*100/len(seq) <= qc:
                                if side=='R':
                                    j=i
                                    while (j<len(scaffold_dic[scaffold_num])):
                                        if j > i:
                                            location-=scaffold_dic[scaffold_num][j-1]['gap']
                                            location-=len(contigSeq[scaffold_dic[scaffold_num][j]['contig']])
                                            std=(std**2+scaffold_dic[scaffold_num][j-1]['std']**2)**0.5
                                        if j==len(scaffold_dic[scaffold_num])-1:
                                            nextGap=float('inf')
                                        else:
                                            nextGap=scaffold_dic[scaffold_num][j]['gap']
                                        if location > nextGap+len(seq)+3*(std**2+scaffold_dic[scaffold_num][j]['std']**2)**0.5:
                                            j+=1
                                        elif location >= -3*std:
                                            scaffold_dic[scaffold_num][j]['reads'].append({'name':name, 'seq':seq, 'phred':phred, 'location':location, 'std':round(std), 'position':'L'})
                                            j+=1
                                        else:
                                            break
                                if side=='L':
                                    j=i
                                    while (j>=1):
                                        if j < i:
                                            location-=scaffold_dic[scaffold_num][j]['gap']
                                            location-=len(contigSeq[scaffold_dic[scaffold_num][j]['contig']])
                                            std=(std**2+scaffold_dic[scaffold_num][j]['std']**2)**0.5
                                        if j==1:
                                            nextGap=float('inf')
                                        else:
                                            nextGap=scaffold_dic[scaffold_num][j-1]['gap']
                                        if location > nextGap+len(seq)+3*(std**2+scaffold_dic[scaffold_num][j-1]['std']**2)**0.5:
                                            j-=1
                                        elif location >= -3*std:
                                            scaffold_dic[scaffold_num][j-1]['reads'].append({'name':name, 'seq':seq, 'phred':phred, 'location':scaffold_dic[scaffold_num][j-1]['gap']-location+len(seq), 'std':round(std), 'position':'R'})
                                            j-=1
                                        else:
                                            break

def runPhrap(gapDic,num,reads,gapSize,gapStd,seqL,seqR):
    reads.sort(key = lambda x:x['location'])
    readLen=sum([len(read['seq']) for read in reads])/len(reads)
    readNum=max(gapSize+3*gapStd+200,200)*rc/readLen
    step=max(round(len(reads)/min(readNum,rn)),1)
    reads=reads[0::step]

    with open('phraps/reads%d.fasta'%num,'w') as fout:
        for r in reads:
            fout.write('>%s\n'%r['name'])
            fout.write('%s\n'%r['seq'])
    os.system('cd phraps; ../phrap -penalty -20 -gap_init -10 -gap_ext -5 -minmatch 10 -bandwidth 0 -minscore 20 -forcelevel 0 -bypasslevel 0 -revise_greedy -force_high -node_space 2 -confirm_penalty -10 -confirm_score 30 reads%d.fasta -ace > /dev/null 2>&1'%num)

    contigs=[]
    for seq_record in SeqIO.parse('phraps/reads%d.fasta.contigs'%num, 'fasta'):
        contigs.append(str(seq_record.seq).upper())
        contigs.append(str(seq_record.seq.reverse_complement()).upper())
    # os.system('cd phraps; rm reads%d.fasta*'%num)

    endLenL=min(100,len(seqL))
    endLenR=min(100,len(seqR))

    contigL=[]
    contigs2=[]
    for contig in contigs:
        alignments=aligner.align(contig, seqL)
        try:
            alignment=alignments[0]
        except IndexError as e:
            continue
        cigar_ref = alignment.aligned[0]
        cigar_query = alignment.aligned[1]
        hoL=min(cigar_ref[0][0], cigar_query[0][0])
        hoR=min(len(seqL)-cigar_query[-1][1], len(contig)-cigar_ref[-1][1])
        if alignment.score>=20 and hoL<=5 and (len(seqL)-cigar_query[-1][1])<=5<=(len(contig)-cigar_ref[-1][1]):
            contigL.append(seqL[len(seqL)-endLenL:]+contig[cigar_ref[-1][1]+len(seqL)-cigar_query[-1][1]:])
        else:
            contigs2.append(contig)

    contigA=[]
    contigR=[]
    for contig in contigL:
        if contig=='':
            continue
        alignments=aligner.align(contig, seqR)
        try:
            alignment=alignments[0]
        except IndexError as e:
            continue
        cigar_ref = alignment.aligned[0]
        cigar_query = alignment.aligned[1]
        hoL=min(cigar_ref[0][0], cigar_query[0][0])
        hoR=min(len(seqR)-cigar_query[-1][1], len(contig)-cigar_ref[-1][1])
        if alignment.score>=20 and hoR<=5 and cigar_query[0][0]<=5<=cigar_ref[0][0]:
            if contig[:cigar_ref[0][0]-cigar_query[0][0]] != '':
                contigA.append(contig[:cigar_ref[0][0]-cigar_query[0][0]]+seqR[:endLenR])

    for contig in contigs2:
        alignments=aligner.align(contig, seqR)
        try:
            alignment=alignments[0]
        except IndexError as e:
            continue
        cigar_ref = alignment.aligned[0]
        cigar_query = alignment.aligned[1]
        hoL=min(cigar_ref[0][0], cigar_query[0][0])
        hoR=min(len(seqR)-cigar_query[-1][1], len(contig)-cigar_ref[-1][1])
        if alignment.score>=20 and hoR<=5 and cigar_query[0][0]<=5<=cigar_ref[0][0]:
            contigR.append(contig[:cigar_ref[0][0]-cigar_query[0][0]]+seqR[:endLenR])
     
    if len(contigA)>0:
        diff=float('inf')
        for contig in contigA:
            if abs(len(contig)-gapSize-200)<diff:
                diff=abs(len(contig)-gapSize-200)
                gapSeq=contig
        gapDic[num]={'A':gapSeq, 'endLenL':endLenL, 'endLenR':endLenR}
    else:
        if len(contigL)>0:
            gapSeqL=contigL[np.array([len(c) for c in contigL]).argmax()]
        else:
            gapSeqL=''
        if len(contigR)>0:
            gapSeqR=contigR[np.array([len(c) for c in contigR]).argmax()]
        else:
            gapSeqR=''
        gapDic[num]={'L':gapSeqL[endLenL:], 'R':gapSeqR[:len(gapSeqR)-endLenR]}


def generateScaffold(scaffold_dic, scaffold_num, alignment_score = 1, mismatch = 2, gap_cost = 5, minscore=20):
    scafSeq = ''
    gapNum = 0 
    closedNum = 0
    old_contigLength_list = []
    new_contigLength_list = []

    for k in range(len(scaffold_dic[scaffold_num])):
        Gap = scaffold_dic[scaffold_num][k]
        if k==0:
            Contig = scaffold_dic[scaffold_num][k+1]
            currentContig = contigSeq[Contig['contig']]
        elif k==len(scaffold_dic[scaffold_num])-1:
            scafSeq+=currentContig
            new_contigLength_list.append(len(currentContig))
            old_contigLength_list.append(len(contigSeq[Gap['contig']]))
        else:
            gapNum+=1
            Contig = scaffold_dic[scaffold_num][k+1]
            old_contigLength_list.append(len(contigSeq[Gap['contig']]))
            gapSize = Gap['gap']
            # gapStd = Gap['std']
            gapStd = 100
            if 'A' in Gap['gapSequence']:
                closedNum+=1
                endLenL = Gap['gapSequence']['endLenL']
                endLenR = Gap['gapSequence']['endLenR']
                currentContig = currentContig[:-endLenL]+Gap['gapSequence']['A']+contigSeq[Contig['contig']][endLenR:]
            else:
                currentContig = currentContig+Gap['gapSequence']['L']
                nextContig = Gap['gapSequence']['R']+contigSeq[Contig['contig']]
                overlapSize = -(gapSize-len(Gap['gapSequence']['L'])-len(Gap['gapSequence']['R']))
                truncate = overlapSize+3*gapStd+100
                if truncate > 0 and truncate <= mT:
                    seq_L = currentContig[-truncate:]
                    seq_R = nextContig[:truncate]
                    alignments = MultiAlignment_func_python.MultiOptimal_merge_alignment(seq_L,seq_R,alignment_score,mismatch,gap_cost,3)
                    # minHangingOutLen=float('inf')
                    score=0
                    overlapDiff=float('inf')
                    for aligned in alignments:
                        # hangingOutLen=min(aligned[2]-1, aligned[4]-1)
                        # hangingOutLen+=min(len(seq_L)-aligned[3], len(seq_R)-aligned[5])
                        hangingOutL=min(aligned[2]-1, aligned[4]-1)
                        hangingOutR=min(len(seq_L)-aligned[3], len(seq_R)-aligned[5])
                        overlapLen = aligned[4]-1+len(aligned[1])+len(seq_L)-aligned[3]
                        if aligned[0]>=minscore and max(hangingOutL,hangingOutR)<=HO and abs(overlapLen-overlapSize)<overlapDiff:
                            score=aligned[0]
                            overlap_seq=aligned[1]
                            start_i=aligned[2]
                            end_i=aligned[3]
                            start_j=aligned[4]
                            end_j=aligned[5]
                            overlapDiff=abs(overlapLen-overlapSize)
                            # minHangingOutLen=hangingOutLen
                    if score >= minscore:
                        closedNum+=1
                        currentContig = currentContig[:-truncate]+seq_L[:start_i-1]+overlap_seq+seq_R[end_j:]+nextContig[truncate:]
                        continue
                if gapSize+3*gapStd<0:
                    currentContig = currentContig[:len(currentContig)-len(Gap['gapSequence']['L'])]
                    scafSeq += currentContig
                    scafSeq += 'N'
                    new_contigLength_list.append(len(currentContig))
                    currentContig = nextContig[len(Gap['gapSequence']['R']):]
                else:
                    new_contigLength_list.append(len(currentContig))
                    scafSeq+=currentContig
                    currentContig = nextContig
                    newGapSize = gapSize-len(Gap['gapSequence']['L'])-len(Gap['gapSequence']['R'])
                    if newGapSize<=0:
                        newGapSize = 1
                    scafSeq+='N'*newGapSize

    scaffold_dic[scaffold_num]={'scafSeq':scafSeq, 'gapNum':gapNum, 'closedNum':closedNum, 'old_contigLength_list':old_contigLength_list,'new_contigLength_list':new_contigLength_list}


def writeSequence(f, sequence):
    length = len(sequence)
    line=int(length/100)
    for k in range(line):
        f.write(sequence[100*k:100*(k+1)]+'\n')
    if length>line*100:
        f.write(sequence[100*line:]+'\n')


def computeN50(length_list):
    length_list.sort()
    total_length=sum(length_list)
    temp=0
    for length in length_list:
        temp += length
        if temp >= (total_length/2):
            N50 = length
            break
    return(N50)


def writeResults(scaffold_dic):
    scafFile = open(output, 'w')
    resultFile = open('statistics.txt', 'w')
    totalGapNum = 0
    totalClosedNum = 0
    old_contigLength_list = []
    new_contigLength_list = []
    scaffoldLength_list = []

    for scaffold_num in scaffold_dic:
        totalGapNum+=scaffold_dic[scaffold_num]['gapNum']
        totalClosedNum+=scaffold_dic[scaffold_num]['closedNum']
        old_contigLength_list+=scaffold_dic[scaffold_num]['old_contigLength_list']
        new_contigLength_list+=scaffold_dic[scaffold_num]['new_contigLength_list']

        scafSeq=scaffold_dic[scaffold_num]['scafSeq']
        scaffoldLength_list.append(len(scafSeq))
        scafFile.write('>scaffold%d\n'%scaffold_num)
        writeSequence(scafFile, scafSeq)

    oldContigN50 = computeN50(old_contigLength_list)
    newContigN50 = computeN50(new_contigLength_list)
    oldTotalContigLen = sum(old_contigLength_list)
    newTotalContigLen = sum(new_contigLength_list)
    scaffoldN50 = computeN50(scaffoldLength_list)

    resultFile.write("Total gap number: "+"{:,}\n".format(totalGapNum))
    resultFile.write("Closed gap number: "+"{:,}\n".format(totalClosedNum))
    resultFile.write("Total contig length from "+"{:,}".format(oldTotalContigLen)+" to "+"{:,}\n".format(newTotalContigLen))
    resultFile.write("Contig N50 from "+"{:,}".format(oldContigN50)+" to "+"{:,}\n".format(newContigN50))
    resultFile.write("\nScaffold N50: "+"{:,}\n".format(scaffoldN50))

    scafFile.close()
    resultFile.close()


parser = argparse.ArgumentParser(description="Local assembly using Phrap.") 

parser.add_argument('-o', default='output_genome.fasta', help='Output file saving gap-closed genome. [output_genome.fasta]')
parser.add_argument('-t', type=int, default=1, help='Number of threads. [1]')
parser.add_argument('-qc', type=float, default=100, help='Maximum expected erroneous bases in the read used for local assembly. [100]')
parser.add_argument('-rc', type=int, default=100, help='Maximum coverage of reads used for local assembly. [100]')
parser.add_argument('-rn', type=int, default=2000, help='Maximum number of reads used for local assembly. [2000]')
parser.add_argument('-mT', type=int, default=1000, help='Maximum truncated length for alignment merging adjacent contigs. [1000]')
parser.add_argument('-mA', type=int, default=1, help='Matching score in alignment merging adjacent contigs. [1]')
parser.add_argument('-mM', type=int, default=2, help='Mismatch penalty in alignment merging adjacent contigs. [2]')
parser.add_argument('-mG', type=int, default=3, help='Gap cost in alignment merging adjacent contigs. [3]')
parser.add_argument('-mS', type=int, default=20, help='Minimum alignment score to merge adjacent contigs. [20]')
parser.add_argument('-HO', type=int, default=5, help='Maximum admissible hanging-out length in alignment merging adjacent contigs. [5]')

args = parser.parse_args()

thread = args.t
output = args.o
qc = args.qc
rc = args.rc
rn = args.rn
mT = args.mT
mA = args.mA
mM = args.mM
mG = args.mG
mS = args.mS
HO = args.HO

aligner = Align.PairwiseAligner(mode = 'local', match_score=1, mismatch_score=-2, open_gap_score=-5, extend_gap_score=-1)

contigSeq, contigInd=importContigs()

scaffold_dic=importEvidence()

importReads(scaffold_dic)

os.system('rm -rf phraps; mkdir phraps')

m=Manager()
gapDic=m.dict()

pool=Pool(thread)
num=0
for scaffold_num in scaffold_dic:
    for k in range(1,len(scaffold_dic[scaffold_num])-1):
        num+=1
        reads=scaffold_dic[scaffold_num][k]['reads']
        gapSize=scaffold_dic[scaffold_num][k]['gap']
        gapStd=scaffold_dic[scaffold_num][k]['std']
        seqL=contigSeq[scaffold_dic[scaffold_num][k]['contig']][-1000:]
        seqR=contigSeq[scaffold_dic[scaffold_num][k+1]['contig']][:1000]
        pool.apply_async(func=runPhrap,args=(gapDic,num,reads,gapSize,gapStd,seqL,seqR))
pool.close()
pool.join()

num=0
for scaffold_num in scaffold_dic:
    for k in range(1,len(scaffold_dic[scaffold_num])-1):
        num+=1
        Gap=scaffold_dic[scaffold_num][k]
        scaffold_dic[scaffold_num][k]={'contig':Gap['contig'], 'gap':Gap['gap'], 'std':Gap['std'], 'gapSequence':gapDic[num]}

scaffold_dic=m.dict(scaffold_dic)

pool=Pool(thread)
for scaffold_num in scaffold_dic:
    pool.apply_async(func=generateScaffold,args=(scaffold_dic, scaffold_num, mA, mM, mG, mS))
pool.close()
pool.join()

writeResults(scaffold_dic)

