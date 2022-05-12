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
import re


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
                endSeqL=contigSeq[contigName][-Length:]
                reads=[{'name':'endSeqL', 'seq':endSeqL, 'phred':[38]*len(endSeqL), 'location':0, 'std':0, 'position':'L'}]
                scaffold_dic[scaffold_num].append({'contig':contigName, 'orientation':orientation, 'gap':int(record[4]), 'std':int(record[5]), 'reads':reads})
                endSeqR=contigSeq[contigName][:Length]
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
                        if len(seq)>0:
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


def setupRegressionModel(reads, gapStd, scale, alignment_score = 1, mismatch = 20, gap_cost = 30, minscore=20, hangingOut=0):

    row=0
    indptr=[0]
    indices=[]
    data=[]
    response=[]
    weight=[]

    size=len(reads)
    for i in range(size-1):
        for j in range(i+1,size):
            l1=len(reads[i]['seq'])
            l2=len(reads[j]['seq'])
            if reads[i]['position']==reads[j]['position']:
                std=int((reads[i]['std']**2+reads[j]['std']**2)**0.5)
                read_to_read_up=l2-minscore
                read_to_read_low=minscore-l1-std*scale
            else:
                std=int((reads[i]['std']**2+reads[j]['std']**2+gapStd**2)**0.5)
                read_to_read_up=l2-minscore+std*scale
                read_to_read_low=minscore-l1-std*scale
            if read_to_read_low<reads[j]['location']-reads[i]['location']<read_to_read_up:
                alignments = MultiAlignment_func_python.MultiOptimal_merge_alignment(reads[i]['seq'],reads[j]['seq'],alignment_score,mismatch,gap_cost,3)
                minDiff = float('inf')
                score=0
                for aligned in alignments:
                    diff=((l2-aligned[5])-(l1-aligned[3]))-(reads[j]['location']-reads[i]['location'])
                    if abs(diff)<minDiff and min(l2-aligned[5],l1-aligned[3])<=hangingOut and min(aligned[4]-1,aligned[2]-1)<=hangingOut:
                        score=aligned[0]
                        start_i=aligned[2]
                        end_i=aligned[3]
                        start_j=aligned[4]
                        end_j=aligned[5]
                        minDiff=abs(diff)
                if score>=minscore and minDiff<=std*scale:
                    row+=1
                    indices.append(j)
                    data.append(1)
                    indices.append(i)
                    data.append(-1)
                    indptr.append(indptr[-1]+2)
                    response+=[(l2-end_j)-(l1-end_i)]
                    weight.append(1/std)

    design=sp.csr_matrix((data, indices, indptr), shape=(row, size))
    response=np.matrix(response).T
    if args.w:
        weight=np.array(weight)
    else:
        weight=np.ones(len(weight))

    return(design, response, weight)


#对csr_matrix格式的矩阵删除行
def deleteRowsCsr(mat, indices):
    indices = list(indices)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask]


#####Huber M-estimate of reads coordinates 
def IRLS(X,Y,W,reads,thr1=2,thr2=3):
    if Y.shape[0]==0:
        return([], [])

    X0=X
    Y0=Y

    W=sp.diags(W)
    X=W.dot(X)
    Y=W.dot(Y) #initial weight multiplies each alignment

    t=X.T
    A=t.dot(X)
    y=sp.csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = sp.linalg.lgmres(A, b, atol=1e-05) #compute ordinary least squares estiamte as initial estimate
    residual=abs((X.dot(sp.csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    threshold=sp.csr_matrix(np.ones(len(residual))*thr1).dot(W).toarray()[0] #threshold thr1 of Huber's weight function
    old_estimate=estimate
    n=0
    
    while n<100:
        index=np.where(residual>threshold)[0]
        reweight=np.ones(len(residual))
        reweight[index]=threshold[index]/residual[index] # update weights for alignments with residuals greater than thr1
        reweight=sp.diags(reweight)
        t=X.T
        A=t.dot(reweight).dot(X)
        y=sp.csr_matrix(Y)
        b=t.dot(reweight).dot(y).todense()
        estimate, exitCode = sp.linalg.lgmres(A, b, atol=1e-05) #compute weighted least squares estimate
        residual=abs((X.dot(sp.csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        diff=max(abs(estimate-old_estimate))
        if diff<=2: # convergence condition of estimates
            break
        else:
            old_estimate=estimate
            n+=1

#remove alignments with residuals greater than thr2
    residual0=abs((X0.dot(sp.csr_matrix(estimate).T)-sp.csr_matrix(Y0)).todense()).T.getA()[0] 
    G = nx.Graph()
    for i in range(X.shape[0]):
        if residual0[i]<=thr2:
            n1=X.indices[2*i]
            n2=X.indices[2*i+1]
            G.add_edge(n1,n2)

#divide reads into separate connected components
    reads_list=[]
    estimates_list=[]
    for c in nx.connected_components(G):
        sub_index=list(G.subgraph(c).nodes)
        sub_estimates=[estimate[i] for i in sub_index]
        sub_reads=[reads[i] for i in sub_index]
        reads_name=[read['name'] for read in sub_reads]
        if 'endSeqL' in reads_name and 'endSeqR' in reads_name:
            if sub_estimates[reads_name.index('endSeqL')] > sub_estimates[reads_name.index('endSeqR')]:
                continue
        if 'endSeqL' in reads_name:
            ind=reads_name.index('endSeqL')
            sub_reads=[sub_reads[i] for i in range(len(sub_estimates)) if sub_estimates[i]>=sub_estimates[ind]]
            reads_name=[read['name'] for read in sub_reads]
            sub_estimates=[sub_estimates[i] for i in range(len(sub_estimates)) if sub_estimates[i]>=sub_estimates[ind]]
        if 'endSeqR' in reads_name:
            ind=reads_name.index('endSeqR')
            sub_reads=[sub_reads[i] for i in range(len(sub_estimates)) if sub_estimates[i]<=sub_estimates[ind]]
            sub_estimates=[sub_estimates[i] for i in range(len(sub_estimates)) if sub_estimates[i]<=sub_estimates[ind]]
        estimates_list.append(list(map(int,np.round(sub_estimates))))
        reads_list.append(sub_reads)

    return(estimates_list, reads_list)


def generateConsensus(reads, estimates):
    reads_name=[read['name'] for read in reads]
    if 'endSeqL' not in reads_name and 'endSeqR' not in reads_name:
        return('')
    
    end=np.array([estimates[i]-len(reads[i]['seq']) for i in range(len(estimates))])
    length=max(estimates)-min(end)
    assembly=[dict() for i in range(length)]
    start_point=min(end)
    for i in range(len(estimates)):
        distance=end[i]-start_point
        for j in range(len(reads[i]['seq'])):
            base=reads[i]['seq'][j]
            if base in assembly[j+distance]:
                assembly[j+distance][base].append(i)
            else:
                assembly[j+distance][base]=[i]
    
    complement_reads=[]
    while(True):
        error_reads=[]
        essential_reads=[]
        for i in range(len(assembly)):
            statistics=assembly[i]
            count=dict()
            for k in statistics:
                count[k]=len(statistics[k])
            if len(statistics.keys())>1:
                major=max(count.values())
                for k in count:
                    if count[k]<major:
                        error_reads+=statistics[k]
            if sum(count.values())==1:
                essential_reads+=sum(statistics.values(),[])
        if error_reads==list():
            break 
        remove_candidate=collections.Counter(error_reads)
        remove_read=max(remove_candidate, key=remove_candidate.get)
        distance=end[remove_read]-start_point
        for j in range(len(reads[remove_read]['seq'])):
            base=reads[remove_read]['seq'][j]
            assembly[j+distance][base].remove(remove_read)
        if remove_read in essential_reads:
            complement_reads.append(remove_read)

    for i in complement_reads:
        distance=end[i]-start_point
        for j in range(len(reads[i]['seq'])):
            base=reads[i]['seq'][j]
            assembly[j+distance][base].append(i)

    for i in range(len(assembly)):
        statistics=assembly[i]
        count=dict()
        for k in statistics:
            count[k]=len(statistics[k])
        if count==dict():
            assembly[i]='N'
        else:
            base=max(count, key=count.get)
            assembly[i]=base

    consensus=''.join(assembly)
    position=''

    if 'endSeqR' in reads_name:
        position+='R'
        ind=reads_name.index('endSeqR')
        delta=reads[ind]['location']-estimates[ind]
    if 'endSeqL' in reads_name:
        position+='L'
        ind=reads_name.index('endSeqL')
        consensus=consensus[estimates[ind]-min(end)-len(reads[ind]['seq']):]
        delta=0-estimates[ind]
    estimates=[e+delta for e in estimates]

    if consensus != '' and position != '':
        return({'seq':consensus, 'position':position, 'reads':reads, 'estimates':estimates})
    else:
        return('')


def integrateConsensus(estimates_list, reads_list):
    gapConsensus={'L':'', 'R':''}
    for i in range(len(estimates_list)):
        estimates=estimates_list[i]
        reads=reads_list[i]
        if estimates != []:
            consensus=generateConsensus(reads, estimates)
            if consensus!='':
                if consensus['position']=='RL':
                    return({'RL':consensus['seq'], 'reads':consensus['reads'], 'estimates':consensus['estimates']})
                elif consensus['position']=='L':
                    gapConsensus['L']=consensus['seq']
                    gapConsensus['reads_L']=consensus['reads']
                    gapConsensus['estimates_L']=consensus['estimates']
                elif consensus['position']=='R':
                    gapConsensus['R']=consensus['seq']
                    gapConsensus['reads_R']=consensus['reads']
                    gapConsensus['estimates_R']=consensus['estimates']
    return(gapConsensus)


def updateConsensus(consensus, reads, estimates, endLenL, endLenR, trim, evaluation=False, adjust=20, minScore=20):

    updatedConsensus=[]
    for base in consensus:
        updatedConsensus.append({'':[],'A':[],'G':[],'C':[],'T':[],'insertion':{'':0}})

    end=np.array([estimates[i]-len(reads[i]['seq']) for i in range(len(estimates))])
    for k in range(len(reads)):
        if reads[k]['name']=='endSeqL' or reads[k]['name']=='endSeqR':
            continue
        query=reads[k]['seq']
        phred=reads[k]['phred']
        startPoint=max(estimates[k]-min(end)-len(query)-adjust, 0)
        endPoint=max(estimates[k]-min(end)+adjust, 0)
        ref=consensus[startPoint:endPoint]

        if ref == '':
            continue
        alignments = aligner.align(query, ref)
        alignment = alignments[0]
        if alignment.score < minScore:
            continue
        cigar_query = alignment.aligned[0]
        cigar_ref = alignment.aligned[1]
        for i in range(len(cigar_ref)):
            L=0
            R=cigar_ref[i][1]-cigar_ref[i][0]
            if i==0 and startPoint+cigar_ref[i][0]>=trim:
                L=trim
            if i==len(cigar_ref)-1 and startPoint+cigar_ref[i][1]<=len(consensus)-trim:
                R=R-trim
            for j in range(L,R):
                base=query[cigar_query[i][0]+j]
                if base in updatedConsensus[startPoint+cigar_ref[i][0]+j]:
                    updatedConsensus[startPoint+cigar_ref[i][0]+j][base].append(1-10**(-phred[cigar_query[i][0]+j]/10))
                if j>0:
                    updatedConsensus[startPoint+cigar_ref[i][0]+j]['insertion']['']+=1

        for i in range(len(cigar_ref)-1):
            left=cigar_ref[i][1]
            right=cigar_ref[i+1][0]
            if right-left > 0:
                for j in range(startPoint+left, startPoint+right):
                    updatedConsensus[j][''].append(1)
                    updatedConsensus[j]['insertion']['']+=1
                updatedConsensus[startPoint+right]['insertion']['']+=1
            else:
                insertion=query[cigar_query[i][1]:cigar_query[i+1][0]]
                if insertion in updatedConsensus[startPoint+left]['insertion']:
                    updatedConsensus[startPoint+left]['insertion'][insertion]+=1
                else:
                    updatedConsensus[startPoint+left]['insertion'][insertion]=1

    phred = ''
    for k in range(len(updatedConsensus)):
        insertion_dic=updatedConsensus[k].pop('insertion')
        insertion=max(insertion_dic, key=insertion_dic.get)

        A=G=C=T=D=0
        for prob in updatedConsensus[k]['A']:
            A+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        for prob in updatedConsensus[k]['G']:
            G+=math.log(0.99995*prob);A+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        for prob in updatedConsensus[k]['C']:
            C+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);A+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        for prob in updatedConsensus[k]['T']:
            T+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);A+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        for prob in updatedConsensus[k]['']:
            D+=math.log(1-0.00005);A+=math.log(0.00005);G+=math.log(0.00005);C+=math.log(0.00005);T+=math.log(0.00005)

        PA=1/(1+np.exp((T-A)*((T-A)<=708))+np.exp((G-A)*((G-A)<=708))+np.exp((C-A)*((C-A)<=708))+np.exp((D-A)*((D-A)<=708)))*((T-A)<=708)*((G-A)<=708)*((C-A)<=708)*((D-A)<=708)
        PG=1/(1+np.exp((A-G)*((A-G)<=708))+np.exp((T-G)*((T-G)<=708))+np.exp((C-G)*((C-G)<=708))+np.exp((D-G)*((D-G)<=708)))*((A-G)<=708)*((T-G)<=708)*((C-G)<=708)*((D-G)<=708)
        PC=1/(1+np.exp((A-C)*((A-C)<=708))+np.exp((G-C)*((G-C)<=708))+np.exp((T-C)*((T-C)<=708))+np.exp((D-C)*((D-C)<=708)))*((A-C)<=708)*((G-C)<=708)*((T-C)<=708)*((D-C)<=708)
        PT=1/(1+np.exp((A-T)*((A-T)<=708))+np.exp((G-T)*((G-T)<=708))+np.exp((C-T)*((C-T)<=708))+np.exp((D-T)*((D-T)<=708)))*((A-T)<=708)*((G-T)<=708)*((C-T)<=708)*((D-T)<=708)
        PD=1/(1+np.exp((A-D)*((A-D)<=708))+np.exp((G-D)*((G-D)<=708))+np.exp((C-D)*((C-D)<=708))+np.exp((T-D)*((T-D)<=708)))*((A-D)<=708)*((G-D)<=708)*((C-D)<=708)*((T-D)<=708)

        if PA==max(PA,PG,PC,PT,PD):
            site='A'
        elif PG==max(PA,PG,PC,PT,PD):
            site='G'
        elif PC==max(PA,PG,PC,PT,PD):
            site='C'
        elif PT==max(PA,PG,PC,PT,PD):
            site='T'
        elif PD==max(PA,PG,PC,PT,PD):
            site=''

        if k < endLenL or k >= len(updatedConsensus)-endLenR:
            site = consensus[k]
        if site != '':
            phred += chr(round(-10*math.log10(max(1-max(PA,PG,PC,PT,PD), 0.00006)))+33)

        if A==G==C==T==D==0:
            site=consensus[k]
            insertion=''
        if k < endLenL or k > len(updatedConsensus)-endLenR:
            insertion = ''
        if evaluation:
            insertion = ''

        updatedConsensus[k]=insertion+site

    updatedConsensus=('').join(updatedConsensus)

    return(updatedConsensus, phred)


def closeGap(gapDic,num,reads,gapSize,gapStd,scale):
    reads.sort(key = lambda x:x['location'])
    readLen=sum([len(read['seq']) for read in reads])/len(reads)
    readNum=max(gapSize+3*gapStd+200,200)*rc/readLen
    step=max(int(len(reads)/readNum),1)
    reads=reads[0::step]
    design, response, weight=setupRegressionModel(reads, gapStd, scale, ma, mm, gc, ms, ho)
    estimates_list, reads_list=IRLS(design, response, weight, reads, thr1, thr2)
    gapConsensus=integrateConsensus(estimates_list, reads_list)

    if 'RL' in gapConsensus:
        consensus=gapConsensus['RL']
        reads=gapConsensus['reads']
        estimates=gapConsensus['estimates']
        reads_name=[read['name'] for read in reads]
        endLenL=len(reads[reads_name.index('endSeqL')]['seq'])
        endLenR=len(reads[reads_name.index('endSeqR')]['seq'])
        gapSeq, gapPhred=updateConsensus(consensus, reads, estimates, endLenL, endLenR, 15, False, 20)
        gapSeq, gapPhred=updateConsensus(gapSeq, reads, estimates, endLenL, endLenR, 2, True, 30)
        gapSequence={'RL':gapSeq, 'phred':gapPhred, 'endLenL':endLenL, 'endLenR':endLenR}
    else:
        gapSequence={}
        if gapConsensus['L'] != '':
            consensus=gapConsensus['L']
            reads=gapConsensus['reads_L']
            estimates=gapConsensus['estimates_L']
            reads_name=[read['name'] for read in reads]
            endLenL=len(reads[reads_name.index('endSeqL')]['seq'])
            gapSeq, gapPhred=updateConsensus(consensus, reads, estimates, endLenL, 0, 15, False, 20)
            gapSeq, gapPhred=updateConsensus(gapSeq, reads, estimates, endLenL, 0, 2, True, 30)
            gapSequence['L']=gapSeq[endLenL:]
            gapSequence['phred_L']=gapPhred[endLenL:]
        else:
            gapSequence['L']=''
            gapSequence['phred_L']=''
        if gapConsensus['R'] != '':
            consensus=gapConsensus['R']
            reads=gapConsensus['reads_R']
            estimates=gapConsensus['estimates_R']
            reads_name=[read['name'] for read in reads]
            endLenR=len(reads[reads_name.index('endSeqR')]['seq'])
            gapSeq, gapPhred=updateConsensus(consensus, reads, estimates, 0, endLenR, 15, False, 20)
            gapSeq, gapPhred=updateConsensus(gapSeq, reads, estimates, 0, endLenR, 2, True, 30)
            gapSequence['R']=gapSeq[:-endLenR]
            gapSequence['phred_R']=gapPhred[:-endLenR]
        else:
            gapSequence['R']=''
            gapSequence['phred_R']=''

    gapDic[num]=gapSequence


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


def writeResults(scaffold_dic, alignment_score = 1, mismatch = 2, gap_cost = 5, minscore=20):
    scafFile = open(output, 'w')
    gapFile = open('gapSequence.fastq', 'w')
    fillFile = open('evidence.fill', 'w')
    resultFile = open('statistics.txt', 'w')
    totalGapNum = 0
    totalClosedNum = 0
    MergedNum = 0
    filledNum = 0
    old_contigLength_list = []
    new_contigLength_list = []
    scaffoldLength_list = []
    for scaffold_num in scaffold_dic:
        scafSeq = ''
        scafFile.write('>scaffold%d\n'%scaffold_num)
        fillFile.write('~~~~~\nscaffold%d\n'%scaffold_num)
        for k in range(len(scaffold_dic[scaffold_num])):
            Gap = scaffold_dic[scaffold_num][k]
            if k==0:
                Contig = scaffold_dic[scaffold_num][k+1]
                fillFile.write('^'.ljust(20)+Contig['contig'].ljust(20)+'0'.ljust(10)+str(len(Gap['gapSequence']['R'])).ljust(10)+'0'.ljust(10)+'0\n')
                currentContig = Gap['gapSequence']['R']+contigSeq[Contig['contig']]
            elif k==len(scaffold_dic[scaffold_num])-1:
                fillFile.write(Gap['contig'].ljust(20)+'$'.ljust(20)+str(len(Gap['gapSequence']['L'])).ljust(10)+'0'.ljust(10)+'0'.ljust(10)+'0\n')
                currentContig += Gap['gapSequence']['L']
                scafSeq+=currentContig
                new_contigLength_list.append(len(currentContig))
                old_contigLength_list.append(len(contigSeq[Gap['contig']]))
            else:
                totalGapNum+=1
                Contig = scaffold_dic[scaffold_num][k+1]
                old_contigLength_list.append(len(contigSeq[Gap['contig']]))
                gapSize = Gap['gap']
                # gapStd = Gap['std']
                gapStd = 100
                if 'RL' in Gap['gapSequence']:
                    endLenL = Gap['gapSequence']['endLenL']
                    endLenR = Gap['gapSequence']['endLenR']
                    currentContig = currentContig[:-endLenL]+Gap['gapSequence']['RL']+contigSeq[Contig['contig']][endLenR:]
                    totalClosedNum+=1
                    filledLen=len(Gap['gapSequence']['RL'])-endLenL-endLenR
                    if filledLen >= 0:
                        filledNum+=1
                        if filledLen > 0:
                            gapFile.write('@%s_%s_filled\n'%(Gap['contig'], Contig['contig']))
                            writeSequence(gapFile, Gap['gapSequence']['RL'][endLenL:-endLenR])
                            gapFile.write('+\n')
                            writeSequence(gapFile, Gap['gapSequence']['phred'][endLenL:-endLenR])
                    else:
                        MergedNum+=1
                    fillFile.write(Gap['contig'].ljust(20)+Contig['contig'].ljust(20)+str(filledLen).ljust(10)+'0'.ljust(10)+'1'.ljust(10)+'0\n')
                else:
                    currentContig += Gap['gapSequence']['L']
                    nextContig = Gap['gapSequence']['R']+contigSeq[Contig['contig']]
                    overlapSize = -(gapSize-len(Gap['gapSequence']['L'])-len(Gap['gapSequence']['R']))
                    truncate = overlapSize+3*gapStd+100
                    if truncate > 0:
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
                            if max(hangingOutL,hangingOutR)<=HO and abs(overlapLen-overlapSize)<overlapDiff:
                                score=aligned[0]
                                overlap_seq=aligned[1]
                                start_i=aligned[2]
                                end_i=aligned[3]
                                start_j=aligned[4]
                                end_j=aligned[5]
                                overlapDiff=abs(overlapLen-overlapSize)
                                # minHangingOutLen=hangingOutLen
                        if score >= minscore:
                            lenL=len(currentContig)-len(Gap['gapSequence']['L'])
                            lenR=len(nextContig)-len(Gap['gapSequence']['R'])
                            currentContig = currentContig[:-truncate]+seq_L[:start_i-1]+overlap_seq+seq_R[end_j:]+nextContig[truncate:]
                            gapSeq=currentContig[lenL:-lenR]
                            if (len(currentContig)-lenL-lenR)>=0:
                                filledNum+=1
                                if len(gapSeq)>0:
                                    lenGapL=min(len(gapSeq),len(Gap['gapSequence']['L'])-(len(seq_L)-end_i))
                                    lenGapR=len(Gap['gapSequence']['phred_R'])-(len(gapSeq)-lenGapL)
                                    gapPhred=Gap['gapSequence']['phred_L'][:lenGapL]+Gap['gapSequence']['phred_R'][lenGapR:]
                                    gapFile.write('@%s_%s_filled\n'%(Gap['contig'], Contig['contig']))
                                    writeSequence(gapFile, gapSeq)
                                    gapFile.write('+\n')
                                    writeSequence(gapFile, gapPhred)
                            else:
                                MergedNum+=1
                            fillFile.write(Gap['contig'].ljust(20)+Contig['contig'].ljust(20)+str(len(currentContig)-lenL-lenR).ljust(10)+'0'.ljust(10)+'1'.ljust(10)+'0\n')
                            totalClosedNum+=1
                            continue
                    if gapSize+3*gapStd<0:
                        currentContig = currentContig[:len(currentContig)-len(Gap['gapSequence']['L'])]
                        scafSeq += currentContig
                        scafSeq += 'N'
                        new_contigLength_list.append(len(currentContig))
                        currentContig = nextContig[len(Gap['gapSequence']['R']):]
                        fillFile.write(Gap['contig'].ljust(20)+Contig['contig'].ljust(20)+'0'.ljust(10)+'0'.ljust(10)+'0'.ljust(10)+str(gapSize)+'\n')
                    else:
                        if len(Gap['gapSequence']['L'])>0:
                            gapFile.write('@%s_%s_left\n'%(Gap['contig'], Contig['contig']))
                            writeSequence(gapFile, Gap['gapSequence']['L'])
                            gapFile.write('+\n')
                            writeSequence(gapFile, Gap['gapSequence']['phred_L'])
                        new_contigLength_list.append(len(currentContig))
                        scafSeq+=currentContig
                        currentContig = nextContig
                        if len(Gap['gapSequence']['R'])>0:
                            gapFile.write('@%s_%s_right\n'%(Gap['contig'], Contig['contig']))
                            writeSequence(gapFile, Gap['gapSequence']['R'])
                            gapFile.write('+\n')
                            writeSequence(gapFile, Gap['gapSequence']['phred_R'])
                        newGapSize = gapSize-len(Gap['gapSequence']['L'])-len(Gap['gapSequence']['R'])
                        fillFile.write(Gap['contig'].ljust(20)+Contig['contig'].ljust(20)+str(len(Gap['gapSequence']['L'])).ljust(10)+str(len(Gap['gapSequence']['R'])).ljust(10)+'0'.ljust(10)+'%d\n'%newGapSize)
                        if newGapSize<=0:
                            newGapSize = 1
                        scafSeq+='N'*newGapSize
        writeSequence(scafFile, scafSeq)
        scaffoldLength_list.append(len(scafSeq))

    oldContigN50 = computeN50(old_contigLength_list)
    newContigN50 = computeN50(new_contigLength_list)
    oldTotalContigLen = sum(old_contigLength_list)
    newTotalContigLen = sum(new_contigLength_list)
    scaffoldN50 = computeN50(scaffoldLength_list)

    resultFile.write("Total gap number: "+"{:,}\n".format(totalGapNum))
    resultFile.write("Closed gap number: "+"{:,}\n".format(totalClosedNum))
    resultFile.write("\tMerged gap number: "+"{:,}\n".format(MergedNum))
    resultFile.write("\tFilled gap number: "+"{:,}\n\n".format(filledNum))
    resultFile.write("Total contig length from "+"{:,}".format(oldTotalContigLen)+" to "+"{:,}\n".format(newTotalContigLen))
    resultFile.write("Contig N50 from "+"{:,}".format(oldContigN50)+" to "+"{:,}\n".format(newContigN50))
    resultFile.write("\nScaffold N50: "+"{:,}\n".format(scaffoldN50))

    scafFile.close()
    gapFile.close()
    fillFile.close()
    resultFile.close()



#读入参数：库长标准差、比对参数、线程数、WLS的残差取舍参数、consensus的base占优比例参数
parser = argparse.ArgumentParser(description="Do local OLC by iterative regression.") 

parser.add_argument('-o', default='output_genome.fasta', help='Output file saving gap-closed genome. [output_genome.fasta]')
parser.add_argument('-t', type=int, default=1, help='Number of threads. [1]')
parser.add_argument('-qc', type=float, default=100, help='Maximum expected erroneous bases in the read used for local assembly. [100]')
parser.add_argument('-l', type=int, default=100, help='Length of the contig end sequence cut out for local assembly. [100]')
parser.add_argument('-rc', type=int, default=100, help='Coverage of reads used for local assembly. [100]')
parser.add_argument('-S', type=float, default=1, help='Scale for standard deviations of priori distance between reads. [1]')
parser.add_argument('-ma', type=int, default=1, help='Matching score in reads pairwise alignment. [1]')
parser.add_argument('-mm', type=int, default=20, help='Mismatch penalty in reads pairwise alignment. [20]')
parser.add_argument('-gc', type=int, default=30, help='Gap cost in reads pairwise alignment. [30]')
parser.add_argument('-ms', type=int, default=20, help='Minimum score to output in reads pairwise alignment. [20]')
parser.add_argument('-ho', type=int, default=0, help='Maximum admissible hanging-out length in reads pairwise alignment. [0]')
parser.add_argument('-w', action="store_true", help='Assign initial weights for pairwise alignments in robust regression OLC. [null]')
parser.add_argument('-thr1', type=float, default=2, help='Tuning constant of weight function in IRLS algorithm. [2]')
parser.add_argument('-thr2', type=float, default=3, help='Excluding samples with residuals greater than thr2 after IRLS algorithm. [3]')
parser.add_argument('-mA', type=int, default=1, help='Matching score in alignment merging adjacent contigs. [1]')
parser.add_argument('-mM', type=int, default=2, help='Mismatch penalty in alignment merging adjacent contigs. [2]')
parser.add_argument('-mG', type=int, default=3, help='Gap cost in alignment merging adjacent contigs. [3]')
parser.add_argument('-mS', type=int, default=20, help='Minimum alignment score to merge adjacent contigs. [20]')
parser.add_argument('-HO', type=int, default=5, help='Maximum admissible hanging-out length in alignment merging adjacent contigs. [5]')

args = parser.parse_args()

thread = args.t
Length = args.l
rc = args.rc
qc = args.qc
S = args.S
thr1 = args.thr1
thr2 = args.thr2
ma = args.ma
mm = args.mm
gc = args.gc
ms = args.ms
ho = args.ho
HO = args.HO
mA = args.mA
mM = args.mM
mG = args.mG
mS = args.mS
output = args.o

aligner = Align.PairwiseAligner(mode = 'local', match_score=1, mismatch_score=-2, open_gap_score=-5, extend_gap_score=-1)

contigSeq, contigInd=importContigs()

scaffold_dic=importEvidence()

importReads(scaffold_dic)

m=Manager()
gapDic=m.dict()

pool=Pool(thread)
num=0
for scaffold_num in scaffold_dic:
    for k in range(len(scaffold_dic[scaffold_num])):
        num+=1
        reads=scaffold_dic[scaffold_num][k]['reads']
        gapSize=scaffold_dic[scaffold_num][k]['gap']
        gapStd=scaffold_dic[scaffold_num][k]['std']
        pool.apply_async(func=closeGap,args=(gapDic,num,reads,gapSize,gapStd,S))
pool.close()
pool.join()

num=0
for scaffold_num in scaffold_dic:
    for k in range(len(scaffold_dic[scaffold_num])):
        num+=1
        scaffold_dic[scaffold_num][k]['gapSequence']=gapDic[num]

writeResults(scaffold_dic, mA, mM, mG, mS)
