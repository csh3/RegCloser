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
                    if abs(diff)<minDiff and aligned[0]>=minscore and min(l2-aligned[5],l1-aligned[3])<=hangingOut and min(aligned[4]-1,aligned[2]-1)<=hangingOut:
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


#####Huber M-estimate, followed by trimmed least squares
def IRLS(X,Y,W,reads,r1=2,r2=3):
    if Y.shape[0]==0:
        return([], [])

    X0=X
    Y0=Y
    W0=W

    W=sp.diags(W0)
    X=W.dot(X0)
    Y=W.dot(Y0) #initial weight multiplies each alignment

    t=X.T
    A=t.dot(X)
    y=sp.csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = sp.linalg.lgmres(A, b, atol=1e-05) #compute ordinary least squares estiamte as initial estimate
    residual=abs((X.dot(sp.csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    threshold=sp.csr_matrix(np.ones(len(residual))*r1).dot(W).toarray()[0] #threshold r1 of Huber's weight function
    old_estimate=estimate
    n=0
    
    while n<100:
        index=np.where(residual>threshold)[0]
        reweight=np.ones(len(residual))
        reweight[index]=threshold[index]/residual[index] # update weights for alignments with residuals greater than r1
        reweight=sp.diags(reweight)
        A=t.dot(reweight).dot(X)
        b=t.dot(reweight).dot(y).todense()
        estimate, exitCode = sp.linalg.lgmres(A, b, atol=1e-05) #compute weighted least squares estimate
        residual=abs((X.dot(sp.csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        diff=max(abs(estimate-old_estimate))
        if diff<=2: # convergence condition of estimates
            break
        else:
            old_estimate=estimate
            n+=1

#remove alignments with residuals greater than r2 and calculate least squares solution
    residual0=abs((X0.dot(sp.csr_matrix(estimate).T)-sp.csr_matrix(Y0)).todense()).T.getA()[0]
    index=np.where(residual0>r2)[0]
    X0=deleteRowsCsr(X0,index)
    Y0=np.delete(Y0,index,0)
    W0=np.delete(W0,index,0)

    W=sp.diags(W0)
    X=W.dot(X0)
    Y=W.dot(Y0)
    t=X.T
    A=t.dot(X)
    y=sp.csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = sp.linalg.lgmres(A, b, atol=1e-05)

#divide reads into separate connected components
    G = nx.Graph()
    for i in range(X.shape[0]):
        n1=X.indices[2*i]
        n2=X.indices[2*i+1]
        G.add_edge(n1,n2)

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

    if ('endSeqL' not in reads_name) and ('endSeqR' in reads_name):
        estimates=[-(estimates[i]-len(reads[i]['seq'])) for i in range(len(estimates))]
        for i in range(len(reads)):
            reads[i]['seq']=reads[i]['seq'][::-1]
            reads[i]['phred']=reads[i]['phred'][::-1]

    ind=np.array(estimates).argsort()
    reads=[reads[i] for i in ind]
    estimates=[estimates[i] for i in ind]
    reads_name=[r['name'] for r in reads]

    endLenL=endLenR=0
    if 'endSeqR' in reads_name:
        ind_start=reads_name.index('endSeqR')
        endLenR=len(reads[ind_start]['seq'])
    if 'endSeqL' in reads_name:
        ind_start=reads_name.index('endSeqL')
        endLenL=len(reads[ind_start]['seq'])

    reads=reads[ind_start:]
    estimates=estimates[ind_start:]
    estimates=[e-(estimates[0]-len(reads[0]['seq'])) for e in estimates]

    consensus=[{c:[0]} for c in reads[0]['seq']]
    phreds=[{reads[0]['seq'][i]:[1-10**(-reads[0]['phred'][i]/10)]} for i in range(len(reads[0]['phred']))]
    readsName=[reads[0]['name']]
    readNum=1
    adjust=10

    query_backup_list=[reads[0]['seq']]
    estimate_backup_list=[estimates[0]]
    coordinate_backup_list=[len(consensus)]

    for i in range(1,len(reads)):
        query=reads[i]['seq']
        phred=reads[i]['phred']
        consensusSeq=''
        for count in consensus:
            base=max(count, key=lambda k:len(count[k]))
            if base=='-':
                bkp=count.pop('-')
                base=max(count, key=lambda k:len(count[k]))
                count['-']=bkp
            consensusSeq+=base

        estimate=estimates[i]
        score1=0
        for k in range(len(query_backup_list)):
            k+=1
            query_backup=query_backup_list[-k]
            estimate_backup=estimate_backup_list[-k]
            coordinate_backup=coordinate_backup_list[-k]
            alignments = MultiAlignment_func_python.MultiOptimal_merge_alignment(query_backup,query,1,2,3,3)
            score=0
            for aligned in alignments:
                if aligned[0]>=20 and len(aligned[1])>=20 and abs((len(query)-aligned[5])-(len(query_backup)-aligned[3])-(estimate-estimate_backup))<=adjust:
                    score=aligned[0]
                    diff=(len(query)-aligned[5])-(len(query_backup)-aligned[3])
            if score>=20:
                coordinate=coordinate_backup+diff
                score1=score
                break
        if score1<20:
            continue

        startPoint=max(coordinate-len(query)-adjust, 0)
        endPoint=max(coordinate+adjust, 0)
        ref=consensusSeq[startPoint:endPoint]
        if ref=='':
            continue
        alignments = MultiAlignment_func_python.MultiOptimal_merge_alignment(ref,query,1,2,3,3)
        score=0
        minDiff=float('inf')
        for aligned in alignments:
            if aligned[0]>=20 and len(aligned[1])>=20 and abs(startPoint+aligned[3]+len(query)-aligned[5]-coordinate)<=min(adjust, minDiff):
                score=aligned[0]
                start1=aligned[2]
                end1=aligned[3]
                start2=aligned[4]
                end2=aligned[5]
                minDiff=abs(startPoint+end1+len(query)-end2-coordinate)
                # hoL=min(aligned[2]-1, aligned[4]-1)
                # hoR=min(len(ref)-aligned[3], len(query)-aligned[5])
        if score>=20:
            if len(query_backup_list)>=10:
                query_backup_list.pop(0)
                estimate_backup_list.pop(0)
                coordinate_backup_list.pop(0)
            query_backup_list.append(query)
            estimate_backup_list.append(estimate)
            coordinate_backup_list.append(startPoint+end1+len(query)-end2)

            alignments = aligner.align(ref[start1-1:end1],query[start2-1:end2])
            alignment = alignments[0]
            cigar_ref = alignment.aligned[0]
            cigar_query = alignment.aligned[1]

            for j in range(len(cigar_ref)):
                c_r=cigar_ref[j]
                c_q=cigar_query[j]
                for k in range(c_r[1]-c_r[0]):
                    base=query[start2-1+c_q[0]+k]
                    prob=1-10**(-phred[start2-1+c_q[0]+k]/10)
                    if base in consensus[startPoint+start1-1+c_r[0]+k]:
                        consensus[startPoint+start1-1+c_r[0]+k][base].append(readNum)
                        phreds[startPoint+start1-1+c_r[0]+k][base].append(prob)
                    else:
                        consensus[startPoint+start1-1+c_r[0]+k][base]=[readNum]
                        phreds[startPoint+start1-1+c_r[0]+k][base]=[prob]

            hoC=len(consensusSeq)-startPoint-end1+len(ref[start1-1:end1])-cigar_ref[-1][1]
            hoQ=len(query)-end2+len(query[start2-1:end2])-cigar_query[-1][1]
            if hoC>=hoQ:
                for j in range(hoQ):
                    base=query[len(query)-hoQ+j]
                    prob=1-10**(-phred[len(query)-hoQ+j]/10)
                    if base in consensus[len(consensusSeq)-hoC+j]:
                        consensus[len(consensusSeq)-hoC+j][base].append(readNum)
                        phreds[len(consensusSeq)-hoC+j][base].append(prob)
                    else:
                        consensus[len(consensusSeq)-hoC+j][base]=[readNum]
                        phreds[len(consensusSeq)-hoC+j][base]=[prob]
            else:
                for j in range(hoC):
                    base=query[len(query)-hoQ+j]
                    prob=1-10**(-phred[len(query)-hoQ+j]/10)
                    if base in consensus[len(consensusSeq)-hoC+j]:
                        consensus[len(consensusSeq)-hoC+j][base].append(readNum)
                        phreds[len(consensusSeq)-hoC+j][base].append(prob)
                    else:
                        consensus[len(consensusSeq)-hoC+j][base]=[readNum]
                        phreds[len(consensusSeq)-hoC+j][base]=[prob]
                for j in range(hoC,hoQ):
                    base=query[len(query)-hoQ+j]
                    prob=1-10**(-phred[len(query)-hoQ+j]/10)
                    consensus.append({base:[readNum]})
                    phreds.append({base:[prob]})

            hoC=startPoint+start1-1+cigar_ref[0][0]
            hoQ=start2-1+cigar_query[0][0]
            if hoQ<=hoC:
                for j in range(hoC-hoQ,hoC):
                    base=query[j-(hoC-hoQ)]
                    prob=1-10**(-phred[j-(hoC-hoQ)]/10)
                    if base in consensus[j]:
                        consensus[j][base].append(readNum)
                        phreds[j][base].append(prob)
                    else:
                        consensus[j][base]=[readNum]
                        phreds[j][base]=[prob]
            else:
                for j in range(hoC):
                    base=query[hoQ-hoC+j]
                    prob=1-10**(-phred[hoQ-hoC+j]/10)
                    if base in consensus[j]:
                        consensus[j][base].append(readNum)
                        phreds[j][base].append(prob)
                    else:
                        consensus[j][base]=[readNum]
                        phreds[j][base]=[prob]
                        
            for j in range(len(cigar_ref)-1):
                left=cigar_ref[j][1]
                right=cigar_ref[j+1][0]
                if right-left>0:
                    for k in range(left,right):
                        if '-' in consensus[startPoint+start1-1+k]:
                            consensus[startPoint+start1-1+k]['-'].append(readNum)
                            phreds[startPoint+start1-1+k]['-'].append(1)
                        else:
                            consensus[startPoint+start1-1+k]['-']=[readNum]
                            phreds[startPoint+start1-1+k]['-']=[1]

            for j in range(len(cigar_ref)-1):
                if cigar_ref[j][1]==cigar_ref[j+1][0]:
                    for k in range(cigar_query[j+1][0]-cigar_query[j][1]):
                        site=startPoint+start1-1+cigar_ref[j][1]+k
                        base=query[start2-1+cigar_query[j][1]+k]
                        prob=1-10**(-phred[start2-1+cigar_query[j][1]+k]/10)
                        insertC={base:[readNum]}
                        insertP={base:[prob]}
                        for l in range(readNum):
                            if (l in sum(consensus[site-1].values(),[])) and (l in sum(consensus[site].values(),[])):
                                if '-' in insertC:
                                    insertC['-'].append(l)
                                    insertP['-'].append(1)
                                else:
                                    insertC['-']=[l]
                                    insertP['-']=[1]
                        consensus.insert(site,insertC)
                        phreds.insert(site,insertP)

            readsName.append(reads[i]['name'])
            readNum+=1

    start=end=0
    for i in range(len(consensus)):
        if 0 in sum(consensus[i].values(),[]):
            end=i
    trim=end-start+1
    trimL=trimR=0
    if ('endSeqL' in readsName) and ('endSeqR' not in readsName):
        position='L'
        trimL=trim
    if ('endSeqL' not in readsName) and ('endSeqR' in readsName):
        position='R'
        trimR=trim
        consensus=consensus[::-1]
        phreds=phreds[::-1]
    if ('endSeqL' in readsName) and ('endSeqR' in readsName):
        position='RL'
        trimL=trim
        end=0
        ind=readsName.index('endSeqR')
        for i in range(len(consensus)):
            if ind in sum(consensus[i].values(),[]):
                end=i
        trim=len(consensus)-(end+1)
        consensus=consensus[:len(consensus)-trim]
        phreds=phreds[:len(phreds)-trim]
        start=0
        for i in range(len(consensus)):
            if ind not in sum(consensus[i].values(),[]):
                start=i
        trimR=end-start

    while(True):
        error_reads=[]
        essential_reads=[]
        for count in consensus:
            if count=={}:
                continue
            base=max(count, key=lambda k:len(count[k]))
            for k in count:
                if k!=base:
                    error_reads+=count[k]
            if len(sum(count.values(),[]))==1:
                for k in count:
                    essential_reads+=count[k]
        if error_reads==[]:
            break
        else:
            count_error_reads=collections.Counter(error_reads)
            for essential_read in essential_reads:
                if essential_read in count_error_reads:
                    count_error_reads.pop(essential_read)
            if count_error_reads=={}:
                break
            else:
                remove_read=max(count_error_reads, key=count_error_reads.get)
                if count_error_reads[remove_read]<=5:
                    break
                else:
                    for i in range(len(consensus)):
                        for k in consensus[i]:
                            if remove_read in consensus[i][k]:
                                ind=consensus[i][k].index(remove_read)
                                consensus[i][k].pop(ind)
                                phreds[i][k].pop(ind) 

    phred = ''
    for i in range(len(phreds)):
        A=G=C=T=D=0
        if 'A' in phreds[i]:
            for prob in phreds[i]['A']:
                A+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        if 'G' in phreds[i]:
            for prob in phreds[i]['G']:
                G+=math.log(0.99995*prob);A+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        if 'C' in phreds[i]:
            for prob in phreds[i]['C']:
                C+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);A+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        if 'T' in phreds[i]:
            for prob in phreds[i]['T']:
                T+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);A+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        if '-' in phreds[i]:
            for prob in phreds[i]['-']:
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

        consensus[i]=site
        if site != '':
            phreds[i] = chr(round(-10*math.log10(max(1-max(PA,PG,PC,PT,PD), 0.00006)))+33)
        else:
            phreds[i] = ''

    if position=='RL':
        filledSeq=consensus[trimL:len(consensus)-trimR]
        filledSeq=('').join(filledSeq)
        filledPhred=phreds[trimL:len(phreds)-trimR]
        filledPhred=('').join(filledPhred)
    if position=='L':
        consensus=consensus[trimL:]
        phreds=phreds[trimL:]
    if position=='R':
        consensus=consensus[:len(consensus)-trimR]
        phreds=phreds[:len(phreds)-trimR]
    consensus=('').join(consensus)
    phreds=('').join(phreds)

    if position=='RL':
        return({'position':position, 'seq':consensus, 'phred':phreds, 'endLenL':endLenL, 'endLenR':endLenR, 'filledSeq':filledSeq, 'filledPhred':filledPhred})
    else:
        return({'position':position, 'seq':consensus, 'phred':phreds})


def closeGap(gapDic,num,reads,gapSize,gapStd,scale):
    reads.sort(key = lambda x:x['location'])
    reads_name_old=[r['name'] for r in reads]
    if 'endSeqL' in reads_name_old:
        endSeqL=reads[reads_name_old.index('endSeqL')]
    if 'endSeqR' in reads_name_old:
        endSeqR=reads[reads_name_old.index('endSeqR')]

    readLen=sum([len(read['seq']) for read in reads])/len(reads)
    readNum=max(gapSize+3*gapStd+200,200)*rc/readLen
    step=max(round(len(reads)/min(readNum,rn)),1)
    reads=reads[0::step]
    reads_name=[r['name'] for r in reads]
    if 'endSeqL' in reads_name_old and 'endSeqL' not in reads_name:
        reads.append(endSeqL)
    if 'endSeqR' in reads_name_old and 'endSeqR' not in reads_name:
        reads.append(endSeqR)
    reads.sort(key = lambda x:x['location'])

    design, response, weight=setupRegressionModel(reads, gapStd, scale, ma, mm, gc, ms, ho)
    estimates_list, reads_list=IRLS(design, response, weight, reads, r1, r2)

    gapSequence={'L':'', 'phred_L':'', 'R':'', 'phred_R':''}
    for i in range(len(estimates_list)):
        estimates=estimates_list[i]
        reads=reads_list[i]
        if estimates != []:
            consensus=generateConsensus(reads, estimates)
            if consensus!='':
                if consensus['position']=='RL':
                    gapSequence={'RL':consensus['seq'], 'phred':consensus['phred'], 'endLenL':consensus['endLenL'], 'endLenR':consensus['endLenR'], 'filledSeq':consensus['filledSeq'], 'filledPhred':consensus['filledPhred']}
                elif consensus['position']=='L':
                    gapSequence['L']=consensus['seq']
                    gapSequence['phred_L']=consensus['phred']
                elif consensus['position']=='R':
                    gapSequence['R']=consensus['seq']
                    gapSequence['phred_R']=consensus['phred']

    gapDic[num]=gapSequence


def generateScaffold(scaffold_dic, scaffold_num, alignment_score = 1, mismatch = 2, gap_cost = 5, minscore=20):
    scafSeq = ''
    gapNum = 0 
    mergedNum = 0
    filledNum = 0
    closedNum = 0
    fillEvidence = []
    fillSeq = []
    old_contigLength_list = []
    new_contigLength_list = []
    for k in range(len(scaffold_dic[scaffold_num])):
        Gap = scaffold_dic[scaffold_num][k]
        if k==0:
            Contig = scaffold_dic[scaffold_num][k+1]
            fillEvidence.append(['^',Contig['contig'],'0',str(len(Gap['gapSequence']['R'])),'0','0'])
            currentContig = Gap['gapSequence']['R']+contigSeq[Contig['contig']]
        elif k==len(scaffold_dic[scaffold_num])-1:
            fillEvidence.append([Gap['contig'],'$',str(len(Gap['gapSequence']['L'])),'0','0','0'])
            currentContig += Gap['gapSequence']['L']
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
            if 'RL' in Gap['gapSequence']:
                endLenL = Gap['gapSequence']['endLenL']
                endLenR = Gap['gapSequence']['endLenR']
                currentContig = currentContig[:-endLenL]+Gap['gapSequence']['RL']+contigSeq[Contig['contig']][endLenR:]
                closedNum+=1
                filledLen=len(Gap['gapSequence']['filledSeq'])
                if filledLen > 0:
                    filledNum+=1
                    fillSeq.append(['@%s_%s_filled'%(Gap['contig'], Contig['contig']),Gap['gapSequence']['filledSeq'],'+',Gap['gapSequence']['filledPhred']])
                else:
                    filledLen=len(Gap['gapSequence']['RL'])-endLenL-endLenR
                    if filledLen>=0:
                        filledLen=0
                        filledNum+=1
                    else:
                        mergedNum+=1
                fillEvidence.append([Gap['contig'],Contig['contig'],str(filledLen),'0','1','0'])
            else:
                currentContig += Gap['gapSequence']['L']
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
                        lenL=len(currentContig)-len(Gap['gapSequence']['L'])
                        lenR=len(nextContig)-len(Gap['gapSequence']['R'])
                        currentContig = currentContig[:-truncate]+seq_L[:start_i-1]+overlap_seq+seq_R[end_j:]+nextContig[truncate:]
                        gapSeq=currentContig[lenL:-lenR]
                        if (len(currentContig)-lenL-lenR)>=0:
                            filledNum+=1
                            if len(gapSeq)>0:
                                lenGapL=min(len(gapSeq),max(len(Gap['gapSequence']['L'])-(len(seq_L)-end_i),0))
                                lenGapR=len(Gap['gapSequence']['R'])-(len(gapSeq)-lenGapL)
                                gapPhred=Gap['gapSequence']['phred_L'][:lenGapL]+Gap['gapSequence']['phred_R'][lenGapR:]
                                fillSeq.append(['@%s_%s_filled'%(Gap['contig'], Contig['contig']),gapSeq,'+',gapPhred])
                        else:
                            mergedNum+=1
                        fillEvidence.append([Gap['contig'],Contig['contig'],str(len(currentContig)-lenL-lenR),'0','1','0'])    
                        closedNum+=1
                        continue
                if gapSize+3*gapStd<0:
                    currentContig = currentContig[:len(currentContig)-len(Gap['gapSequence']['L'])]
                    scafSeq += currentContig
                    scafSeq += 'N'
                    new_contigLength_list.append(len(currentContig))
                    currentContig = nextContig[len(Gap['gapSequence']['R']):]
                    fillEvidence.append([Gap['contig'],Contig['contig'],'0','0','0',str(gapSize)])
                else:
                    if len(Gap['gapSequence']['L'])>0:
                        fillSeq.append(['@%s_%s_left'%(Gap['contig'], Contig['contig']),Gap['gapSequence']['L'],'+',Gap['gapSequence']['phred_L']])
                    new_contigLength_list.append(len(currentContig))
                    scafSeq+=currentContig
                    currentContig = nextContig
                    if len(Gap['gapSequence']['R'])>0:
                        fillSeq.append(['@%s_%s_right'%(Gap['contig'], Contig['contig']),Gap['gapSequence']['R'],'+',Gap['gapSequence']['phred_R']])
                    newGapSize = gapSize-len(Gap['gapSequence']['L'])-len(Gap['gapSequence']['R'])
                    fillEvidence.append([Gap['contig'],Contig['contig'],str(len(Gap['gapSequence']['L'])),str(len(Gap['gapSequence']['R'])),'0',str(newGapSize)])
                    if newGapSize<=0:
                        newGapSize = 1
                    scafSeq+='N'*newGapSize

    scaffold_dic[scaffold_num]={'scafSeq':scafSeq,'gapNum':gapNum,'mergedNum':mergedNum,'filledNum':filledNum,'closedNum':closedNum,'fillEvidence':fillEvidence,'fillSeq':fillSeq,'old_contigLength_list':old_contigLength_list,'new_contigLength_list':new_contigLength_list}


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
    gapFile = open('gapSequence.fastq', 'w')
    fillFile = open('evidence.fill', 'w')
    resultFile = open('statistics.txt', 'w')
    totalGapNum = 0
    totalClosedNum = 0
    totalMergedNum = 0
    totalFilledNum = 0
    old_contigLength_list = []
    new_contigLength_list = []
    scaffoldLength_list = []
    
    for scaffold_num in scaffold_dic:
        totalGapNum+=scaffold_dic[scaffold_num]['gapNum']
        totalClosedNum+=scaffold_dic[scaffold_num]['closedNum']
        totalMergedNum+=scaffold_dic[scaffold_num]['mergedNum']
        totalFilledNum+=scaffold_dic[scaffold_num]['filledNum']
        old_contigLength_list+=scaffold_dic[scaffold_num]['old_contigLength_list']
        new_contigLength_list+=scaffold_dic[scaffold_num]['new_contigLength_list']
        
        scafSeq=scaffold_dic[scaffold_num]['scafSeq']
        scaffoldLength_list.append(len(scafSeq))
        scafFile.write('>scaffold%d\n'%scaffold_num)
        writeSequence(scafFile, scafSeq)
        
        fillSeq=scaffold_dic[scaffold_num]['fillSeq']
        for record in fillSeq:
            gapFile.write(record[0].strip()+'\n')
            writeSequence(gapFile, record[1])
            gapFile.write(record[2].strip()+'\n')
            writeSequence(gapFile, record[3])

        fillEvidence=scaffold_dic[scaffold_num]['fillEvidence']
        fillFile.write('~~~~~\nscaffold%d\n'%scaffold_num)
        for record in fillEvidence:
            fillFile.write(record[0].ljust(20)+record[1].ljust(20)+record[2].ljust(10)+record[3].ljust(10)+record[4].ljust(10)+record[5]+'\n')
        
    oldContigN50 = computeN50(old_contigLength_list)
    newContigN50 = computeN50(new_contigLength_list)
    oldTotalContigLen = sum(old_contigLength_list)
    newTotalContigLen = sum(new_contigLength_list)
    scaffoldN50 = computeN50(scaffoldLength_list)

    resultFile.write("Total gap number: "+"{:,}\n".format(totalGapNum))
    resultFile.write("Closed gap number: "+"{:,}\n".format(totalClosedNum))
    resultFile.write("\tMerged gap number: "+"{:,}\n".format(totalMergedNum))
    resultFile.write("\tFilled gap number: "+"{:,}\n\n".format(totalFilledNum))
    resultFile.write("Total contig length from "+"{:,}".format(oldTotalContigLen)+" to "+"{:,}\n".format(newTotalContigLen))
    resultFile.write("Contig N50 from "+"{:,}".format(oldContigN50)+" to "+"{:,}\n".format(newContigN50))
    resultFile.write("\nScaffold N50: "+"{:,}\n".format(scaffoldN50))

    scafFile.close()
    gapFile.close()
    fillFile.close()
    resultFile.close()



#读入参数：库长标准差、比对参数、线程数、WLS的残差取舍参数、consensus的base占优比例参数
parser = argparse.ArgumentParser(description="Local assembly using the robust regression OLC approach.") 

parser.add_argument('-o', default='output_genome.fasta', help='Output file saving gap-closed genome. [output_genome.fasta]')
parser.add_argument('-t', type=int, default=1, help='Number of threads. [1]')
parser.add_argument('-qc', type=float, default=100, help='Maximum expected erroneous bases in the read used for local assembly. [100]')
parser.add_argument('-l', type=int, default=100, help='Length of the contig end sequence cut out for local assembly. [100]')
parser.add_argument('-rc', type=int, default=100, help='Maximum coverage of reads used for local assembly. [100]')
parser.add_argument('-rn', type=int, default=2000, help='Maximum number of reads used for local assembly. [2000]')
parser.add_argument('-S', type=float, default=1, help='Scale for standard deviations of priori distance between reads. [1]')
parser.add_argument('-ma', type=int, default=1, help='Matching score in reads pairwise alignment. [1]')
parser.add_argument('-mm', type=int, default=20, help='Mismatch penalty in reads pairwise alignment. [20]')
parser.add_argument('-gc', type=int, default=30, help='Gap cost in reads pairwise alignment. [30]')
parser.add_argument('-ms', type=int, default=20, help='Minimum score to output in reads pairwise alignment. [20]')
parser.add_argument('-ho', type=int, default=0, help='Maximum admissible hanging-out length in reads pairwise alignment. [0]')
parser.add_argument('-w', action="store_true", help='Assign initial weights for pairwise alignments in robust regression OLC. [null]')
parser.add_argument('-r1', type=float, default=2, help='Tuning constant of weight function in IRLS algorithm. [2]')
parser.add_argument('-r2', type=float, default=3, help='Excluding samples with residuals greater than this value after IRLS algorithm. [3]')
parser.add_argument('-mT', type=int, default=1000, help='Maximum truncated length for alignment merging adjacent contigs. [1000]')
parser.add_argument('-mA', type=int, default=1, help='Matching score in alignment merging adjacent contigs. [1]')
parser.add_argument('-mM', type=int, default=2, help='Mismatch penalty in alignment merging adjacent contigs. [2]')
parser.add_argument('-mG', type=int, default=3, help='Gap cost in alignment merging adjacent contigs. [3]')
parser.add_argument('-mS', type=int, default=20, help='Minimum alignment score to merge adjacent contigs. [20]')
parser.add_argument('-HO', type=int, default=5, help='Maximum admissible hanging-out length in alignment merging adjacent contigs. [5]')

args = parser.parse_args()

thread = args.t
Length = args.l
rc = args.rc
rn = args.rn
qc = args.qc
S = args.S
r1 = args.r1
r2 = args.r2
ma = args.ma
mm = args.mm
gc = args.gc
ms = args.ms
ho = args.ho
HO = args.HO
mT = args.mT
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
        Gap=scaffold_dic[scaffold_num][k]
        scaffold_dic[scaffold_num][k]={'contig':Gap['contig'], 'gap':Gap['gap'], 'std':Gap['std'], 'gapSequence':gapDic[num]}

scaffold_dic=m.dict(scaffold_dic)

pool=Pool(thread)
for scaffold_num in scaffold_dic:
    pool.apply_async(func=generateScaffold,args=(scaffold_dic, scaffold_num, mA, mM, mG, mS))
pool.close()
pool.join()

writeResults(scaffold_dic)
