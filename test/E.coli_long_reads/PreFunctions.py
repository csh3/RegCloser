from Bio.Seq import Seq
from Bio import SeqIO, Align
import numpy as np
from scipy.sparse import linalg,csr_matrix,csc_matrix,diags
import networkx as nx


def importLongReads(filename):
    ReadSeq={}
    ReadInd=[]
    
    for seq_record in SeqIO.parse(filename, 'fasta'):
        read_id = seq_record.id.split('.')[1]
        ReadSeq[read_id+"_F"] = str(seq_record.seq).upper()
        ReadInd.append(read_id+"_F")
        ReadSeq[read_id+"_R"] = str(Seq(seq_record.seq).reverse_complement()).upper()
        ReadInd.append(read_id+"_R")        
    return(ReadSeq,ReadInd)  
    
def STreeDFS(T,node,DFlag,Direction):
    sorted_nei = sorted(T.neighbors(node),key=lambda x:abs(T.edges[x,node]['weight']),reverse=True)
    for nei in sorted_nei:
        if Direction[nei]!=0:
            continue
        Direction[nei]=int(Direction[node]*DFlag[nei,node])
        #print(Direction[node],Direction[nei])
        STreeDFS(T,nei,DFlag,Direction)
    return Direction

def TreeDFS(T,node,DFlag,Direction):
    for nei in T.neighbors(node):
        if Direction[nei]!=0:
            continue
        Direction[nei]=int(Direction[node]*DFlag[nei,node])
        #print(Direction[node],Direction[nei])
        TreeDFS(T,nei,DFlag,Direction)
    return Direction
    
def sigmoid(x):
    return 1/(1+np.exp(-x))

def DirectMST(direct_dict,reads):
    DG = nx.Graph()
    ReadsCount = len(reads)
    print(ReadsCount)
    ReadsDirection = [0 for i in range(ReadsCount)]
    DFlag = np.zeros((ReadsCount,ReadsCount))
    J_mat = np.zeros((ReadsCount,ReadsCount))
    for n1 in range(ReadsCount):
        for n2 in range(n1,ReadsCount):
            c1 = reads[n1]
            c2 = reads[n2]
            # if c1 in trouble_reads or c2 in trouble_reads:
            # 	continue
            key = str(min(c1,c2))+'-'+str(max(c1,c2))
            d_list = []
            if key in direct_dict:	
                d_list = direct_dict[key]
                # print(key,d_list)				
            if len(d_list)<1:
                continue
            d_sum = sum(d_list)
            d_abssum = sum(np.abs(d_list))
            posi = (d_sum+d_abssum)/2
            nega = (d_abssum-d_sum)/2
            J_mat[n1,n2] = d_sum #np.tanh(d_sum)#sigmoid(d_sum)#np.sign(d_sum)*np.sqrt(abs(d_sum))
            J_mat[n2,n1] = d_sum #np.tanh(d_sum)#sigmoid(d_sum)#np.sign(d_sum)*np.sqrt(abs(d_sum))#d_sum
            if max(posi,nega)>0:
                DG.add_edge(n1,n2,weight=-max(posi,nega)) #np.tanh
                DFlag[n1,n2] = np.sign(d_sum)
                DFlag[n2,n1] = DFlag[n1,n2]
                # print(c1,c2,DFlag[n1,n2],"weight",-max(posi,nega))
            else:
                print("Why?",key,d_list,d_sum,d_abssum,posi,nega)
    MST = nx.minimum_spanning_tree(DG)
    s_node = 0 #sorted(G.degree(),key = lambda x:x[1],reverse=True)[0][0]
    ReadsDirection[s_node]=1
    Direction = STreeDFS(MST,s_node,DFlag,ReadsDirection)
    return (Direction,J_mat)

def setupRegressionModel_direct(reads,D, edge_dict):    
    row=0
    indptr=[0]
    indices=[]
    data=[]
    response=[]

    size=len(reads)
    response= []
    for i in range(size-1):
        for j in range(i+1,size):
            A = reads[i]
            B = reads[j]
            key = str(A)+"-"+str(B)
            if key in edge_dict:
                for overlap in edge_dict[key]:##总是reads[j]['start']-reads[i]['start']，且以i的方向为正方向
                    if overlap['B_Orientation']==1 and D[i]!=D[j]: #异向overlap
                        if D[i]==1: #reads[i] --> reads[j] <--
                            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']
                        elif D[i]==-1: #reads[i] <-- reads[j] --> 
                            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']#overlap['B_length']-overlap['B_end']-overlap['A_length']+overlap['A_end']
                            diff = -diff
                        else:
                            print("Error!")
                    elif overlap['B_Orientation']==0 and D[i]==D[j]: #同向overlap                        
                        if D[i]==1: #reads[i] --> reads[j] -->    
                            diff = overlap['A_start']-overlap['B_start']
                        elif D[i]==-1:#reads[i] <-- reads[j] <--
                            diff = overlap['A_start']-overlap['B_start']#overlap['A_end']-overlap['A_length']+overlap['B_end']
                            diff = -diff
                        else:
                            print("Error!")
                        # print(A,B,D[i],D[j],diff)
                    else:#else的情况说明该overlap的方向与估计的方向关系不符，应该被滤除
                        print("Inconsistent-direction overlap",overlap,D[i],D[j])
                        continue
                    row+=1
                    indices.append(j) 
                    data.append(1)
                    indices.append(i)
                    data.append(-1)                   
                    indptr.append(indptr[-1]+2)
                    response+=[diff]                             
                    
    design=csr_matrix((data, indices, indptr), shape=(row, size))
    response=np.matrix(response).T

    return(design, response)

def deleteRowsCsr(mat, indices):
    indices = list(indices)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask]

#####Huber M-estimate of reads coordinates 
def IRLS(X,Y,reads,perc=90,thr2=50):
    if Y.shape[0]==0:
        return([], [])

    residual_list = []
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute ordinary least squares estiamte
    residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    residual_list.append(residual)
    max_residual=max(residual)
    print("Initial max_residual:", round(max_residual,5),"Regression Size:",Y.shape[0])
    # print("Quantile 50 55 60 ... 95:",[np.percentile(residual,5*x+50) for x in range(10)],np.percentile(residual,99))
    thr1 = np.percentile(residual,perc)
    threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] #threshold of Huber's weight function
    old_estimate=estimate
    n=0
    diff = 100
    
    while n<100 and max_residual>thr2:
        index=np.where(residual>threshold)[0]
        reweight=np.ones(len(residual))
        print("Quantile 50 55 60 ... 95,99:",[np.percentile(residual,50+5*x) for x in range(1,9)],np.percentile(residual,99))
    
        thr1 = np.percentile(residual,perc) #set thr1 to 90/95 percentile of all residuals
        threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] 
        reweight[index]=threshold[index]/residual[index] # update weights for alignments with residuals greater than thr1
        reweight=diags(reweight)
        t=X.T
        A=t.dot(reweight).dot(X)
        y=csr_matrix(Y)
        b=t.dot(reweight).dot(y).todense()
        estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute weighted least squares estimate
        residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        residual_list.append(residual)
        max_residual=max(residual)
        #print("max_residual:", round(max_residual,5))
        diff=max(abs(estimate-old_estimate))
        if diff<100: # convergence condition of estimates
            break
        else:
            old_estimate=estimate
            n+=1
    print("IRLS Stop at the %d rounds with max diff %d, max residual %f"%(n,diff,round(max_residual,5)))
    
    outlier_overlap=[]
    for i in range(len(residual)):        
        n1=X.indices[2*i]
        n2=X.indices[2*i+1]
        if residual[i]>thr2:
            outlier_overlap.append((reads[n1],reads[n2],residual[i]))
    
    index=np.where(residual>thr2)[0]
    X=deleteRowsCsr(X,index)
    Y=np.delete(Y,index,0)
    
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05)
    residual_2=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0] ## residuals after step 2
    #max_residual=max(residual)

    #remove alignments with residuals greater than thr2
    G = nx.Graph()
    for i in range(X.shape[0]):
        n1=X.indices[2*i]
        n2=X.indices[2*i+1]
        G.add_edge(n1,n2)

#divide reads into separate connected components
    reads_list=[]
    estimates_list=[]
    index_list=[]
    for c in nx.connected_components(G):
        # sub_index=list(G.subgraph(c).nodes)
        # sub_estimates=[estimate[i] for i in sub_index]
        # sub_reads=[reads[i] for i in sub_index]
        # estimates_list.append(list(map(int,np.round(sub_estimates))))
        # reads_list.append(sub_reads)
        # index_list.append(sub_index)
        sub_index=list(G.subgraph(c).nodes)
        sub_reads=[reads[i] for i in sub_index]
        sub_estimates=[estimate[i] for i in sub_index]
        sortIndex=np.argsort(sub_estimates)
        sub_index=[sub_index[k] for k in sortIndex]
        sub_reads=[sub_reads[k] for k in sortIndex]
        sub_estimates=[sub_estimates[k] for k in sortIndex]
        estimates_list.append(list(map(int,np.round(sub_estimates))))
        reads_list.append(sub_reads)
        index_list.append(sub_index)
    return(estimates_list, reads_list, index_list, outlier_overlap, residual_list, residual_2)

def MAD(list):
    m = np.median(list)
    return np.median([abs(x-m) for x in list])

#####Huber M-estimate of reads coordinates 
def IRLS_MAD(X,Y,reads,thr2=50,iter=500):
    if Y.shape[0]==0:
        return([], [])

    residual_list = []
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute ordinary least squares estiamte
    residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    residual_list.append(residual)
    max_residual=max(residual)
    print("Initial max_residual:", round(max_residual,5),"Regression Size:",Y.shape[0])
    # print("Quantile 50 55 60 ... 95:",[np.percentile(residual,5*x+50) for x in range(10)],np.percentile(residual,99))
    thr1 = 1.4826*MAD(residual)#np.percentile(residual,perc)
    threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] #threshold of Huber's weight function
    old_estimate=estimate
    n=0
    diff = 100
    
    while n<iter and max_residual>thr2:
        index=np.where(residual>threshold)[0]
        reweight=np.ones(len(residual))
        print("Quantile 50 55 60 ... 95,99:",[np.percentile(residual,50+5*x) for x in range(1,9)],np.percentile(residual,99))
    
        thr1 = 1.4826*MAD(residual)#np.percentile(residual,perc) #set thr1 to 90/95 percentile of all residuals
        threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] 
        reweight[index]=threshold[index]/residual[index] # update weights for alignments with residuals greater than thr1
        reweight=diags(reweight)
        t=X.T
        A=t.dot(reweight).dot(X)
        y=csr_matrix(Y)
        b=t.dot(reweight).dot(y).todense()
        estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute weighted least squares estimate
        residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        residual_list.append(residual)
        max_residual=max(residual)
        #print("max_residual:", round(max_residual,5))
        diff=max(abs(estimate-old_estimate))
        if diff<100: # convergence condition of estimates
            break
        else:
            old_estimate=estimate
            n+=1
    print("IRLS Stop at the %d rounds with max diff %d, max residual %f"%(n,diff,round(max_residual,5)))
    
    outlier_overlap=[]
    for i in range(len(residual)):        
        n1=X.indices[2*i]
        n2=X.indices[2*i+1]
        if residual[i]>thr2:
            outlier_overlap.append((reads[n1],reads[n2],residual[i]))
    
    index=np.where(residual>thr2)[0]
    X=deleteRowsCsr(X,index)
    Y=np.delete(Y,index,0)
    
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05)
    residual_2=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0] ## residuals after step 2
    #max_residual=max(residual)

    #remove alignments with residuals greater than thr2
    G = nx.Graph()
    for i in range(X.shape[0]):
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
        estimates_list.append(list(map(int,np.round(sub_estimates))))
        reads_list.append(sub_reads)

    return(estimates_list, reads_list, outlier_overlap, residual_list, residual_2)

def IRLS_Huber(X,Y,reads,thr1=100,thr2=200,iter=100):
    if Y.shape[0]==0:
        return([], [])

    residual_list = []
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute ordinary least squares estiamte
    residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    residual_list.append(residual)
    max_residual=max(residual)
    print("Initial max_residual:", round(max_residual,5),"Regression Size:",Y.shape[0])
    # print("Quantile 50 55 60 ... 95:",[np.percentile(residual,5*x+50) for x in range(10)],np.percentile(residual,99))
    threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] #threshold of Huber's weight function
    old_estimate=estimate
    n=0
    diff = 100
    
    while n<iter and max_residual>thr2:
        index=np.where(residual>threshold)[0]
        reweight=np.ones(len(residual))
        print("Quantile 50 55 60 ... 95,99:",[np.percentile(residual,50+5*x) for x in range(1,9)],np.percentile(residual,99))
    
        # thr1 = np.percentile(residual,perc) #set thr1 to 90/95 percentile of all residuals
        threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] 
        reweight[index]=threshold[index]/abs(residual[index]) # update weights for alignments with residuals greater than thr1
        reweight=diags(reweight)
        t=X.T
        A=t.dot(reweight).dot(X)
        y=csr_matrix(Y)
        b=t.dot(reweight).dot(y).todense()
        estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute weighted least squares estimate
        residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        residual_list.append(residual)
        max_residual=max(residual)
        #print("max_residual:", round(max_residual,5))
        diff=np.percentile(abs(estimate-old_estimate),95)#max(abs(estimate-old_estimate))
        if diff<10: # convergence condition of estimates
            break
        else:
            old_estimate=estimate
            n+=1
    print("IRLS Stop at the %d rounds with max diff %d, max residual %f"%(n,diff,round(max_residual,5)))
    
    outlier_overlap=[]
    for i in range(len(residual)):        
        n1=X.indices[2*i]
        n2=X.indices[2*i+1]
        if residual[i]>thr2:
            outlier_overlap.append((reads[n1],reads[n2],residual[i]))
    
    index=np.where(residual>thr2)[0]
    X=deleteRowsCsr(X,index)
    Y=np.delete(Y,index,0)
    
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05)
    residual_2=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0] ## residuals after step 2
    #max_residual=max(residual)

    #remove alignments with residuals greater than thr2
    G = nx.Graph()
    for i in range(X.shape[0]):
        n1=X.indices[2*i]
        n2=X.indices[2*i+1]
        G.add_edge(n1,n2)

#divide reads into separate connected components
    reads_list=[]
    estimates_list=[]
    index_list = []
    for c in nx.connected_components(G):
        # sub_index=list(G.subgraph(c).nodes)
        # sub_estimates=[estimate[i] for i in sub_index]
        # sub_reads=[reads[i] for i in sub_index]
        # estimates_list.append(list(map(int,np.round(sub_estimates))))
        # reads_list.append(sub_reads)
        sub_index=list(G.subgraph(c).nodes)
        sub_reads=[reads[i] for i in sub_index]
        sub_estimates=[estimate[i] for i in sub_index]
        sortIndex=np.argsort(sub_estimates)
        sub_index=[sub_index[k] for k in sortIndex]
        sub_reads=[sub_reads[k] for k in sortIndex]
        sub_estimates=[sub_estimates[k] for k in sortIndex]
        estimates_list.append(list(map(int,np.round(sub_estimates))))
        reads_list.append(sub_reads)
        index_list.append(sub_index)

    return(estimates_list, reads_list, index_list, outlier_overlap, residual_list, residual_2)