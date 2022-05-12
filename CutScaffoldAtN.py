# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

import argparse
import re
from Bio import SeqIO
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('-t', type=int, default=1, help='Number of threads. [1]')
parser.add_argument("-f", type=int, default=0, help="Contigs shorter than this value will be removed from the draft genome before gap closing. [0]")
parser.add_argument('-rs', action="store_true", help='Re-scaffold using SSPACE. [null]')
args = parser.parse_args()

def writeResults(f, number, length, contig):
	f.write(">contig%d\n"%number)
	line=int(length/100)
	for k in range(line):
		f.write(contig[100*k:100*(k+1)]+'\n')
	if length>line*100:
		f.write(contig[100*line:]+'\n')

if args.rs:
	fout = open('initial_contig.fa','w')
	fout_length = open('contigLength.txt','w')
	filtered_seq = open('filtered.fasta','w')

	contigNum = 0
	totalContigNum = 0
	for seq_record in SeqIO.parse('draft_genome.fasta', 'fasta'):
		scaf_content=str(seq_record.seq)
		scaf_break_arr = re.split('N+', scaf_content)
		for new_contig in scaf_break_arr:
			length=len(new_contig)
			totalContigNum+=1
			if length>=args.f:
				contigNum += 1
				fout_length.write('contig%d\t%d\n'%(contigNum, length))
				writeResults(fout, contigNum, length, new_contig)
			else:
				writeResults(filtered_seq, totalContigNum, length, new_contig)

	fout.close()
	fout_length.close()
	filtered_seq.close()

else:
	contigNum = 0
	scafNum = 0
	evidence = {}
	for seq_record in SeqIO.parse('draft_genome.fasta', 'fasta'):
		scaf_content = str(seq_record.seq)	
		scaf_break_arr = re.split('(N+|n+)', scaf_content)
		contig_count = int(np.ceil(len(scaf_break_arr)/2))
		scafNum += 1
		evidence[scafNum]=[{'len':len(scaf_content),'gap':0}]
		for i in range(contig_count):
			new_contig = scaf_break_arr[2*i]
			if i==contig_count-1:
				scaf_N_count = 0
			else:
				scaf_N_count = len(scaf_break_arr[2*i+1])
			length=len(new_contig)
			contigNum+=1
			evidence[scafNum]+=[{'contig':contigNum,'length':length,'seq':new_contig,'gap':scaf_N_count}]

	fout = open('initial_contig.fa','w')
	fout_length = open('contigLength.txt','w')
	fout_evidence = open('evidence_ini.txt','w')
	filtered_seq = open('filtered.fasta','w')

	for i in range(1,scafNum+1):
		j = len(evidence[i])-1
		while j>0:
			ctg = evidence[i][j]['contig']
			if evidence[i][j]['length']<args.f:
				writeResults(filtered_seq,ctg,evidence[i][j]['length'],evidence[i][j]['seq'])
				evidence[i][j-1]['gap']+=evidence[i][j]['length']+evidence[i][j]['gap']
				del evidence[i][j]
			j=j-1

	ctg_id=0
	for i in range(1,scafNum+1):
		if len(evidence[i])<1:
			continue
		scaf_len = evidence[i][0]['len']
		if len(evidence[i])>1:
			fout_evidence.write("~~~~~\n%s\t%d\tN\n"%(i,scaf_len))
		for j in range(1,len(evidence[i])):
			ctg = evidence[i][j]['contig']
			ctg_id +=1
			fout_length.write('contig%d\t%d\n'%(ctg_id, evidence[i][j]['length']))
			fout_evidence.write('%d\tf\t%d\t0\t%d\t%d\t%d\t0\n'%(ctg_id, evidence[i][j]['length'], evidence[i][j]['gap'], evidence[i][j]['gap'], evidence[i][j]['gap']))
			writeResults(fout, ctg_id, evidence[i][j]['length'], evidence[i][j]['seq'])

	fout.close()
	fout_length.close()
	fout_evidence.close()
	filtered_seq.close()