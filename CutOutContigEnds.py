# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

from Bio import SeqIO
import sys
import argparse

descript="This program cut out contig ends for mapping reads.\n"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-l', required=True, help='Reads library, required.')
parser.add_argument('-tl', type=int, required=True, help='Cut out length of the contig end, required.')

args = parser.parse_args()

library = args.l
trimLength = args.tl

with open('reference_%s.fasta'%library, 'w') as fout:
    for seq_record in SeqIO.parse('initial_contig.fa', 'fasta'):
        contigName = seq_record.id
        contigSeq = str(seq_record.seq)
        length = len(contigSeq)
        if length <= 2*trimLength:
            fout.write('>'+contigName+'_A\n')
            line=int(length/100)
            for k in range(line):
                fout.write(contigSeq[100*k:100*(k+1)]+'\n')
            if length > line*100:
                fout.write(contigSeq[100*line:]+'\n')
        else:
            fout.write('>'+contigName+'_L\n')
            leftSeq = contigSeq[:trimLength]
            line=int(trimLength/100)
            for k in range(line):
                fout.write(leftSeq[100*k:100*(k+1)]+'\n')
            if trimLength>line*100:
                fout.write(leftSeq[100*line:]+'\n')
            fout.write('>'+contigName+'_R\n')
            rightSeq = contigSeq[-trimLength:]
            line=int(trimLength/100)
            for k in range(line):
                fout.write(rightSeq[100*k:100*(k+1)]+'\n')
            if trimLength>line*100:
                fout.write(rightSeq[100*line:]+'\n')
