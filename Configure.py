# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

'''
需要预先配置好prerequisite文件，格式如下：
reads_directory:(the directory in which you save reads)
library_name	filename1	filename2	library-mean	library-std		direction	LocalAssembly/ReEstimateGapSize
(lib_name 		file_name_1	file_name_2	    mean	       std		      FR/RF 	  1/2/3)
'''

import argparse
import os

parser = argparse.ArgumentParser(description="Configure scripts for pipeline.")

parser.add_argument('-p', type=str, default='prerequisite', help='A formatted file specifying the code path, reads directory, and library information. [Prerequisite]')
parser.add_argument('-t', type=int, default=1, help='Number of threads. [1]')
parser.add_argument('-ml', type=int, default=3,help="Contig end with insert size + ml * std bases were cut out for mapping reads. [3]")
parser.add_argument('-mk', type=int, default=19, help='Minimum seed length in BWA mapping. [19]')
parser.add_argument('-mT', type=int, default=30, help='Minimum score to output in BWA mapping. [30]')
parser.add_argument('-c', type=int, default=5, help='Maximum number of soft-clipped bases on either end of a mapped read. [5]')
parser.add_argument('-mq', type=int, default=60, help='Mapped Reads with mapping quality greater than this value will be identified as anchored reads. [60]')
parser.add_argument('-nr', action="store_true", help='Not re-map anchored reads to the whole draft genome. [null]')
parser.add_argument('-ra', type=float, default=1.8, help='Consecutive bases with coverage higher then ra * mode coverage will be marked as high coverage regions. [1.8]')
parser.add_argument('-k', type=int, default=5, help='Minimum number of links to compute scaffold in SSPACE. [5]')
parser.add_argument('-a', type=float, default=0.7, help='Maximum link ratio between two best contig pairs in SSPACE. [0.7]')


args = parser.parse_args()

if not os.path.exists("scripts"):
	os.system("mkdir scripts")
else:
	os.system("rm -rf scripts; mkdir scripts")

library_info={}
with open(args.p,'r') as fin: #'prerequisite'
	longest_ins = -1000
	longest_lib = ' '
	path=fin.readline().strip().split(':')[1].strip()
	if path.endswith('/'):
		path=path[:-1]
	reads_directory=fin.readline().strip().split(':')[1].strip()
	if reads_directory.endswith('/'):
		reads_directory=reads_directory[:-1]
	for line in fin:
		record = line.strip().split()
		if len(record) == 7:
			library_info[record[0]]={'filename1':record[1],'filename2':record[2], 'mean':int(record[3]), 'std':int(record[4]),'direction':record[5],'flag':record[6]}
			library_info[record[0]]['trimlength'] = int(record[3])+args.ml*int(record[4])+200
			if int(record[3])>longest_ins:
				longest_ins = int(record[3])
				longest_lib = record[0]

#Mapping.sh
with open('scripts/Mapping.sh','w+') as fout:
	fout.write('ln -s initial_contig.fa reference_all.fasta\n\n')
	for library in library_info:
		fout.write('echo Cutting out contig ends for lib%s\n'%library)
		fout.write('python %s/CutOutContigEnds.py -l %s -tl %d &\n'%(path, library, library_info[library]['trimlength']))
	fout.write('\nwait\n\n')
	fout.write('echo bwa indexing\n')
	fout.write("bwa index reference_all.fasta 2>Mapping.log &\n")
	for library in library_info:
		fout.write("bwa index reference_%s.fasta 2>>Mapping.log &\n"%library)
	fout.write('\nwait\n\n')
	for library in library_info:
		fout.write('echo bwa mapping %s to contig ends\n'%library_info[library]['filename1'])
		fout.write("bwa mem -a -k %d -T %d -t %d reference_%s.fasta %s/%s > %s_R1.sam 2>>Mapping.log\n"%(args.mk, args.mT, args.t, library, reads_directory, library_info[library]['filename1'], library))
		fout.write('echo bwa mapping %s to contig ends\n'%library_info[library]['filename2'])
		fout.write("bwa mem -a -k %d -T %d -t %d reference_%s.fasta %s/%s > %s_R2.sam 2>>Mapping.log\n"%(args.mk, args.mT, args.t, library, reads_directory, library_info[library]['filename2'], library))
		fout.write('\n')
	for library in library_info:
		fout.write('echo Filtering reads mapped to contig ends from lib%s\n'%library)
		fout.write("python %s/FilterMappedReads.py -c %d -l %s -d %s -f1 mapped_%s -f2 mapped_%s &\n"%(path, args.c, library, library_info[library]['direction'], library_info[library]['filename1'], library_info[library]['filename2']))
	fout.write('\nwait\n\n')
	for library in library_info:
		fout.write('echo bwa re-mapping mapped_%s to whole draft genome\n'%library_info[library]['filename1'])
		fout.write("bwa mem -k %d -T %d -t %d reference_all.fasta mapped_%s > %s_remap_R1.sam 2>>Mapping.log\n"%(args.mk, args.mT, args.t, library_info[library]['filename1'], library))
		fout.write('echo bwa re-mapping mapped_%s to whole draft genome\n'%library_info[library]['filename2'])
		fout.write("bwa mem -k %d -T %d -t %d reference_all.fasta mapped_%s > %s_remap_R2.sam 2>>Mapping.log\n"%(args.mk, args.mT, args.t, library_info[library]['filename2'], library))
		fout.write('\n')
	fout.write('rm reference* mapped*\n')

#HighCoverage
with open('scripts/HighCoverage.sh','w+') as fout:
	command = 'python %s/HighCoverage.py -t %d -c %d -tl %d -ra %f -hs %s_R1.sam %s_R2.sam'%(path, args.t, args.c, library_info[longest_lib]['trimlength'], args.ra, longest_lib, longest_lib)
	fout.write('echo identifying high coverage regions\n')
	fout.write(command+'\n')

#CollectReads
with open('scripts/CollectReads.sh','w+') as fout:
	fout.write('rm -rf FILEs\n')
	fout.write('mkdir FILEs\n\n')
	for library in library_info:
		fout.write('echo Collecting reads for local assembly and making tab files from lib%s\n'%library)
		if args.nr:
			fout.write('python %s/CollectReads.py -nr -c %d -l %s -d %s -m %d -s %d -f %s -tl %d -ml %d -mq %d &\n'%(path, args.c, library, library_info[library]['direction'], library_info[library]['mean'], library_info[library]['std'], library_info[library]['flag'], library_info[library]['trimlength'], args.ml, args.mq))
		else:
			fout.write('python %s/CollectReads.py -c %d -l %s -d %s -m %d -s %d -f %s -tl %d -ml %d -mq %d &\n'%(path, args.c, library, library_info[library]['direction'], library_info[library]['mean'], library_info[library]['std'], library_info[library]['flag'], library_info[library]['trimlength'], args.ml, args.mq))
	fout.write('\nwait\n')


#Scaffold
with open('scripts/libraries.txt','w+') as fout:
	for library in library_info:
		fout.write('lib%s\tTAB\ttab%s.txt\t%d\t%f\t%s\n'%(library,library,library_info[library]['mean'],2*library_info[library]['std']/library_info[library]['mean'], library_info[library]['direction']))
with open('scripts/SSPACE.sh','w+') as fout:
	fout.write('perl %s/sspace/SSPACE_Standard_v3.0.pl -l libraries.txt -s initial_contig.fa -k %d -a %f -T %d -b sspace_k%d\n'%(path, args.k, args.a, args.t, args.k))
	fout.write('python %s/ParseEvidence.py sspace_k%d/sspace_k%d.final.evidence evidence.txt\n'%(path, args.k, args.k))
	# fout.write('rm tab*')