# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

import argparse
import os
import datetime

start = datetime.datetime.now()

descript="RegCloser is a genome gap-closing tool based on a robust regression OLC approach, using paired sequence reads.\n\n Contact: cao.shenghao@foxmail.com"

parser = argparse.ArgumentParser(description=descript)

parser.add_argument('-p', type=str, default='prerequisite', help='A formatted file specifying the code path, reads directory, and library information. [Prerequisite]')
parser.add_argument('-g', required=True, help='Draft genome, required.')
parser.add_argument('-d', required=True, help='Working directory saving intermediate and output files, required.')
parser.add_argument('-o', default='output_genome.fasta', help='Output file saving gap-closed genome. [output_genome.fasta]')
parser.add_argument('-t', type=int, default=1, help='Number of threads. [1]')
parser.add_argument('-s', default='Start', help='Starting module. [Start]')
parser.add_argument('-e', default='End', help='Ending module. [End]')
parser.add_argument('-rs', action="store_true", help='Re-scaffold using SSPACE. [null]')
parser.add_argument("-f", type=int, default=0, help="Contigs shorter than this value will be removed from the draft genome before gap closing. [0]")
parser.add_argument('-ml',type=int, default=3,help="Contig end with insert size + ml * std bases were cut out for mapping reads. [3]")
parser.add_argument('-mk', type=int, default=19, help='Minimum seed length in BWA mapping. [19]')
parser.add_argument('-mT', type=int, default=30, help='Minimum score to output in BWA mapping. [30]')
parser.add_argument('-c', type=int, default=5, help='Maximum number of soft-clipped bases on either end of a mapped read. [5]')
parser.add_argument('-mq', type=int, default=60, help='Mapped Reads with mapping quality greater than this value will be identified as anchored reads. [60]')
parser.add_argument('-nr', action="store_true", help='Not re-map anchored reads to the whole draft genome. [null]')
parser.add_argument('-hf', action="store_true", help="Filter out anchored reads falling in the high coverage regions. [null]")
parser.add_argument('-ra', type=float, default=1.8, help='Consecutive bases with coverage higher then ra * mode coverage will be marked as high coverage regions. [1.8]')
parser.add_argument('-k', type=int, default=5, help='Minimum number of links to compute scaffold in SSPACE. [5]')
parser.add_argument('-a', type=float, default=0.7, help='Maximum link ratio between two best contig pairs in SSPACE. [0.7]')
parser.add_argument('-sd', type=int, default=100, help="Default standard deviation of gap size. [100]")
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
parser.add_argument('-r1', type=float, default=2, help='Tuning constant of weight function in IRLS algorithm. [2]')
parser.add_argument('-r2', type=float, default=3, help='Excluding samples with residuals greater than this value after IRLS algorithm. [3]')
parser.add_argument('-mA', type=int, default=1, help='Matching score in alignment merging adjacent contigs. [1]')
parser.add_argument('-mM', type=int, default=2, help='Mismatch penalty in alignment merging adjacent contigs. [2]')
parser.add_argument('-mG', type=int, default=3, help='Gap cost in alignment merging adjacent contigs. [3]')
parser.add_argument('-mS', type=int, default=20, help='Minimum alignment score to merge adjacent contigs. [20]')
parser.add_argument('-HO', type=int, default=5, help='Maximum admissible hanging-out length in alignment merging adjacent contigs. [5]')

args = parser.parse_args()

print("\n--------------------------------------\n")
print("Running RegCloser: -d %s -g %s -t %d -s %s -e %s\n"%(args.d, args.g, args.t, args.s, args.e))

with open(args.p,'r') as fin:
	path=fin.readline().strip().split(':')[1].strip()
	if path.endswith('/'):
		path=path[:-1]

if args.nr:
	os.system("python %s/Configure.py -nr -t %d -c %d -ml %d -k %d -a %f -ra %f -p %s -mk %d -mT %d -mq %d"%(path, args.t, args.c, args.ml, args.k, args.a, args.ra, args.p, args.mk, args.mT, args.mq))
else:
	os.system("python %s/Configure.py -t %d -c %d -ml %d -k %d -a %f -ra %f -p %s -mk %d -mT %d -mq %d"%(path, args.t, args.c, args.ml, args.k, args.a, args.ra, args.p, args.mk, args.mT, args.mq))

progress_running = 0

if(args.s=="Start"):
	progress_running = 1
	if os.path.exists(args.d):
		os.system("mv %s %s_old"%(args.d,args.d))
	os.system("mkdir %s"%args.d)

#InitialContig
if(args.s=="InitialContig" or progress_running == 1):
	progress_running = 1
	print("\n--------------------------------------\nInitialContig running.\n")
	os.system("cp %s %s/draft_genome.fasta"%(args.g, args.d))
	if args.rs:
		os.system("cd %s; python %s/CutScaffoldAtN.py -rs -t %d -f %d"%(args.d, path, args.t, args.f))
	else:
		os.system("cd %s; python %s/CutScaffoldAtN.py -t %d -f %d"%(args.d, path, args.t, args.f))
if(args.e=="InitialContig"):
	progress_running = 0

#Mapping
if(args.s=="Mapping" or progress_running == 1):
	progress_running = 1
	print("\n--------------------------------------\nMapping running.\n")
	os.system("cp scripts/Mapping.sh %s"%args.d)
	os.system("cd %s; bash ./Mapping.sh"%args.d)
if(args.e=="Mapping"):
	progress_running = 0

#HighCoverage
if(args.s=="HighCoverage" or progress_running == 1):
	progress_running = 1
	if args.hf:
		print("\n--------------------------------------\nHighCoverage running.\n")
		os.system("cp scripts/HighCoverage.sh %s"%args.d)
		os.system("cd %s; bash ./HighCoverage.sh"%args.d)
	else:
		os.system("cd %s; touch highCoverage.bed"%args.d)
if(args.e=="HighCoverage"):
	progress_running = 0

#CollectReads
if(args.s=="CollectReads" or progress_running == 1):
	progress_running = 1
	print("\n--------------------------------------\nCollectReads running.\n")
	os.system("cp scripts/CollectReads.sh %s"%args.d)
	os.system("cd %s; bash ./CollectReads.sh"%args.d)
if(args.e=="CollectReads"):
	progress_running = 0

#Scaffold
if args.rs and (args.s=="Scaffold" or progress_running == 1):
	progress_running = 1
	print("\n--------------------------------------\nScaffold running.\n")
	os.system("cp scripts/SSPACE.sh scripts/libraries.txt %s"%args.d)
	os.system("cd %s; bash ./SSPACE.sh >> Scaffold.log 2>&1"%args.d)
if (not args.rs) and (args.s=="Scaffold" or progress_running == 1):
	progress_running = 1
	os.system("cd %s; mv evidence_ini.txt evidence.txt"%args.d)	
if(args.e=="Scaffold"):
	progress_running = 0

#ReEstimateGapSize
if(args.s=="ReEstimateGapSize" or progress_running == 1):
	progress_running = 1
	print("\n--------------------------------------\nReEstimateGapSize running.\n")
	os.system("cd %s; python %s/ReEstimateGapSize.py -p ../%s -sd %d > ReEstimateGapSize.log 2>&1"%(args.d,path,args.p, args.sd))
if(args.e=="ReEstimateGapSize"):
	progress_running = 0

#LocalAssembly
if(args.s=="LocalAssembly" or progress_running == 1):
	progress_running = 1
	print("\n--------------------------------------\nLocalAssembly running.\n")
	#os.system("cp Banded_Alignment.so %s"%args.d)
	if args.w:
		os.system("cd %s; python %s/LocalAssembly.py -w -t %d -l %d -rc %d -S %f -qc %f -r1 %f -r2 %f -ma %d -mm %d -gc %d -ms %d -ho %d -mA %d -mM %d -mG %d -mS %d -HO %d -o %s > LocalAssembly.log 2>&1"%(args.d, path, args.t, args.l, args.rc, args.S, args.qc, args.r1, args.r2, args.ma, args.mm, args.gc, args.ms, args.ho, args.mA, args.mM, args.mG, args.mS, args.HO, args.o))
	else:
		os.system("cd %s; python %s/LocalAssembly.py -t %d -l %d -rc %d -S %f -qc %f -r1 %f -r2 %f -ma %d -mm %d -gc %d -ms %d -ho %d -mA %d -mM %d -mG %d -mS %d -HO %d -o %s > LocalAssembly.log 2>&1"%(args.d, path, args.t, args.l, args.rc, args.S, args.qc, args.r1, args.r2, args.ma, args.mm, args.gc, args.ms, args.ho, args.mA, args.mM, args.mG, args.mS, args.HO, args.o))
if(args.e=="LocalAssembly"):
	progress_running = 0

end = datetime.datetime.now()
print("\nExecution time:", end-start, "\n")
