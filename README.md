Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

# RegCloser

## 1. Introduction
RegCloser is a genome gap-closing tool using paired sequence reads. It closes gaps via local assembly of reads in the gap regions. The local assembly applies a robust regression OLC approach. It first performs pairwise alignment of reads guided by the priori distance information from insert size. Then it generates a robust layout by minimizing the global Huber loss function, followed by a trimmed regression. RegCloser works well in resolving tandem repeats.

## 2. Installation
You can download the software package by the command:

```
git clone https://github.com/csh3/RegCloser.git
```

or click the **Download Zip** button and decompress the software package.

## 3. Dependencies
The current version requires:

1. `Python3` with the following modules: 
`os`, `sys`, `re`, `argparse`, `biopython`, `numpy`, `math`, `networkx`, `scipy`, `collections`, `datetime`, `multiprocessing`

2. `BWA (version 0.7.17)`
3. `MultiAlignment_func_python.so`

The firt two can be installed through the [Bioconda](https://bioconda.github.io/) channel. The third one has been included in the software package or you can compile the source code `MultiAlignment_func_python.cpp` on your machine using the following command.

```
g++ -fPIC MultiAlignment_func_python.cpp -o MultiAlignment_func_python.so -shared -I/home/miniconda3/include/python3.6m 
# Here /home/miniconda3/include/python3.6m is a directory storing the head file Python.h
```

## 4. Usage
### 4.1. Pipeline
The main program is `RunPipeline.py`, and the pipeline consists of the following 7 modules. You can start or end with any one module by the option `-s` or `-e`.

```
1. InitialContig	Break the draft genome into contigs
2. Mapping 		Map sequence reads to the draft genome using BWA and identify anchored reads
3. HighDepth		Identify high depth regions in the contig ends (optional)
4. CollectReads		Collect reads in the gap regions for local assembly and make tab files of linking information between contigs 	
5. Scaffold		Generate new scaffolds from initial contigs using SSPACE_Standard_v3.0 (optional)
6. ReEstimateGapSize	Re-estimate gap sizes and their standard deviations in the scaffolds
7. LocalAssembly 	Assemble gap sequences via the robust regression OLC approach
```

We recommend to use *RegCloser* in an iterative way that you can take the output genome as the input of the next iteration, and perform several times until no more gaps to be filled.

### 4.2. Prerequisite file

A prerequisite file is needed to specify the path of software package and the directory storing the sequence reads, as well as the information of different libraries. An example is illustrated below.

```
Code_path: /home/RegCloser
Reads_directory: /home/E.coli/reads
frag		frag_1.fastq		frag_2.fastq		300     20     FR    Y
shortjump	shortjump_1.fastq	shortjump_2.fastq	3600    298    RF    Y
```
The first line specifies the code path, and the second line specifies the reads directory. From the third line, each line describes one reads library and contains 7 columns, separated by spaces or tabs. Each column is explained in more detail below.

```
Column 1: 	Name of the library
		Each library should be designated a different name. 
Column 2 & 3: 	Fastq files for both ends
		For each paired reads, one of the reads should be in the first file, and the other one in the second file. The paired reads are required to be on the same line.
Column 4:	Average insert size between paired reads
Column 5: 	Standard deviation of insert size between paired reads
Column 6: 	Orientation of paired reads, FR or RF
		F stands for --> orientation, and R for <-- orientation. Paired-end libraries are FR, and mate-pair libraries are RF.
Column 7: 	Whether the libary is used for local assembly or making tab files, 1, 2, or 3
		1 stands for only local assembly, 2 stands for only making tab files, and 3 for both.
```

### 4.3. Basic usage
Run with 40 threads

```
python RunPipeline.py -g draft_genome.fasta -d iter-1 -t 40
```

Re-run the *LocalAssembly* module

```
python RunPipeline.py -g draft_genome.fasta -d iter-1 -t 40 -s LocalAssembly
```

Iterate over the result of *RegCloser*

```
python RunPipeline.py -g iter-1/output_genome.fasta -d iter-2 -t 40
```

### 4.4. Output files
The intermediate and output files are saved under the directory specified by the option `-d`. *RegCloser* outputs 4 result files. They are described in details below.

**output_genome.fasta** saves the output genome sequence with gaps closed by *RegCloser*. You can specify the filename using option `-o` according to your preference.

**gapSequence.fastq** saves the assembled sequences in the gap regions and their Phred quality scores (ASCII_BASE 33). The identifier of each sequence records the gap it comes from. For example, *@contig1\_contig2\_filled* means the sequence filled the gap between *contig1* and contig2; *@contig2\_contig3\_left* means the sequence extended from the left boundary of the gap between *contig2* and *contig3*; *@contig3\_contig4\_right* means the sequence extended from the right boundary of the gap between *contig3* and *contig4*.
 
**evidence.fill** records the information of the gaps in the output genome. Each line describes one gap and contains 6 columns.
 
```
Column 1: 	Left contig of a gap
Column 2:   	Right contig of a gap
Column 3 & 4:	Length of sequences extending from the left and right boundaries of a gap
		If the gap was filled, column 4 is 0, and column 3 is the length of the filling sequence.
		If column 3 is a negative value, it means the two adjacent contigs flanking the gap were merged, and column 3 tells the overlap length.
Column 5: 	Status of a gap, 0 or 1
		If the gap was filled, the flag was set to 1, otherwise 0.
Column 6: 	Current gap size 
		If the gap was filled, the value is 0.
```
 
 **statistics.txt** records the statistics of the genome sequence after gap closing. It includes `closed gap number`, classified into `merged gap number` and `filled gap number`. It also includes `total contig length`, `contig N50`, and `scaffold N50`.

### 4.5. Command line options

|      Option     | Type  | Description |
| :----------: | :-------------: | :------------------ |
|    <code><b>-p</b></code>  | <i><font size=2>STR</font></i> | A formatted file specifying the code path, reads directory, and library information. [Prerequisite]|
|    <code><b>-g</b></code>  | <i><font size=2>STR</font></i> | Draft genome, required. |
|    <code><b>-d</b></code>  | <i><font size=2>STR</font></i> | Working directory saving intermediate and output files, required.|
|    <code><b>-o</b></code>  | <i><font size=2>STR</font></i> | Output file saving gap-closed genome. [output_genome.fasta] |
|    <code><b>-t</b></code>  | <i><font size=2>INT</font></i> | Number of threads. [1] |
|    <code><b>-s</b></code>  | <i><font size=2>STR</font></i> | Starting module. [Start] |
|    <code><b>-e</b></code>  | <i><font size=2>STR</font></i> | Ending module. [End] |
|    <code><b>-rs</b></code>  | | Re-scaffold using SSPACE. [null] |
|    <code><b>-f</b></code>  | <i><font size=2>INT</font></i> | Contigs shorter than this value will be removed from the draft genome before gap closing. [0] |
|    <code><b>-ml</b></code>  | <i><font size=2>INT</font></i> | Contig end with insert size + ml * std bases were cut out for mapping anchored reads. [3] |
|    <code><b>-mk</b></code>  | <i><font size=2>INT</font></i> | Minimum seed length in BWA mapping. [19] |
|    <code><b>-mT</b></code>  | <i><font size=2>INT</font></i> | Minimum score to output in BWA mapping. [30] |
|    <code><b>-c</b></code>  | <i><font size=2>INT</font></i> | Maximum number of soft-clipped bases on either end of a mapped read. [5] |
|    <code><b>-mq</b></code>  | <i><font size=2>INT</font></i> | Mapped Reads with mapping quality greater than this value will be identified as anchored reads. [60] |
|    <code><b>-nr</b></code>  | | Not re-map anchored reads to the whole draft genome to exclude multi-mapped reads. [null] |
|    <code><b>-hf</b></code>  |  | Filter out anchored reads falling in the high coverage regions. [null] |
|    <code><b>-ra</b></code>  | <i><font size=2>FLOAT</font></i> | Consecutive bases with coverage higher then ra * mode coverage will be marked as high coverage regions. [1.8] |
|    <code><b>-k</b></code>  | <i><font size=2>INT</font></i> | Minimum number of links to compute scaffold in SSPACE. [5] |
|    <code><b>-a</b></code>  | <i><font size=2>FLOAT</font></i> | Maximum link ratio between two best contig pairs in SSPACE. [0.7] |
|    <code><b>-sd</b></code>  | <i><font size=2>INT</font></i> | Default standard deviation of gap size. [100] |
|    <code><b>-qc</b></code>  | <i><font size=2>FLOAT</font></i> | Maximum expected erroneous bases in the read used for local assembly. [100] |
|    <code><b>-l</b></code>  | <i><font size=2>INT</font></i> | Length of the contig end sequence cut out for local assembly. [100] |
|    <code><b>-rc</b></code>  | <i><font size=2>INT</font></i> | Coverage of reads used for local assembly. [100] |
|    <code><b>-S</b></code>  | <i><font size=2>FLOAT</font></i> | Scale for standard deviations of priori distance between reads. [1] |
|    <code><b>-ma</b></code>  | <i><font size=2>INT</font></i> | Matching score in reads pairwise alignment. [1] |
|    <code><b>-mm</b></code>  | <i><font size=2>INT</font></i> | Mismatch penalty in reads pairwise alignment. [20] |
|    <code><b>-gc</b></code>  | <i><font size=2>INT</font></i> | Gap cost in reads pairwise alignment. [30] |
|    <code><b>-ms</b></code>  | <i><font size=2>INT</font></i> | Minimum score to output in reads pairwise alignment. [20] |
|    <code><b>-ho</b></code>  | <i><font size=2>INT</font></i> | Maximum admissible hanging-out length in reads pairwise alignment. [0] |
|    <code><b>-w</b></code>  | | Assign initial weights for pairwise alignments in robust regression OLC. [null] |
|    <code><b>-r1 </b></code>  | <i><font size=2>FLOAT</font></i> | Tuning constant of weight function in IRLS algorithm. [2] |
|    <code><b>-r2 </b></code>  | <i><font size=2>FLOAT</font></i> | Excluding samples with residuals greater than this value after IRLS algorithm. [3] |
|    <code><b>-mA</b></code>  | <i><font size=2>INT</font></i> | Matching score in alignment merging adjacent contigs. [1] |
|    <code><b>-mM</b></code>  | <i><font size=2>INT</font></i> | Mismatch penalty in alignment merging adjacent contigs. [2] |
|    <code><b>-mG</b></code>  | <i><font size=2>INT</font></i> | Gap cost in alignment merging adjacent contigs. [3] |
|    <code><b>-mS</b></code>  | <i><font size=2>INT</font></i> | Minimum alignment score to merge adjacent contigs. [20] |
|    <code><b>-HO</b></code>  | <i><font size=2>INT</font></i> | Maximum admissible hanging-out length in alignment merging adjacent contigs. [5] |


## 5. Current version

The version of the current release is v1.0.


## 6. Contact

Please contact <cao.shenghao@foxmail.com> for any questions.

## 7. License

**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License**

For details, please read `RegCloser/License`.
