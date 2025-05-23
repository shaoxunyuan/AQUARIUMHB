

## Table of Contents

1. [Introduction](#1-introduction)  
2. [Quick Start with Example Data](#2-quick-start-with-example-data)  
3. [Data Preprocessing](#3-data-preprocessing)  
4. [Identifying Circular RNAs Using CIRI-full](#4-identifying-circular-rnas-using-ciri-full)  
5. [Visualizing and Estimating Isoform Abundance Using CIRI-vis](#5-visualizing-and-estimating-isoform-abundance-using-ciri-vis)  
6. [Reconstruction of Partial-Length CircRNAs](#6-reconstruction-of-partial-length-circrnas)
7. [Author Information](#7-author-information)

## 1. Introduction

AQUARIUMHB is an integrated computational pipeline designed for comprehensive analysis of circular RNAs (circRNAs) in human blood samples from RNA-seq data. 

## 2. Quick Start with Example Data

Using PRJNA429023 as example data

The PRJNA429023 dataset contains 12 RNA-seq samples derived from blood:

- 6 healthy control samples  
- 6 intellectual disability patient samples  

For detailed sample information, please visit:  [GSE108887](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108887)

1. Obtain the SRR accession list file (`SRR_Acc_List.txt`) from:   [NCBI SRA Study](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA429023&o=acc_s%3Aa)  

2. Download all SRA files using ```prefetch```:  
   
   ```bash
   prefetch SRR_Acc_List.txt
   ```

3. Batch convert with compression using ```fasterq-dump```
   
   ```bash
   for sra in SRR{6450118..6450129}.sra; do
   
       fasterq-dump --split-3 "$sra" 
   
   done
   ```

## 3. Align Paired-End Reads to Reference Genome

```bash
# Path of reference genome and annotation, index files for reference genome must made first using `bwa index`

fa=Homo_sapiens.GRCh38.dna_sm.chromosomes.fa

gtf=Homo_sapiens.GRCh38.94.chr.gtf

# BioSample ID. Use SRR6450118 as an example:

BioSample=SRR6450118

dir_detect=${BioSample}/full

mkdir -p ${dir_detect}

# Choose multi-threading based on the specific situation.
THREAD_COUNT=8 

bwa mem -T 10 ${fa} ${BioSample}_1.fastq ${BioSample}_2.fastq -o ${BioSample}/full/align.sam -t $THREAD_COUNT

```

## 4. Identify circular RNAs using [CIRI-full](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-full.html#) pipeline

```bash
# Choose multi-threading based on the specific situation.
THREAD_COUNT=8 

# Use SRR6450118 as an example
BioSample=SRR6450118

# CIRI2: Circular RNA identification based on multiple seed matching
perl CIRI2.pl --in ${dir_detect}/align.sam --out ${dir_detect}/ciri.report --ref_file ${fa} --anno ${gtf} --thread_num $THREAD_COUNT

# CIRI-AS is a detection tool for circRNA internal components and alternative splicing events.
perl CIRI-AS.pl --sam ${dir_detect}/align.sam --ciri ${dir_detect}/ciri.report --out ${dir_detect}/as --ref_file ${fa} --anno ${gtf} --output_all yes

# The CIRI-full Pipeline module is an automatic pipeline for detecting and reconstructing circRNAs.
java -jar $CIRI-full.jar RO1 -1 ${BioSample}_1.fastq.gz -2 ${BioSample}_2.fastq.gz -o ${dir_detect}/full -t $THREAD_COUNT

bwa mem -T 19 -t $THREAD_COUNT ${fa} ${dir_detect}/full_ro1.fq -o ${dir_detect}/full_ro1.sam -t $THREAD_COUNT

java -jar ${scriptdir}/CIRI-full.jar RO2 -r ${fa} -s ${dir_detect}/full_ro1.sam -l 150 -o ${dir_detect}/full

java -jar ${scriptdir}/CIRI-full.jar Merge -r ${fa} -a ${gtf} -c ${dir_detect}/ciri.report -as ${dir_detect}/as_jav.list -ro ${dir_detect}/full_ro2_info.list -o ${dir_detect}/full

```

Output files in ```dir_detect```  for subsequent analysis, circRNA information at BSJ levels

```
SRR6450118/full/ciri.report
SRR6450119/full/ciri.report
SRR6450120/full/ciri.report
SRR6450121/full/ciri.report
SRR6450122/full/ciri.report
SRR6450123/full/ciri.report
SRR6450124/full/ciri.report
SRR6450125/full/ciri.report
SRR6450126/full/ciri.report
SRR6450127/full/ciri.report
SRR6450128/full/ciri.report
SRR6450129/full/ciri.report
```

## 5. Visualize and estimate abundance of isoforms using [CIRI-vis](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-vis.html)

```bash

dir_vis=${BioSample}/vis

mkdir -p ${dir_vis}

java -jar CIRI-vis.jar -i ${dir_detect}/full_merge_circRNA_detail.anno -l ${dir_detect}/as_library_length.list -d ${dir_vis} -r ${fa} -min 1
```

Output files in ```dir_vis```  for subsequent analysis, circRNA information at isoform levels

```
SRR6450118/vis/stout.list
SRR6450119/vis/stout.list
SRR6450120/vis/stout.list
SRR6450121/vis/stout.list
SRR6450122/vis/stout.list
SRR6450123/vis/stout.list
SRR6450124/vis/stout.list
SRR6450125/vis/stout.list
SRR6450126/vis/stout.list
SRR6450127/vis/stout.list
SRR6450128/vis/stout.list
SRR6450129/vis/stout.list
```

## 6. Make a Full-length Reference Set

```r
library(AQUARIUMHB)

MakeReferenceIsoform(datapathfile = "DataPathFile.txt",outputfile = "ReferenceIsoformFinal.txt")

```

```datapathfile```：The DataPathFile.txt contains path information of input files, with two columns: SampleID and SamplePath.

| SampleID   | SamplePath                    |
|------------|-------------------------------|
| SRR6450118 | PRJNA429023/SRR6450118/       |
| SRR6450119 | PRJNA429023/SRR6450119/       |
| SRR6450120 | PRJNA429023/SRR6450120/       |
| SRR6450121 | PRJNA429023/SRR6450121/       |
| SRR6450122 | PRJNA429023/SRR6450122/       |
| SRR6450123 | PRJNA429023/SRR6450123/       |
| SRR6450124 | PRJNA429023/SRR6450124/       |
| SRR6450125 | PRJNA429023/SRR6450125/       |
| SRR6450126 | PRJNA429023/SRR6450126/       |
| SRR6450127 | PRJNA429023/SRR6450127/       |
| SRR6450128 | PRJNA429023/SRR6450128/       |
| SRR6450129 | PRJNA429023/SRR6450129/       |


```outputfile```：The filename of the output file.




## 7. Author information

* **Author**: Shaoxun Yuan  
* **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China  
* **Contact**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)  


```