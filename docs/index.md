

## Table of Contents

1. [Introduction](#1-introduction)  
2. [Quick Start with Example Data](#2-quick-start-with-example-data)  
3. [Data Preprocessing](#3-data-preprocessing)  
4. [Identify circular RNAs using CIRI-full](#4-identify-circular-rnas-using-ciri-full)  
5. [Visualize and estimate abundance of isoforms using CIRI-vis](#5-visualize-and-estimate-abundance-of-isoforms-using-ciri-vis)  
6. [Author Information](#6-author-information)

## 1. Introduction

AQUARIUMHB is an integrated computational pipeline designed for comprehensive analysis of circular RNAs (circRNAs) in human blood samples from RNA-seq data. 

## 2. Quick Start with Example Data

Using PRJNA429023 as example data

The PRJNA429023 dataset contains 12 RNA-seq samples derived from blood:

- 6 healthy control samples  
- 6 intellectual disability patient samples  

For detailed sample information, please visit:  [GSE108887](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108887)

1. Obtain the SRR accession list file (`SRR_Acc_List.txt`) from:   [NCBI SRA Study](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA429023&o=acc_s%3Aa)  

2. Download all SRA files using:  
   
   ```bash
   prefetch SRR_Acc_List.txt
   ```

3. Batch convert with compression
   
   ```bash
   for sra in SRR{6450118..6450129}.sra; do
   
       fasterq-dump --split-3 "$sra" 
   
   done
   ```

## 3. Align Paired-End Reads to Reference Genome

```bash
# Path of reference genome and annotation, index files for reference genome must  made first using `bwa index`

fa=Homo_sapiens.GRCh38.dna_sm.chromosomes.fa

gtf=Homo_sapiens.GRCh38.94.chr.gtf

BioSample=SRR6450118

bwa mem -T 10 ${fa} ${BioSample}_1.fastq ${BioSample}_2.fastq -o ${BioSample}/full/align.sam"

```

## 4. Identify circular RNAs using [CIRI-full](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-full.html#) pipeline

```bash
# Choose multi-threading based on the specific situation.
THREAD_COUNT=8 

# Use SRR6450118 as an example
BioSample=SRR6450118

dir_detect=${BioSample}/full

mkdir -p ${dir_detect}

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

Output files in ```dir_detect```  for subsequent analysis:

```ciri.report``` , circRNA information at BSJ levels


## 5. Visualize and estimate abundance of isoforms using [CIRI-vis](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-vis.html)

```bash

dir_vis=${BioSample}/vis

mkdir -p ${dir_vis}

java -jar CIRI-vis.jar -i ${dir_detect}/full_merge_circRNA_detail.anno -l ${dir_detect}/as_library_length.list -d ${dir_vis} -r ${fa} -min 1
```

Output files in ```dir_vis```  for subsequent analysis:

```stout.list``` , circRNA information at isoform levels

## 6. Author information

* **Author**: Shaoxun Yuan  
* **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China  
* **Contact**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)  
