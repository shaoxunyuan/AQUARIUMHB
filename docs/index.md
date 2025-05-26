

## Table of Contents

1. [Introduction](#1-introduction)  

2. [Quick Start with Example Data](#2-quick-start-with-example-data)  

3. [Alignment to Reference Genome](#3-alignment-to-reference-genome)  

4. [Identifying Circular RNAs Using CIRI-full](#4-identifying-circular-rnas-using-ciri-full)  

5. [Visualizing and Estimating Isoform Abundance Using CIRI-vis](#5-visualizing-and-estimating-isoform-abundance-using-ciri-vis)  

6. [Constructing a Full-length Reference Set](#6-constructing-a-full-length-reference-set)

7. [Reconstruction of partial length circRNA isoforms](#7-reconstruction-of-partial-length-circrna-isoforms)

8. [Quantifications of both linear and cicurlar isoforms](#8-quantifications-of-both-linear-and-cicurlar-isoforms)

9. [References](#9-references)

10. [Author Information](#10-author-information)

## 1. Introduction

AQUARIUMHB is a comprehesive toolkit for identify, annotate, quantify, and analyze human blood circular RNAs from RNA-seq data. 

## 2. Quick Start with Example Data

Using PRJNA429023 as example data

The PRJNA429023 dataset contains 12 RNA-seq samples derived from blood:

* 6 healthy control samples  
* 6 intellectual disability patient samples  

For detailed sample information, please visit:  [GSE108887](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108887)

- Obtain the SRR accession list file (`SRR_Acc_List.txt`) from:   [NCBI SRA Study](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA429023&o=acc_s%3Aa)  

- Download all SRA files using ```prefetch```:  

```bash
prefetch ./SRR_Acc_List.txt
```

- Batch convert with compression using ```fasterq-dump```

```bash
for sra in SRR{6450118..6450129}.sra; do

    fasterq-dump --split-3 "$sra" 

done
```

## 3. Align Paired-End Reads to Reference Genome

```bash
# Path of reference genome and annotation, index files for reference genome must made first using `bwa index`

fa=./Homo_sapiens.GRCh38.dna_sm.chromosomes.fa

gtf=./Homo_sapiens.GRCh38.94.chr.gtf

# Choose multi-threading based on the specific situation.

THREAD_COUNT=8 

# BioSample ID. Use SRR6450118 as an example:

BioSample=SRR6450118

dir_detect=${BioSample}/full

mkdir -p ${dir_detect}

bwa mem -T 10 ${fa} ${BioSample}_1.fastq ${BioSample}_2.fastq -o ${dir_detect}/align.sam -t ${THREAD_COUNT}

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

java -jar CIRI-full.jar RO2 -r ${fa} -s ${dir_detect}/full_ro1.sam -l 150 -o ${dir_detect}/full

java -jar CIRI-full.jar Merge -r ${fa} -a ${gtf} -c ${dir_detect}/ciri.report -as ${dir_detect}/as_jav.list -ro ${dir_detect}/full_ro2_info.list -o ${dir_detect}/full

```

Output files in ```dir_detect```  for subsequent analysis, circRNA information at BSJ levels

| SRR6450118/full/ciri.report |
| SRR6450119/full/ciri.report |
| ......                      |
| SRR6450129/full/ciri.report |

## 5. Visualize and estimate abundance of isoforms using [CIRI-vis](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-vis.html)

```bash
dir_vis=${BioSample}/vis

mkdir -p ${dir_vis}

java -jar CIRI-vis.jar -i ${dir_detect}/full_merge_circRNA_detail.anno -l ${dir_detect}/as_library_length.list -d ${dir_vis} -r ${fa} -min 1
```

Output files in ```dir_vis```  for subsequent analysis, circRNA information at isoform levels

| SRR6450118/vis/stout.list |
| SRR6450119/vis/stout.list |
| ......                    |
| SRR6450129/vis/stout.list |

## 6. Constructing a Full-length Reference Set

```r
R

library(AQUARIUMHB)

MakeReferenceIsoform(datapathfile = "PRJNA429023/DataPathFile.txt",outputfile = "ReferenceIsoformFinal.txt")
```

- `datapathfile`：input file contains path information of input files, with two columns needed: `SampleID` and `SamplePath`.

| SampleID   | SamplePath              |
| ---------- | ----------------------- |
| SRR6450118 | PRJNA429023/SRR6450118/ |
| SRR6450119 | PRJNA429023/SRR6450119/ |
| ......     | ......                  |
| SRR6450129 | PRJNA429023/SRR6450129/ |

- `outputfile`：output file contains all information of full-length transcripts from FLcircAS, IsoCirc, and blood samples. Example is as follows:

| chr | bsj                     | start     | end       | isoformID                                                       | strand | exon_count | exon_length | exon_total_length | ReferenceSource                                          |
| --- | ----------------------- | --------- | --------- | --------------------------------------------------------------- | ------ | ---------- | ----------- | ----------------- | -------------------------------------------------------- |
| 10  | 10:100036449\|100036604 | 100036604 | 100036449 | chr10\|100036449\|100036604\|+                                  | +      | 1          | 156         | 156               | FLcircAS_Liver                                           |
| 3   | 3:48548532\|48550234    | 48548532  | 48550234  | chr3\|48548532,48549864,48550118\|48548537,48549960,48550234\|- | -      | 3          | 6,97,117    | 220               | Full                                                     |
| ... | ...                     | ...       | ...       | ...                                                             | ...    | ...        | ...         | ...               | ...                                                      |
| 2   | 2:231497095\|231503204  | 231497095 | 231503204 | chr2\|231497095,231503081\|231497252,231503204\|-               | -      | 2          | 158,124     | 282               | Full,FLcircAS_HeLa,IsoCirc_SkeletalMuscle,IsoCirc_Testis |

`ReferenceSource` Sources of full-length isoforms. For example:


    `FLcircAS_Liver`: represents isoform from FLcircAS dataset, 

    `Full`: represents isoform from blood samples according pipeline,

    `Full,FLcircAS_HeLa,IsoCirc_SkeletalMuscle,IsoCirc_Testis`: represents isoform from FLcircAS, IsoCirc and blood samples.

## 7. Reconstruction of partial length circRNA isoforms

Due to reasons such as short sequencing read lengths or insufficient sequencing depth, the circRNAs identified by CIRT - full may fall into the following three types:

- **Full**: Those that contain the full - length information within the transcript.

- **Break**: Those with partial internal transcript information missing.

- **Only**: Those that contain only the information of the bsj locus.

AQUARIUMHB uses the Full - length Reference Set constructed in the previous step to reconstruct the Break and Only types of circRNAs, and then generates a gtf file for subsequent quantitative analysis for `salmon`.

- **circRNA_full.gtf**

```r
circRNA_full.gtf(datapathfile = "PRJNA429023/DataPathFile.txt",referencefile = "ReferenceIsoformFinal.txt")
```

Ouput gtf files in `quant` directory for each sample:  

| Sample ID  | Output File Path                  |
| ---------- | --------------------------------- |
| SRR6450118 | SRR6450118/quant/circRNA_full.gtf |
| SRR6450119 | SRR6450119/quant/circRNA_full.gtf |
| ...        | ...                               |
| SRR6450129 | SRR6450129/quant/circRNA_full.gtf |

- **circRNA_break.gtf**

```r
circRNA_break.gtf(datapathfile = "PRJNA429023/DataPathFile.txt",referencefile = "ReferenceIsoformFinal.txt")
```

Ouput gtf files in `quant` directory for each sample:

| Sample ID  | Output File Path                   |
| ---------- | ---------------------------------- |
| SRR6450118 | SRR6450118/quant/circRNA_break.gtf |
| SRR6450119 | SRR6450119/quant/circRNA_break.gtf |
| ...        | ...                                |
| SRR6450129 | SRR6450129/quant/circRNA_break.gtf |

- **circRNA_only.gtf**

```r
circRNA_only.gtf(datapathfile = "PRJNA429023/DataPathFile.txt",referencefile = "ReferenceIsoformFinal.txt")
```

Ouput gtf files in `quant` directory for each sample:

| Sample ID  | Output File Path                  |
| ---------- | --------------------------------- |
| SRR6450118 | SRR6450118/quant/circRNA_only.gtf |
| SRR6450119 | SRR6450119/quant/circRNA_only.gtf |
| ...        | ...                               |
| SRR6450129 | SRR6450129/quant/circRNA_only.gtf |

## 8. Quantifications of both linear and cicurlar isoforms

To quantify both linear and circular isoforms simultaneously, users need to set the following variables:

- path to the reference genome (`fa`)

- path to the reference genome annotation (`gtf)

- sample name for quantification (`BioSample`)

- paths to the paired-end sequencing FASTQ files (`fastq1` and `fastq2`)

- output path for the quantification results (`dir_quant`)

- number of threads (`THREAD_COUNT`)

```bash

## Path of reference genome and annotation

fa=./Homo_sapiens.GRCh38.dna_sm.chromosomes.fa

gtf=./Homo_sapiens.GRCh38.94.chr.gtf

THREAD_COUNT=8 

## Select sample for quantification

BioSample=SRR6450118

fastq1=./${BioSample}_1.fastq.gz

fastq2=./${BioSample}_2.fastq.gz

dir_quant=./${BioSample}/quant


```

To pseudo-linearize circular isoforms, the `make.adapt.sh` script is required to add adapter sequences to the circular isoform sequence file circRNA_raw.fa. 

`make.adapt.sh` can be downloaded from [AQUARIUMHB/tree/main/inst/scripts](https://github.com/shaoxunyuan/AQUARIUMHB/tree/main/inst/scripts#)



```bash
# linear / cicular reference fasta and annotation

gffread $gtf -g $fa -ME -w ${dir_quant}/ref_linear.fa

cat ./${dir_quant}/circRNA_full.gtf ./${dir_quant}/circRNA_break.gtf ./${dir_quant}/circRNA_only.gtf> ./${dir_quant}/circRNA_final.gtf

gffread ./${dir_quant}/circRNA_final.gtf -g $fa -ME -w ./${dir_quant}/circRNA_raw.fa

readlen=$(zcat "$fastq1" | awk 'NR==2 {print length($0)}')  

sh ./make.adapt.sh $readlen ./${dir_quant}/circRNA_raw.fa ./${dir_quant}/circRNA_final.fa

cat ./${dir_quant}/ref_linear.fa ./${dir_quant}/circRNA_final.fa > ./${dir_quant}/final.fa

cat $gtf ./${dir_quant}/circRNA_final.gtf > ./${dir_quant}/final.gtf
```

```bash
# salmon index

salmon index --kmerLen 31 --transcripts ./${dir_quant}/final.fa --index ./${dir_quant}/index_final --keepDuplicates
```

```bash
# salmon quant

salmon quant --index ./${dir_quant}/index_final --libType IU --geneMap ./${dir_quant}/final.gtf --threads $THREAD_COUNT --seqBias --gcBias --validateMappings \

                       --mates1 $fastq1 --mates2 $fastq2  \

                       --output ./${dir_quant}/profile_results
```

Ouput quantification result files in quant directory for each sample:

| Sample ID  | Output File Path                          |
| ---------- | ----------------------------------------- |
| SRR6450118 | SRR6450118/quant/profile_results/quant.sf |
| SRR6450119 | SRR6450119/quant/profile_results/quant.sf |
| ...        | ...                                       |
| SRR6450129 | SRR6450129/quant/profile_results/quant.sf |

Quantification result for `SRR6450118` in `SRR6450118/quant/profile_results/quant.sf` : 

| Name                                                                                                 | Length | EffectiveLength | TPM      | NumReads |
| ---------------------------------------------------------------------------------------------------- | ------ | --------------- | -------- | -------- |
| ENST00000456328                                                                                      | 1657   | 1382.022        | 0.553826 | 8.485    |
| ENST00000450305                                                                                      | 632    | 391.000         | 0.000000 | 0.000    |
| ENST00000488147                                                                                      | 1351   | 1330.162        | 1.968974 | 29.034   |
| ......                                                                                               | ...... | .......         | ......   | .......  |
| chr18\|79695225,79695977,79697860,79701813,79704767\|79695308,79696070,79697988,79702284,79704917\|+ | 1031   | 885.066         | 1.349396 | 13.240   |
| chr18\|79728907,79736355\|79729069,79736521\|+                                                       | 431    | 276.218         | 1.349472 | 4.132    |
| chrY\|11532940,11730889\|11532945,11730985\|+                                                        | 204    | 35.000          | 0.000000 | 0.000    |

## 9. References

[CIRI-full](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-full.html#) 

[AQUARIUM_HB: a bioinformatics pipeline for human blood circular RNA analysis](https://www.researchsquare.com/article/rs-5657706/v1#)

## 10. Author information

* **Author**: Shaoxun Yuan  
* **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China  
* **Contact**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)  


