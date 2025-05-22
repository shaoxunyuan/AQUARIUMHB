

## Table of Contents

1. [Introduction](#1-introduction)  
2. [Required](#2-required)  
3. [Quick Start with Example Data](#3-quick-start-with-example-data)  
4. [Data Preprocessing](#4-data-preprocessing)  

## 1. Introduction

AQUARIUMHB is an integrated computational pipeline designed for comprehensive analysis of circular RNAs (circRNAs) in human blood samples from RNA-seq data. 

## 3. Quick Start with Example Data

### PRJNA429023 Dataset Demonstration

PRJNA429023 dataset 包含12个来源于血液样本的RNA-seq测序数据

- 6 healthy control samples  
- 6 intellectual disability patient samples  

Please visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108887  查看样本数据具体信息。

### PRJNA429023 Dataset Download

1. Obtain the SRR accession list file (SRR_Acc_List.txt) from:  
   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA429023&o=acc_s%3Aa  

2. Download all SRA files using:  
   
   ```bash
   prefetch SRR_Acc_List.txt
   ```
   
   
* **Author**: Shaoxun Yuan  
* **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China  
* **Contact**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)  


