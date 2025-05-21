# AQUARIUMHB: A Comprehensive Toolkit for Human Blood Circular RNA Analysis

* **Author**: Shaoxun Yuan  
* **Affiliation**: School of Artificial Intelligence and Information Technology, Nanjing University of Chinese Medicine, China  
* **Contact**: [yuanshaoxun@njucm.edu.cn](mailto:yuanshaoxun@njucm.edu.cn)  

## Table of Contents  
1. [Introduction](#1-introduction)  
2. [Installation](#2-installation)  
3. [Quick Start with Example Data](#3-quick-start-with-example-data)  
4. [Data Preprocessing](#4-data-preprocessing)  

## 1. Introduction  

AQUARIUMHB is an integrated computational pipeline designed for comprehensive analysis of circular RNAs (circRNAs) in human blood samples from RNA-seq data. The toolkit provides complete solutions for:  

- circRNA identification  
- Functional annotation  
- Quantitative analysis  
- Differential expression analysis  
- Visualization of blood-specific circRNAs  

## 2. Installation  

[Installation instructions will be provided here]  

## 3. Quick Start with Example Data  

### PRJNA429023 Dataset Demonstration  
The workflow will be demonstrated using the PRJNA429023 dataset containing:  
- 6 healthy control samples  
- 6 intellectual disability patient samples  

## 4. Data Preprocessing  

### Step 1: Download SRA Data  
1. Obtain the SRR accession list from:  
   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA429023&o=acc_s%3Aa  
2. Download all SRA files using:  
   ```bash
   prefetch SRR_Acc_List.txt
