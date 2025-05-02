# Fire limits soil microbial dispersal and differentially impacts bacterial and fungal communities

**Jacob R. Hopkins, Giuseppina Vizzari, Alison E. Bennett, Antonino Malacrinò**

## Abstract

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. All raw data is available at NCBI SRA as specified in `01_samples/01_study_metadata.txt`. 

Our pipeline included:
* `nf-core/ampliseq` v2.10.0 [https://github.com/nf-core/ampliseq/](https://github.com/nf-core/ampliseq/)
* `MAFFT` [https://academic.oup.com/nar/article/30/14/3059/2904316](https://academic.oup.com/nar/article/30/14/3059/2904316)
* `FastTree` [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)
* `R` v4.4.1 [https://www.R-project.org/](https://www.R-project.org/)

# Code

### **1.** [Data download and curation](/03_organize_data)
Data and code for downloading and organizing raw data.

### **2.** [Raw data processing](/04_ampliseq)
Code for processing raw data using `nf-core/ampliseq`.

### **3.** [Data analysis - step 1](/05_analysis1)
Code for producing three major datasets for data analysis.

### **4.** [Data analysis - step 2](/06_analysis2)
Data and code for reproducing all data analysis presented in the manuscript. Given that this meta-analysis requires high computational resources, we decided to release the processed data together with the code to perform the analyses.
