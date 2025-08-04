# Fire limits soil microbial dispersal and differentially impacts bacterial and fungal communities

**Jacob R. Hopkins, Giuseppina Vizzari, Alison E. Bennett, Antonino Malacrinò**

Published in *Global Change Biology*: [10.1111/gcb.70386](https://doi.org/10.1111/gcb.70386)

Preprint: [10.1101/2025.05.02.651892](https://doi.org/10.1101/2025.05.02.651892)

## Abstract

Fire is a globally pervasive force reshaping ecosystems, yet its influence on the ecological processes structuring soil microbiomes remains poorly understood. Using a meta-analysis of >2,600 amplicon sequencing samples across 19 global studies, we tested whether fire alters soil microbiome assembly processes, diversity, and ecological selection for pyrophilic specialists. Contrary to prevailing assumptions, we found that fire did not significantly shift ecological selection processes in bacteria or fungi but instead constrained dispersal, particularly reducing dispersal in bacterial and fungal communities, and increasing ecological drift in fungi. Despite limited evidence for ecological selection, fire consistently filtered for specialist taxa, increasing their relative abundance across microbial communities. Fire also reduced fungal diversity and evenness, while bacterial communities exhibited greater dominance and loss of rare taxa. These findings support the idea that fire promotes microbial post-fire niche specialization while disrupting dispersal pathways. Our results indicate that increasing fire frequency and severity under climate change may homogenize soil microbial communities, reduce microbial resilience, and constrain ecosystem recovery.

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
