# CML VUS screening

This code is organized to analyze the results of the CML VUS screening experiment.

### SuppleCode1: epegRNA abundance library design  
This code is for designing the library for pegRNA abundance-based CML VUS screening. It uses the [DeepPrime-FT (Yu et al., Cell, 2023)](https://doi.org/10.1016/j.cell.2023.03.034) implemented in the Python package [GenET](https://github.com/Goosang-Yu/genet) to design all possible pegRNAs that can create each variant and predict their prime editing efficiency. The optimal pegRNAs were selected to construct the library.


### SuppleCode2: epegRNA abundance screening analysis  
This code was used to preprocess and analyze the data obtained from pegRNA abundance-based CML VUS screening. After preprocessing the raw FASTQ files uploaded to SRA to keep only the parts to be analyzed, read counts for each barcode sequence were counted, and UMI clustering was performed to organize the read counts for each variant. The variant-specific read counts obtained were analyzed using the [MAGeCK](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4) pipeline.

### SuppleCode3: SynPrime library design  
This code is for designing the library for SynPrime-based CML VUS screening.


### SuppleCode4: SypPrime screening analysis  
This code was used to preprocess and analyze the data obtained from SynPrime-based CML VUS screening. Alignment was performed using [CRISPResso2](https://github.com/pinellolab/CRISPResso2), and subsequent analyses were conducted based on the frequency table created after alignment. The code also analyzed the patterns of reads matching the pegRNAs designed to create each variant and the remaining reads. Additionally, it includes analysis code used for further cell-level validation of hit variants identified from the screening results.

### SuppleCode5: In vivo validation analysis  
This code analyzes the experimental results validating the CML TKI resistant variants identified by SynPrime screening in an in vivo xenograft model.


# Raw sequencing data
The raw NGS files used in this analysis are all uploaded to SRA (BioProject; [PRJNA1048659](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1048659)). The accession numbers of the FASTQ files used in each pipeline can be found in the README of each SupplementaryCode.

# Environments
These codes were tested in Ubuntu 22.04 LTS or CentOS 7 environments.


# Requirements
- Python >= 3.7
- biopython
- pandas
- genet
- moepy
- scipy
- numpy
- 
