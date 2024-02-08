**LongStrain for longitudinal metagenomic data**
===========================================

This code was used for the paper "An integrated strain-level analytic pipeline utilizing longitudinal metagenomic data".  
Please report to *boyanzhou1992@gmail.com* if you have any questions.  

## Requirements

1.  Kraken2 (tested with 2.0.8)  
https://ccb.jhu.edu/software/kraken2/
2.  Bowtie2 (tested with 0.16.0.1)  
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
3.  samtools (tested with 1.9)  
http://www.htslib.org/doc/
4.  This code was written and tested on python 3.6.5, and requires the following packages:
    - numpy (tested with 1.19.5)
    - pysam (tested with 0.16.0.1)
    - pandas (tested with 1.1.5)
    
    **Note:** Kraken2, Bowtie2, and samtools need to be added to the PATH Environment Variable. 
    If you encounter issues, please try to run in an environment with these packages.

## Install

1. Install the required software and packages

2. Download all files and deposit them under a directory.

## Usage
The detailed introduction can be found in the WIKI: https://github.com/BoyanZhou/LongStrain/wiki


* ```-l PATH/all_subjects.log``` is the absolute path of log file.

**Warnning:** LongStrain works on longitudinal or concurrent samples, Does NOT work on a single sample.


**Please contact boyanzhou1992@gmail.com for any questions or bug reporting.**
