**LongStrain for longitudinal metagenomic data**
===========================================

This code was used for the paper ""


## Requirements

1.  Kraken2 (tested with 2.0.8)
2.  Bowtie2
3.  This code was written and tested on python 3.6.5, and requires the following packages:
    - numpy (tested with )
    - pysam (tested with )
    - pandas (tested with )

    If you encounter issues, please try to run in an environment with these packages.

## Install

1. Download the files from 

```
   
```

## Usage

There are following steps for analysis of longitudinal metagenomic data using LongStrain.

### step1: fastq_preprocess
Process the longitudinal samples of one subject at each time.
```python longstrain.py -m fastq_preprocess -i PATH/subject1 -o PATH/all_subject_results -l PATH/all_subjects.log -s subject1 -p subject_t1,subject_t2 -f subject_t1_R1.fastq.gz,subject_t1_R2.fastq.gz:subject_t2.fastq.gz --kraken_database PATH_to_kraken_database/NCBI_standard```

```-i PATH/subject1``` is the directory of all longitudinal fq files from subject1;
```-o PATH/all_subject_results``` is the directory to store results from all subjects;
```-l PATH/all_subjects.log``` is the absolute path of log file;
```-s subject1``` is the name of an input subject, one subject each time;
```-p subject_t1,subject_t2``` are prefix of longitudinal samples, separated by ",";
```-f subject_t1_R1.fastq.gz,subject_t1_R2.fastq.gz:subject_t2.fastq.gz```  are corresponding longitudinal samples, fastqs from different samples separated by ':', paired fastqs from one sample separated by ',';
```--kraken_database PATH_to_kraken_database/NCBI_standard``` is the path to the directory of standard database of Kraken.

### step2: relative_abundance_aggregation


### step3: build_species_reference

### step4: reads_assignment


### step5: longitudinal_analysis



**NOTE:** LongStrain WILL NOT work on a single sample.


**Please contact boyanzhou1992@gmail.com for any questions or bug reporting.**
