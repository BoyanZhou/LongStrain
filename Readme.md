**LongStrain for longitudinal metagenomic data**
===========================================

This code was used for the paper "An integrated strain-level analytic pipeline utilizing longitudinal metagenomic data".


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

2. 

```
   
```

## Usage

There are following steps for analysis of longitudinal metagenomic data using LongStrain.

### step1: fastq_preprocess
Process the longitudinal samples of one subject at each time. Repeat this step until all samples are processed.
```python longstrain.py -m fastq_preprocess -i PATH/subject1 -o PATH/all_subject_results -l PATH/all_subject_results/all_subjects.log -s subject1 -p subject_t1,subject_t2 -f subject_t1_R1.fastq.gz,subject_t1_R2.fastq.gz:subject_t2.fastq.gz --kraken_database PATH_to_kraken_database/NCBI_standard```

* ```-m fastq_preprocess``` is the processing mode for step1;
* ```-i PATH/subject1``` is the directory of all longitudinal fq files from subject1;
* ```-o PATH/all_subject_results``` is the directory to store results from all subjects;
* ```-l PATH/all_subjects.log``` is the absolute path of log file;
* ```-s subject1``` is the name of an input subject, one subject each time;
* ```-p subject_t1,subject_t2``` are prefix of longitudinal samples, separated by ",";
* ```-f subject_t1_R1.fastq.gz,subject_t1_R2.fastq.gz:subject_t2.fastq.gz```  are corresponding longitudinal samples, fastqs from different samples separated by ':', paired fastqs from one sample separated by ',';
* ```--kraken_database PATH_to_kraken_database/NCBI_standard``` is the path to the directory of standard database of Kraken.

### step2: relative_abundance_aggregation
This step is optional. 
```python longstrain.py -m relative_abundance_aggregation -i PATH/all_subject_results -o PATH/all_subject_results -l PATH/all_subject_results/all_subjects.log --subject_list subject1,subject2,subject3 --output_combined_file PATH/all_subject_results/combined_relative_abundance_file.txt```

* ```-m relative_abundance_aggregation``` is the processing mode for step2;
* ```-i PATH/all_subject_results``` is the output directory of the step1 and the input directory of the step2;
* ```-o PATH/all_subject_results``` is the output directory, which should be same with the input directory;
* ```-l PATH/all_subjects.log``` is the absolute path of log file;
* ```--subject_list subject1,subject2,subject3``` is the name of all input subject, separated by ",";
* ```--output_combined_file PATH/all_subject_results/combined_relative_abundance_file.txt``` is the absolute path of aggregated relative abundance of all species.

### step3: build_species_reference



### step4: reads_assignment


### step5: longitudinal_analysis



**NOTE:** LongStrain WILL NOT work on a single sample.


**Please contact boyanzhou1992@gmail.com for any questions or bug reporting.**
