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

2. Download all files and deposit them under a directory.

## Usage

There are following steps for analysis of longitudinal metagenomic data using LongStrain.

### step1: fastq_preprocess
Process the longitudinal samples of one subject at each time. Repeat this step until all samples are processed.  
```python longstrain.py -m fastq_preprocess -i PATH/subject1 -o PATH/all_subject_results -l PATH/all_subjects.log -s subject1 -p subject_t1,subject_t2 -f subject_t1_R1.fastq.gz,subject_t1_R2.fastq.gz:subject_t2.fastq.gz --kraken_database PATH_to_kraken_database/NCBI_standard```

* ```-m fastq_preprocess``` is the processing mode for step1;
* ```-i PATH/subject1``` is the directory of all longitudinal fq files from subject1;
* ```-o PATH/all_subject_results``` is the directory to store results from all subjects;
* ```-l PATH/all_subjects.log``` is the absolute path of log file;
* ```-s subject1``` is the name of an input subject, one subject each time;
* ```-p subject_t1,subject_t2``` are prefix of longitudinal samples, separated by ",";
* ```-f subject_t1_R1.fastq.gz,subject_t1_R2.fastq.gz:subject_t2.fastq.gz```  are corresponding longitudinal samples, fastqs from different samples separated by ':', paired fastqs from one sample separated by ',';
* ```--kraken_database PATH_to_kraken_database/NCBI_standard``` is the path to the directory of standard database of Kraken.

### step2: relative_abundance_aggregation
This step is optional. This processing gets the relative abundances of species and lists them in the order of decreasing. Generally, we can choose the top 30 species or species with relative abundance >1% in enough subjects. If you already have target species list, just skip this step.  
```python longstrain.py -m relative_abundance_aggregation -i PATH/all_subject_results -o PATH/all_subject_results -l PATH/all_subjects.log --subject_list subject1,subject2,subject3 --output_combined_file PATH/all_subject_results/combined_relative_abundance_file.txt```

* ```-m relative_abundance_aggregation``` is the processing mode for step2;
* ```-i PATH/all_subject_results``` is the output directory of the step1 and the input directory of the step2;
* ```-o PATH/all_subject_results``` is the output directory, which should be same with the input directory;
* ```-l PATH/all_subjects.log``` is the absolute path of log file;
* ```--subject_list subject1,subject2,subject3``` is the name of all input subject, separated by ",";
* ```--output_combined_file PATH/all_subject_results/combined_relative_abundance_file.txt``` is the absolute path of aggregated relative abundance of all species.

### step3: build_species_reference  
After determining the target species, we need to build the reference database for mapping. For each species, this step only needs to be done once. In other word, for a given species, if its reference database has been built, we do not need to do it again.    
```python longstrain.py -m build_species_reference --reference_built_path /Path_to_LongStrain_database --target_species_list Bifidobacterium_longum,Akkermansia_muciniphila -l PATH/all_subjects.log --kraken_database PATH_to_kraken_database/NCBI_standard```

* ```-m build_species_reference``` is the processing mode for step3;
* ```--reference_built_path``` is the absolute path of directory where the species's reference database is built;
* ```--target_species_list Bifidobacterium_longum,Akkermansia_muciniphila``` is the list of target species separated by ",", in which blank should be replaced by "_";
* ```-l PATH/all_subjects.log``` is the absolute path of log file;
* ```--kraken_database PATH_to_kraken_database/NCBI_standard``` is the path to the directory of standard database of Kraken.

### step4: reads_assignment
Assign reads marked by Kraken2 to target species under the directory of previous output.  
```python longstrain.py -m reads_assignment -i PATH/all_subject_results --taxon_species_json Path_to/taxon_species_name.json --subject_list subject1,subject2,subject3 --reference_built_path /Path_to_LongStrain_database --target_species_list Bifidobacterium_longum,Akkermansia_muciniphila -l PATH/all_subjects.log```

* ```-m reads_assignment``` is the processing mode for step4;
* ```-i PATH/all_subject_results``` is the input and output directory of the step4;
* ```--taxon_species_json Path_to/taxon_species_name.json``` is the absolute path of the database of taxon ID, this file is downloaded with the code and under the same directory with longstrain.py;
* ```--subject_list subject1,subject2,subject3``` is the name of all input subject, separated by ",";
* ```--reference_built_path /Path_to_LongStrain_database``` is the absolute path of directory where the species's reference database is built in step3;
* ```--target_species_list Bifidobacterium_longum,Akkermansia_muciniphila``` is the list of target species separated by ",", in which blank should be replaced by "_";
* ```-l PATH/all_subjects.log``` is the absolute path of log file.

### step5: longitudinal_analysis
The output of LongStrain are strains' variants in VCF format and proportions of the primary and secondary strain at each time point.  
```python longstrain.py -m longitudinal_analysis -i PATH/all_subject_results -o PATH/all_species_summary --taxon_species_json Path_to/taxon_species_name.json --reference_built_path /Path_to_LongStrain_database --subject_list subject1,subject2,subject3 --target_species_list Bifidobacterium_longum,Akkermansia_muciniphila -l PATH/all_subjects.log```

* ```-m reads_assignment``` is the processing mode for step4;
* ```-i PATH/all_subject_results``` is the input of the step5, which is also the input and output of the step5;
* ```-o PATH/all_species_summary``` is the output directory of summary result of target species;
* ```--taxon_species_json Path_to/taxon_species_name.json``` is the absolute path of the database of taxon ID, this file is downloaded with the code and under the same directory with longstrain.py;
* ```--reference_built_path /Path_to_LongStrain_database``` is the absolute path of directory where the species's reference database is built in step3;
* ```--subject_list subject1,subject2,subject3``` is the name of all input subject, separated by ",";
* ```--target_species_list Bifidobacterium_longum,Akkermansia_muciniphila``` is the list of target species separated by ",", in which blank should be replaced by "_";
* ```-l PATH/all_subjects.log``` is the absolute path of log file.

**Warnning:** LongStrain works on longitudinal or concurrent samples, Does NOT work on a single sample.


**Please contact boyanzhou1992@gmail.com for any questions or bug reporting.**
