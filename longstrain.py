#!/usr/bin/python3
# -*- coding: utf-8 -*-


""" Copyright (c) 2019  Boyan Zhou

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.


:Authors: Boyan Zhou
:Contact: boyanzhou1992@gmail.com
:Date: Aug 2021
:Type: tool
:Input: NGS of longitudinal metagenomic data
:Output: vcf and estimations of proportions of identified strains

------------------------------------------------------------------------------------------------------------------------
************************
* Version of software: *
************************
Kraken2         Version: 2.0.8 (test on 2.0.8)
bowtie2         Version: 2.3.5.1
samtools        Version: 1.9

------------------------------------------------------------------------------------------------------------------------
**************
* Update log *
**************
Date:   2021/06/23
By:     Boyan Zhou
Change:
"""


import longstrain
import logging
import time
import sys
import os
import json


def get_subject_sample_list(subject_sample_file_path):
    """
    Get subject list and sample list from a given txt file, like below:
    subject1:sample1_1,sample1_2,sample1_3
    subject2:sample2_1,sample2_2,sample2_3
    :param subject_sample_file_path: the full path of info file
    :return:
    """
    with open(subject_sample_file_path) as subject_sample_f:
        subject_list = []
        sample_list = []
        for line in subject_sample_f:
            cols = line.split(":")
            subject_list.append(cols[0].strip())
            sample_list.append([i.strip() for i in cols[1].split(",") if len(i.strip()) > 0])
    return subject_list, sample_list


def main():
    start_time = time.time()
    # get parsed arguments
    options = longstrain.argument_parse.arg_parsed()

    ################
    # set log file #
    ################
    logging.basicConfig(filename=options.logfile, format='%(asctime)s\t%(levelname)s\t%(name)s: %(message)s',
                        level=logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    logging.getLogger().addHandler(handler)

    my_logger = logging.getLogger("main")
    my_logger.info("Started with the command: " + " ".join(sys.argv) + "\n")

    # *************************
    # step1: fastq_preprocess *
    # *************************
    if options.model == "fastq_preprocess":
        subject_name = options.subject
        sample_name_list = [i.strip() for i in options.prefix.split(",")]
        fq_name_list = [i.split(",") for i in options.fastq.split(":")]
        # check whether the number of samples matches fqs
        if len(sample_name_list) != len(fq_name_list):
            print(f"Error! The number of samples does not match fastqs. Terminated!")
            my_logger.info(f"Error! The number of samples does not match fastqs. Terminated!")
            sys.exit()

        # get the abs path of all fqs and check the existence of fqs, likes:
        # [["path/x.fq"], ["path/x_R1.fq", "path/x_R2.fq"]]
        all_sample_fq_path_list = []
        for fq_names in fq_name_list:
            # fq_names = ["x.fq"] or ["x_R1.fq", "x_R2.fq"]
            one_sample_fq_path_list = [os.path.join(options.inputdir, i.strip()) for i in fq_names]
            for one_sample_fq_path in one_sample_fq_path_list:
                if not os.path.exists(one_sample_fq_path):
                    print(f"Error! {one_sample_fq_path} does not exist. Terminated!")
                    my_logger.info(f"Error! {one_sample_fq_path} does not exist. Terminated!")
                    sys.exit()
            all_sample_fq_path_list.append(one_sample_fq_path_list)

        # create folders under output directory, "output_path/"
        if not os.path.exists(options.outputdir):
            os.system(f"mkdir {options.outputdir}")
        # subject folder, "output_path/subject1"
        subject_folder = os.path.join(options.outputdir, subject_name)
        if not os.path.exists(subject_folder):
            os.system(f"mkdir {subject_folder}")
        # sample folder, ["output_path/subject1/sample1", "output_path/subject1/sample2"]
        sample_output_folder_list = [os.path.join(subject_folder, sample_name) for sample_name in sample_name_list]
        for sample_output_folder in sample_output_folder_list:
            os.system(f"mkdir {sample_output_folder}")     # output folder
        longstrain.kraken_process.kraken_process_total(sample_name_list, all_sample_fq_path_list,
                                                       sample_output_folder_list, options.kraken_database, my_logger)

    # ***************************************
    # step2: relative_abundance_aggregation *
    # ***************************************
    elif options.model == "relative_abundance_aggregation":
        # the input path should be the output path of "fastq_preprocess"
        # ******************************************
        # get subject_name_list, samples_name_list *
        # ******************************************
        # check whether these subjects have been processed by Kraken, existence of sample1_report.txt
        if options.subject_sample_file:
            # get subject_name_list, samples_name_list from text info file
            if not os.path.exists(options.subject_sample_file):
                print(f"Error! Can not find {options.subject_sample_file}. Terminated!")
                my_logger.info(f"Error! Can not find {options.subject_sample_file}. Terminated!")
                sys.exit()
            subject_name_list, samples_name_list = get_subject_sample_list(options.subject_sample_file)
        else:
            # find samples in subject folders "output_path/subject1"
            subject_name_list = options.subject_list.split(",")
            samples_name_list = []
            for subject_name in subject_name_list:
                subject_folder = os.path.join(options.inputdir, subject_name)
                samples_name_list.append(sorted([i for i in os.listdir(os.path.join(subject_folder)) if
                                                 os.path.isdir(os.path.join(subject_folder, i))]))

        # **********************
        # get input parameters *
        # **********************
        subject_name_process_list = []  # subject1, subject2, subject3
        samples_name_process_list = []  # [subject1_sample1, subject1_sample2, subject2_sample1]
        samples_kraken_report_path_process_list = []
        # for each subject and each sample, check the existence of
        for subject_name, samples_name in zip(subject_name_list, samples_name_list):
            # subject_name1, [sample1_name, sample2_name]
            samples_name_with_report = []   # list of samples which has "sample_name_report.txt"
            samples_kraken_report_path = []
            samples_name_without_report = []
            for sample_name in samples_name:
                kraken_report_path = os.path.join(options.inputdir, subject_name, sample_name,
                                                  f"{sample_name}_report.txt")
                if os.path.exists(kraken_report_path):
                    samples_name_with_report.append(sample_name)
                    samples_kraken_report_path.append(kraken_report_path)
                else:
                    samples_name_without_report.append(sample_name)

            if len(samples_name_without_report) > 0:
                # record samples without kraken report
                my_logger.info(
                    f"Warning! For {subject_name}, can not find the kraken report of {samples_name_without_report}.")
            if len(samples_name_with_report) > 0:
                # record samples with kraken report
                subject_name_process_list.append(subject_name)
                samples_name_process_list.extend([f"{subject_name}_{i}" for i in samples_name_with_report])
                samples_kraken_report_path_process_list.extend(samples_kraken_report_path)
            else:
                my_logger.info(f"Warning! There is no kraken report of {subject_name}.")
        if len(subject_name_process_list) == 0:
            print(f"Error! Can not find kraken reports for all subject. Terminated!")
            my_logger.info(f"Error! Can not find kraken reports for all subject. Terminated!")
            sys.exit()
        else:
            my_logger.info(f"Combine Kraken reports for the following samples {samples_name_process_list}")
            longstrain.kraken_report_aggregation.combine_kraken_result_to_RA(samples_name_process_list,
                                                                             samples_kraken_report_path_process_list,
                                                                             options.output_combined_file)

    # ********************************
    # step3: build_species_reference *
    # ********************************
    elif options.model == "build_species_reference":
        target_species_list = []
        if options.target_species_list:
            target_species_list = [i.replace("_", " ") for i in options.target_species_list.split(",") if len(i) > 0]
        elif options.target_species_file:
            target_species_list = get_subject_sample_list(options.target_species_file)
        else:
            print(f"Error! Please specify the target species by --target_species_list or --target_species_file. "
                  f"Terminated")
            my_logger.info(f"Error! Please specify the target species by --target_species_list or "
                           f"--target_species_file. Terminated")

        # assembly_summary_path = "/gpfs/data/lilab/home/zhoub03/software/kraken2/NCBI_standard/
        # library/bacteria/assembly_summary.txt", we use bacteria, archaea and viral
        # assembly_summary_path = os.path.join(options.kraken_database, "library", "bacteria", "assembly_summary.txt")
        longstrain.reference_build.ref_build(options.reference_built_path, target_species_list, options.kraken_database,
                                             my_logger)

    elif options.model == "reads_assignment" or options.model == "longitudinal_analysis":
        # the input path should be the output path of "fastq_preprocess"
        # ******************************************
        # get subject_name_list, samples_name_list *
        # ******************************************
        # check whether these subjects have been processed by Kraken, existence of sample1_report.txt
        if options.subject_sample_file:
            # get subject_name_list, samples_name_list from text info file
            if not os.path.exists(options.subject_sample_file):
                print(f"Error! Can not find {options.subject_sample_file}. Terminated!")
                my_logger.info(f"Error! Can not find {options.subject_sample_file}. Terminated!")
                sys.exit()
            subject_name_list, samples_name_list = get_subject_sample_list(options.subject_sample_file)
        else:
            if options.subject_list:
                # if given subject_list, find samples in subject folders "input_path/subject1"
                subject_name_list = options.subject_list.split(",")
            else:
                # if not given subject_list, process all subject_names under the input directory
                subject_name_list = sorted(
                    [i for i in os.listdir(options.inputdir) if os.path.isdir(os.path.join(options.inputdir, i))])

            samples_name_list = []
            for subject_name in subject_name_list:
                subject_folder = os.path.join(options.inputdir, subject_name)
                samples_name_list.append(sorted([i for i in os.listdir(os.path.join(subject_folder)) if
                                                 os.path.isdir(os.path.join(subject_folder, i))]))
                print(f"Get {','.join(samples_name_list[-1])} under the folder {subject_folder}.")
                my_logger.info(f"Get {','.join(samples_name_list[-1])} under the folder {subject_folder}.")

        # ********************
        # get target species *
        # ********************
        target_species_list = []
        if options.target_species_list:
            target_species_list = [i.replace("_", " ") for i in options.target_species_list.split(",") if len(i) > 0]
        elif options.target_species_file:
            target_species_list = get_subject_sample_list(options.target_species_file)
        else:
            print(f"Error! Please specify the target species by --target_species_list or --target_species_file. "
                  f"Terminated")
            my_logger.info(f"Error! Please specify the target species by --target_species_list or "
                           f"--target_species_file. Terminated")

        # ******************************
        # get species names and taxids *
        # ******************************
        # get species comparison table
        taxid_species_taxid, taxid_species_name, species_name_taxid = json.load(
            open(options.taxon_species_json, encoding="utf-8"))
        my_logger.info(f"All target species: {target_species_list}")
        species_taxids = []  # get its species taxids
        species_names_with_taxids = []  # some species does not have taxids will be skipped
        for species_name in target_species_list:
            try:
                species_taxids.append(species_name_taxid[species_name])
                species_names_with_taxids.append(species_name)
            except KeyError:
                my_logger.info(f"Warning! {species_name} does not have species taxid!")
                continue
        if len(species_names_with_taxids) == 0:
            print(f"Error! No species in the input has species taxid! Can't assign reads.")
            my_logger.info(f"Error! No species in the input has species taxid! Can't assign reads.")
            sys.exit()
        # species_names_with_taxids  # only species have taxids included in the following analysis
        # species_taxids

        # *************************
        # step4: reads_assignment *
        # *************************
        if options.model == "reads_assignment":
            my_logger.info("Assigning Reads ... ...")
            # already get subject_name_list, samples_name_list
            for subject_name, samples_name in zip(subject_name_list, samples_name_list):
                # samples_name is a list, ["sample_t0", "sample_t1", "sample_t2"]
                for sample_name in samples_name:
                    # **************************************
                    # get fq path under sample's directory *
                    # **************************************
                    sample_dir = os.path.join(options.inputdir, subject_name, sample_name)
                    assigned_reads_dir = os.path.join(sample_dir, "assigned_reads")
                    if not os.path.exists(sample_dir):
                        my_logger.info(f"Warning! {sample_dir} does not exist! Skip the sample {sample_name}.")
                        continue
                    # make directory for assigned reads
                    if not os.path.exists(assigned_reads_dir):
                        os.system(f"mkdir {assigned_reads_dir}")
                        my_logger.info(f"mkdir {assigned_reads_dir}")

                    if f"{sample_name}_#.fq" in os.listdir(sample_dir):
                        # single end fq
                        single_fq_path = os.path.join(sample_dir, f"{sample_name}_#.fq")
                        longstrain.assign_reads.single_reads_assign(sample_name, single_fq_path,
                                                                    species_names_with_taxids, species_taxids,
                                                                    taxid_species_taxid, assigned_reads_dir)
                    elif f"{sample_name}__1.fq" in os.listdir(sample_dir) and f"{sample_name}__2.fq" in os.listdir(sample_dir):
                        # paired end fqs
                        paired_fq_path = [os.path.join(sample_dir, f"{sample_name}__1.fq"),
                                          os.path.join(sample_dir, f"{sample_name}__2.fq")]
                        longstrain.assign_reads.paired_reads_assign(sample_name, paired_fq_path,
                                                                    species_names_with_taxids, species_taxids,
                                                                    taxid_species_taxid, assigned_reads_dir)

        # ******************************
        # step5: longitudinal_analysis *
        # ******************************
        else:
            # options.model == "longitudinal_analysis"
            # ******************************
            # for each species, processing *
            # ******************************
            my_logger.info("LongStrain analyzing ... ...")
            for species_name, species_taxid in zip(species_names_with_taxids, species_taxids):
                species_name_joint = species_name.replace(" ", "_")     # Akkermansia_muciniphila
                # create a species dir in output path
                species_output_path = os.path.join(options.outputdir, species_name_joint)
                if not os.path.exists(species_output_path):
                    os.system(f"mkdir {species_output_path}")

                # get reference path
                species_reference_fas_path = os.path.join(options.reference_built_path, species_name_joint,
                                                          f"{species_name_joint}.fas")
                species_reference_bowtie_path = os.path.join(options.reference_built_path, species_name_joint,
                                                             species_name_joint)

                # ******************************
                # for each subject, processing *
                # ******************************
                for subject_name, samples_name in zip(subject_name_list, samples_name_list):
                    # samples_name is a list, ["sample_t0", "sample_t1", "sample_t2"]
                    sample_name_used_list = []
                    fqs_path_list = []
                    sample_name_with_fq = []
                    for sample_name in samples_name:
                        # **************************************
                        # get fq path under sample's directory *
                        # **************************************
                        sample_dir = os.path.join(options.inputdir, subject_name, sample_name)
                        assigned_reads_dir = os.path.join(sample_dir, "assigned_reads")
                        single_fq_name = f"{sample_name}_{species_name_joint}_{species_taxid}.fq"
                        paired_fq_name = [f"{sample_name}_{species_name_joint}_{species_taxid}_1.fq",
                                          f"{sample_name}_{species_name_joint}_{species_taxid}_2.fq"]
                        if single_fq_name in os.listdir(assigned_reads_dir):
                            sample_name_used_list.append(sample_name)
                            fqs_path_list.append([os.path.join(assigned_reads_dir, single_fq_name)])
                            sample_name_with_fq.append(sample_name)
                        elif paired_fq_name[0] in os.listdir(assigned_reads_dir) and paired_fq_name[1] in \
                                os.listdir(assigned_reads_dir):
                            sample_name_used_list.append(sample_name)
                            fqs_path_list.append([os.path.join(assigned_reads_dir, i) for i in paired_fq_name])
                            sample_name_with_fq.append(sample_name)

                    if len(fqs_path_list) > 1:
                        longstrain.longitudinal_microbiome.fqs_species_process(subject_name, sample_name_with_fq,
                                                                               fqs_path_list,
                                                                               species_reference_fas_path,
                                                                               species_reference_bowtie_path,
                                                                               output_path=species_output_path)
                    else:
                        my_logger.info(f"Warning! {subject_name} does not have enough samples (at least two samples) "
                                       f"for {species_name}, skip.")


if __name__ == "__main__":
    sys.exit(main())
