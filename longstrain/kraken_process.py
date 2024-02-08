"""
This script processes raw reads data by Kraken
Must include "--use-mpa-style --report-zero-counts" in parameters

Update log:
Date: 2021/06/25
Author: Boyan Zhou
1. make it a module in LongStrain
2. version of kraken/2.0.8
"""

import os
import sys


def kraken_process_total(sample_name_list, all_sample_fq_path_list, sample_output_folder_list, kraken_database_path,
                         logger, whether_use_mpa_style=False):
    """
    Get processed fq and relative abundance report by Kraken
    :param sample_name_list: ["sample1", "sample2"]
    :param all_sample_fq_path_list: [["path/x.fq"], ["path/x_R1.fq", "path/x_R2.fq"]]
    :param sample_output_folder_list: ["output_path/subject1/sample1", "output_path/subject1/sample2"]
    :param kraken_database_path: "/gpfs/data/lilab/home/zhoub03/software/kraken2/NCBI_standard"
    :param logger: a logging object
    :return: sample_output.txt, sample_report.txt, for single,
    """

    use_mpa_style = ""
    if whether_use_mpa_style:
        use_mpa_style = "--use-mpa-style "
    for sample_index in range(len(sample_name_list)):
        sample_name = sample_name_list[sample_index]
        sample_fq_path_list = all_sample_fq_path_list[sample_index]
        sample_output_folder = sample_output_folder_list[sample_index]

        # paired-end or single-end
        if len(sample_fq_path_list) == 1:
            # single
            command_kraken2 = f"kraken2 --db {kraken_database_path} --report {sample_output_folder}/{sample_name}_" \
                              f"report.txt {use_mpa_style}--use-names --report-zero-counts --classified-out " \
                              f"{sample_output_folder}/{sample_name}_#.fq {sample_fq_path_list[0]} " \
                              f"--output {sample_output_folder}/{sample_name}_output.txt"
        elif len(sample_fq_path_list) == 2:
            # paied
            command_kraken2 = f"kraken2 --db {kraken_database_path} --report {sample_output_folder}/{sample_name}_" \
                              f"report.txt {use_mpa_style}--use-names --report-zero-counts --paired --classified-out " \
                              f"{sample_output_folder}/{sample_name}_#.fq {sample_fq_path_list[0]} " \
                              f"{sample_fq_path_list[1]} --output {sample_output_folder}/{sample_name}_output.txt"
        else:
            logger.info(f"Error! {sample_name} has incorrect number of fqs: {sample_fq_path_list} Terminated.")
            sys.exit()

        logger.info(f"Processing {sample_name} by Kraken: {command_kraken2}")
        os.system(command_kraken2)
