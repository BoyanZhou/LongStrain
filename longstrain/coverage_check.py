"""
This script aims to check the coverage of a microbial genome by raw fq reads
module add samtools/1.9

Update log:
Date: 2021/06/20
By: Boyan Zhou
1. make it can deal with paired and single fqs at the same time
"""

import os
import pysam


def coverage_single(species_name, species_ref_fas, one_fq_file, output_path, effective_rate=0.8, read_length=100,
                    threshold=5, whether_paired=True):
    """
    Calculate the coverage depth for a single fq pair, given species and reference
    :param species_name:
    :param species_ref_fas: full path of species' reference fas
    :param one_fq_file: sample_1.fq (just one file)
    :param output_path:
    :param effective_rate: proportion of effective reads
    :param read_length: average
    :param threshold: threshold of whether report 1, default is 5X
    :param whether_paired: whether the paired fq file exists (sample_1.fq)
    :return:
    """
    # check fai of species_ref_fas
    species_ref_fas_fai = species_ref_fas + ".fai"
    if not os.path.exists(species_ref_fas_fai):
        os.system(f"samtools faidx {species_ref_fas}")

    # get pysam fasta, example:
    # species_pysam_ref = pysam.FastaFile("/gpfs/data/lilab/home/zhoub03/software/
    # my_strain3/Bacteroides_fragilis/Bacteroides_fragilis.fas")
    species_pysam_ref = pysam.FastaFile(species_ref_fas)
    species_ref_total_length = sum(species_pysam_ref.lengths)
    print(f"The total length of {species_name} is {species_ref_total_length}.")

    # raw coverage rate
    os.system(f"wc -l {one_fq_file} > {os.path.join(output_path, 'row_count_temp.txt')}")
    # os.system(f"wc -l /gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/TEDDY_T1D/1/case/
    # SRR7557199/SRR7557199_Bifidobacterium_breve_1685_1.fq > row_count_temp.txt")
    with open(os.path.join(output_path, 'row_count_temp.txt'), "r") as count_temp:
        row_count = count_temp.readline().split(" ")[0]
    if row_count.isnumeric():
        row_count = int(int(row_count)/4)
        raw_depth = row_count * effective_rate * read_length/species_ref_total_length
        if whether_paired:
            raw_depth = raw_depth * 2
    else:
        raw_depth = 0
    return [raw_depth >= threshold, raw_depth]


def coverage_longitudinal_sample(species_name, species_ref_fas, fq_files, output_path, effective_rate=0.8,
                                 read_length=100, threshold=5):
    """
    Check whether longitudinal samples pass the examination two samples > 5X or total > 10X
    :param species_name:
    :param species_ref_fas:
    :param fq_files: [["x.fq"], ["y_1.fq", "y_2.fq"]]， full path
    :param output_path:
    :param effective_rate:
    :param read_length:
    :param threshold:
    :return:
    """
    longitudinal_samples_coverage_bool = []
    longitudinal_samples_coverage = []
    for fq_file in fq_files:
        # fq_file = ["x.fq"] or ["y_1.fq", "y_2.fq"]; paired or single fq
        # sample_record = [bool, raw_depth]
        if len(fq_file) == 1:
            whether_paired = False
        elif len(fq_file) == 2:
            whether_paired = True
        else:
            print("Warning! " + str(fq_file) + " is neither paired nor single reads! Skip it.")
            continue
        sample_record = coverage_single(species_name, species_ref_fas, fq_file[0], output_path, effective_rate,
                                        read_length, threshold, whether_paired)
        longitudinal_samples_coverage_bool.append(sample_record[0])
        longitudinal_samples_coverage.append(sample_record[1])

    print(longitudinal_samples_coverage_bool)
    print(longitudinal_samples_coverage)
    if longitudinal_samples_coverage_bool.count(True) > 1 or sum(longitudinal_samples_coverage) > 10:
        print(f"{species_name} passes the coverage check.")
        return True
    else:
        print(f"{species_name} does not pass the coverage check.")
        return False
