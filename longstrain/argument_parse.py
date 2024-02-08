#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = ("Boyan Zhou (boyanzhou1992@gmail.com)")
__version__ = '1.1'
__date__ = '01 June 2021'


import argparse
import os


def arg_parsed():
    parser = argparse.ArgumentParser(prog="LongStrain",
                                     description=
                                     f"DESCRIPTION\n"
                                     f" LongStrain version {__version__} ({__date__}): \n"
                                     f" Strain Identification For Longitudinal Microbiome Data\n\n"
                                     f"AUTHORS: {__author__}\n\n"
                                     f"COMMON COMMANDS\n\n"
                                     f"We assume here that Kraken2 and Bowtie2 are in the system path \n"
                                     f"\n========== LongStrain Process ================= \n\n"
                                     f""
                                     f"",
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )

    parser.add_argument("-m", dest="model", required=True, type=str,
                        choices=["fastq_preprocess", "relative_abundance_aggregation", "build_species_reference",
                                 "reads_assignment", "longitudinal_analysis"],
                        help="choice of working model, must be one these five options.")
    parser.add_argument('-i', dest="inputdir", type=str, help="Directory of input files")
    parser.add_argument('-o', dest="outputdir", type=str, help=f"Directory of output files")
    parser.add_argument('-s', "--subject", dest="subject", type=str, help="Name of an input subject")
    parser.add_argument('-p', "--prefix", dest="prefix", type=str,
                        help="Prefix of longitudinal samples, separated by ','", metavar="Example: sample_t0,sample_t1")

    parser.add_argument('-f', "--fastq", dest="fastq", type=str,
                        help="Fastqs corresponding to longitudinal samples\n"
                             "fastqs from different samples separated by ':'\n"
                             "paired fastqs from one sample separated by ','",
                        metavar="Example: sample_t0.fq.gz:sample_t1_R1.fq.gz,sample_t1_R2.fq.gz")
    parser.add_argument('-l', "--logfile", dest="logfile", required=True, type=str, help="Absolute path of log file",
                        metavar="Example: PATH/example.log")
    parser.add_argument("--target_species_list", dest="target_species_list", type=str,
                        help="List of target species, ' ' must be replaced by '_', different species are separated by ','",
                        metavar="Example: Bifidobacterium_longum,Bacteroides_sp._A1C1,Lactococcus_lactis,Bifidobacterium_breve")

    parser.add_argument("--target_species_file", type=str, dest="target_species_file",
                        help="Path of text file storing target species, can have multiple lines, different species are separated by ','",
                        metavar="Example: Bifidobacterium longum,Bacteroides sp. A1C1,Lactococcus lactis,Bifidobacterium breve")

    parser.add_argument('-r', "--reference", type=str, help="Reference Genome: fasta file with Bowtie2 index",
                        metavar="Example: PATH/Bifidobacterium_longum.fas")
    # parser.add_argument("--reference", required=True, type=str, help="Reference Genome: fasta file with index",  metavar="PATH/reference.fas")

    # option group "Fastq Preprocess"
    group1 = parser.add_argument_group("Fastq Preprocess", "Process the raw reads by Kraken")
    arg = group1.add_argument
    arg('--kraken_database', dest="kraken_database", type=str, help="The path of Kraken database")

    # option group "Relative Abundance Aggregation"
    group2 = parser.add_argument_group("Relative Abundance Aggregation", "Combine Relative Abundance Of All Sample")
    arg = group2.add_argument
    arg('--subject_list', dest="subject_list", type=str, help="Names of all target subjects, separated by ','")
    arg('--subject_sample_file', dest="subject_sample_file", type=str, help="Text file storing subjects and samples for relative abundance aggregation, "
                                                                            "each line store a subject in the format 'subject_name:sample1,sample2'.")
    arg('--output_combined_file', dest="output_combined_file", type=str, help="Full path of output of relative abundance of files. Required in this step.")

    # option group "build_species_database"
    group3 = parser.add_argument_group("Build Species Reference", "Build Genome Reference Database For Target Species")
    arg = group3.add_argument
    arg('--reference_built_path', dest="reference_built_path", type=str, help="Full path where LongStrain reference database should be built.")
    arg('--GTDB_dict', dest="GTDB_dict", type=str, help="Full path of GTDB_dict_json, use this option only if you use the GTDB of Kraken2.")
    arg('--GTDB_rep_genome_dir', dest="GTDB_rep_genome_dir", type=str, help="The directory where the representative genomes of GTDB are uncompressed, use this option only if you use the GTDB of Kraken2.")

    # option group "reads_assignment"
    group4 = parser.add_argument_group("Reads Assignment", "Assign reads produced by Kraken to target species")
    arg = group4.add_argument
    arg('--NCBI_taxon_species_json', dest="NCBI_taxon_species_json", type=str, help="Full path of json file storing taxon and species name in NCBI database.\n"
                                                                                    "Only use NCBI_taxon_species_json or GTDB_taxon_species_json")
    arg('--GTDB_taxon_species_json', dest="GTDB_taxon_species_json", type=str, help="Full path of json file storing taxon and species name in GTDB database.\n"
                                                                                    "Only use NCBI_taxon_species_json or GTDB_taxon_species_json")
    # option group "longitudinal_analysis"
    # group5 = parser.add_argument_group("Longitudinal Analysis", "Analyze longitidinal assigned reads by previous steps")
    # arg = group5.add_argument
    # arg('--taxon_species_json', dest="taxon_species_json", type=str, help="Full path of json file storing taxon and species name.")

    args = parser.parse_args()
    return args


def file_exists(file):
    # type for checking file exists
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError(f"{file} does not exist.")
    return file


