"""
This script aims to assign reads to target species

Update log:
Date: 2021/06/19
Author: Boyan Zhou
1. complete the whole script pipeline, new function: single_reads_assign for single end reads
"""

import os
import json
import logging.handlers
import sys


def paired_reads_assign(prefix, paired_fq_path, species_name_list, species_taxid_list, taxid_species_taxid_dict,
                        output_path, db_name):
    """
    Assign paired reads of Kraken2 results that fall into species_taxid we need
    :param prefix: prefix of sample and output
    :param paired_fq_path: list of length 2, ["path/prefix__1.fq", "path/prefix__2.fq"]
    :param species_name_list:
    :param species_taxid_list: the same order as species_name_list
    :param taxid_species_taxid_dict:
    :param output_path: collect taxid of sub-species level to species level
    :param db_name: must be NCBI or GTDB
    :return:
    """
    # check exist of paired fqs
    if os.path.exists(paired_fq_path[0]) and os.path.exists(paired_fq_path[1]):
        fq1 = open(paired_fq_path[0], "r")
        fq2 = open(paired_fq_path[1], "r")
    else:
        print(f"{prefix}__1.fq does not exist! Skip this sample.")
        return

    species_taxid_out = {}
    for species_taxid, species_name in zip(species_taxid_list, species_name_list):
        species_name = "_".join(species_name.split(" "))
        print(prefix, species_taxid, species_name)
        # taxid_fq.update({species_taxid: (f"{srr}_{species_taxid}_1.fq", f"{srr}_{species_taxid}_2.fq")})
        species_taxid_out.update({species_taxid: [
            open(os.path.join(output_path, f"{prefix}_{species_name}_{species_taxid}_1.fq"), "w"),
            open(os.path.join(output_path, f"{prefix}_{species_name}_{species_taxid}_2.fq"), "w")]})

    while True:
        line1 = fq1.readline()
        if not line1:
            break
        line2 = fq2.readline()
        taxid1 = line1.rstrip().split("|")[-1]
        taxid2 = line2.rstrip().split("|")[-1]

        if db_name == "GTDB":
            if taxid1 == taxid2:
                try:
                    species_taxid_out[taxid1][0].write(line1)
                    species_taxid_out[taxid1][1].write(line2)
                    for i in range(3):
                        species_taxid_out[taxid1][0].write(fq1.readline())
                        species_taxid_out[taxid1][1].write(fq2.readline())
                except KeyError:
                    # taxid1 not in the species_taxid_list
                    for i in range(3):
                        fq1.readline()
                        fq2.readline()
            else:
                for i in range(3):
                    fq1.readline()
                    fq2.readline()

        else:
            # db_name == "NCBI"
            try:
                # reads from fq1 and fq2 belong to the same species, else deprecated
                read1_species_taxid = taxid_species_taxid_dict[taxid1]
                read2_species_taxid = taxid_species_taxid_dict[taxid2]
                if read1_species_taxid == read2_species_taxid:
                    species_taxid_out[read1_species_taxid][0].write(line1)
                    species_taxid_out[read1_species_taxid][1].write(line2)
                    for i in range(3):
                        species_taxid_out[read1_species_taxid][0].write(fq1.readline())
                        species_taxid_out[read1_species_taxid][1].write(fq2.readline())
                else:
                    for i in range(3):
                        fq1.readline()
                        fq2.readline()
            except KeyError:
                for i in range(3):
                    fq1.readline()
                    fq2.readline()

    for fq_files in species_taxid_out.values():
        fq_files[0].close()
        fq_files[1].close()
    fq1.close()
    fq2.close()


def single_reads_assign(prefix, single_fq_path, species_name_list, species_taxid_list, taxid_species_taxid_dict,
                        output_path, db_name):
    """
    Assign paired reads of Kraken2 results that fall into species_taxid we need
    :param prefix: prefix of sample and output
    :param single_fq_path: "path/prefix_#.fq"
    :param species_name_list:
    :param species_taxid_list: the same order as species_name_list
    :param taxid_species_taxid_dict:
    :param output_path: collect taxid of sub-species level to species level
    :param db_name: must be NCBI or GTDB
    :return:
    """
    if os.path.exists(single_fq_path):
        fq1 = open(single_fq_path, "r")
    else:
        print(f"{single_fq_path} does not exist! Skip this sample.")
        return

    species_taxid_out = {}
    for species_taxid, species_name in zip(species_taxid_list, species_name_list):
        species_name = "_".join(species_name.split(" "))
        print(species_taxid, species_name)
        # taxid_fq.update({species_taxid: (f"{srr}_{species_taxid}_1.fq", f"{srr}_{species_taxid}_2.fq")})
        species_taxid_out.update(
            {species_taxid: open(os.path.join(output_path, f"{prefix}_{species_name}_{species_taxid}.fq"), "w")})

    while True:
        line1 = fq1.readline()
        if not line1:
            break
        taxid1 = line1.rstrip().split("|")[-1]
        if db_name == "GTDB":
            try:
                species_taxid_out[taxid1].write(line1)
                for i in range(3):
                    species_taxid_out[taxid1].write(fq1.readline())
            except KeyError:
                # read1_species_taxid not in species_taxid_list
                for i in range(3):
                    fq1.readline()
        else:
            # db_name == "NCBI"
            try:
                # reads from fq1 belong to species_taxid_list
                read1_species_taxid = taxid_species_taxid_dict[taxid1]
                species_taxid_out[read1_species_taxid].write(line1)
                for i in range(3):
                    species_taxid_out[read1_species_taxid].write(fq1.readline())
            except KeyError:
                # read1_species_taxid not in species_taxid_list
                for i in range(3):
                    fq1.readline()

    for fq_file in species_taxid_out.values():
        fq_file.close()
    fq1.close()
