"""
### Construct match table between taxid, species name, and representative reference name in GTDB ###
Warning: This is only for GTDB, not for NCBI Ref Seq

### Software dependency ###
python/3.6.5
bowtie2/2.3.5.1
samtools/1.9

### Update log ###
Date: 2023/04/28
Author: Boyan Zhou
1. Initiation
"""


import pandas as pd
import os
import json


def match_taxid_species_name(kraken2_report_path):
    """
    build dict for taxid and species name in Kraken2-GTDB, using a Kraken2 report (original format)
    :param kraken2_report_path: kraken2_report_path = "C:/Users/boyan/OneDrive - NYU Langone Health/research/SNV_CNV_microbiome/paper/manuscript/submission_mSystems_resubmit/database_update/Bacteroides_fragilis_repeat0_t0_report.txt"
    :return: {'14325': 'Bacteroides fragilis', '14326': 'Bacteroides fragilis_A'}
    """
    taxid_species_name_dict = {}
    with open(kraken2_report_path, "r") as report_f:
        for line in report_f:
            # cols = ['78.88', '210114', '210114', 'S', '14325', '              Bacteroides fragilis']
            cols = line.strip().split("\t")
            # only species level information
            if cols[3] == "S":
                taxid_species_name_dict.update({cols[4]: cols[5].strip()})
    return taxid_species_name_dict


def match_species_name_representative_genome(sp_clusters_path):
    """
    Only collect species level information, build dict for species name and representative genome id
    Currently use r89 of GTDB, uncompress the representative genomes to get two folders "archaea" and "bacteria"
    GB_GCA_000007185.1_genomic.fna.gz, or RS_GCF_000591055.1_genomic.fna.gz, UBA7939_genomic.fna.gz
    :param sp_clusters_path: sp_clusters_path = "C:/Users/boyan/OneDrive - NYU Langone Health/research/SNV_CNV_microbiome/paper/manuscript/submission_mSystems_resubmit/database_update/sp_clusters_r89.tsv"
    :return: {"Staphylococcus aureus": "RS_GCF_001027105.1"}
    """
    species_name_rep_genome_dict = {}
    with open(sp_clusters_path, "r") as sp_clusters_f:
        for line in sp_clusters_f:
            cols = line.strip().split("\t")
            # cols[:2] = ["RS_GCF_001027105.1",	"s__Staphylococcus aureus"]
            if not cols[1].startswith("s__"):
                # skip the line not species level
                print(cols[:2])
                continue
            species_name = cols[1][3:]
            rep_genome_id = cols[0]     # like "GB_GCA_002400405.1"
            species_name_rep_genome_dict.update({species_name: rep_genome_id})
    # check the number of representative genomes under release 89 of GTDB
    # find /gpfs/data/lilab/home/zhoub03/GTDB/release89/archaea -type f | wc -l     # 1248
    # find /gpfs/data/lilab/home/zhoub03/GTDB/release89/bacteria -type f | wc -l    # 23458
    # The number of rep genomes summarized from sp_clusters_r89.tsv is 24706, exactly match with above
    # add "_genomic.fna.gz" where searching for the genome
    return species_name_rep_genome_dict


""" create json of GTDB """
kraken2_report_path = "C:/Users/boyan/OneDrive - NYU Langone Health/research/SNV_CNV_microbiome/paper/manuscript/submission_mSystems_resubmit/database_update/Bacteroides_fragilis_repeat0_t0_report.txt"
sp_clusters_path = "C:/Users/boyan/OneDrive - NYU Langone Health/research/SNV_CNV_microbiome/paper/manuscript/submission_mSystems_resubmit/database_update/sp_clusters_r89.tsv"
GTDB_dict_json_path = "C:/Users/boyan/OneDrive - NYU Langone Health/research/SNV_CNV_microbiome/paper/manuscript/submission_mSystems_resubmit/database_update/gtdb_dict.json"

# match taxid and species_name
taxid_species_name_dict = match_taxid_species_name(kraken2_report_path)
# taxid_species_name_dict_json = json.dumps(taxid_species_name_dict)

# match species_name and representative_genome
species_name_rep_genome_dict = match_species_name_representative_genome(sp_clusters_path)
# species_name_rep_genome_dict_json = json.dumps(species_name_rep_genome_dict)

with open(GTDB_dict_json_path, 'w') as gtdb_json:
    json.dump({"dict1": taxid_species_name_dict, "dict2": species_name_rep_genome_dict}, gtdb_json)


