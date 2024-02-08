"""
Build a bowtie reference for a given species from GTDB database
Software dependency
module add python/cpu/3.6.5
module add bowtie2/2.3.5.1
module add samtools/1.9

Update log:
Date: 2023/04/29
Author: Boyan Zhou
1. Complete the whole function
Date: 2023/08/16
1. deal with "_" in species name in GTDB, like "Neisseria meningitidis_B", but input species name "_" replace by " "
"""


import os
import json


# get assembly info of NCBI database
# assembly_summary =
# "/gpfs/data/lilab/home/zhoub03/software/kraken2/NCBI_standard/library/bacteria/assembly_summary.txt"
# assembly = pd.read_table(assembly_summary, skiprows=[0])  # skip the first row of annotation
# assembly_complete_genome = assembly[assembly["assembly_level"] == "Complete Genome"]

def ref_build_gtdb(output_parent_path, target_species_list, GTDB_dict_json_path, GTDB_rep_genome_uncompress_dir, logger):
    """
    :param output_parent_path: path containing species folder, "output_parent_path/species_name"
    :param target_species_list: list of target species, ["Bifidobacterium longum", "Bacteroides uniformis"]
    :param GTDB_dict_json_path: path of GTDB_dict_json
    :param GTDB_rep_genome_uncompress_dir: parent directory where the rep genomes of GTDB were released
    GTDB_rep_genome_uncompress_dir = "/gpfs/data/lilab/home/zhoub03/GTDB/release89"
    :param logger: a logging object
    :return: build a database from a list of target species
    """
    with open(GTDB_dict_json_path, "r") as json_f:
        data = json.load(json_f)
        # taxid_species_name_dict = data["dict1"]
        species_name_genome_dict = data["dict2"]
    # deal with species name like "Neisseria meningitidis_B", but no "_" in input species name
    species_name_genome_dict_add = {species_name.replace("_", " "): species_name_genome_dict[species_name] for
                                    species_name in species_name_genome_dict.keys() if "_" in species_name}
    species_name_genome_dict.update(species_name_genome_dict_add)
    if not os.path.exists(output_parent_path):
        os.system(f"mkdir -p {output_parent_path}")
    gtdb_archaea_dir = os.path.join(GTDB_rep_genome_uncompress_dir, "archaea")
    gtdb_bacteria_dir = os.path.join(GTDB_rep_genome_uncompress_dir, "bacteria")
    gtdb_archaea_genome_list = os.listdir(gtdb_archaea_dir)
    gtdb_bacteria_genome_list = os.listdir(gtdb_bacteria_dir)

    # add bowtie2 to environment variable
    # for each species in the list
    for species_name in target_species_list:
        # # step1: species_name = "Staphylococcus aureus", like "GB_GCA_002400405.1" + "_genomic.fna.gz"
        rep_genome_name = species_name_genome_dict[species_name] + "_genomic.fna.gz"
        if rep_genome_name in gtdb_bacteria_genome_list:
            rep_genome_path = os.path.join(gtdb_bacteria_dir, rep_genome_name)      # if is bacteria
        elif rep_genome_name in gtdb_archaea_genome_list:
            rep_genome_path = os.path.join(gtdb_archaea_dir, rep_genome_name)       # if is archaea
        else:
            logger.info(f"Warning! The rep genome {rep_genome_name} of {species_name} is not in the current GTDB!")
            continue

        # # step2: bowtie2 build
        species_name_joint = species_name.replace(" ", "_")
        output_species_dir = os.path.join(output_parent_path, species_name_joint)
        species_bowtie2_build(species_name_joint, rep_genome_path, output_species_dir, logger)


def species_bowtie2_build(species_name_joint, rep_genome_path, output_species_dir, logger):
    """
    build reference using bowtie2
    :param species_name_joint:
    :param rep_genome_path:
    :param output_species_dir:
    :return:
    """
    if not os.path.exists(output_species_dir):
        os.system(f"mkdir -p {output_species_dir}")
    os.chdir(output_species_dir)
    bowtie_build_file = [i for i in os.listdir(output_species_dir) if i.endswith(".bt2")]
    if len(bowtie_build_file) > 0:
        logger.info(f"There have been some bowtie2 build files under the target folder {output_species_dir}.\n"
                    f"Skip building reference for {species_name_joint}")
    else:
        os.system(f"cp {rep_genome_path} ./")
        reference_gz = os.path.split(rep_genome_path)[1]
        os.system(f"gunzip -c {reference_gz} > {species_name_joint}.fas")
        # os.system(f"gunzip -c {rep_genome_path} > {species_name_joint}.fas")
        logger.info(f"Build bowtie2 reference index for {species_name_joint} using {rep_genome_path} under {output_species_dir}")

        """ check whether bowtie2 has built """
        files_under_build_dir = os.listdir(output_species_dir)
        if f"{species_name_joint}.1.bt2" in files_under_build_dir:
            logger.info(f"Warning! The bowtie reference under {output_species_dir} has been built. Skip rebuilding it!")
        else:
            bowtie2_build_command = f"bowtie2-build {reference_gz} {species_name_joint}"
            logger.info(bowtie2_build_command)
            os.system(bowtie2_build_command)


"""
# GTDB_dict_json_path = "C:/Users/boyan/OneDrive - NYU Langone Health/research/SNV_CNV_microbiome/paper/manuscript/submission_mSystems_resubmit/database_update/gtdb_dict.json"
GTDB_dict_json_path = "/gpfs/data/lilab/home/zhoub03/GTDB/gtdb_r89_dict.json"
GTDB_released_dir = "/gpfs/data/lilab/home/zhoub03/GTDB/release89"      # two folders "archaea" and "bacteria" under it

with open(GTDB_dict_json_path, "r") as json_f:
    data = json.load(json_f)
    taxid_species_name_dict = data["dict1"]
    species_name_rep_genome_dict = data["dict2"]
"""
