"""
This script performs strain analysis using our method, for species assigned fq data
It follows the script: multiple_species_reads_assignment.py
module add python/cpu/3.6.5
module add bowtie2/2.3.5.1
module add samtools/1.9

Update log:
Date: 2021/04/21
By: Boyan Zhou
1. make the ID range more clear

Date: 2021/04/26
By: Boyan Zhou
1. add logging function
2. make it can accept ID list， not only ID range

"""

import os
import longitudinal_microbiome
import coverage_check
import json
import sys
import logging
import logging.handlers
import datetime


# This is the whole function for analysis
def all_target_species_for_t1d(parent_path, species_name_list, parent_output_path, database_parent_path, data_set, logger,
                               process_ids):
    """
    longitudinal analysis function for all T1D samples
    :param species_name_list: list of target species
    :param parent_path: parent directory of all id-directory
    :param parent_output_path: parent
    :param database_parent_path:
    :param process_ids: case_control ids that need to be processed, int list, [] means process all IDs
    :param data_set: "T1D" or "T1D_IA_overlap" or "IA"
    :param logger:
    :return:
    """
    # get all target species name, list format, standard names
    # target_species_list = tdfc.get_target_species_name(tdfc.target_species_name_file)

    # get srr number by case-control ID (and case is True, control is False)
    # ****************
    # get needed IDs *
    # ****************
    info_related = tdfc.InfoRelated()
    all_ids = info_related.get_case_control_id_t1d(data_set)    # all IDs in given data set, int list
    if len(process_ids) != 0:
        # it would be an int list
        process_ids_set = set(process_ids)
        all_ids = list(process_ids_set.intersection(set(all_ids)))
    all_ids.sort()
    logger.info(f"Process IDs : {all_ids}")

    for case_control_id in all_ids:
        logger.info(f"Processing {case_control_id}")
        # make result directory of this ID
        if not os.path.exists(os.path.join(parent_output_path, str(case_control_id))):
            os.system(f"mkdir {os.path.join(parent_output_path, str(case_control_id))}")

        # get the directory path of each id, id is int
        id_directory_path = os.path.join(parent_path, str(case_control_id))
        if not os.path.exists(id_directory_path):
            continue

        # ****************
        # get needed SRR *
        # ****************
        case_srr_numbers = info_related.get_srr_by_t1d_case_control_id(case_control_id, data_set, True)
        control_srr_numbers = info_related.get_srr_by_t1d_case_control_id(case_control_id, data_set, False)

        if len(case_srr_numbers) > 0:
            logger.info(f"For {case_control_id} Case, SRR exists.")
            input_parent_path = os.path.join(id_directory_path, "case")
            output_path = os.path.join(parent_output_path, str(case_control_id), "case")
            os.system(f"mkdir {output_path}")
            strain_analysis_for_one_id(case_srr_numbers, species_name_list, input_parent_path, output_path,
                                       case_control_id, database_parent_path, logger)
        else:
            logger.warning(f"Warning! For {case_control_id} Case, SRR does not exists.")

        if len(control_srr_numbers) > 0:
            logger.info(f"For {case_control_id} Control, SRR exists.")
            input_parent_path = os.path.join(id_directory_path, "control")
            output_path = os.path.join(parent_output_path, str(case_control_id), "control")
            os.system(f"mkdir {output_path}")
            strain_analysis_for_one_id(control_srr_numbers, species_name_list, input_parent_path, output_path,
                                       case_control_id, database_parent_path, logger)
        else:
            logger.warning(f"Warning! For {case_control_id} Control, SRR does not exists.")


def strain_analysis_for_longitudinal_samples(fqs_abs_path, species_name_list, input_parent_path, output_path, case_control_id, database_parent_path, logger):
    """
    :param fqs_abs_path: list of fqs' path, may be paired or not, [["x.fq"], ["y_1.fq", "y_2.fq"]]
    :param species_name_list: list of target species
    :param input_parent_path: path/id/case or path/id/control
    :param output_path: path/id/
    :param case_control_id: number
    :param database_parent_path: default /gpfs/data/lilab/home/zhoub03/software/my_strain3
    :return:
    """
    taxid_species_taxid, taxid_species_name, species_name_taxid = json.load(open(
        "/gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/TEDDY_T1D/taxon_species_name.json", encoding="utf-8"))

    #########################################################
    # for each target species, analyze longitudinal samples #
    #########################################################
    # species_name_list = tdfc.get_target_species_name(tdfc.target_species_name_file)
    for species_name in species_name_list:
        logger.info(f"Processing {species_name} ... ...")
        species_taxid = species_name_taxid[species_name]
        species_name_joint = species_name.replace(" ", "_")

        ########################################################
        # get fqs of one species from all longitudinal samples #
        ########################################################
        species_bowtie_ref = os.path.join(database_parent_path, species_name_joint, species_name_joint)
        refs_genome = species_bowtie_ref + ".fas"

        logger.info(species_bowtie_ref)
        logger.info(refs_genome)
        if os.path.exists(refs_genome):
            # if the reference file exists, coverage check
            whether_pass_check = coverage_check.coverage_longitudinal_sample(species_name, refs_genome,
                                                                             [i[0] for i in fqs_abs_path], output_path)
            if whether_pass_check:
                longitudinal_microbiome.fqs_species_process(f"our_method_{case_control_id}_{species_name_joint}",
                                                            fqs_abs_path, refs_genome, species_bowtie_ref, output_path)
        else:
            logger.info(f"The reference of {species_name} is not built, skip!")
            continue

        # record effective SRR number
        with open(os.path.join(output_path, "effective_SRR.txt"), "w") as effective_srr_record:
            effective_srr_record.write(
                str(case_control_id) + "\t" + species_name + "\t" + "\t".join(effective_srr) + "\n")


if __name__ == "__main__":
    my_logger = logging.getLogger('mylogger')
    my_logger.setLevel(logging.DEBUG)
    rf_handler = logging.handlers.TimedRotatingFileHandler(sys.argv[-1], when='midnight', interval=1, backupCount=7,
                                                           atTime=datetime.time(0, 0, 0, 0))
    rf_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    my_logger.addHandler(rf_handler)

    t1d_parent_path = "/gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/TEDDY_T1D"
    output_parent_path = "/gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/TEDDY_T1D/longitudinal_strains_results"
    species_name_list = tdfc.get_target_species_name(tdfc.target_species_name_file)
    database_parent_path = tdfc.our_method_database_path    # LongStrain database, "strain3" here

    my_logger.info("The command line is : " + " ".join(sys.argv))
    my_logger.info("The total parent folder is : " + t1d_parent_path)
    my_logger.info("The output path is : " + output_parent_path)
    my_logger.info("The target species list is : " + "\t".join(species_name_list))
    my_logger.info("The database used is : " + database_parent_path)

    process_mode = sys.argv[1]
    data_set = sys.argv[2]
    if process_mode == "all":
        # python longitudinal_strain_analysis_for_T1D.py all T1D_IA_overlap all TEDDY1_1.log
        # process all IDs
        all_target_species_for_t1d(t1d_parent_path, species_name_list, output_parent_path, database_parent_path,
                                   data_set, my_logger, [])
    elif process_mode == "id_range":
        # python longitudinal_strain_analysis_for_T1D.py id_range T1D_IA_overlap 105,120 TEDDY1_1.log
        # given a ID range, like "0, 20"
        process_ids_range = sys.argv[3]
        start_id, end_id = process_ids_range.split(",")
        process_ids = [i for i in range(int(start_id), int(end_id))]
        all_target_species_for_t1d(t1d_parent_path, species_name_list, output_parent_path, database_parent_path,
                                   data_set, my_logger, process_ids)
    elif process_mode == "id_list":
        # python longitudinal_strain_analysis_for_T1D.py id_list T1D_IA_overlap 1,5,7,10 TEDDY1_1.log
        # given a ID list, like "1, 5, 10"
        process_ids = sys.argv[3]
        process_ids = [int(i) for i in process_ids.split(",")]
        all_target_species_for_t1d(t1d_parent_path, species_name_list, output_parent_path, database_parent_path,
                                   data_set, my_logger, process_ids)
    else:
        my_logger.error(f"{process_mode} is not an eligible mode!")
