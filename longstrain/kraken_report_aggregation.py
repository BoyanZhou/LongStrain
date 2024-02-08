"""
This script is used to transform Kraken2 report to relative abundance format of MetaPhlAn2
Usage example: python3 Kraken2_report_to_RA.py input_path out_path
if kraken2 output all taxa, all reports will have the same length

Update log:
Date: 2021/06/25
1.  Complete the whole function
Date: 2022/02/22
1.  Add the parameter of host_species_to_remove to remove the possible host counts in RA calculation,
    default "Homo sapiens"
2.  Add the function to output
Date: 2022/03/11

"""

import numpy as np


class Kraken2Report:
    number2taxon = {1: "k", 2: "p", 3: "c", 4: "o", 5: "f", 6: "g", 7: "s"}
    rank = ["d", "k", "p", "c", "o", "f", "g", "s"]

    def __init__(self, filename, sample_name, host_species_to_remove=""):
        """
        Get all relative abundance information to self.RA_record
        :param filename: full path of the Kraken2 report
        :param sample_name:
        :param host_species_to_remove: usually the host, "Homo sapiens"
        """
        # usually format of filename is name_report.txt
        self.count_record = dict()
        self.sample = sample_name
        # the result looks like {"d":{"taxon1": 111, "taxon2": 222}, "k":{}, "p":{},
        # "c":{}, "o":{}, "f":{}, "g":{},
        # "s":{"d__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo sapiens" : 99}}
        # --------------------------------------------------------------------------------------------------------------
        # get counts record of each taxon at each rank
        for i in Kraken2Report.rank:
            self.count_record.update({i: {}})
        with open(filename, "r") as file:
            for line in file:
                # kraken2 report use d represent kingdom? why
                cols = line.strip().split("\t")
                taxon_level = cols[0].split("|")[-1][0]     # in ["d", "k", "p", "c", "o", "f", "g", "s"]
                self.count_record[taxon_level].update({cols[0]: int(cols[1])})
        # --------------------------------------------------------------------------------------------------------------
        # calculate relative abundance
        self.library_size = sum(list(self.count_record["d"].values()))
        # remove possible host counts, subtract from library size
        if len(host_species_to_remove) > 0:
            for full_taxon in self.count_record["s"].keys():
                # d__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo sapiens
                if full_taxon.endswith(host_species_to_remove):
                    self.library_size = self.library_size - self.count_record["s"][full_taxon]
        if self.library_size == 0:
            self.library_size = 1   # avoid zero

        self.RA_record = {}     # attention!!! higher rank of host is not removed
        for key0, value0 in self.count_record.items():
            self.RA_record.update({key0: {}})
            for key1, value1 in value0.items():
                self.RA_record[key0].update({key1: value1/self.library_size})
        # pop the host species. attention!!! higher rank of host is not removed
        if len(host_species_to_remove) > 0:
            for full_taxon in self.count_record["s"].keys():
                if full_taxon.endswith(host_species_to_remove):
                    self.RA_record["s"].pop(full_taxon)

    def output_RA(self, out_file_name):
        with open(out_file_name, "w") as out:
            out.write("#SampleID\tKraken2_Analysis\n")
            out.write("#" + "\t".join([j + ":" + str(sum(self.count_record[j].values())) for j in Kraken2Report.rank])
                      + "\n")
            for j in Kraken2Report.rank:
                for taxon, count in self.count_record[j].items():
                    out.write(taxon + "\t" + str(count/self.library_size) + "\n")


class CombinedReportsRA:
    """
    # combine results of many kraken2 reports
    # these reports must have same taxa (orders can be different)
    # for Kraken2, it should include "--use-mpa-style --report-zero-counts" in parameters
    """
    def __init__(self, kraken2report):
        """
        :param kraken2report: an object of Kraken2Report
        """
        self.sample_names = [kraken2report.sample]
        self.result_dict = {}
        # iterate the first sample as a template
        for key0, value0 in kraken2report.RA_record.items():
            self.result_dict.update({key0: {}})
            for key1, value1 in value0.items():
                self.result_dict[key0].update({key1: [value1]})

    def add_report(self, report_added):
        """
        Add other Kraken2report to CombinedReportsRA
        :param report_added:
        :return:
        """
        # check whether reports have the same number of "s"
        if len(self.result_dict["s"]) == len(report_added.RA_record["s"]):
            self.sample_names.append(report_added.sample)
            for key0, value0 in report_added.RA_record.items():
                for key1, value1 in value0.items():
                    self.result_dict[key0][key1].append(value1)
        else:
            print(f"Warnning! Added report does not the same number of species! Need to check it!")

    def report_combined_RA(self, relative_abundance_output_path, taxon_rank="s", mean_RA_threshold=0.0001):
        """
        :param relative_abundance_output_path:
        :param taxon_rank: ["d", "k", "p", "c", "o", "f", "g", "s"]
        :param mean_RA_threshold: not output RA lower than this threshold
        :return:
        """
        RA_target_rank = list(self.result_dict[taxon_rank].items())
        # get relative abundance result of target taxon rank, order by sum of RA, from large to small
        RA_target_rank.sort(key=lambda x: sum(x[1]), reverse=True)
        with open(relative_abundance_output_path, "w") as out_f:
            # the header of RA file
            out_f.write("Taxon\t" + "\t".join(self.sample_names) + "\n")
            for taxon_RA in RA_target_rank:
                # [('d', [-3, -6, -1]), ('c', [-1, -2, 3]), ('a', [1, 2, 4]), ('b', [2, 4, 1])]
                average_RA = np.mean(np.array(taxon_RA[1]))
                if average_RA < mean_RA_threshold:
                    # not output RA lower than the threshold
                    break
                out_f.write(taxon_RA[0] + "\t" + "\t".join(np.array(taxon_RA[1]).astype("str")) + "\n")

##########################################################################
##########################################################################
# #   Main Function: Generate combined Kraken2 report in RA or Count   # #
##########################################################################
##########################################################################


def combine_kraken_result_to_RA(samples_name_process_list, samples_kraken_report_path_process_list,
                                combined_relative_abundance_path, host_species_to_remove="Homo sapiens"):
    """
    Combine the Kraken reports of samples to get relative abundance file.
    Must include "--use-mpa-style --report-zero-counts" in Kraken parameters
    :param samples_name_process_list: [subject1_sample1, subject1_sample2, subject2_sample1]
    :param samples_kraken_report_path_process_list: full path of corresponding kraken reports
    :param combined_relative_abundance_path: full output path of combined relative_abundance
    :param host_species_to_remove: default is "Homo sapiens"
    :return:
    """
    # get the report of the first subject
    sample0_RA_report = Kraken2Report(samples_kraken_report_path_process_list[0], samples_name_process_list[0],
                                      host_species_to_remove)
    combined_reports = CombinedReportsRA(sample0_RA_report)
    # and add others one by one
    if len(samples_name_process_list) > 1:
        for i in range(1, len(samples_name_process_list)):
            print(f"Combine Kraken reports. Processing {samples_name_process_list[i]} ... ... ")
            combined_reports.add_report(Kraken2Report(samples_kraken_report_path_process_list[i],
                                                      samples_name_process_list[i], host_species_to_remove))

    # output combined RA above the threshold
    combined_reports.report_combined_RA(combined_relative_abundance_path)


def combine_Kraken_counts_report(combined_raw_counts_output_path, input_reports_abs_path_list, sample_names_list):
    """
    The order of taxa is determined by the first sample in the list. All reports must have the same length of taxa.
    Collect the counts of all taxa of each sample's report in the dict: all_sample_combined_counts_dict
    Then output the counts of all samples in the order of taxa order in the first sample
    :param combined_raw_counts_output_path: path of combined output reports
    :param input_reports_abs_path_list: by default the length should be > 1, reports need to be combined
    :param sample_names_list: should have the same len with input_reports_abs_path_list
    :return:
    """
    if len(input_reports_abs_path_list) < 2:
        print("The number of reports that are needed to be combined is smaller that 2.\n"
              "We do not need to combine it.")
    elif len(input_reports_abs_path_list) != len(sample_names_list):
        print("Error! The number of reports does not equal to the number of sample names.")
    else:
        # store the order of all taxa, and build the original dict
        taxa_list = []
        all_sample_combined_counts_dict = {}    # {taxon: [sample1_count, sample2_count]}
        with open(input_reports_abs_path_list[0], "r") as first_report:
            for line in first_report:
                cols = line.split("\t")
                taxa_list.append(cols[0])
                all_sample_combined_counts_dict.update({cols[0]: []})

        sample_exist_list = []  # samples that have all same taxon
        for report_abs_path, sample_name in zip(input_reports_abs_path_list, sample_names_list):
            #####################################
            # get the count info of each report #
            #####################################
            print(f"Processing {report_abs_path} ... ...")
            count_dict_current_sample = {}
            with open(report_abs_path, "r") as report_f:
                for line in report_f:
                    cols = line.strip().split("\t")
                    count_dict_current_sample.update({cols[0]: cols[1]})
            # if the taxa match exactly
            if len(count_dict_current_sample) == len(taxa_list):
                sample_exist_list.append(sample_name)
                for taxon in taxa_list:
                    all_sample_combined_counts_dict[taxon].append(count_dict_current_sample[taxon])
            else:
                ###############################################
                # omit the sample does not exactly same taxon #
                ###############################################
                print(f"{sample_name} does not have same taxon. Omit it!")
                continue

        # output the results of all sample
        with open(combined_raw_counts_output_path, "w") as counts_out_f:
            counts_out_f.write("Taxon\t" + "\t".join(sample_exist_list) + "\n")
            # each ordered taxa
            for taxon_name in taxa_list:
                counts_out_f.write(taxon_name + "\t" + "\t".join(all_sample_combined_counts_dict[taxon_name]) + "\n")
