# this program identifies two strains from longitudinal microbiome sample from one individual
"""
module add bowtie2/2.3.5.1
module add samtools/1.9

Update log:
Date: 2021/08/01
Author: Boyan Zhou
1. modify strain_initial(), can initial with variant sites
Date: 2022/03/02
1. Fix a small bug caused by np.column_stack(tuple)
Date: 2022/03/28
modify fqs_species_process: not rm temporary bam file and not run LongStrain algorithm
"""

import pysam
import numpy as np
import sys
import longstrain.haplotype as haplotype
import json
import os


class TaxonSpecies:
    def __init__(self):
        taxid_species_taxid, taxid_species_name, species_name_taxid = json.load(open(
            "/gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/TEDDY_T1D/taxon_species_name.json", encoding="utf-8"))
        self.taxid_species_taxid = taxid_species_taxid
        self.taxid_species_name = taxid_species_name
        self.species_name_taxid = species_name_taxid

    def species_taxid_to_species_name(self, input_id):
        if input_id in self.taxid_species_name:
            return self.taxid_species_name[input_id]
        else:
            print(f"ID-{input_id} is not in the database!")
            return None


# get the inter_nodes of chromosome by setting length of bin
def get_nodes(chr_len, bin_len):
    bin_left = [i*bin_len for i in range(chr_len//bin_len + 1)]
    bin_right = [i * bin_len for i in range(1, chr_len // bin_len + 1)] + [chr_len]
    boundary = [[i, j] for i, j in zip(bin_left, bin_right)]
    return boundary


def get_effective_site(bin_results_ref_num, bin_results_alt1_num, bin_results_alt2_num, bin_results_depth):
    """
    # Get effective site index
    :param bin_results_ref_num:
    :param bin_results_alt1_num:
    :param bin_results_alt2_num:
    :param bin_results_depth:
    :return: effective_site, variant_site (bool vector)
    """
    # At least 2 time points with heterozygote
    # effective_site_index: at least two samples with (0.1 < heter_proportion < 0.9)
    ref_num_sum = np.sum(bin_results_ref_num, axis=1)
    alt1_num_sum = np.sum(bin_results_alt1_num, axis=1)
    alt2_num_sum = np.sum(bin_results_alt2_num, axis=1)
    heter_proportion = bin_results_ref_num / (bin_results_depth + 0.0001)

    ######################
    # get effective site #
    ######################
    # mixted ref and alt at more than one time point
    effective_site_index1 = np.sum((heter_proportion > 0.1001) * (heter_proportion < 0.8999),
                                   axis=1) > 1  # vector of d1 like array([False,  True,  True])
    # all ref and all alt at some time points, 0.1 < total allele frequency < 0.9
    allele_frequency = ref_num_sum / (np.sum(bin_results_depth, axis=1) + 0.0001)
    effective_site_index2 = (allele_frequency > 0.1001) * (allele_frequency < 0.8999)
    # return bool vector
    effective_site = (effective_site_index1 + effective_site_index2) * (alt1_num_sum > 1) * (alt2_num_sum < 2)

    ####################
    # get variant site #
    ####################
    # "variant site" old version
    # variant_site_index = np.sum((heter_proportion < 0.8999) * (bin_results_depth > 1), axis=1) > 1
    variant_site_index1 = np.sum((heter_proportion < 0.8999) * (bin_results_depth > 2), axis=1) > 1
    variant_site_index2 = (allele_frequency < 0.8999) * (np.sum(bin_results_depth, axis=1) > 5)
    variant_site = (variant_site_index1 + variant_site_index2) * (alt2_num_sum < 2)
    # make sure effective_site not overlap with variant_site
    variant_site = variant_site * (~effective_site)
    return effective_site, variant_site


def get_nonzero_quantile_bins(depth_array, q=0.8, bin_num=5):
    """
    return nonzero quantile, several bins of depth array
    :param depth_array:
    :param q: quantile
    :param bin_num: number of bins
    :return: default 5 rows with three cols: start, end and depth (all sample combined)
    """
    # depth_array has three cols: start, end and depth (all sample combined)
    depth_nonzero = depth_array[depth_array[:, 2] != 0, ]
    bin_ordered_index = np.argsort(depth_nonzero[:, 2])
    if len(bin_ordered_index) == 0:
        # if there is no nonzero bin
        return None
    else:
        q_index = int(np.rint(q*len(bin_ordered_index))) - 1
        if bin_num % 2 == 1:
            index_start = max(0, q_index-(bin_num-1) // 2)
            index_end = min(len(bin_ordered_index) - 1, q_index + (bin_num - 1) // 2)
        else:
            index_start = max(0, q_index-(bin_num-2) // 2)
            index_end = min(len(bin_ordered_index) - 1, q_index + bin_num // 2)
        return depth_nonzero[bin_ordered_index[index_start: index_end + 1], ]


# subject has multiple microbiome samples at several time points
class Subject:
    def __init__(self, subject_id, samples_name_list, reference):
        self.ID = subject_id
        self.reference = reference  # fasta file
        self.pysam_ref = pysam.FastaFile(self.reference)    # imported fasta file
        self.sample_names = samples_name_list               # sample names
        self.samples = []                                   # pysam files, 1-to-1 with self.sample_names
        self.reads_in_bin_by_chromosome = {}
        self.strain_RA = np.array([])
        # strain_RA_initial need to be initialized in regions with adequate coverage
        self.strain_RA_initial = np.array([])
        self.fqs = []
        self.haplotype_filename = ""

    # add fqs, [fq1, fq2]
    def add_fqs(self, paired_fqs):
        self.fqs.append(paired_fqs)

    # processing of fqs to bams
    def fqs_to_bams(self, bowtie_reference, thread=4, q_threhold=20):
        bam_list = []
        for sample_name, fq_pair in zip(self.sample_names, self.fqs):
            if os.path.exists(f"{sample_name}_q{q_threhold}_sorted_rmdup.bam"):
                bam_list.append(f"{sample_name}_q{q_threhold}_sorted_rmdup.bam")
            else:
                if len(fq_pair) == 2:
                    # pair end fqs
                    bowtie2_map_command = f"bowtie2 -p {thread} -x {bowtie_reference} -1 {fq_pair[0]} -2 {fq_pair[1]} -S {sample_name}.sam"
                else:
                    # single end fq
                    bowtie2_map_command = f"bowtie2 -p {thread} -x {bowtie_reference} -U {fq_pair[0]} -S {sample_name}.sam"
                print(bowtie2_map_command)
                os.system(bowtie2_map_command)
                os.system(f"samtools view -bS -q {q_threhold} {sample_name}.sam > {sample_name}_q{q_threhold}.bam")
                os.system(f"samtools sort {sample_name}_q{q_threhold}.bam -o {sample_name}_q{q_threhold}_sorted.bam")
                os.system(f"samtools rmdup {sample_name}_q{q_threhold}_sorted.bam {sample_name}_q{q_threhold}_sorted_rmdup.bam")
                os.system(f"samtools index {sample_name}_q{q_threhold}_sorted_rmdup.bam")
                # os.system(f"rm {sample_name}.sam {sample_name}_q{q_threhold}.bam {sample_name}_q{q_threhold}_sorted.bam")
                bam_list.append(f"{sample_name}_q{q_threhold}_sorted_rmdup.bam")
        return bam_list

    # add microbiome samples to the Subject according to time points
    def add_bams(self, bam_list):
        for sample_bam in bam_list:
            # self.sample_names.append(sample_name)
            self.samples.append(pysam.AlignmentFile(sample_bam, "rb", reference_filename=self.reference))

    # horizontal genome coverage summary: covered by at least several reads, must after add_bams
    def horizontal_genome_coverage(self):
        chr_name_list = list(self.samples[0].references)
        chr_len_list = list(self.samples[0].lengths)
        horizontal_coverage_f = open(f"{self.ID}_horizontal_genome_coverage.txt", "w")

        # write the header of horizontal_genome_coverage file
        file_header = ["Contig", "Length"]
        for sample_name in self.sample_names:
            file_header.extend([f"{sample_name}_1X", f"{sample_name}_5X", f"{sample_name}_10X"])
        horizontal_coverage_f.write("\t".join(file_header) + "\n")
        n_at_depth_all = []     # each element is a chr
        # record coverage rate in each chromosome
        for chr_name, chr_len in zip(chr_name_list, chr_len_list):
            chr_record = [chr_name, str(chr_len)]
            n_at_depth_chr = []
            for pysam_obj in self.samples:
                # tuple of four bases count, row number is 4, col number is total site
                count_tuple = pysam_obj.count_coverage(chr_name, start=None, stop=None, region=None, quality_threshold=15, read_callback='all')
                count_array = np.array(count_tuple).T
                depth_all_sites = np.sum(count_array, axis=1)   # series, length is the chr length
                n_at_depth = [np.sum(depth_all_sites > 0), np.sum(depth_all_sites > 5), np.sum(depth_all_sites > 10)]
                n_at_depth_chr.extend(n_at_depth)
                # record coverage rate at different depth level for this sample
                chr_record.extend([str(i/chr_len) for i in n_at_depth])
            horizontal_coverage_f.write("\t".join(chr_record) + "\n")
            n_at_depth_all.append(n_at_depth_chr)   # add count record of each chr
        # total coverage rate across all chromosome
        n_at_depth_sum = np.sum(np.array(n_at_depth_all), axis=0)   # count of samples at all chr, Series
        chr_len_sum = sum(chr_len_list)
        chr_sum_record = ["All_Contigs", str(chr_len_sum)]
        chr_sum_record.extend(list((n_at_depth_sum/chr_len_sum).astype(str)))
        horizontal_coverage_f.write("\t".join(chr_sum_record) + "\n")
        horizontal_coverage_f.close()

    # coverage summarize by bins per chromosome
    def coverage_summary(self, bin_len=10000):
        pysam_ref = pysam.FastaFile(self.reference)  # import fasta file
        # for each chromosome of pysam_ref
        for chromosome, chromosome_len in zip(pysam_ref.references, pysam_ref.lengths):
            reads_bin_chromosome = []   # row: bins, col: start, end, reads of all sample
            # for each bin, 10000bp example
            for boundary in get_nodes(chromosome_len, bin_len):
                # get reads of all samples at this bin
                boundary.extend(
                    [pysam_sample.count(contig=chromosome, start=boundary[0], stop=boundary[1]) for pysam_sample in
                     self.samples])
                reads_bin_chromosome.append(boundary)
            # update reads in bins
            self.reads_in_bin_by_chromosome.update({chromosome: np.array(reads_bin_chromosome)})

    def get_pos_read_names(self, effective_pos_read_names_dict, variant_pos_read_names_dict, target_bin, chromosome):
        """
        # Get a dict storing {pos: {"read_names", "0-1", "sample index"}} from one target bin, "read_names", "0-1",
        "sample index" are three numpy array
        # and a dict storing {pos: {0:"A", 1:"G"}}
        :param effective_pos_read_names_dict: {pos: {"read_names", "0-1", "sample index"}} or empty {}
        :param variant_pos_read_names_dict: {pos: {"read_names", "0-1", "sample index"}} or empty {}
        :param target_bin: [start, end]
        :param chromosome: chromosome name
        :return: three dicts, and np.array format
        """
        fasta_bin = self.pysam_ref.fetch(reference=chromosome, start=target_bin[0], end=target_bin[1], region=None)
        fasta_bin = fasta_bin.upper()
        fasta_len = target_bin[1] - target_bin[0]

        # build array to store results of each pos, sample
        bin_results_depth = np.zeros([fasta_len, len(self.samples)], dtype=int)
        bin_results_ref_num = np.zeros([fasta_len, len(self.samples)], dtype=int)
        bin_results_alt1 = np.empty([fasta_len, len(self.samples)], dtype=(str, 100))
        bin_results_alt1_num = np.zeros([fasta_len, len(self.samples)], dtype=int)
        bin_results_alt2 = np.empty([fasta_len, len(self.samples)], dtype=(str, 100))
        bin_results_alt2_num = np.zeros([fasta_len, len(self.samples)], dtype=int)

        # for each sample at different time points, collect variation information
        for sample_i in range(len(self.sample_names)):
            sample_pileup = self.samples[sample_i].pileup(chromosome, target_bin[0], target_bin[1],
                                                          truncate=True, stepper="samtools",
                                                          fastafile=self.pysam_ref)
            # pc_index = -1
            for pc in sample_pileup:
                # add pileup information to pre-defined arrays
                # pc_ref_pos = pc.reference_pos
                pc_index = pc.reference_pos - target_bin[0]
                try:
                    piled_bases = pc.get_query_sequences(mark_matches=True, mark_ends=False, add_indels=True)
                except AssertionError:
                    continue

                # pos_in_ref = boundary[0] + pc_index     # record position in ref
                bin_results_depth[pc_index, sample_i] = len(piled_bases)  # record depth
                # bin_results[pc_index, sample_i, 2] = fasta_bin[pc_index]        # record ref
                # ref_num, alt1, alt1_num , alt2, alt2_num
                bin_results_ref_num[pc_index, sample_i], bin_results_alt1[pc_index, sample_i], bin_results_alt1_num[
                    pc_index, sample_i], bin_results_alt2[pc_index, sample_i], bin_results_alt2_num[
                    pc_index, sample_i] = Subject._parse_piled_bases(piled_bases)

        # get effective site and variant site in the target bin, two bool vector
        effective_site_index, variant_site_index = \
            get_effective_site(bin_results_ref_num, bin_results_alt1_num, bin_results_alt2_num, bin_results_depth)
        total_pos_index = effective_site_index + variant_site_index     # one dim bool np.array

        pos_genotype_base_dict = {}  # dict of geno to base at each site, like {pos: {0:"A", 1:"T"}}
        # if there are effective sites
        if np.sum(total_pos_index) > 0:
            ############################################################################################################
            """ collect information of effective sites and variant sites,  {pos: {0:"A", 1:"G"}} """
            effective_pos = np.where(effective_site_index)[0] + target_bin[0]
            variant_pos = np.where(variant_site_index)[0] + target_bin[0]
            total_pos = np.where(total_pos_index)[0] + target_bin[0]  # pos in ref
            # print("effective_pos is")
            # print(effective_pos)
            # print("variant_pos is")
            # print(variant_pos)
            # print("total_pos is")
            # print(total_pos)
            for pos_i in total_pos:
                pos_i_original = pos_i - target_bin[0]
                alt_base = Subject._get_alt_base(bin_results_alt1[pos_i_original])
                if alt_base is None:
                    continue
                pos_genotype_base_dict.update({pos_i: {0: fasta_bin[pos_i_original], 1: alt_base}})
                # build dict for effective and variant respectively
                if pos_i in effective_pos:
                    effective_pos_read_names_dict.update(
                        {pos_i: [np.array([], dtype=str), np.array([], dtype=int), np.array([], dtype=int)]})
                elif pos_i in variant_pos:
                    variant_pos_read_names_dict.update(
                        {pos_i: [np.array([], dtype=str), np.array([], dtype=int), np.array([], dtype=int)]})
            # update effective_pos and variant_pos
            effective_pos = list(effective_pos_read_names_dict.keys())
            variant_pos = list(variant_pos_read_names_dict.keys())
            ############################################################################################################
            """ store read_names, "0-1" and sample index, which combine all samples's reads together """
            for index, sample_i in enumerate(self.samples):
                # get the pileup of target bin
                sample_pileup = sample_i.pileup(chromosome, target_bin[0], target_bin[1], truncate=True,
                                                stepper="samtools", fastafile=self.pysam_ref)
                # pc_index = -1
                for pc in sample_pileup:
                    if pc.reference_pos in total_pos:
                        if len(pc.get_query_names()) == 0:  # no reads
                            continue
                        # if pc.reference_pos not in pos_read_names_dict: continue

                        # store 0-1: whether it is a mutation
                        try:
                            piled_bases = np.array(pc.get_query_sequences(mark_matches=True, mark_ends=False, add_indels=True))
                        except AssertionError:
                            continue
                        whether_mutation = np.full(len(piled_bases), 1)
                        whether_mutation[piled_bases == "."] = 0
                        whether_mutation[piled_bases == ","] = 0
                        if pc.reference_pos in effective_pos:
                            # store reads names
                            effective_pos_read_names_dict[pc.reference_pos][0] = np.append(
                                effective_pos_read_names_dict[pc.reference_pos][0], np.array(pc.get_query_names()))
                            effective_pos_read_names_dict[pc.reference_pos][1] = np.append(
                                effective_pos_read_names_dict[pc.reference_pos][1], whether_mutation)
                            # store sample index
                            effective_pos_read_names_dict[pc.reference_pos][2] = np.append(
                                effective_pos_read_names_dict[pc.reference_pos][2], np.full(len(piled_bases), index))

                        elif pc.reference_pos in variant_pos:
                            # store reads names
                            variant_pos_read_names_dict[pc.reference_pos][0] = np.append(
                                variant_pos_read_names_dict[pc.reference_pos][0], np.array(pc.get_query_names()))
                            variant_pos_read_names_dict[pc.reference_pos][1] = np.append(
                                variant_pos_read_names_dict[pc.reference_pos][1], whether_mutation)
                            # store sample index
                            variant_pos_read_names_dict[pc.reference_pos][2] = np.append(
                                variant_pos_read_names_dict[pc.reference_pos][2], np.full(len(piled_bases), index))
        return effective_pos_read_names_dict, variant_pos_read_names_dict, pos_genotype_base_dict

    @staticmethod
    def _parse_piled_bases(piled_bases):
        # build a dict to record base types
        base_type = {".": 0, "A": 0, "G": 0, "C": 0, "T": 0}
        for base in piled_bases:
            if base == "." or base == ",":
                base_type["."] += 1
            # if base is an alpha
            else:
                base = base.upper()     # turn all bases into upper
                if base in ["A", "G", "C", "T"]:
                    base_type[base] += 1
                elif len(base) > 2 and base[2:].isalnum():
                    # update the indel number in base_type dict
                    base_indel = base[1:]
                    base_type[base_indel] = base_type.setdefault(base_indel, 0) + 1

        ref_num = base_type["."]
        # sort the base_type result
        base_type = sorted(base_type.items(), key=lambda x: x[1], reverse=True)
        if base_type[0][0] == ".":
            alt1, alt1_num = base_type[1]
            alt2, alt2_num = base_type[2]
        else:
            alt1, alt1_num = base_type[0]
            if base_type[1][0] == ".":
                alt2, alt2_num = base_type[2]
            else:
                alt2, alt2_num = base_type[1]
        if alt1_num == 0:
            alt1 = ""
        if alt2_num == 0:
            alt2 = ""
        # return base type and its number
        return ref_num, alt1, alt1_num, alt2, alt2_num

    @staticmethod
    def _get_alt_base(alt_bases):
        # input is an one-dimensional array of alt bases, like array["", "A", "A", ""]
        # return the base with max count
        alt_base_dict = {}
        for base in alt_bases:
            if base != "":
                if base in alt_base_dict:
                    alt_base_dict[base] += 1
                else:
                    alt_base_dict[base] = 1
        if alt_base_dict == {}:
            return None
        return max(alt_base_dict.items(), key=lambda x: x[1])[0]

    @staticmethod
    def _output_haplotypes(haplotypes_list, pos_genotype_base, posterior_probability_dict, chromosome):
        """
        :param haplotypes_list: three haplotype generate by haplotype.py
        :param pos_genotype_base: {pos:{0:"A", 1:"G"}
        :param posterior_probability_dict:
        :param chromosome: chromosome name
        :return: a dict, {pos: output_record}
        """
        # like {pos1: genotype}, {183742: 0}
        pos_effective = []
        pos_depth_haplotypes = []   # store dicts of pos_depth of three haplotypes
        for h in haplotypes_list:
            pos_effective.extend(list(h.pos_genotype.keys()))
            pos_depth_haplotypes.append(h.get_pos_depth())

        pos_effective = list(set(pos_effective))
        pos_effective.sort()

        output_dict = {}
        for pos in pos_effective:
            if pos not in pos_genotype_base:
                print(pos)
                continue
            # output pos, ref, alt, likelihood
            one_line = [chromosome, str(pos+1), pos_genotype_base[pos][0], pos_genotype_base[pos][1],
                        "\t".join(posterior_probability_dict[pos].astype(str))]
            # output genotype and depth
            for i in range(3):
                if pos in haplotypes_list[i].pos_genotype:
                    one_line.append(str(haplotypes_list[i].pos_genotype[pos]))

                else:
                    one_line.append(".")
                if pos in pos_depth_haplotypes[i]:
                    one_line.append(str(pos_depth_haplotypes[i][pos]))
                else:
                    one_line.append(".")
            output_dict.update({pos: "\t".join(one_line) + "\n"})
        return output_dict

    @staticmethod
    def _output_shared_variant_site(chromosome, pos_geno_read_number, pos_genotype_base):
        """
        # Generate dict of record of variant site sharing reads with effective sites
        :param chromosome: like NZ_AP012323.1
        :param pos_geno_read_number: {pos: ([1, 0, "."], [12, 11, 0])}
        :param pos_genotype_base: {pos: {0:"A", 1:"G"}}
        :return: one example: {pos: "NZ_AP012323.1 2140    C   T   nan nan nan nan nan nan 0   1763    1   468 1   53"}
        """
        variant_record_dict = {}
        for pos, geno_read_number in pos_geno_read_number.items():
            variant_record = chromosome + "\t" + str(pos+1) + "\t"
            genotype_base = pos_genotype_base[pos]
            variant_record += genotype_base[0] + "\t" + genotype_base[1] + "\t" + "-9\t-9\t-9\t-9\t-9\t-9\t"
            for geno, count_num in zip(geno_read_number[0], geno_read_number[1]):
                variant_record += str(geno) + "\t" + str(count_num) + "\t"
            variant_record_dict.update({pos: variant_record + "\n"})
        return variant_record_dict

    @staticmethod
    def _find_indel_surrounding_pos(sorted_effective_pos, pos_genotype_base, surrounding_range=10):
        """
        Identify pos of INDEL and
        index (bool) of surrounding pos in sorted_effective_pos including indel
        :param sorted_effective_pos: an array, sorted pos of effective sites, array([1, 3, 6])
        :param pos_genotype_base: {pos: {0:"A", 1:"+G"}}
        :param surrounding_range: definition of surrounding pos, default +- 10 bp
        :return: np.array of indel pos; bool array of indel_surrounding in sorted_effective_pos
        """
        indel_surrounding_index_bool = np.full(len(sorted_effective_pos), False)
        indel_pos_list = [i for i in sorted_effective_pos if pos_genotype_base[i][1][0] in ["-", "+"]]
        if len(indel_pos_list) > 0:
            # get all intervals of indel surrounding
            total_indel_range = []      # [[10, 30], [50, 75], [100, 120]]
            current_range = [indel_pos_list[0]-surrounding_range, indel_pos_list[0]+surrounding_range]
            for indel_pos in indel_pos_list:
                if (indel_pos-10) <= current_range[1]:
                    current_range[1] = indel_pos + 10
                else:
                    total_indel_range.append(current_range)
                    current_range = [indel_pos-10, indel_pos+10]
            total_indel_range.append(current_range)

            for current_range in total_indel_range:
                # bool index of pos in current_range
                in_range_index = (sorted_effective_pos >= current_range[0]) * (sorted_effective_pos <= current_range[1])
                indel_surrounding_index_bool = indel_surrounding_index_bool + in_range_index
        return np.array(indel_pos_list, dtype=int), indel_surrounding_index_bool

    def strain_initial(self, chromosome):
        """
        Initialize the strain counts matrix
        :param chromosome:
        :return: return the np.array of estimated proportion, M by 3 matrix
        """
        print("Strain initial bin ... ...")
        reads_bin_chromosome = self.reads_in_bin_by_chromosome[chromosome]
        # row: bins, col: start, end, reads of all sample; numpy array, col number is sample number +2
        # combine all sample to three columns and get target bins (high coverage)
        print(reads_bin_chromosome)
        reads_bin_all_sample = np.c_[
            reads_bin_chromosome[:, 0:2], np.sum(reads_bin_chromosome[:, 2:reads_bin_chromosome.shape[1]], 1)]
        # target bins are bins' coverage at quantile 0.8, by default 5 bins
        # three cols: start, end and depth (all sample combined)
        target_bins = get_nonzero_quantile_bins(reads_bin_all_sample, q=0.8, bin_num=5)
        if target_bins is None:
            return None

        # build an dict to store {pos: {"read_names", "0-1", "sample index"}}
        # collect all sites information in all target bins (default 5)
        pos_read_names_effective = {}
        pos_read_names_variant = {}
        for target_bin in target_bins:
            pos_read_names_effective, pos_read_names_variant, pos_genotype_base = \
                Subject.get_pos_read_names(self, pos_read_names_effective,
                                           pos_read_names_variant, target_bin, chromosome)

        if pos_read_names_effective == {}:
            # all are variant sites; get the coverage of target bins in all samples
            # find the index of target bins in reads_bin_chromosome
            all_target_bin_start, all_bin_index, target_bin_index = np.intersect1d(reads_bin_chromosome[:, 0],
                                                                                   target_bins[:, 0],
                                                                                   return_indices=True)
            # assigned all reads count to the dominant strain
            sample_reads_count_in_target_bins = np.sum(reads_bin_chromosome[all_bin_index, 2:], axis=0)
            # initial counts of secondary strain, can not < 1
            secondary_initial_counts = np.maximum(1, (sample_reads_count_in_target_bins/20).astype(int))
            reads_count_by_strains = np.column_stack((sample_reads_count_in_target_bins,
                                                     secondary_initial_counts,
                                                     np.ones(np.shape(sample_reads_count_in_target_bins)[0], int)))
            return reads_count_by_strains

        else:
            haplotypes, reads_count_by_strains, posterior_probability_dict = haplotype.haplotype_identification(
                pos_read_names_effective, len(self.sample_names))
            """
            # record strain proportions
            with open(f"{self.ID}_strain_proportion_record.txt", "w") as record_f:
                for i in reads_count_by_strains:
                    record_f.write("\t".join([str(j) for j in i]) + "\n")
            # return the np.array of estimated proportion, M by 3 matrix
            """
            return reads_count_by_strains

    def generate_vcf(self, output_prefix):
        if len(self.haplotype_filename):
            print(f"No haplotype file is generated!")
        else:
            vcf_filename = f"{output_prefix}_primary_secondary.vcf"

    def strain_identification(self, output_prefix, bin_len=10000):
        """
        # Estimate strain RA and strain's variations
        :param output_prefix:
        :param bin_len:
        :return:
        """
        ########################
        # initialize strain RA #
        ########################
        """ initialize strain RA in regions with adequate coverage from the longest chr to shortest chr """
        # once initialized, not need to do in other chromosome
        reads_count_by_strain_initial = None
        chromosome_len_dict = {}
        # assume three strains
        reads_count_by_strains = np.zeros((len(self.samples), 3))

        for chromosome, chromosome_len in zip(self.pysam_ref.references, self.pysam_ref.lengths):
            chromosome_len_dict.update({chromosome: chromosome_len})
        # get sorted chromosome and length
        chromosome_items = list(chromosome_len_dict.items())
        chromosome_items.sort(key=lambda x: x[1], reverse=True)
        for chromosome, chromosome_len in chromosome_items:
            reads_count_by_strain_initial = Subject.strain_initial(self, chromosome)
            if reads_count_by_strain_initial is not None:
                reads_count_by_strains = reads_count_by_strain_initial.copy()
                break
        if reads_count_by_strain_initial is not None:
            print("Initial estimation succeeds.")
        else:
            print("Fail! No region with adequate coverage for Initial estimation.")
            return
        # get reads_count_by_strain_initial and reads_count_by_strains if pass
        # output strain counts
        # np.savetxt(f"{output_prefix}_reads_count_by_strain_initial.txt", reads_count_by_strain_initial, fmt="%d", delimiter="\t")

        output_file = open(f"{output_prefix}_haplotypes.txt", "w")
        self.haplotype_filename = f"{output_prefix}_haplotypes.txt"
        ###############################
        # screen bins and chromosomes #
        ###############################
        """ snv calling for all chromosome """
        for chromosome, chromosome_len in zip(self.pysam_ref.references, self.pysam_ref.lengths):
            print(f"processing chromosome: {chromosome}")
            # build a list of dict with length=3 to store depth of each sample of haplotypes (for estimating lambda)
            haplotypes_sample_depths = [{}, {}, {}]

            # build dicts to store haplotype results, {pos: output_record}
            chr_haplotype_output = {}

            # build dicts to store bins without effective sites (isolated variant sites)
            chr_pos_read_names_variant = {}
            chr_pos_genotype = {}

            for boundary in get_nodes(chromosome_len, bin_len):
                # boundary is like [10000, 20000]
                # get a dict containing information that bin
                pos_read_names_effective, pos_read_names_variant, pos_genotype_base = \
                    Subject.get_pos_read_names(self, {}, {}, boundary, chromosome)

                # if no result, continue
                if len(pos_read_names_effective) == 0 and len(pos_read_names_variant) == 0:
                    continue

                # get INDEL pos in effective sites
                effective_indel_pos_array, effective_indel_surrounding_index_bool = self._find_indel_surrounding_pos(
                    np.array(sorted(list(pos_read_names_effective.keys())), dtype=int), pos_genotype_base)
                effective_pos_number_remove_indel_surrounding = np.sum(~effective_indel_surrounding_index_bool)

                # save the reads count of strains, M by 3 matrix
                current_reads_count_by_strains = reads_count_by_strains.copy()

                if pos_read_names_effective == {}:
                    #####################
                    # no effective site #
                    #####################
                    # only variant site, store them in chr dict
                    chr_pos_read_names_variant.update(pos_read_names_variant)
                    chr_pos_genotype.update(pos_genotype_base)

                else:
                    ########################
                    # have effective sites #
                    ########################
                    # get strain results of effective sites in this bin
                    haplotypes, reads_count_by_strains, posterior_probability_dict = haplotype.haplotype_identification(
                        pos_read_names_effective, len(self.sample_names), reads_count_by_strains)
                    """ store the variation of effective sites """
                    chr_haplotype_output.update(Subject._output_haplotypes(haplotypes, pos_genotype_base,
                                                                           posterior_probability_dict, chromosome))

                    # update sample depth of haplotypes (for estimating lambda)
                    for h_index, haplotype_a in enumerate(haplotypes):
                        haplotypes_sample_depths[h_index].update(haplotype_a.get_samples_pos_depth())

                    """ store the variation of variant sites """
                    if pos_read_names_variant != {}:
                        # pos_variant_site_reserved is list of pos, like [100, 108, 111, 134]
                        # pos_geno_read_number is a dict {pos: ([1, 0, "."], [12, 11, 0])}
                        pos_variant_site_reserved, pos_geno_read_number = \
                            haplotype.incorporate_variant_sites_into_haplotype(haplotypes, pos_read_names_variant)
                        # update isolate variants to total dict of the whole chromosome
                        if len(pos_variant_site_reserved) > 0:
                            for pos_a in pos_variant_site_reserved:
                                chr_pos_read_names_variant.update({pos_a: pos_read_names_variant[pos_a]})
                                chr_pos_genotype.update({pos_a: pos_genotype_base[pos_a]})
                        # update variants around effective sites to chr_haplotype_output
                        if pos_geno_read_number != {}:
                            chr_haplotype_output.update(
                                Subject._output_shared_variant_site(chromosome, pos_geno_read_number,
                                                                    pos_genotype_base))

                if effective_pos_number_remove_indel_surrounding <= 0.1 * len(pos_read_names_variant):
                    #################################
                    # update reads_count_by_strains #
                    #################################
                    # low proportion of effective sites, variant site should be counted
                    assigned_reads_in_this_bin = np.sum(reads_count_by_strains - current_reads_count_by_strains,
                                                        axis=1)  # len is sample number M
                    # get the total number of reads in this bin from all bin record
                    reads_in_bin_in_this_chromosome = self.reads_in_bin_by_chromosome[chromosome]
                    reads_in_this_bin = reads_in_bin_in_this_chromosome[reads_in_bin_in_this_chromosome[:, 0] == boundary[0], 2:][0]  # len is sample number M
                    # add unassigned reads to dominant strain's column
                    dominant_index = np.argmax(np.sum(reads_count_by_strains), axis=0)  # 0 or 1 or 2
                    reads_count_by_strains[:, dominant_index] += reads_in_this_bin - assigned_reads_in_this_bin
                    ###################################
                    # update haplotypes_sample_depths #
                    ###################################

            ###################################################
            # process isolate variants sites of the whole chr #
            ###################################################
            # get lambda_matrix
            lambda_matrix = []
            for haplotype_sample_depths in haplotypes_sample_depths:
                # sample lambda for a haploytype
                if haplotype_sample_depths == {}:
                    # if no site is recorded, haplotype_lambda is dim one, length of M, numpy array
                    haplotype_lambda = np.zeros(len(self.sample_names))
                else:
                    # if there are sites, get average depths
                    haplotype_lambda = np.mean(np.array(list(haplotype_sample_depths.values())), axis=0)
                lambda_matrix.append(list(haplotype_lambda))
            lambda_matrix = np.array(lambda_matrix)     # 3 by M matrix
            print("lambda matrix is")
            print(lambda_matrix)
            # calculate the results of isolate variants sites
            isolate_variants_result_dict = haplotype.infer_isolated_variant_sites(chr_pos_read_names_variant, chr_pos_genotype, chromosome, lambda_matrix)
            chr_haplotype_output.update(isolate_variants_result_dict)

            # output all result of this chromosome
            for pos in sorted(chr_haplotype_output.keys()):
                output_file.write(chr_haplotype_output[pos])

        # print initial reads count for all chromosome
        with open(f"{output_prefix}_reads_count_by_strain_initial.txt", "w") as count_result_f:
            count_result_f.write(f"Sample\tPrimary\tSecondary\tNoise\n")
            for sample_name, sample_count in zip(self.sample_names, reads_count_by_strain_initial):
                count_result_f.write(sample_name + "\t" + "\t".join(sample_count.astype(dtype=str)) + "\n")
        # print cumulative reads count for all chromosome
        output_read_count_array = (reads_count_by_strains - reads_count_by_strain_initial).astype(dtype=str)
        with open(f"{output_prefix}_reads_count_by_strain_final.txt", "w") as count_result_f:
            count_result_f.write(f"Sample\tPrimary\tSecondary\tNoise\n")
            for sample_name, sample_count in zip(self.sample_names, output_read_count_array):
                count_result_f.write(sample_name + "\t" + "\t".join(sample_count.astype(dtype=str)) + "\n")

        # print cumulative proportions for all chromosome
        with open(f"{output_prefix}_proportion_by_strain_initial.txt", "w") as proportion_result_f:
            proportion_result_f.write(f"Sample\tPrimary\tSecondary\tNoise\n")
            for sample_name, sample_count in zip(self.sample_names, reads_count_by_strain_initial):
                sample_count = sample_count.astype(dtype=int)
                sample_prop = sample_count / np.sum(sample_count)
                proportion_result_f.write(sample_name + "\t" + "\t".join(sample_prop.astype(dtype=str)) + "\n")

        # print cumulative proportions for all chromosome
        with open(f"{output_prefix}_proportion_by_strain_final.txt", "w") as proportion_result_f:
            proportion_result_f.write(f"Sample\tPrimary\tSecondary\tNoise\n")
            for sample_name, sample_count in zip(self.sample_names, output_read_count_array):
                sample_count = sample_count.astype(dtype=int)
                sample_prop = sample_count/np.sum(sample_count)
                proportion_result_f.write(sample_name + "\t" + "\t".join(sample_prop.astype(dtype=str)) + "\n")

        output_file.close()


def fqs_species_process(subject_name, sample_name_list, fqs_path_list, species_ref, bowtie_ref, output_path):
    """
    Process fqs with known species name, one species a time
    :param subject_name: Prefix of the subject
    :param sample_name_list: ["X-t0", "X-t1"]
    :param fqs_path_list: fqs is like [[path/X_1.fq,path/X_2.fq],[path/X.fq]]
    :param species_ref: fas reference of one species
    :param bowtie_ref: corresponding bowtie ref of that fas
    :param output_path:
    :return:
    """
    # ### change to output path ###
    if not os.path.exists(output_path):
        os.system(f"mkdir -p {output_path}")
    os.chdir(output_path)
    print(f"subject name is: {subject_name}\nfastq list is: {fqs_path_list}\nmapping reference is: {species_ref}\n"
          f"output path is: {output_path}")
    subject1 = Subject(subject_name, sample_name_list, species_ref)
    # add fastqs to Subject class
    for i in fqs_path_list:
        subject1.add_fqs(i)
    # get bam file
    bam_list = subject1.fqs_to_bams(bowtie_ref)
    subject1.add_bams(bam_list)
    subject1.coverage_summary()
    # main function of strain identification
    subject1.strain_identification(subject_name)
    # generate vcf file
    # subject1.


def bams_species_process(subject_name, sample_name_list, bams_path_list, species_ref, output_path):
    """
    Process fqs with known species name, one species a time
    :param subject_name: Prefix of the subject
    :param bams_path_list: abs path of bams is like [path/X.bam, path/Y.bam]
    :param species_ref: fas reference of one species
    :param output_path:
    :return:
    """
    # ### change to output path ###
    if not os.path.exists(output_path):
        os.system(f"mkdir -p {output_path}")
    os.chdir(output_path)
    print(f"subject name is: {subject_name}\nfastq list is: {bams_path_list}\nmapping reference is: {species_ref}\n"
          f"output path is: {output_path}")
    subject1 = Subject(subject_name, sample_name_list, species_ref)

    # add virtual fastqs to Subject class
    for i in range(len(bams_path_list)):
        subject1.add_fqs(f"virtual_{i}.fq")

    # get bam file
    # bam_list = subject1.fqs_to_bams(bowtie_ref)
    subject1.add_bams(bams_path_list)
    subject1.coverage_summary()
    # main function of strain identification
    subject1.strain_identification(subject_name)


if __name__ == "__main__":
    # for one test: python3 longitudinal_microbiome.py
    # /gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/TEDDY_T1D/1 1 1685
    # target_path, case_control_id, species_id = sys.argv[1:4]
    # os.chdir(target_path)
    # species_process(species_id, case_control_id)

    """ for fqs to strain"""
    model = sys.argv[1]
    if model == "fq2strain":
        subject_name, sample_names, fqs, species_ref, bowtie_ref, output_path = sys.argv[2:8]
        fqs = fqs.split(":")
        # output_path = "/gpfs/data/lilab/home/zhoub03/Blaser_data/LongStrain_test"
        sample_name_list = sample_names.split(",")
        fqs_species_process(subject_name, sample_name_list, [i.split(",") for i in fqs], species_ref, bowtie_ref, output_path)
        command1 = "python3 longitudinal_microbiome.py fq2strain Bifidobacterium_breve_simu1 " \
                   "Bifidobacterium_breve_simu1_t0_1.fq,Bifidobacterium_breve_simu1_t0_2.fq:" \
                   "Bifidobacterium_breve_simu1_t1_1.fq,Bifidobacterium_breve_simu1_t1_2.fq:" \
                   "Bifidobacterium_breve_simu1_t2_1.fq,Bifidobacterium_breve_simu1_t2_2.fq " \
                   "/gpfs/data/lilab/home/zhoub03/software/my_strain/Bifidobacterium_breve/Bifidobacterium_breve.fas " \
                   "/gpfs/data/lilab/home/zhoub03/software/my_strain/Bifidobacterium_breve/Bifidobacterium_breve"

    elif model == "bam2strain":
        subject_name, sample_names, bams_path_list, species_ref, output_path = sys.argv[2:7]
        bams_path_list = bams_path_list.split(",")
        sample_name_list = sample_names.split(",")
        bams_species_process(subject_name, sample_name_list, bams_path_list, species_ref, output_path)

    else:
        print("Model error! Model must be fq2strain or bam2strain.")
