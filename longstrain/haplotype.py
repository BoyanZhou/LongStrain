# construct microbiome haplotype

import numpy as np
import sys


class Haplotype:
    """
    each haplotype represents a strain, default 3 Haplotypes in one sample
    """
    def __init__(self, index, sample_size, reads_count_by_sample=np.array([])):
        self.sample_size = sample_size                  # number of sample contained in this haplotype
        self.number = index                             # index of the haplotype (0, 1, or 2)
        self.sample_read_names = [[]] * sample_size     # reads name contained by each sample, not used now
        self.read_names_sample_index = {}               # indicate which sample the read belongs to
        self.read_names_positions = {}              # positions interrogated in that reads, like {"read1": [12, 38, 59]}
        if reads_count_by_sample.shape == (0, ):
            self.reads_count_by_sample = np.ones(sample_size, dtype=int)
        else:
            self.reads_count_by_sample = reads_count_by_sample  # vector, length is sample size
        self.pos_genotype = {}                          # like {pos1: genotype}, {183742: 0}
        # self.reads_number = 0

    def haplotype_extend(self, geno_index, data_dict, pos, geno):
        """
        Add data at pos to existing Haplotype
        :param geno_index: bool value, index of reads interrogated by geno
        :param data_dict: the entire data in the block with many pos
        :param pos: the target pos
        :param geno: 0 or 1
        :return:
        """
        # add reads of a new locus to the haplotype
        data_list = data_dict[pos]
        # add a new locus to pos_genotype
        self.pos_genotype.update({pos: geno})
        for index in np.where(geno_index)[0]:
            # for each reads
            read_name = data_list[0][index]
            read_sample_index = data_list[2][index]
            self.read_names_sample_index.update({read_name: read_sample_index})
            if read_name in self.read_names_positions:
                self.read_names_positions[read_name].append(pos)
            else:
                self.read_names_positions.update({read_name: [pos]})

            if read_name not in self.sample_read_names[read_sample_index]:
                self.sample_read_names[read_sample_index].append(read_name)
                self.reads_count_by_sample[read_sample_index] += 1          # update counts by sample
        # self.reads_number = len(self.read_names_sample_index.keys())

    def haplotype_subtract(self, geno_index, data_dict, pos):
        # subtract one read from the haplotype
        data_list = data_dict[pos]
        read_name = data_list[0][geno_index]
        read_sample_index = data_list[2][geno_index]

        if read_name in self.read_names_sample_index:
            del self.read_names_sample_index[read_name]
        if read_name in self.sample_read_names[read_sample_index]:
            self.sample_read_names[read_sample_index].remove(read_name)
            self.reads_count_by_sample[read_sample_index] -= 1       # update counts by sample
        read_positions = self.read_names_positions.pop(read_name)
        read_positions_geno = [self.pos_genotype[i] for i in read_positions]
        # print(f"subtract {read_positions_geno}")
        # self.reads_number = len(self.read_names_sample_index.keys())
        # return a list containing positions
        return read_positions, read_positions_geno

    def haplotype_add(self, geno_index, data_dict, pos, read_positions, read_positions_geno):
        # add one read to the haplotype, pos is not included in read_positions
        # geno_index is an int, data_list has three elements: read names, geno(0 or 1), sample index
        data_list = data_dict[pos]
        read_name = data_list[0][geno_index]
        read_sample_index = data_list[2][geno_index]

        read_positions.append(pos)
        read_positions_geno.append(data_list[1][geno_index])
        self.read_names_sample_index.update({read_name: read_sample_index})
        if read_name not in self.sample_read_names[read_sample_index]:
            self.sample_read_names[read_sample_index].append(read_name)
            self.reads_count_by_sample[read_sample_index] += 1       # update counts by sample
        # read_positions is a list
        if read_name in self.read_names_positions:
            self.read_names_positions[read_name].extend(read_positions)
        else:
            self.read_names_positions.update({read_name: read_positions})
        # self.reads_number = len(self.read_names_sample_index.keys())

        for p, geno in zip(read_positions, read_positions_geno):
            if p not in self.pos_genotype:
                # print(read_name, p, geno)
                self.pos_genotype.update({p: geno})

    def get_pos_depth(self):
        """
        # Get the coverage depth at each site
        :return: total depth of this haplotype at pos
        """
        pos_depth_dict = {}
        for pos_list in self.read_names_positions.values():
            for pos in pos_list:
                if pos in pos_depth_dict:
                    pos_depth_dict[pos] += 1
                else:
                    pos_depth_dict[pos] = 1
        return pos_depth_dict

    def get_samples_pos_depth(self):
        """
        # Get the coverage depth of samples at each site
        :return: a dict, depth of each sample {pos: [dp1, dp2, dp3, dp4]}, length is sample size
        """
        pos_sample_depth_dict = {}
        for read_name in self.read_names_positions.keys():
            read_sample_index = self.read_names_sample_index[read_name]
            for pos in self.read_names_positions[read_name]:
                # get pos in pos-list of that read
                if pos in pos_sample_depth_dict:
                    pos_sample_depth_dict[pos][read_sample_index] += 1
                else:
                    # pos has not been record
                    pos_depth_list = [0] * self.sample_size
                    pos_depth_list[read_sample_index] += 1
                    pos_sample_depth_dict.update({pos: pos_depth_list})
        return pos_sample_depth_dict


def whether_skip(likelihood_array, skip_threshold):
    # one dimension likelihood array
    order_index = np.argsort(likelihood_array)
    genotype_to_haplotype3 = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0, 1]])
    geno_difference = genotype_to_haplotype3[order_index[-1]] - genotype_to_haplotype3[order_index[-2]]
    if np.sum(geno_difference == 0) == 2:
        # if there two same genotypes
        max_likelihood = likelihood_array[order_index[-1]] + likelihood_array[order_index[-2]]
    else:
        max_likelihood = likelihood_array[order_index[-1]]
    if max_likelihood < skip_threshold:
        return True
    else:
        return False


def calculate_maximum_likelihood_matrix(value_matrix, bool_matrix):
    # this function calculate the sum of maximum non-zero likelihood of each row
    # value_matrix likes array([[-1.63657945, -0.24074084, -3.94713272], [-1.63657945, -0.24074084, -3.94713272],
    # [-1.63657945, -0.24074084, -3.94713272], [-1.63657945, -0.24074084, -3.94713272]])
    # bool_matrix likes array([[False, False,  True], [False, False,  True],
    # [ True,  True, False], [ True,  True, False]])
    matrix_with_true_index = value_matrix * bool_matrix
    row_index_more_than_one_true = np.sum(bool_matrix, axis=1) > 1
    likelihood_sum = np.sum(matrix_with_true_index) - np.sum(np.min(matrix_with_true_index[row_index_more_than_one_true], axis=1))
    return likelihood_sum


def haplotype_identification(pos_reads_mutation_dict, sample_size, reads_count_by_strains=np.array([]),
                             skip_threshold=0.99):
    """
    Main function of identifying haplotypes
    :param pos_reads_mutation_dict: is like {pos: [read_names, "0-1" genotype, sample index]}
    :param sample_size: number of samples in this individual
    :param reads_count_by_strains:
    :param skip_threshold:
    :return: haplotypes=[Haplotype0, Haplotype1, Haplotype2], reads_count_by_strains=M by 3 matrix storing reads counts,
    posterior_probability_dict= {pos: probability_six_combination}
    """
    genotype_to_haplotype = np.array([[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]])
    genotype_to_haplotype3 = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0, 1]])
    error_rate = np.log(0.005)
    # print("haplotype identification start ... ...")
    posterior_probability_dict = {}

    # initial bin
    if reads_count_by_strains.shape == (0, ):
        bin_initial = True
    else:
        bin_initial = False

    # matrix for storing counts of each strains in each sample, samples * strains
    if bin_initial:
        # if this is haplotype initialization, assign 0 to the count
        reads_count_by_strains = np.ones([sample_size, 3], dtype=int)
        # initialize the list storing haplotypes
        haplotypes = [Haplotype(i, sample_size) for i in range(3)]
    else:
        # if haplotypes have been initialized, assign cumulative counts to haplotypes
        # print(reads_count_by_strains)
        haplotypes = [Haplotype(i, sample_size, reads_count_by_strains[:, i]) for i in range(3)]
    pos_last = -1000   # update after each position, give it a negative value as start
    haplotype_last = np.array([], dtype=int)    # strain index of reads (0, 1 or 2 ) at last position

    for pos in sorted(pos_reads_mutation_dict.keys()):
        # print(pos)
        # need to estimate which strains these new reads belong to, update for each position
        strains_proportions = reads_count_by_strains / np.tile(np.sum(reads_count_by_strains, axis=1), (3, 1)).T
        if pos == -2:
            print(strains_proportions)
        strains_proportions = np.log(strains_proportions)

        """ first pos in dict; start of the chain """
        if bin_initial:
            index_0 = pos_reads_mutation_dict[pos][1] == 0  # which genotype is 0
            index_1 = pos_reads_mutation_dict[pos][1] == 1
            haplotype_last = np.full(len(index_0), 0)
            if np.sum(index_0) >= np.sum(index_1):
                haplotypes[0].haplotype_extend(index_0, pos_reads_mutation_dict, pos, 0)
                haplotypes[1].haplotype_extend(index_1, pos_reads_mutation_dict, pos, 1)
                haplotype_last[index_1] = 1
            else:
                haplotypes[1].haplotype_extend(index_0, pos_reads_mutation_dict, pos, 0)
                haplotypes[0].haplotype_extend(index_1, pos_reads_mutation_dict, pos, 0)
                haplotype_last[index_0] = 1
            pos_last = pos
            # update reads count matrix
            reads_count_by_strains[:, 0] = haplotypes[0].reads_count_by_sample.copy()
            reads_count_by_strains[:, 1] = haplotypes[1].reads_count_by_sample.copy()

        else:
            """ not the first pos, the chain has started """
            # compare read names of last position with present position to find shared reads
            # three outcomes: shared names, index in pos_last, index in pos_current

            if pos_last != -1000:
                read_shared = np.intersect1d(pos_reads_mutation_dict[pos_last][0], pos_reads_mutation_dict[pos][0],
                                             assume_unique=False, return_indices=True)
            else:
                read_shared = ([],)

            if len(read_shared[0]) == 0:
                # there is no shared read between two positions
                # get sample index, row of reads_strains_proportions = reads number
                reads_strains_proportions = strains_proportions[pos_reads_mutation_dict[pos][2]]

                # reads_genotype_matrix stores genotype for each reads
                reads_genotype_matrix = np.tile(pos_reads_mutation_dict[pos][1], (2, 1))
                reads_genotype_matrix[0, :] = 1 - reads_genotype_matrix[0, :]
                likelihood_six = []

                for reads_to_strains in np.dot(genotype_to_haplotype, reads_genotype_matrix):
                    likelihood_six.append(np.sum(reads_strains_proportions * np.eye(3, dtype=int)[reads_to_strains]))
                likelihood_max_index = np.array(likelihood_six).argmax()

                # corresponding relationship between geno and strain; e.g. [0, 2]
                geno_strain = genotype_to_haplotype[likelihood_max_index]
                probability_six = np.exp(likelihood_six)
                probability_six = probability_six/np.sum(probability_six)
                # if the maximum probability < 0.99, skip this position and not record this position
                if probability_six[likelihood_max_index] < skip_threshold and not bin_initial:
                    # print(pos)
                    # print(probability_six)
                    continue
                posterior_probability_dict.update({pos: probability_six})

                if pos == -2:
                    print(read_shared)
                    print("posterior probability is:")
                    print(probability_six)
                    print("pos is:")
                    print(pos)
                    print(pos_reads_mutation_dict[pos])
                    print("last pos is:")
                    print(pos_last)
                    print(pos_reads_mutation_dict[pos_last])
                    print("haplotype last is:")
                    print(haplotype_last)
                    sys.exit()

                # update pos information
                haplotypes[geno_strain[0]].haplotype_extend(pos_reads_mutation_dict[pos][1] == 0,
                                                            pos_reads_mutation_dict, pos, 0)
                haplotypes[geno_strain[1]].haplotype_extend(pos_reads_mutation_dict[pos][1] == 1,
                                                            pos_reads_mutation_dict, pos, 1)
                pos_last = pos
                reads_count_by_strains[:, geno_strain[0]] = haplotypes[geno_strain[0]].reads_count_by_sample.copy()
                reads_count_by_strains[:, geno_strain[1]] = haplotypes[geno_strain[1]].reads_count_by_sample.copy()
                haplotype_last = pos_reads_mutation_dict[pos][1].copy()
                haplotype_last[pos_reads_mutation_dict[pos][1] == 0] = geno_strain[0]
                haplotype_last[pos_reads_mutation_dict[pos][1] == 1] = geno_strain[1]

            else:
                # there are shared reads
                """ for six strains-genotype combinations, get the one with max likelihood """
                likelihood_six_combination = []
                geno_actual = pos_reads_mutation_dict[pos][1][read_shared[2]].copy()
                for combination in genotype_to_haplotype3:
                    if pos == -2:
                        print(combination)
                    # calculate the likelihood of two part
                    likelihood_reads_consistency = 0
                    likelihood_strain_in_sample = 0
                    # haplotype_last = np.array([1, 1, 0, 0, 2])
                    # combination = genotype_to_haplotype3[0], e.g. [0, 1, 1]
                    # geno_actual = np.array([1, 1, 0, 0, 0])
                    geno_predicted = combination[haplotype_last[read_shared[1]]]

                    if np.sum(geno_predicted != geno_actual) > 0:
                        # reads inconsistency exists
                        inconsistent_index = read_shared[2][geno_actual != geno_predicted]  # index, not bool
                        # likelihood of strain convert
                        combination_stack = np.tile(combination, [len(inconsistent_index), 1])
                        strain_available_a = combination_stack == np.tile(
                            pos_reads_mutation_dict[pos][1][inconsistent_index], [3, 1]).T
                        # the former is combination_stack, the later is stacked true genotypes
                        # [[0, 1, 1],       [[1, 1, 1],
                        #  [0, 1, 1],        [1, 1, 1],
                        #  [0, 1, 1],        [0, 0, 0],
                        #  [0, 1, 1],]       [0, 0, 0],]

                        # get set of strains of last position, set them to unavailable
                        strain_available = np.array([True, True, True])
                        strain_available[np.array(list(set(haplotype_last)))] = False
                        strain_available = np.tile(strain_available, [len(inconsistent_index), 1])

                        # get reads index of possible strain convert, three columns, bool
                        reads_converted = strain_available_a * strain_available
                        reads_converted_index = np.sum(reads_converted, axis=1) > 0

                        # if strain converted is possible
                        if np.sum(reads_converted_index) > 0:
                            strain_converted_likelihood = strains_proportions[pos_reads_mutation_dict[pos][2][
                                inconsistent_index[reads_converted_index]]] * reads_converted[reads_converted_index]
                            likelihood_reads_consistency += np.sum(
                                np.maximum(error_rate, np.sum(strain_converted_likelihood, axis=1)))

                        # add the likelihood of sequencing error part
                        if np.sum(reads_converted_index) < reads_converted.shape[0]:
                            likelihood_reads_consistency += error_rate * (
                                        reads_converted.shape[0] - np.sum(reads_converted_index))

                    else:
                        # no reads inconsistency, calculate log_likelihood of shared reads (given strain)
                        reads_strains_proportions = strains_proportions[
                            pos_reads_mutation_dict[pos_last][2][read_shared[1]]]
                        likelihood_strain_in_sample += np.sum(
                            reads_strains_proportions * np.eye(3, dtype=int)[haplotype_last[read_shared[1]]])
                    if pos == -2:
                        print(likelihood_reads_consistency, likelihood_strain_in_sample)

                    if len(read_shared[2]) < len(pos_reads_mutation_dict[pos][0]):
                        # print(pos)
                        # there are unshared new reads
                        # likelihood of unshared reads in new position
                        unshared_index = np.full(len(pos_reads_mutation_dict[pos][0]), True)
                        unshared_index[read_shared[2]] = False
                        # combination = np.array(array([0, 1, 1]))
                        # print(unshared_index)
                        combination_stack = np.tile(combination, [np.sum(unshared_index), 1])
                        # print(combination_stack)
                        # print(np.tile(pos_reads_mutation_dict[pos][1][unshared_index], [3, 1]).T)
                        strains_index_unshared_reads = combination_stack == np.tile(
                            pos_reads_mutation_dict[pos][1][unshared_index], [3, 1]).T  # bool matrix, three columns

                        # strains_index_unshared_reads is matrix of three columns
                        likelihood_unshared_reads_maximum = calculate_maximum_likelihood_matrix(strains_proportions[pos_reads_mutation_dict[pos][2][unshared_index]], strains_index_unshared_reads)

                        # get the maximum likelihood of each read, then sum
                        if pos == -2:
                            print(strains_proportions[pos_reads_mutation_dict[pos][2][unshared_index]])
                            print(strains_index_unshared_reads)
                            print(likelihood_unshared_reads_maximum)
                        likelihood_strain_in_sample += likelihood_unshared_reads_maximum

                    # total likelihood of this combination
                    if pos == -2:
                        print(likelihood_reads_consistency, likelihood_strain_in_sample)
                    likelihood_six_combination.append(likelihood_reads_consistency + likelihood_strain_in_sample)

                max_combination_index = np.array(likelihood_six_combination).argmax()
                combination_max = genotype_to_haplotype3[max_combination_index]
                probability_six_combination = np.exp(likelihood_six_combination)
                probability_six_combination = probability_six_combination/np.sum(probability_six_combination)
                # whether skip the position
                if whether_skip(probability_six_combination, skip_threshold) and not bin_initial:
                    continue

                posterior_probability_dict.update({pos: probability_six_combination})

                if pos == -2:
                    print(read_shared)
                    print("posterior probability is:")
                    print(probability_six_combination)
                    print("pos is:")
                    print(pos)
                    print(pos_reads_mutation_dict[pos])
                    print("last pos is:")
                    print(pos_last)
                    print(pos_reads_mutation_dict[pos_last])
                    print("haplotype last is:")
                    print(haplotype_last)
                    sys.exit()

                """ update information under combination with max likelihood: combination_max """
                geno_predicted = combination_max[haplotype_last[read_shared[1]]]
                haplotype_last_temp = pos_reads_mutation_dict[pos][1].copy()    # just prototype
                # for reads unshared part, if there are unshared reads
                if len(read_shared[2]) < len(pos_reads_mutation_dict[pos][0]):
                    unshared_index = np.full(len(pos_reads_mutation_dict[pos][0]), True)
                    unshared_index[read_shared[2]] = False
                    combination_stack = np.tile(combination_max, [np.sum(unshared_index), 1])
                    # row number of strains_index_unshared_reads is number of unshared reads; three columns; bool
                    strains_index_unshared_reads = combination_stack == np.tile(
                        pos_reads_mutation_dict[pos][1][unshared_index], [3, 1]).T
                    # the index of maximum in each row is the strain of that read
                    counts_unshared_reads = reads_count_by_strains[pos_reads_mutation_dict[pos][2][
                        unshared_index]] * strains_index_unshared_reads
                    strains_index_unshared_reads = np.argmax(counts_unshared_reads, axis=1)
                    # update haplotype_last_temp which is in the same order of current reads
                    haplotype_last_temp[unshared_index] = strains_index_unshared_reads
                    # update three haplotypes
                    for haplotype_index in range(3):
                        if np.sum(strains_index_unshared_reads == haplotype_index) > 0:
                            temp_index = np.full(len(pos_reads_mutation_dict[pos][1]), False)
                            temp_index[
                                np.where(unshared_index)[0][strains_index_unshared_reads == haplotype_index]] = True
                            haplotypes[haplotype_index].haplotype_extend(temp_index, pos_reads_mutation_dict, pos,
                                                                         combination_max[haplotype_index])

                # for reads shared part
                # strains_shared_reads = haplotype_last[read_shared[1]]   # has its own order
                haplotype_last_temp[read_shared[2]] = haplotype_last[read_shared[1]]

                """ genotype inconsistent part"""
                # if genotype inconsistent exists, modify haplotype_last_temp and haplotype record
                if np.sum(geno_predicted != geno_actual) > 0:
                    # inconsistent_index has own order; index, not bool, index in the current position
                    inconsistent_index = read_shared[2][geno_actual != geno_predicted]
                    # likelihood of strain convert
                    combination_stack = np.tile(combination_max, [len(inconsistent_index), 1])
                    strain_available_a = combination_stack == np.tile(
                        pos_reads_mutation_dict[pos][1][inconsistent_index], [3, 1]).T
                    # get set of strains of last position, set them to unavailable
                    strain_available = np.array([True, True, True])
                    strain_available[np.array(list(set(haplotype_last)))] = False
                    strain_available = np.tile(strain_available, [len(inconsistent_index), 1])

                    # get reads index of possible strain convert, three columns, bool
                    # row number is len(inconsistent_index)
                    reads_converted = strain_available_a * strain_available
                    reads_converted_index = np.sum(reads_converted, axis=1) > 0

                    # if strain converted is possible
                    if np.sum(reads_converted_index) > 0:
                        strain_converted_likelihood = reads_count_by_strains[pos_reads_mutation_dict[pos][2][
                            inconsistent_index[reads_converted_index]]] * reads_converted[reads_converted_index]
                        # update haplotype_last
                        haplotype_last_temp[inconsistent_index[reads_converted_index]] = np.argmax(
                            strain_converted_likelihood, axis=1)
                        for strain_temp, index_a in zip(np.argmax(strain_converted_likelihood, axis=1),
                                                        np.where(reads_converted_index)[0]):
                            # index_b is the single index of read in shared reads
                            index_b = np.where(geno_actual != geno_predicted)[0][index_a]
                            read_positions, read_positions_geno = haplotypes[
                                haplotype_last[read_shared[1][index_b]]].haplotype_subtract(read_shared[2][index_b],
                                                                                            pos_reads_mutation_dict,
                                                                                            pos)
                            haplotypes[strain_temp].haplotype_add(read_shared[2][index_b], pos_reads_mutation_dict, pos,
                                                                  read_positions, read_positions_geno)

                """ genotype consistent part"""
                # if genotype consistent exists, modify haplotype record, but not modify haplotype_last_temp
                if np.sum(geno_predicted == geno_actual) > 0:
                    # consistent_index has own order; index, not bool, index in the current position
                    consistent_index = read_shared[2][geno_actual == geno_predicted]
                    for c_index in consistent_index:
                        haplotypes[haplotype_last_temp[c_index]].haplotype_add(c_index, pos_reads_mutation_dict, pos,
                                                                               [], [])
                pos_last = pos
                haplotype_last = haplotype_last_temp
                reads_count_by_strains[:, 0] = haplotypes[0].reads_count_by_sample.copy()
                reads_count_by_strains[:, 1] = haplotypes[1].reads_count_by_sample.copy()
                reads_count_by_strains[:, 2] = haplotypes[2].reads_count_by_sample.copy()

    return haplotypes, reads_count_by_strains, posterior_probability_dict


def incorporate_variant_sites_into_haplotype(haplotypes, pos_reads_mutation_dict_variant):
    """
    :param haplotypes: haplotypes generated by haplotype_identification (list if length 3), class Haplotype
    :param pos_reads_mutation_dict_variant: {pos: {"read_names", "0-1", "sample index"}}, np.array
    :return: variant_site_reserved, pos_geno_read_number
    """

    variant_site_reserved = []      # store the pos of potential deleted cnv
    # store geno number at each pos covered by reads of haplotypes, {pos:([0 , 1 , .], [10, 11, 0])}
    pos_geno_read_number = {}
    # build a list list of arrays storing read names for each haplotype, default length is 3
    read_names_haplotypes = [np.array(list(i.read_names_sample_index.keys())) for i in haplotypes]

    # process each pos of variant sites
    for pos in sorted(pos_reads_mutation_dict_variant.keys()):
        pos_h_geno = []             # len=3, like [0 , 1 , .]
        pos_h_number = []           # len=3, like [10, 11, 0], record number of bases of haloptypes
        # get read names at a variant site
        variant_read_names_array = pos_reads_mutation_dict_variant[pos][0]
        for h_index, read_names_haplotype in enumerate(read_names_haplotypes):
            ###############################################################
            # check whether reads of pos overlap with reads of haplotypes #
            ###############################################################
            read_names_shared, read_names_variant_index, read_names_haplotype_index = np.intersect1d(
                variant_read_names_array, read_names_haplotype, assume_unique=False, return_indices=True)

            # if no shared reads with this haplotype with index h_index
            if len(read_names_shared) == 0:
                pos_h_geno.append(".")
                pos_h_number.append(0)
            # if there are shared reads with this haplotype with index h_index
            else:
                genotypes_shared = pos_reads_mutation_dict_variant[pos][1][read_names_variant_index]
                geno_0_num = np.sum(genotypes_shared == 0)
                geno_1_num = np.sum(genotypes_shared == 1)
                if geno_0_num > geno_1_num:
                    if geno_0_num > 9*geno_1_num or geno_1_num < 2:
                        pos_h_geno.append(0)
                        pos_h_number.append(geno_0_num)
                    else:
                        pos_h_geno.append(".")
                        pos_h_number.append(len(genotypes_shared))
                elif geno_1_num > geno_0_num:
                    if geno_1_num > 9*geno_0_num or geno_0_num < 2:
                        pos_h_geno.append(1)
                        pos_h_number.append(geno_1_num)
                    else:
                        pos_h_geno.append(".")
                        pos_h_number.append(len(genotypes_shared))
                else:
                    pos_h_geno.append(".")
                    pos_h_number.append(len(genotypes_shared))

        #####################
        # summarize the pos #
        #####################
        if sum(pos_h_number) == 0:
            # if no shared reads, store this pos as variant site for further analysis
            variant_site_reserved.append(pos)
        else:
            if set(pos_h_geno) != {"."}:
                # if at least a haplotype pass the filter
                pos_geno_read_number.update({pos: (pos_h_geno, pos_h_number)})
    return variant_site_reserved, pos_geno_read_number


def geno_depth_refine(read_geno_sample_index, pos):
    """
    # Refine the {pos: {["read_names", "0-1", "sample index"}]} information
    :param read_geno_sample_index:
    :param pos:
    :return:
    """
    geno_0_index = read_geno_sample_index[1] == 0
    geno_1_index = read_geno_sample_index[1] == 1
    geno_01_number = np.array([np.sum(geno_0_index), np.sum(geno_1_index)])
    geno_max_index = np.argmax(geno_01_number)
    if geno_01_number[geno_max_index] >= 0.9*len(read_geno_sample_index[1]) or geno_01_number[1-geno_max_index] <= 1:
        index_temp = [geno_0_index, geno_1_index][geno_max_index]
        return {pos: {"geno": [0, 1][geno_max_index], "sample_index": read_geno_sample_index[2][index_temp], "read_name": read_geno_sample_index[0][index_temp]}}
    else:
        return None


def infer_isolated_variant_sites(chr_pos_read_names_variant, chr_pos_genotype, chromosome, lambda_matrix, threshold_prob=0.99):
    """
    # infer isolated variant sites for the whole chromosomel
    :param chr_pos_read_names_variant: {pos: {["read_names", "0-1", "sample index"}]}, np.array
    :param chr_pos_genotype: {pos: {0:"A", 1:"G"}}
    :param lambda_matrix: 3 by M numpy matrix
    :param threshold_prob
    :return: a dict of record, {pos: "NZ_AP012323.1 2140    C   T   nan nan nan nan nan nan 0   1763    1   468 1   53"}
    """
    sample_size = lambda_matrix.shape[1]
    likelihood_ratio_threshold = np.log(threshold_prob/(1-threshold_prob))
    # build the dict storing refined results, {pos: {"geno": , "sample_index": , "read_name": }}
    pos_geno_depth_dict = {}
    #########################################
    # step1: divide pos into potential zone #
    #########################################
    pos_zones = []
    current_pos_list = []
    last_pos = -1
    for pos in sorted(chr_pos_read_names_variant.keys()):
        refined_record = geno_depth_refine(chr_pos_read_names_variant[pos], pos)

        # whether refined_record is effective
        if refined_record:
            pos_geno_depth_dict.update(refined_record)
            if last_pos == -1:
                # start of zones
                current_pos_list.append(pos)
                last_pos = pos
            else:
                # not a start
                read_names_last_pos = pos_geno_depth_dict[last_pos]["read_name"]
                read_names_current_pos = refined_record[pos]["read_name"]
                reads_shared = np.intersect1d(read_names_last_pos, read_names_current_pos, assume_unique=False, return_indices=True)
                if len(reads_shared[0]) > 0:
                    # there are shared reads
                    current_pos_list.append(pos)
                    last_pos = pos
                else:
                    # no shared reads
                    pos_zones.append(current_pos_list)
                    current_pos_list = [pos]
                    last_pos = pos
    # add the last zone
    if len(current_pos_list) > 0:
        pos_zones.append(current_pos_list)

    ###########################################
    # step2: infer which haplotype belongs to #
    ###########################################
    variant_output_dict = {}
    lambda_matrix = np.maximum(lambda_matrix, 0.001)        # avoid zero
    lambda_sum = np.sum(lambda_matrix, axis=1)              # length of 3
    lambda2_sum_minus_lambda0_sum = lambda_sum[2] - lambda_sum[0]
    lambda2_sum_minus_lambda1_sum = lambda_sum[2] - lambda_sum[1]
    lambda_matrix_log = np.log(lambda_matrix)               # 3 by M
    lambda0_log_minus_lambda2_log = lambda_matrix_log[0] - lambda_matrix_log[2]
    lambda1_log_minus_lambda2_log = lambda_matrix_log[1] - lambda_matrix_log[2]
    for pos_list in pos_zones:
        if len(pos_list) == 0:
            continue
        samples_depth = []
        for pos in pos_list:
            # get the samples depth at pos
            samples_depth.append(list(np.bincount(pos_geno_depth_dict[pos]["sample_index"], weights=None, minlength=sample_size)))
        # get mean depth of poses for each sample in this region, one dim array, length M
        samples_depth_mean = np.mean(np.array(samples_depth), axis=0)
        # log likelihood ratio, h2 = 0
        h_likelihood_ratio_dict = dict({0: 0, 1: 0, 2: 0})
        h_likelihood_ratio_dict[0] = np.sum(samples_depth_mean * lambda0_log_minus_lambda2_log) + lambda2_sum_minus_lambda0_sum
        h_likelihood_ratio_dict[1] = np.sum(samples_depth_mean * lambda1_log_minus_lambda2_log) + lambda2_sum_minus_lambda1_sum
        # sort haplotypes likelihood ratio
        h_likelihood_ratio = list(h_likelihood_ratio_dict.items())
        h_likelihood_ratio.sort(key=lambda x: x[1], reverse=True)
        # ratio of the largest and second largest
        ratio_difference = h_likelihood_ratio[0][1] - h_likelihood_ratio[1][1]
        max_h = h_likelihood_ratio[0][0]
        if ratio_difference >= likelihood_ratio_threshold:
            for pos in pos_list:
                one_line_record = chromosome + "\t" + str(pos + 1) + "\t"
                geno_dict_temp = chr_pos_genotype[pos]
                one_line_record += geno_dict_temp[0] + "\t" + geno_dict_temp[1] + "\t" + str(h_likelihood_ratio_dict[0]) + "\t" + str(h_likelihood_ratio_dict[1]) + "\t" + "0" + "\t"
                one_line_record += "-9\t-9\t-9\t"
                geno_depth_temp = [".", ".", ".", ".", ".", "."]
                geno_depth_temp[max_h*2] = str(pos_geno_depth_dict[pos]["geno"])
                geno_depth_temp[max_h * 2 + 1] = str(len(pos_geno_depth_dict[pos]["sample_index"]))
                one_line_record += "\t".join(geno_depth_temp) + "\n"
                variant_output_dict.update({pos: one_line_record})

    return variant_output_dict


""" assume test data have four samples """
test_data = {3: [np.array(["read1", "read2", "read3", "read4", "read5", "read6", "read7", "read8", "read9", "read10"]),
                 np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 1]), np.array([0, 0, 0, 1, 1, 1, 2, 3, 3, 3])],
             10: [np.array(["read1", "read2", "read3", "read4", "read5", "read6", "read7", "read11", "read12", "read8",
                            "read9", "read10"]),
                  np.array([0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0]),
                  np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3])],
             20: [np.array(["read13", "read14", "read15", "read16", "read17", "read18", "read19", "read20", "read21",
                            "read8", "read9", "read10"]),
                  np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1]),
                  np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3])]}

# test_haplotypes, test_reads_count_by_strains = haplotype_identification(test_data, 4)

