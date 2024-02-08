"""
This script convert primary and secondary strains into vcf format
"""


def haplotype2vcf(abs_path_haplotype, abs_path_vcf):
    """

    :param abs_path_haplotype:
    :param abs_path_vcf:
    :return:
    """
    abs_path_haplotype_f = open(abs_path_haplotype, "r")
    abs_path_vcf_f = open(abs_path_vcf, "w")

    for line in abs_path_haplotype_f:
        cols = line.split("\t")




    abs_path_haplotype_f.close()
    abs_path_vcf_f.close()
