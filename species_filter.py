#!/usr/bin/python3
# -*- coding: utf-8 -*-


""" Copyright (c) 2019  Boyan Zhou

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.


:Authors: Boyan Zhou
:Contact: boyanzhou1992@gmail.com
:Date: Oct 2023
:Type: tool
:Input: Combined Relative Abundance By LongStrain Step2
:Output: Target Species List After Filtering According to the criteria provided by users

------------------------------------------------------------------------------------------------------------------------
**************
* Update log *
**************
Date:   2023/12/10
By:     Boyan Zhou
"""


import logging
import time
import sys
import numpy as np
import pandas as pd
import json
import argparse
import os


def combine_species_names(species_name_list):
    """
    combine full names as a single string separated by ","
    :param species_name_list: (['d__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola vulgatus',
       'd__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella copri',
       'd__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides uniformis',
       'd__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides stercoris',
       'd__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Faecalibacterium|s__Faecalibacterium prausnitzii'],
      dtype='object', name='Taxon')
    :return: 'Phocaeicola_vulgatus,Prevotella_copri,Bacteroides_uniformis,Bacteroides_stercoris,Faecalibacterium_prausnitzii'
    """
    simplified_species_name_list = []
    for full_species_name in species_name_list:
        # full_species_name = "d__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Faecalibacterium|s__Faecalibacterium prausnitzii"
        simplified_species_name_list.append(full_species_name.split("|")[-1].split("__")[-1])
    simplified_species_name_list = [i.replace(" ", "_") for i in simplified_species_name_list]
    return ",".join(simplified_species_name_list)


def main():
    parser = argparse.ArgumentParser(
        prog="species filter for LongStrain\n",
        description=f"DESCRIPTION\n\n"
                    f"Filter species from the output of Step2 of LongStrain (relative_abundance_aggregation)\n"
                    f"Usually, if you know the target species you are interested in, you can skip this step\n"
                    f"and provide the species name list in the subsequent steps of LongStrain.\n"
                    f"However, if you have no specific target species, you can filter the species\n"
                    f"according to the rank of mean relative abundance or \n"
                    f"the threshold of relative abundance in a specified proportion of samples.")

    # Adding command-line arguments
    parser.add_argument('-i', '--input_file', required=True, type=str,
                        help='Path to the input file: aggregated relative abundance file')
    parser.add_argument('-o', '--output_file', required=True, type=str,
                        help='Path to the output file: contain the concatenated filtered species name')
    # Filter1: by rank of mean relative abundance
    parser.add_argument('-n', '--top_n', type=int, help='Select top N species in the order of mean relative abundance')

    # Filter2: by threshold of mean relative abundance
    parser.add_argument('-m', '--mean_RA', type=float,
                        help='the threshold of mean relative abundance to filter species, e.g. 0.01')

    # Filter3: by threshold of relative abundance and prevalence; Creating a mutually exclusive group
    group1 = parser.add_mutually_exclusive_group(required=False)
    # Adding the mutually exclusive arguments to the group
    group1.add_argument('-r', '--relative_abundance', type=float,
                        help='the threshold of relative abundance to filter species, e.g. 0.05;\n'
                             'should be used together with -p')
    group1.add_argument('-p', '--prevalence', type=float,
                        help='the proportion of samples pass the threshold of relative abundance, e.g. 0.3;\n'
                             'should be used together with -r')

    # Parsing the command-line arguments
    args = parser.parse_args()

    # read the file into numpy, test_filename = "/gpfs/data/lilab/home/zhoub03/Blaser_data/MIME_result/Kraken2_result_summary/MIME_RA_combined_standard_plusPF.txt"
    ra_df = pd.read_table(args.input_file, header=0, index_col=0)
    # Calculate row means
    ra_means = ra_df.mean(axis=1)
    # Get the indices that would sort the row means in decreasing order
    sorted_indices = np.argsort(ra_means)[::-1]
    # Reorder the DataFrame based on the sorted indices
    ra_df_sorted = ra_df.iloc[sorted_indices]

    # output file path of extracted species names
    output_joint_species_name_f = open(args.output_file, "w")

    if args.top_n:
        # 1. if we use the first filter criterion, output top_n species name
        if ra_df_sorted.shape[0] >= args.top_n:
            print(f"Extracting the names of top {args.top_n} species by mean relative abundance.")
            output_joint_species_name_f.write(combine_species_names(ra_df_sorted.index[:args.top_n]))

        else:
            print(f"The total number of species is {ra_df_sorted.shape[0]}, extracting their species names.")
            output_joint_species_name_f.write(combine_species_names(ra_df_sorted.index))
    elif args.mean_RA:
        # 2. if we use the second filter criterion, mean relative abundance threshold
        if np.sum(ra_df_sorted.mean(axis=1) > args.mean_RA) > 0:
            print(f"Extracting the names of species passing the threshold of mean relative abundance {args.mean_RA}.")
            output_joint_species_name_f.write(combine_species_names(ra_df_sorted.index[ra_df_sorted.mean(axis=1) > args.mean_RA]))
        else:
            print(f"No species passes the threshold of mean relative abundance.")
    elif args.relative_abundance:
        # 3. if we use the third filter criterion
        passed_species_bool = (ra_df_sorted > args.relative_abundance).mean(axis=1) >= args.prevalence
        if np.sum(passed_species_bool) > 0:
            print(f"Extracting the names of species passing the threshold of relative abundance "
                  f"{args.relative_abundance} in at least {args.prevalence * 100}% of sample.")
            output_joint_species_name_f.write(
                combine_species_names(ra_df_sorted.index[passed_species_bool]))
        else:
            print(f"No species satisfies the requirement of passing the threshold of relative abundance "
                  f"{args.relative_abundance} in at least {args.prevalence * 100}% of sample.")

    output_joint_species_name_f.close()


if __name__ == "__main__":
    sys.exit(main())
