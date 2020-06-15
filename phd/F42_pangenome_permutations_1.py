import argparse
import pandas as pd
import itertools
import numpy as np
import statistics


parser = argparse.ArgumentParser()

parser.add_argument("--isolate_list", help="isolate list, one isolate per line")
parser.add_argument("--pirate_matrix", help="PIRATE binary presence/absence matrix")
parser.add_argument("--output_matrix", help="output per patient binary presence absence matrix", default="")
parser.add_argument("--threshold", help="core genome threshold in # of isolates gene must be present into be considered core", type=int)

args = parser.parse_args()

# isolate_list = []
# with open(args.isolate_list, 'r') as infile1:
#     for line in infile1:
#         isolate_list.append(line.strip())

patient_st_lists = [["A013_H12_09_07_2003_formatted"], ["A034_H64_14_09_2011_formatted"], ["A037_H11_04_08_2004_formatted", "A037_H140_24_03_2010_formatted",
                    "A037_H25_15_11_2006_formatted"], ["A037_H55_19_01_2011_formatted"], ["A043_H61_13_07_2011_formatted"], ["A058_H09_28_09_2012_formatted"],
                    ["A059_H146_25_03_2012_formatted"], ["A076_H46_02_11_2009_formatted"], ["A077_H138_28_05_2009_formatted", "A077_H36_15_04_2009_formatted"],
                    ["A090_H183_08_02_2016_formatted"], ["A098_H26_26_11_2003_formatted"], ["A132_H67_15_11_2011_formatted"], ["A143_H152_23_10_2013_formatted",
                    "A143_H157_17_09_2014_formatted", "A143_H176_18_11_2015_formatted", "A143_H185_16_03_2016_formatted", "A143_H56_26_01_2011_formatted",
                    "A143_H74_31_10_2012_formatted", "A143_H80_01_12_2004_formatted"], ["A274_H23_12_07_2006_formatted"], ["A274_H82_25_04_2005_formatted"],
                    ["A276_H42_02_09_2009_formatted"], ["A290_H02_10_09_2012_formatted"], ["A295_H24_12_07_2006_formatted"], ["A300_H135_24_05_2006_formatted"],
                    ["A300_H44_04_02_2004_formatted"], ["A312_H174_29_06_2015_formatted"], ["A319_H179_18_01_2016_formatted"], ["A319_H65_21_09_2011_formatted"],
                    ["A323_H173_21_03_2011_formatted", "A323_H60_18_06_2011_formatted"], ["A335_H158_05_11_2014_formatted"], ["A337_H38_10_06_2009_formatted"],
                    ["A360_H04_26_03_2012_formatted"], ["A360_H57_09_02_2011_formatted"], ["A366_H151_08_07_2013_formatted"], ["A366_H58_13_04_2011_formatted"],
                    ["A367_H13_19_11_2012_formatted", "A367_H144_28_09_2011_formatted", "A367_H153_25_11_2013_formatted", "A367_H159_03_12_2014_formatted",
                     "A367_H182_27_01_2016_formatted"], ["A367_H28_21_07_2008_formatted"], ["A367_H51_01_09_2010_formatted"], ["A370_H14_21_11_2012_formatted",
                    "A370_H149_25_03_2013_formatted", "A370_H166_27_05_2015_formatted", "A370_H186_20_04_2016_formatted"], ["A406_H162_20_04_2015_formatted"],
                    ["A407_H184_02_03_2016_formatted"]]

permuted_lists = list(itertools.product(*patient_st_lists))

# Read original PIRATE matrix:
pres_abs_df = pd.read_table(args.pirate_matrix, sep="\t", header=0, index_col=0)

core_genome_counts = []
accessory_genome_counts = []
for plist in permuted_lists:
    core_gene_count = 0
    acc_gene_count = 0
    p_df = pres_abs_df[list(plist)]

    for index, row in p_df.iterrows():
        row_list = list(row)
        if sum(row_list) >= args.threshold:
            core_gene_count += 1
        elif sum(row_list) < args.threshold:
            acc_gene_count += 1

    core_genome_counts.append(core_gene_count)
    accessory_genome_counts.append(acc_gene_count)

print("Number of permutations: ", len(core_genome_counts))
print("Mean # of core genes: ", statistics.mean(core_genome_counts))
print("Median # of core genes: ", statistics.median(core_genome_counts))
print("Q1: ", np.quantile(core_genome_counts, q=0.25))
print("Q3: ", np.quantile(core_genome_counts, q=0.75))
print("IQR: ", np.quantile(core_genome_counts, q=0.75) - np.quantile(core_genome_counts, q=0.25))
print("Min: ", min(core_genome_counts))
print("Max: ", max(core_genome_counts))

print("Mean # of accessory genes: ", statistics.mean(accessory_genome_counts))
print("Median # of accessory genes: ", statistics.median(accessory_genome_counts))
print("Q1: ", np.quantile(accessory_genome_counts, q=0.25))
print("Q3: ", np.quantile(accessory_genome_counts, q=0.75))
print("IQR: ", np.quantile(accessory_genome_counts, q=0.75) - np.quantile(accessory_genome_counts, q=0.25))
print("Min: ", min(accessory_genome_counts))
print("Max: ", max(accessory_genome_counts))

# Compute core/accessory genome sizes for all isolates/patients, counting ALL isolates
all_iso_core_genome_size = 0
all_iso_acc_genome_size = 0
for index, row in pres_abs_df.iterrows():
    row_list = list(row)
    if sum(row_list) >= 51:
        all_iso_core_genome_size += 1
    elif sum(row_list) < 51:
        all_iso_acc_genome_size += 1

print("All isolate core genome size: ", all_iso_core_genome_size)
print("All isolate accessory genome size: ", all_iso_acc_genome_size)
