from sys import argv
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()

# Required arguments:
parser.add_argument("--snps", help="snp alignment to characterize; VCF format from snp-sites")
parser.add_argument("--num_isolates", help="total number of isolates in alignment")
parser.add_argument("--out", help="output file")

args = parser.parse_args()

vcf = pd.read_csv(args.snps, sep='\t', skiprows=3)

# Length of snp alignment:
print("Number of bases in snp alignment: ", len(vcf))

# Get counts of various statistics:
multiallelic_sites = 0
# ambiguous sites includes 'N' and '-':
ambiguous_sites = 0
biallelic_sites = 0
ambiguous_and_multiallelic_sites = 0
ambiguous_only = 0

# Ambiguous counts dict:
amb_dict = {}
for i in range(0,int(args.num_isolates)+1):
    amb_dict[i] = 0

alt_bases = vcf["ALT"].tolist()
for base, index in zip(alt_bases, [i for i in range(0,len(alt_bases))]):
    bases = base.strip().split(',')
    if len(bases) == 1 and "*" not in bases:
        biallelic_sites += 1
    # this next scenario should never happen since snp-sites doesn't output sites that just have an 'N':
    elif len(bases) == 1 and "*" in bases:
        ambiguous_only += 1

    elif len(bases) > 1 and "*" not in bases:
        multiallelic_sites += 1

    elif len(bases) == 2 and "*" in bases:
        ambiguous_sites += 1

        #print(list(vcf.iloc[index]))
        if bases[0] == "*":
            ambiguous_code = 1
            amb_count = 0
            gt_list = list(vcf.iloc[index])[9:]
            for item in gt_list:
                if item == ambiguous_code:
                    amb_count += 1
            amb_dict[amb_count] += 1
            #print("ambiguous code = 1, ambiguous count: ", amb_count)

        elif bases[1] == "*":
            ambiguous_code = 2
            amb_count = 0
            gt_list = list(vcf.iloc[index])[9:]
            for item in gt_list:
                if item == ambiguous_code:
                    amb_count += 1
            amb_dict[amb_count] += 1
            #print("ambiguous code = 2, ambiguous count: ", amb_count)

    elif len(bases) > 2 and "*" in bases:
        ambiguous_and_multiallelic_sites += 1

sum_counts = multiallelic_sites + ambiguous_sites + biallelic_sites + ambiguous_and_multiallelic_sites + ambiguous_only
# Sum counts should be the same as the length of the alignment (i.e. each site should be classified)
print(multiallelic_sites, ambiguous_sites, ambiguous_only, biallelic_sites, ambiguous_and_multiallelic_sites)
print(sum_counts)
print(amb_dict)

to_write = "Number of loci in snp alignment: " + str(len(vcf)) + "\nNumber of biallelic sites:\t" + str(biallelic_sites) + "\nNumber of multiallelic_sites:\t" + str(multiallelic_sites) + "\nNumber of ambigous_sites:\t" + str(ambiguous_sites) + "\nNumber of ambigous only sites:\t" + str(ambiguous_only) + "\nNumber of ambiguous and multiallelic sites:\t" + str(ambiguous_and_multiallelic_sites)
print(to_write)
dict_to_write = ''
for key in amb_dict:
    to_write2 = str(key) + "\t" + str(amb_dict[key]) + '\n'
    dict_to_write = dict_to_write + to_write2

# Write above data to file:
with open(args.out, 'w') as outfile1:
    outfile1.write(to_write + '\n#Ambiguous sites histogram (below):\n' + dict_to_write)

# Plot histogram showing # of isolates with an ambiguos base call (x) vs. number of times it occurs (y):
# pos = np.arange(len(amb_dict.keys()))
# width = 1.0
#
# ax = plt.axes()
# ax.set_xticks(pos + (width / 2))
# ax.set_xticklabels(amb_dict.keys())
#
# plt.bar(amb_dict.keys(), amb_dict.values(), width, color='g')
# plt.show()
