#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@outlook.com>
Date   : 2024-02-01
Purpose: Parse hypermutation genes for LES
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Parse hypermutation genes for LES',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v',
                        '--input_vcf',
                        help='input vcf',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-i',
                        '--isolate',
                        help='isolate name',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output_vcf',
                        help='output vcf; can exist already!',
                        metavar='str',
                        type=str,
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    hyper_mutations = []
    with open(args.input_vcf, "r") as infile1:
        for line in infile1:
            if not line.startswith("#"):
                line_elements = line.strip().split("\t")
                gene = getGene(int(line_elements[1]))

                if gene != "":
            
                    annotation = line_elements[7].split(";")[7].split("|")
                    mutation = f"{annotation[3]}_{annotation[1]}_{gene}_{annotation[10]}"

                    outstring = f"{args.isolate}\t{mutation}\t{annotation[1]}\t{annotation[2]}\t{gene}\t{annotation[3]}\t{annotation[9]}\t{annotation[10]}\t{annotation[12]}\t{annotation[13]}\n"

                    hyper_mutations.append(outstring)
    
    header = "isolate\tmutation\tmutation_type\timpact\tlocus_tag\tdna_mutation\taa_mutation\tdna_pos\taa_pos\n"
    with open(args.output_vcf, "a+") as outfile1:
        outfile1.write("".join(hyper_mutations))

def getGene(coord):

    # Hypermutations coordinates from ~/les-manuscript-2024/snp-calling/hypermutation-analysis/lesb58-hypermutations-genes.bed
    genes_dict = {
        "pfpI":[393232, 393771],
        "mutM":[394870, 395682],
        "ung":[5049266, 5049961],
        "dnaQ":[3893872, 3894612],
        "mfd":[2230847, 2234293],
        "mutS":[1523212, 1525779],
        "sodB":[5227405, 5227986],
        "mutT":[5264456, 5265403],
        "sodM":[5331148, 5331759],
        "radA":[5502419, 5503780],
        "mutL":[5882115, 5884016],
        "mutY":[6130482, 6131549],
        "oxyR":[6347994, 6348926],
        "uvrD":[6467002, 6469188],
        "polA":[6519699, 6522440]
    }

    # for gene in genes_dict:
    #     # gene_ = ""
    #     if coord >= genes_dict[gene][0] and coord <= genes_dict[gene][1]:
    #         gene_ = gene
    #     else:
    #         gene_ = ""

    for gene, gene_range in genes_dict.items():
        if gene_range[0] <= coord <= gene_range[1]:
            return gene
    return ""


# --------------------------------------------------
if __name__ == '__main__':
    main()
