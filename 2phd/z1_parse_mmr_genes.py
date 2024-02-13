#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conradizydorczyk@outlook.com>
Date   : 2023-12-10
Purpose: Parse VCF for mmr mutations
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Parse VCF for mmr mutations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v',
                        '--input_vcf',
                        help='Input vcf file with MMR mutations from bedtools intersect',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-i',
                        '--isolate',
                        help='isolate',
                        metavar='str',
                        type=str,
                        required=True)


    return parser.parse_args()


# --------------------------------------------------
def main():
    """Parse MMR genes"""

    args = get_args()
    
    with open(args.input_vcf, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            gene = getGene(int(line_elements[1]))
            
            annotation = line_elements[7].split(";")[7].split("|")
            
            if annotation[2] != "LOW":
                outstring = f"{args.isolate}\t{annotation[1]}\t{gene}\t{annotation[3]}\t{annotation[9]}\t{annotation[10]}\t{annotation[12]}\t{annotation[13]}"
                print(outstring)



# --------------------------------------------------
# Other functions:

def getGene(coord):
    genes_dict = {
        "pfpI":[399493, 400032],
        "mutM":[401131, 401943],
        "ung":[818003, 818698],
        "dnaQ":[1973470, 1974210],
        "mfd":[3360875, 3364321],
        "mutS":[4054525, 4057092],
        "sodB":[4893697, 4894278],
        "mutT":[4930748, 4931695],
        "sodM":[4997439, 4998050],
        "radA":[5167284, 5168645],
        "mutL":[5549780, 5551681],
        "mutY":[5795954, 5797021],
        "oxyR":[6012047, 6012979],
        "uvrD":[6131088, 6133274],
        "polA":[6183784, 6186525]
    }

    for gene in genes_dict:
        if coord >= genes_dict[gene][0] and coord <= genes_dict[gene][1]:
            gene_ = gene
    
    return(gene_)

# --------------------------------------------------
if __name__ == '__main__':
    main()


    