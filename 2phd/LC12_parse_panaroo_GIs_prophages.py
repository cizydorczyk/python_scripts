#!/usr/bin/env python3
"""
Author : conrad <conrad@localhost>
Date   : 2023-06-27
Purpose: Get LESB58 GIs/prophages & note which panaroo genes fall in them.
"""

import argparse
import pandas as pd
import os.path


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Rock the Casbah',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--input_bed",
                        help="lesb58 bed file with prophages & gis",
                        required=True,
                        type=str,
                        metavar='str')
    
    parser.add_argument("--panaroo_presabs",
                        help="panaroo gene presence absence CSV file",
                        required=True,
                        type=str,
                        metavar='str')
    
    parser.add_argument('--lesb58_gff',
                        help='lesb58 gff file used to run panaroo',
                        type=str
                        )
    
    parser.add_argument('--output_genes',
                        help='output file to write genes in prophages/GIs to',
                        type=str
                        )

    parser.add_argument('--output_families',
                        help='output file to write families in prophages/GIs to',
                        type=str
                        )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    # read and parse lesb58 prophages & GIs
    # input file should be:
    # chrom\tname\tstart\tstop
    lesb58_elements = {}
    with open(args.input_bed, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            element_name = line_elements[1]
            element_start = line_elements[2]
            element_stop = line_elements[3]
            if element_name not in lesb58_elements:
                lesb58_elements[element_name] = [element_start, element_stop]
    
    # read lesb58 gff
    raw_file = open(args.lesb58_gff, "r")
    
    # read in and split off fasta portion:
    lines = raw_file.read().replace(",", "") # copying from panaroo code here (post_run_gff_output.py)
    split = lines.split("##FASTA")

    body = []
    for line in split[0].splitlines():
        if "##" in line:
            continue
        else:
            body.append(line)
    
    # get LESB58 gi/prophage genes:
    les_gi_prophage_genes = []
    for line in body:
        line_elements = line.split("\t")
        gene_start = line_elements[3]
        gene_stop = line_elements[4]
        gene_id = line_elements[8].split(";")[0].split("=")[1]

        for element in lesb58_elements:
            element_start = lesb58_elements[element][0]
            element_stop = lesb58_elements[element][1]

            if gene_start >= element_start and gene_start <= element_stop:
                if gene_stop >= element_start and gene_stop <= element_stop:
                    if gene_id not in les_gi_prophage_genes:
                        les_gi_prophage_genes.append(gene_id)


    # read panaroo presence/absence:
    panaroo_df = pd.read_csv(filepath_or_buffer=args.panaroo_presabs, 
                             sep=",",
                             index_col=0,
                             header=0)

    panaroo_dict = panaroo_df.to_dict(orient='list')
    del panaroo_dict["Non-unique Gene name"], panaroo_dict["Annotation"] # remove non isolate columns

    # 
    gi_prophage_families = []
    panaroo_families = panaroo_df.index
    for gene in panaroo_dict["LESB58"]:
        if pd.notnull(gene):
            genes = gene.strip().split(";")
            if len(genes) == 1:
                if genes[0] in les_gi_prophage_genes:
                    gene_index = panaroo_dict["LESB58"].index(gene)
                    gene_family = panaroo_families[gene_index]
                    gi_prophage_families.append(gene_family)
            elif len(genes) > 1:    
                if any(gene_ in panaroo_dict["LESB58"] for gene_ in genes):
                    gene_index = panaroo_dict["LESB58"].index(gene)
                    gene_family = panaroo_families[gene_index]
                    gi_prophage_families.append(gene_family)
    
    print(f"\nNumber of LESB58 genes in GIs/Prophages:\t{len(les_gi_prophage_genes)}\n"
          f"Number of corresponding Panaroo gene families:\t{len(gi_prophage_families)}")

    with open(args.output_genes, "w") as outfile1:
        outfile1.write("\n".join(les_gi_prophage_genes))

    with open(args.output_families, "w") as outfile2:
        outfile2.write("\n".join(gi_prophage_families))

# --------------------------------------------------
if __name__ == '__main__':
    main()
