#!/usr/bin/env python3
"""
Author : conrad <conrad@localhost>
Date   : 2023-06-23
Purpose: Rock the Casbah
"""

import argparse
import os.path
import pandas as pd
from itertools import chain


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Rock the Casbah',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_rgps',
                        help='ppanggolin plastic_regions.tsv',
                        type=str
                        )
    
    parser.add_argument('--input_spots',
                        help='ppanggolin spots.tsv',
                        type=str
                        )
    
    parser.add_argument('--input_summarize_spots',
                        help='ppanggolin summarize_spots.tsv',
                        type=str
                        )
    
    parser.add_argument('--gffs_path',
                        help='path to folder with ALL gffs used to run ppanggolin',
                        type=str
                        )
    
    parser.add_argument('--panaroo_presabs',
                        help='panaroo gene_presence_absence.csv file',
                        type=str
                        )
    
    parser.add_argument('--output_presabs',
                        help='output binary presence absence file with genes in RGPs noted',
                        type=str
                        )
    

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    # read rgps:
    rgps_dict = {}
    with open(args.input_rgps, "r") as infile1:
        for line in infile1:
            if not line.startswith("region"):
                line_elements = line.strip().split("\t")
                isolate = line_elements[1]
                rgp = line_elements[0]
                rgp_contig = line_elements[2]
                rgp_start = line_elements[3]
                rgp_stop = line_elements[4]

                rgp_entry = [rgp, rgp_contig, rgp_start, rgp_stop]
                
                if isolate not in rgps_dict:
                    rgps_dict[isolate] = [rgp_entry]
                else:
                    rgps_dict[isolate].append(rgp_entry)
    
    # read gffs:
    parsed_gffs = {}
    for isolate in rgps_dict:
        #parsed_gff = {}
        isolate_gff = os.path.join(args.gffs_path, isolate + "_panaroo.gff")
        raw_file = open(isolate_gff, "r")
        
        # read in and split off fasta portion:
        lines = raw_file.read().replace(",", "") # copying from panaroo code here (post_run_gff_output.py)
        split = lines.split("##FASTA")

        # if len(split) > 1:
        #     parsed_gff["fasta"] = split[1]
        # else:
        #     parsed_gff["fasta"] = None
        
        # get header and body:
        header = []
        body = []
        for line in split[0].splitlines():
            if "##" in line:
                header.append(line)
            else:
                body.append(line)
        
        #parsed_gff["header"] = "\n".join(header)
        #parsed_gff["body"] = body

        parsed_gffs[isolate] = body
    
    # get genes in rgps for each isolate:
    rgp_genes_dict = {}
    for isolate in rgps_dict:
        rgp_genes = []
        for rgp in rgps_dict[isolate]:
            rgp_start = int(rgp[2])
            rgp_stop = int(rgp[3])
            rgp_contig = rgp[1]

            #gff_body = parsed_gffs[isolate]["body"]
            for line in parsed_gffs[isolate]:
                line_elements = line.strip().split("\t")
                contig = line_elements[0]
                gene_start = int(line_elements[3])
                gene_stop = int(line_elements[4])
                gene_annot_type = line_elements[2]

                if gene_annot_type == "CDS": # ppanggolin reports only CDS entries
                    if contig == rgp_contig:
                        if gene_start >= rgp_start and gene_start <= rgp_stop:# or (gene_stop >= rgp_start and gene_stop <= rgp_stop):
                            if gene_stop >= rgp_start and gene_stop <= rgp_stop:
                                rgp_genes.append(line_elements[8].split(";")[0].split("=")[1])
        rgp_genes_dict[isolate] = rgp_genes
    
    # read panaroo presence/absence:
    panaroo_df = pd.read_csv(filepath_or_buffer=args.panaroo_presabs, 
                             sep=",",
                             index_col=0,
                             header=0)

    panaroo_dict = panaroo_df.to_dict(orient='list')
    del panaroo_dict["Non-unique Gene name"], panaroo_dict["Annotation"] # remove non isolate columns

    # convert pres/abs data to binary data for in RGP or not
    rgp_binary_dict = {isolate:[] for isolate in panaroo_dict}
    
    for isolate in panaroo_dict:
        for gene in panaroo_dict[isolate]:
            if pd.notnull(gene):
                genes = gene.strip().split(";")
                if len(genes) == 1:
                    if genes[0] in rgp_genes_dict[isolate]:
                        rgp_binary_dict[isolate].append(1)
                    else:
                        rgp_binary_dict[isolate].append(0)
                elif len(genes) > 1:
                    # for gene2 in genes: # uncomment this code to get the same #s of genes as produced by ppanggolin
                    #     if gene2 in rgp_genes_dict[isolate]:
                    #         rgp_binary_dict[isolate].append(1)
                    #         print(isolate, gene2)
                    #     else:
                    #         # shouldn't ever trigger ## why shouldn't it trigger? I don't know why i put this...
                    #         print("Check data", isolate, gene2)    
                    if any(gene_ in rgp_genes_dict[isolate] for gene_ in genes):
                        rgp_binary_dict[isolate].append(1)
                    else:
                        rgp_binary_dict[isolate].append(0)
            elif pd.isnull(gene):
                rgp_binary_dict[isolate].append(0)


    rgp_binary_df = pd.DataFrame(rgp_binary_dict)
    rgp_binary_df = rgp_binary_df.set_index(panaroo_df.index)

    # Write output df:
    rgp_binary_df.to_csv(args.output_presabs, sep="\t", header=True)




# --------------------------------------------------
if __name__ == '__main__':
    main()
