#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-10-02
Purpose: Generate a db from NCBI AMRFinderPlus usable withh SRST2
"""

import argparse
from Bio import SeqIO


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Generate a db from NCBI AMRFinderPlus usable withh SRST2',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-db',
                        '--amr_cds',
                        help='AMR_CDS file from AMRFinderPlus db.',
                        metavar='str',
                        type=str,
                        required=True),
    
    parser.add_argument("-r",
                        '--reference_catalog',
                        help="reference catalog from AMRFinderPlus db",
                        type=str,
                        metavar="file",
                        required=True),
    
    parser.add_argument("-o",
                        "--output_fasta",
                        help="output fasta formatted for SRST2",
                        type=str,
                        metavar="file",
                        required=True),


    return parser.parse_args()


# --------------------------------------------------
def main():
    """Change db fasta names to work with SRST2."""

    args = get_args()

    fasta_seqs = list(SeqIO.parse(args.amr_cds, "fasta"))

    family_dict = {}

    for seq in fasta_seqs:
        seq_header = seq.id.split("|")
        seq_family = seq_header[6]
        if seq_family not in family_dict:
            family_dict[seq_family] = [seq]
        elif seq_family in family_dict:
            family_dict[seq_family].append(seq)
    
    output_seqs = []
    family_counter = 0
    for family in family_dict:
        family_counter += 1
        allele_counter = 0
        for seq in family_dict[family]:
            allele_counter += 1
            seq_header = seq.id.split("|")
            new_id = f'{family_counter}__{family}__{seq_header[5]}__{allele_counter}'
            seq.id = new_id
            output_seqs.append(seq)
    
    # Get only AMR genes; remove virulence/stress genes; remove point/susceptibility subtypes of amr genes
    reference_accessions = []
    bad_types = 0
    good_types = 0
    with open(args.reference_catalog, "r") as infile1:
        for line in infile1:
            if not line.startswith("gene_family"):
                line_elements = line.strip().split("\t")
                line_type = line_elements[3]
                line_subtype = line_elements[4]
                accession = line_elements[8]

                if line_type != "AMR":
                    bad_types+=1
                elif line_type == "AMR":
                    if line_subtype != "AMR":
                        bad_types+=1
                    elif line_subtype == "AMR":
                        if accession == "":
                            accession = line.strip().split("\t")[11]
                            reference_accessions.append(accession)
                            good_types +=1
                        else:
                            reference_accessions.append(accession)
                            good_types +=1
    
    print(len(output_seqs))
    
    output_seqs = [seq for seq in output_seqs if seq.description.split("|")[2] in reference_accessions]
    # for seq in output_seqs:
    #     seq.description = ""
    
    SeqIO.write(output_seqs, args.output_fasta, "fasta")


# --------------------------------------------------
if __name__ == '__main__':
    main()
