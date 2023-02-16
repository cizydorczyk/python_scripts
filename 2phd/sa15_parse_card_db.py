#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-10-03
Purpose: Parse latest version of CARD protein homolog db and produce version appropriate for running with SRST2.
"""

import argparse
from Bio import SeqIO


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Parse latest version of CARD protein homolog db and produce version appropriate for running with SRST2.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-db',
                        '--card_db',
                        help='Card protein homolog nucleotide fasta file.',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-i',
                        '--aro_index',
                        help='card aro_index.tsv file.',
                        required=True)
    
    parser.add_argument('-o',
                        '--output_fasta',
                        help='output fasta file usable with SRST2.',
                        required=True)


    return parser.parse_args()


# --------------------------------------------------
def main():
    """Convert Card protein homolog db into version compatible with SRST2."""

    args = get_args()
    
    fasta_seqs = list(SeqIO.parse(args.card_db, "fasta"))

    family_allele_dict = {}
    shortname_family_dict = {}
    family_count_dict = {}

    with open(args.aro_index, "r") as infile1:
        for line in infile1:
            if not line.startswith("ARO Accession"):
                line_elements = line.strip().split("\t")
                gene_family = "_".join(line_elements[8].split(" "))
                gene_shortname = line_elements[11]
                
                if gene_shortname not in shortname_family_dict:
                    shortname_family_dict[gene_shortname] = gene_family
                
                elif gene_shortname in shortname_family_dict:
                    print("Error: duplicate gene shortname.")
                    break
                
                if gene_family not in family_allele_dict:
                    family_allele_dict[gene_family] = 0
                elif gene_family in family_allele_dict:
                    continue
    
    fam_count = 0
    for family in family_allele_dict:
        fam_count +=1
        family_count_dict[family] = fam_count
        
    for seq in fasta_seqs:
        seq_header = seq.id.split("|")
        seq_shortname = seq_header[5].split(" ")[0]
        seq_family = shortname_family_dict[seq_shortname]
        seq_family_count = family_count_dict[seq_family]
        family_allele_dict[seq_family] += 1

        tr_table = str.maketrans(dict.fromkeys('!@#$%^&*();:,<>/?_=`~', '-'))
        tr_seq_family = seq_family.translate(tr_table)
        tr_seq_shortname = seq_shortname.translate(tr_table)

        new_seq_id = f"{seq_family_count}__{tr_seq_family}__{tr_seq_shortname}__{family_allele_dict[seq_family]}"
        seq.id = new_seq_id
        seq.description = seq.description.translate(tr_table)
    
    SeqIO.write(fasta_seqs, args.output_fasta, "fasta")





    # print(len(gene_shortname_list), len(set(gene_shortname_list)))
    # import collections
    # print([item for item, count in collections.Counter(gene_shortname_list).items() if count > 1])


    # for seq in fasta_seqs:
    #     seq_header = seq.id.split("|")
    #     seq_family = seq_header[5].split(" ")[0]
    #     if "-" in seq_family:
    #         seq_fam_elements = seq_family.split("-")
    #         # print("-".join(seq_fam_elements[0:-1]), seq_fam_elements[-1])
    #         seq_family = "-".join(seq_fam_elements[0:-1])
    #         seq_allele = seq_fam_elements[-1]

    #         if seq_family not in family_dict:
    #             family_dict[seq_family] = [seq]
    #         elif seq_family in family_dict:
    #             family_dict[seq_family].append(seq)

    #     elif "-" not in seq_family:
    #         if seq_family not in family_dict:
    #             family_dict[seq_family] = [seq]
    #         elif seq_family in family_dict:
    #             family_dict[seq_family].append(seq)
        
    


# --------------------------------------------------
if __name__ == '__main__':
    main()
