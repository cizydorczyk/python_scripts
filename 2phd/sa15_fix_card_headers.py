#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-10-04
Purpose: Fix card db headers for srst2.
"""

import argparse
from Bio import SeqIO


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Fix card db headers for srst2.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-db',
                        '--card_db',
                        help='Card protein homolog nucleotide fasta file.',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output_fasta',
                        help='output fasta file usable with SRST2.',
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    fasta_seqs = list(SeqIO.parse(args.card_db, "fasta"))

    tr_table = str.maketrans(dict.fromkeys("!@#$%^&*();:,<>/?_=`~|+.'", "-"))

    seq_counter = 0
    for seq in fasta_seqs:
        seq_counter += 1
        seq_header = seq.id.split("|")
        seq_shortname = seq_header[5].split(" ")[0]
        tr_seq_shortname = seq_shortname.translate(tr_table)

        new_seq_id = f"{seq_counter}__{tr_seq_shortname}__{tr_seq_shortname}__{seq_counter}"
        seq.id = new_seq_id

        seq.description = ""
    
    SeqIO.write(fasta_seqs, args.output_fasta, "fasta")


# --------------------------------------------------
if __name__ == '__main__':
    main()
