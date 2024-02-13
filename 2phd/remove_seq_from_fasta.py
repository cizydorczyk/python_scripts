#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-08-24
Purpose: Remove a sequence by name from a fasta file.
"""

import argparse
from Bio import SeqIO


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Remove a sequence by name from a fasta file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input_fasta',
                        help='Input fasta file from which to remove sequence.',
                        metavar='str',
                        type=str,
                        )
    
    parser.add_argument('-s',
                        '--sequence',
                        help='Name of sequence(s) to remove. Comma separated, no spaces if multiple seq.',
                        metavar='str',
                        type=str,
                        )
    
    parser.add_argument('-o',
                        '--output_fasta',
                        help='Output fasta with sequence(s) removed.',
                        metavar='str',
                        type=str,
                        )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    fasta_seq = list(SeqIO.parse(args.input_fasta, "fasta"))

    print(f'Original number of sequences: {len(fasta_seq)}')

    seq_to_remove = args.sequence.strip().split(",")

    print(f'Removing seqs: {seq_to_remove}')

    seq_to_keep = []

    for seq in fasta_seq:
        if seq.name not in seq_to_remove:
            seq_to_keep.append(seq)
        else:
            continue

    print(f'Number of seq after removal: {len(seq_to_keep)}')
    
    SeqIO.write(seq_to_keep, args.output_fasta, "fasta")


# --------------------------------------------------
if __name__ == '__main__':
    main()
