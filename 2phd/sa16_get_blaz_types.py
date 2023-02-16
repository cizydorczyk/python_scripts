#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-11-02
Purpose: Get blaZ sequence types from nucl sequences.
"""

import argparse
from Bio import SeqIO


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Get blaZ types from nuc sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-i',
                        '--input',
                        help='Input blaZ sequence alignment',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output',
                        help='Output file with blaZ types',
                        metavar='str',
                        type=str,
                        required=True)


    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    blaz_seqs = list(SeqIO.parse(args.input, "fasta"))

    complete_blaz_seqs = []
    incomplete_blaz_seqs = []
    blaz_types = {}

    for seq in blaz_seqs:

        seq_len = len(seq.seq.replace("-", ""))

        codon_1 = str(seq.seq[354:357].upper())
        codon_2 = str(seq.seq[618:621].upper())
        aa_1 = getCodon(codon_1)
        aa_2 = getCodon(codon_2)
        blaz_type = getBlazType(aa_1, aa_2)
    
        blaz_types[seq.id] = f'{seq.id}\t{blaz_type}\t{seq_len}\t{codon_1}\t{aa_1}\t{codon_2}\t{aa_2}'

    header = [f'ID\tBlaz_Type\tSeq_Length\tCodon_1\tAA_1\tCodon_2\tAA_2']
    to_write = [blaz_types[key] for key in blaz_types]
    to_write = header + to_write
    to_write = "\n".join(to_write)
    with open(args.output, "w") as outfile1:
        outfile1.write(to_write)

# --------------------------------------------------
def getCodon(codon):
    aa_dict = {
        "GGG":"G",
        "GGA":"G",
        "GGC":"G",
        "GGT":"G",
        "GAG":"E",
        "GAA":"E",
        "GAC":"D",
        "GAT":"D",
        "GCG":"A",
        "GCA":"A",
        "GCC":"A",
        "GCT":"A",
        "GTG":"V",
        "GTA":"V",
        "GTC":"V",
        "GTT":"V",
        "AGG":"R",
        "AGA":"R",
        "AGC":"S",
        "AGT":"S",
        "AAG":"K",
        "AAA":"K",
        "AAC":"N",
        "AAT":"N",
        "ACG":"T",
        "ACA":"T",
        "ACC":"T",
        "ACT":"T",
        "ATG":"M",
        "ATA":"I",
        "ATC":"I",
        "ATT":"I",
        "CGG":"R",
        "CGA":"R",
        "CGC":"R",
        "CGT":"R",
        "CAG":"Q",
        "CAA":"Q",
        "CAC":"H",
        "CAT":"H",
        "CCG":"P",
        "CCA":"P",
        "CCC":"P",
        "CCT":"P",
        "CTG":"L",
        "CTA":"L",
        "CTC":"L",
        "CTT":"L",
        "TGG":"W",
        "TGA":"*",
        "TGC":"C",
        "TGT":"C",
        "TAG":"*",
        "TAA":"*",
        "TAC":"Y",
        "TAT":"Y",
        "TCG":"S",
        "TCA":"S",
        "TCC":"S",
        "TCT":"S",
        "TTG":"L",
        "TTA":"L",
        "TTC":"F",
        "TTT":"F",
    }
    if len(codon) != 3:
        return("--")
    
    if "-" in codon:
        return("--")
    elif "-" not in codon:
        return(aa_dict[codon])

def getBlazType(aa1, aa2):
    
    blaz_type = ""
    if aa1 == "T" and aa2 == "S":
        blaz_type = "A"
    elif aa1 == "K" and aa2 == "N":
        blaz_type = "B"
    elif aa1 == "T" and aa2 == "N":
        blaz_type = "C"
    elif aa1 == "A" and aa2 == "S":
        blaz_type = "D"
    else:
        blaz_type = "Unknown"
    
    return(blaz_type)

# --------------------------------------------------
if __name__ == '__main__':
    main()
