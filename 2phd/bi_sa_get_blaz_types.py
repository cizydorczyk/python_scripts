#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-10-12
Purpose: Get blaZ types from nuc sequences.
"""

import argparse
from Bio import SeqIO
import os.path


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

    parser.add_argument('-cd',
                        '--complete_dir',
                        help='Directory in which to write complete blaZ sequences',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-id',
                        '--incomplete_dir',
                        help='Directory in which to write incomplete blaZ sequences',
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

        if seq_len == 846:
            complete_blaz_seqs.append(seq)
        else:
            incomplete_blaz_seqs.append(seq)

        codon_1 = str(seq.seq[355:358].upper())
        codon_2 = str(seq.seq[620:623].upper())
        aa_1 = getCodon(codon_1)
        aa_2 = getCodon(codon_2)
        blaz_type = getBlazType(aa_1, aa_2)

        no_gaps_seq = seq.seq.replace("-", "")
        no_gaps_codon_1 = no_gaps_seq[354:357].upper()
        no_gaps_codon_2 = no_gaps_seq[618:621].upper()
        no_gaps_aa_1 = getCodon(no_gaps_codon_1)
        no_gaps_aa_2 = getCodon(no_gaps_codon_2)
        no_gaps_blaz_type = getBlazType(no_gaps_aa_1, no_gaps_aa_2)

        
        
        blaz_types[seq.id] = f'{seq.id}\t{blaz_type}\t{seq_len}\t{codon_1}\t{aa_1}\t{codon_2}\t{aa_2}'
                            # f'\t{no_gaps_blaz_type}\t{no_gaps_codon_1}\t{no_gaps_aa_1}\t{no_gaps_codon_2}'
                            # f'\t{no_gaps_aa_2}'
                            
    
    for seq in complete_blaz_seqs:
        aligned_handle = os.path.join(args.complete_dir, f'{seq.id}-aligned-blaZ.fasta')
        full_handle = os.path.join(args.complete_dir, f'{seq.id}-blaZ.fasta')
        
        SeqIO.write(seq, aligned_handle, "fasta")

        seq.seq = seq.seq.replace("-", "")
        SeqIO.write(seq, full_handle, "fasta")
    
    for seq in incomplete_blaz_seqs:
        aligned_handle = os.path.join(args.incomplete_dir, f'{seq.id}-aligned-blaZ.fasta')
        full_handle = os.path.join(args.incomplete_dir, f'{seq.id}-blaZ.fasta')
        
        SeqIO.write(seq, aligned_handle, "fasta")

        seq.seq = seq.seq.replace("-", "")
        SeqIO.write(seq, full_handle, "fasta")

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
