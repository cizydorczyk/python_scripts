import os
import argparse
from Bio import SeqIO

def RemoveReference(seq_to_remove, fasta, output_fasta):
    print(fasta)

    records = list(SeqIO.parse(fasta, "fasta"))
    records2 = [i for i in records if i.id not in seq_to_remove.split(",")]

    with open(output_fasta, "w") as outfile:
	       SeqIO.write(records2, outfile, "fasta")
