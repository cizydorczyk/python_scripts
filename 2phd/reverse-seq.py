from Bio import SeqIO

fasta_file = SeqIO.read("/home/conrad/pa-nanopore/structural-variant-analysis/A180-P2360-15-03-1989-ONT.fasta", "fasta")

def my_function(x):
  return x[::-1]

fasta_file.seq = my_function(fasta_file.seq)

SeqIO.write(fasta_file, "/home/conrad/pa-nanopore/structural-variant-analysis/P2360-reversed/P2360-reversed.fasta", "fasta")