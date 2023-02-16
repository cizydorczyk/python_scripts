from Bio import SeqIO


sequences = list(SeqIO.parse("/home/conrad/sa-ncfb/sa11-srst2/AMR_CDS.fasta", "fasta"))

for sequence in sequences:
    sequence_id = sequence.id.strip().split("|")[1]
    sequence.id = sequence_id
    sequence.description = ""

SeqIO.write(sequences, "/home/conrad/sa-ncfb/sa11-srst2/AMR_CDS_fixed.fasta", "fasta")