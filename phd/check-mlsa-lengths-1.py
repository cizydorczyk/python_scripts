import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument("input_mlsa", help="input mlsa alignment from mlst pipeline")
parser.add_argument("ts_mlst_output", help='Torsten Seemanns mlst output file')
parser.add_argument("output_mlsa", help="output mlsa with all proper length sequences")
parser.add_argument("mlsa_length", type=int, help="expected length of concatenated mlst genes; shorter/longer sequences will be "
                                        "discarded")
parser.add_argument("log_file", help="log file with info on removed alignments")
parser.add_argument("removed_alignments_fasta", help="fasta with removed alignments")

args = parser.parse_args()

orig_seq_count = 0
partial_alignments = 0
missing_gene_alignments = 0
multiple_genes = 0
to_remove_isolates = []
with open(args.ts_mlst_output, 'r') as infile1:
    for line in infile1:
        orig_seq_count += 1
        line_elements = line.strip().split('\t')
        genes_elements = ','.join(line_elements[3:10])
        if "?" in genes_elements:
            to_remove_isolates.append(line_elements[0].split(".")[0])
            partial_alignments += 1
        if "(-)" in genes_elements:
            to_remove_isolates.append(line_elements[0].split(".")[0])
            missing_gene_alignments += 1
        if len(genes_elements.split(',')) > 7:
            to_remove_isolates.append(line_elements[0].split(".")[0])
            multiple_genes += 1


# Keep only sequence records that are equal to expected length with no gaps:
isolate_sequence_records = []
gapped_alignments = 0
to_remove_isolates_set = set(to_remove_isolates)
removed_seq_records = []

for seq_record in SeqIO.parse(args.input_mlsa, "fasta"):
    if seq_record.id in to_remove_isolates_set:
        removed_seq_records.append(seq_record)
    # if len(seq_record.seq) > args.mlsa_length or len(seq_record.seq) < args.mlsa_length:
        # to_remove_isolates.append(seq_record.id)
        # removed_seq_records.append(seq_record)
        # print("Short/long alignment removed: %s" % seq_record.id)
    # elif "-" in seq_record.seq:
    #     # print("Gapped alignment removed: %s" % seq_record.id)
    #     gapped_alignments += 1
    else:
        isolate_sequence_records.append(seq_record)

print("Original  number of sequences: %s" % orig_seq_count)
print("Number of sequences after filtering: %s" % len(isolate_sequence_records))
# print("Gapped alignments: %s" % gapped_alignments)
print("Number of removed alignments: %s" % len(to_remove_isolates))
print("\tMultiple gene alignments: %s" % multiple_genes)
print("\tPartial alignments: %s" % partial_alignments)
print("\tMissing gene alignments: %s" % missing_gene_alignments)

# Write good alignments to file:
SeqIO.write(isolate_sequence_records, args.output_mlsa, "fasta")

# Write log file:
with open(args.log_file, 'w') as outfile1:
    to_write = []
    to_write.append("Original  number of sequences: %s" % orig_seq_count)
    to_write.append("Number of sequences after filtering: %s" % len(isolate_sequence_records))
    to_write.append("Number of removed alignments: %s" % len(to_remove_isolates))
    to_write.append("Multiple gene alignments: %s" % multiple_genes)
    to_write.append("Partial alignments: %s" % partial_alignments)
    to_write.append("Missing gene alignments: %s" % missing_gene_alignments)
    outfile1.write('\n'.join(to_write))

# Write removed alignments to file:
SeqIO.write(removed_seq_records, args.removed_alignments_fasta, "fasta")
