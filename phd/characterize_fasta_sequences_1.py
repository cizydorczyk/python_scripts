from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--aln", help="alignment to remove reference from")
parser.add_argument("--out", help="output fasta file with no reference sequence")

args = parser.parse_args()

records = list(SeqIO.parse(args.aln, "fasta"))

char_records = []
for record in records:
    seq = str(record.seq)
    length = len(seq)
    chars = set(record.seq)

    chars_dict = {}
    for char in chars:
        chars_dict[char] = 0

    chars_list_entry = ''
    for key in sorted(chars_dict):
        key_count = seq.count(key)
        char_frac = float(round(key_count/length,4))
        chars_dict[key] = float(key_count/length)
        chars_list_entry = chars_list_entry + '\t' + key + '\t' + str(char_frac)

    chars_list_entry = record.id + '\t' + chars_list_entry
    char_records.append(chars_list_entry)

print('\n'.join(char_records))

with open(args.out, 'w') as outfile:
    outfile.write('\n'.join(char_records))
