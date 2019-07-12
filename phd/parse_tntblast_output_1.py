import argparse
import os.path


parser = argparse.ArgumentParser()

parser.add_argument("--tntblast_file", help="raw tntblast output file")
parser.add_argument("--isolate", help="isolate number/id")
parser.add_argument("--output_summary_file", help="output_file")

args = parser.parse_args()

primers = {"A_CladeAspe4-YFR5":"A", "B_prfC-1615spe0-YF1/d2034-YR1":"B", "C_C-SNP1-700spe-YF1-762spe-YR2":"C", "C1_C1-578spe-YF1/898-YR1":"C1", "C1-M27_M27aer-spe-YF1/YR2":"C1M27", "C2_nrdI-534spe2-YF1/678R":"C2", "ST131_R19-YF1/YR1":"ST131"}

### Something?

iso_primers = []
primer_clades = []
amplicon_lengths = []
contigs = []
sequences = []

with open(args.tntblast_file, 'r') as infile1:
    for line in infile1:
        if line.startswith("name ="):
            primer = line.strip().split(' ')[-1]
            iso_primers.append(primer)
            primer_clades.append(primers[primer])
        if line.startswith("amplicon length ="):
            amplicon_length = line.strip().split(' ')[-1]
            amplicon_lengths.append(amplicon_length)
        if line.startswith(">"):
            contig = line.strip()
            contigs.append(contig)
            sequence = next(infile1).strip()
            sequences.append(sequence)

if "A" in primer_clades and "ST131" in primer_clades:
    clade = "A"
elif "B" in primer_clades and "ST131" in primer_clades:
    clade = "B"
elif "C" in primer_clades and "C1" in primer_clades and "ST131" in primer_clades and "C1M27" not in primer_clades:
    clade = "C1"
elif "C" in primer_clades and "C1" in primer_clades and "C1M27" in primer_clades and "ST131" in primer_clades:
    clade = "C1M27"
elif "C" in primer_clades and "C2" in primer_clades and "ST131" in primer_clades:
    clade = "C2"
elif "C" in primer_clades and "ST131" in primer_clades:
    clade = "C0"
elif "ST131" in primer_clades and "A" not in primer_clades and "B" not in primer_clades and "C" not in primer_clades and "C1" not in primer_clades and "C1M27" not in primer_clades and "C2" not in primer_clades:
    clade = "other-st131"
else:
    clade = "non-st131_or_mixed"

output_line = '\t'.join([args.isolate, clade, ','.join(primer_clades), ','.join(iso_primers), ';'.join(amplicon_lengths), ','.join(contigs), ','.join(sequences), args.tntblast_file])

if os.path.isfile(args.output_summary_file):
    with open(args.output_summary_file, 'a') as outfile1:
        outfile1.write("\n" + output_line)

else:
    header = "isolate" + "\t" + "clade" + "\t" + "primer_clades" + "\t" + "primers" + "\t" + "amplicon_lengths" + "\t" + "contig_names" + "\t" + "sequences"
    to_write = header + "\n" + output_line
    with open(args.output_summary_file, 'w') as outfile1:
        outfile1.write(to_write)
