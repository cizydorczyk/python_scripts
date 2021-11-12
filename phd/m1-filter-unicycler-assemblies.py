import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_assembly", help="input assembly", required=True)
parser.add_argument("--output_assembly", help="filtered assembly", required=True)
parser.add_argument("--output_removed_contigs", help="fasta with removed contigs", required=True)
parser.add_argument("--minlen", help="minimum contig length to keep", required=True, type=int)
parser.add_argument("--mindepth", help="minimum contig depth to keep", required=True, type=float)
parser.add_argument("--isolate", help="isolate identifier", required=True)

args = parser.parse_args()

records_to_keep = []
removed_contigs = []
original_num_records = 0
for record in SeqIO.parse(args.input_assembly, "fasta"):
    original_num_records += 1

    record_length = int(record.description.strip().split(" ")[1].split("=")[1])
    record_depth = float(record.description.strip().split(" ")[2].split("=")[1][:-1])

    if record_length > args.minlen and record_depth > args.mindepth:
        records_to_keep.append(record)
    else:
        removed_contigs.append(record)

print("\tIsolate\tOriginal Num Contigs\tFiltered Num Contigs\tRemoved Contigs\n")
print(f"\t{args.isolate}\t{original_num_records}\t{len(records_to_keep)}\t{','.join([record.description for record in removed_contigs])}")

with open(args.output_assembly, 'w') as outfile1:
    SeqIO.write(records_to_keep, outfile1, "fasta")

with open(args.output_removed_contigs, 'w') as outfile2:
    SeqIO.write(removed_contigs, outfile2, "fasta")

### Sample bash call
# (base) conrad@conrad-Precision-Tower-3620:~/m-sm-notebook/de-novo-assemblies/k1-unicycler-assemblies$ for i in $(cat ~/m-sm-notebook/sm_isolate_list.txt); do python ~/python_scripts/phd/m1-filter-unicycler-assemblies.py --input_assembly ~/m-sm-notebook/de-novo-assemblies/k1-unicycler-assemblies/$i/assembly.fasta --output_assembly ~/m-sm-notebook/de-novo-assemblies/k1-unicycler-assemblies/$i/$i"-filtered-assembly.fasta" --output_removed_contigs ~/m-sm-notebook/de-novo-assemblies/k1-unicycler-assemblies/$i/$i"-removed-contigs.fasta" --minlen 