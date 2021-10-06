import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_snps_tab", help="snps tab file from snippy")

args = parser.parse_args()

hypermutation_terms = ["mut", "polymerase", "helicase", "sod", "superoxide", "oxyr", "ung", "uvrd", "dam", "adenine", "epsilon", "pola", "dnaq"]

records_to_keep = []
with open(args.input_snps_tab, "r") as infile1:
    for line in infile1:
        if not line.startswith("CHROM"):
            if any(i in line.lower() for i in hypermutation_terms):
                records_to_keep.append(line.strip())

for record in records_to_keep:
    print(record)
