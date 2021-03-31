import argparse
from Bio import SeqIO

################### NOTE ###################
# The script currently does NOT subsample reads...but it can! Just need to
# uncomment lines 55-60 (inclusive) and 64-72 (inclusive). In its current
# state, it just prints estimated genome sizes, depths, and subsampling factors
# to the terminal.

parser = argparse.ArgumentParser()

parser.add_argument("--phaster_regions", help="phaster phage regions file")
parser.add_argument("--phaster_summary", help="phaster summary file")
parser.add_argument("--isolate", help="isolate name to include in new fasta\
                        headers")
parser.add_argument("--output_fasta", help="fasta with only intact/questionable\
                        phage regions kept")
parser.add_argument("--keep", help="type of phage predictions to keep; options\
                    include 'intact', 'questionable', and 'incomplete'. Picking\
                    a less strict option (e.g. questionable) will include \
                    the stricter phages as well (i.e. intact). Default = \
                    questionable.", default="questionable")

args = parser.parse_args()

# Read phage regions
phage_nodes = []
phage_regions = list(SeqIO.parse(args.phaster_regions, "fasta"))
for i in phage_regions:
    phage_nodes.append(i.description.split("\t")[-1].strip())

# print(len(phage_nodes))
# print(len(phage_regions))

# Read summary file
with open(args.phaster_summary, "r") as infile1:
    summary_records = infile1.readlines()[34:]

# Remove incomplate records
intact_questionable_records = []
for line in summary_records:
    if args.keep == "intact":
        if "intact" in line:
            intact_questionable_records.append(line)
    elif args.keep == "questionable":
        if "intact" in line or "questionable" in line:
            intact_questionable_records.append(line)
    elif args.keep == "incomplete":
        if "intact" in line or "questionable" in line or "incomplete" in line:
            intact_questionable_records.append(line)


# Select only fasta records that match intact/questionable phages
nodes_to_keep = []
for node in phage_nodes:
    if any(node in record for record in intact_questionable_records):
        nodes_to_keep.append(node)
        print(node, "\t--\tKEPT")
    else:
        print(node, "\t--\t\tREMOVED")


# Select phage regions to keep and edit header to include isoalte
regions_to_write = []
for i in phage_regions:
    node_id = i.description.split("\t")[0].strip()
    node = i.description.split("\t")[-1].strip()

    if node in nodes_to_keep:
        i.id = node_id + "|" + args.isolate + "|" + node
        regions_to_write.append(i)

# Write intact/complete regions fasta file
SeqIO.write(regions_to_write, args.output_fasta, "fasta")
