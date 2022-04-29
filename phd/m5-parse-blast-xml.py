import argparse
# import pandas as pd
from Bio.Blast import NCBIXML
from Bio import SeqIO

parser = argparse.ArgumentParser()

# General script options:
parser.add_argument("--blast_file")
parser.add_argument("--genus", default="any")
parser.add_argument("--blast_contigs")
parser.add_argument("--reference")
parser.add_argument("--output_contigs")

## Below options for testing for genus using taxids; not currently working. Code kept for reference.
# parser.add_argument("--taxidlist")
# parser.add_argument("--scorestokeep", help="number of unique, top scoring subject sequences to keep (only the top hsp from each sequence is retained).", type=int)

args = parser.parse_args()

blast_handle = open(args.blast_file) #
raw_blast_output = NCBIXML.parse(blast_handle)

records_to_keep = []
for record in raw_blast_output:
    record_alignment_titles = [alignment.title for alignment in record.alignments]
    if args.genus == "any":
        print("Keeping all records as no genus specified...")
        records_to_keep.append(record.query)
    else:
        if any(args.genus in title for title in record_alignment_titles):
            print(f"Contig {record.query} has at least one hit in the specified genus; keeping contig...")
            records_to_keep.append(record.query)
        else:
            print(f"Contig {record.query} has NO hits in the specified genus; removing contig...")

# Read in contigs & reference fasta files:
blast_contigs_records = list(SeqIO.parse(args.blast_contigs, "fasta"))
reference_records = list(SeqIO.parse(args.reference, "fasta"))

# Filter blast contigs to keep only those in records_to_keep:
blast_contigs_to_keep = [record for record in blast_contigs_records if record.id in records_to_keep]
print(f"Keeping contigs {' '.join([record.id for record in blast_contigs_to_keep])}")

# Combine reference + blast contigs to keep:
for record in blast_contigs_to_keep: # add _unmapped to contig names, since reference already contains contigs named 1, 2, 3, etc...
    record.id = record.id + "_unmapped"

combined_output_records = reference_records + blast_contigs_to_keep

with open(args.output_contigs, "w") as outfile1:
    SeqIO.write(combined_output_records, outfile1, "fasta")

blast_handle.close()

###################################
### Code below for using outfmt 6 (see column names in sseqid_df below for columns to output) & testing for genus membership by taxid.
### Not currently working, since sorting of best hsp per alignment relies on a single best hsp, but there may be more than 1. But the script does not currently handle such a case...
###################################


# # Create set of taxids:
# taxidlist_handle = open(args.taxidlist)
# taxid_set = {line.strip() for line in taxidlist_handle}
# taxidlist_handle.close()

# # Read & parse raw blast output (fmt 6):
# raw_blast_output = pd.read_csv(filepath_or_buffer=args.blast_file, sep="\t", header=None)
# raw_blast_output.columns = ["qseqid", "sseqid", "sacc", "bitscore", "evalue", "pident", "length", "staxids", "sscinames"]

# # For each unique subject sequence, get max bitscore & associated taxid:
# unique_sseqids = raw_blast_output.iloc[:,1].unique()

# sseqid_lists = []

# for unique in unique_sseqids:
#     unique_subset = raw_blast_output[raw_blast_output['sseqid'] == unique]
#     taxid = list(unique_subset['staxids'])[0]
#     evalue = list(unique_subset['evalue'])[0]
#     max_bitscore = max(unique_subset['bitscore'])
#     qseqid = list(unique_subset['qseqid'])[0]

#     to_append = [qseqid, unique, max_bitscore, taxid, evalue]
#     sseqid_lists.append(to_append)

# sseqid_df = pd.DataFrame.from_records(sseqid_lists, columns=['qseqid', 'sseqid', 'bitscore', 'taxid', 'evalue']).sort_values(by=['qseqid', 'bitscore'], ascending=False) # main output of this section
# sseqid_df = sseqid_df.reset_index(drop=True)
# print(sseqid_df)

# # Subset hits based on bitscores; keep the top {args.scorestokeep} subject sequences:
# # sseqid_df_subset = sseqid_df.iloc[0:args.scorestokeep,:]

# # Test if any of taxids in sseqid_df_subset are in genus taxid list:
