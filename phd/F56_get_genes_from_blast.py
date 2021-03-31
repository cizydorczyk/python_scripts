import argparse
from Bio import SeqIO
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--blast_results", help="blast results file; outfmt 7")
parser.add_argument("--gff_file", help="isolate gff file")
parser.add_argument("--extend_len", help="number of bases to extend ORF regions by.\
                    Default = 2000. Requires --extend to be set.", default=2000,
                    type=int)
parser.add_argument("--extend", help="extend ORF regions past their ends? Set \
                    flag to extend.", action="store_true")
parser.add_argument("--reference", help="reference genome from which to extract\
                    region of interest.")
parser.add_argument("--output_fasta", help="output fasta to write seq regions to")

args = parser.parse_args()

hits = []
with open(args.blast_results, "r") as infile1:
    for line in infile1:
        if not line.startswith("#"):
            hit = line.strip().split("\t")[0].split("|")[1]
            if hit not in hits:
                hits.append(hit)

# parse gff
chromosome = ""
gff_entries =  {}
with open(args.gff_file, "r") as infile2:
    for line in infile2:
        if any(hit in line for hit in hits):
            line_elements = line.strip().split("\t")
            lhit = line_elements[8].split(";")[0].split("|")[1]
            if any(hit == lhit for hit in hits):

                lstart = line_elements[3]
                lend = line_elements[4]
                gff_entries[lhit] = [int(lstart), int(lend)]

                # set csome
                if chromosome == "":
                    chromosome = line_elements[0].split("|")[-1]



# if == 2 entries, check if consecutive ORFs
# if consecutive, merge coords
# if not consecutive, keep separate ORFs

# DOES NOT HANDLE >2 CONSECUTIVE ORFs
# DELETES 1st PAIR OF CONSECUTIVE ENTRIES AND THUS WILL MISS 3RD ENTRY THAT
# IS CONSECUTIVE TO THE PREVIOUS 2 IF PRESENT

if len(gff_entries) > 1:
    entry_list = [keys for keys in gff_entries]
    orf_nums = [int(key.split(".")[-1]) for key in gff_entries]

    for entry in entry_list:
        prev_entry = int(entry.split(".")[-1])-1
        next_entry = int(entry.split(".")[-1])+1

        if prev_entry in orf_nums:
            prev_full_entry = ".".join(entry.split(".")[0:-1] + [str(prev_entry)])
            prev_entry_coords = gff_entries[prev_full_entry]
            current_entry_coords = gff_entries[entry]

            new_coords = [min(prev_entry_coords + current_entry_coords), max(prev_entry_coords + current_entry_coords)]
            new_entry = prev_full_entry + "_" + entry

            gff_entries[new_entry] = new_coords

            # remove entries already tested
            gff_entries.pop(entry)
            gff_entries.pop(prev_full_entry)

            # remove from entry_list, as this is what we are iterating over
            entry_list.remove(entry)
            entry_list.remove(prev_full_entry)

            # print("prev entry triggered")

        elif next_entry in orf_nums:
            next_full_entry = ".".join(entry.split(".")[0:-1] + [str(next_entry)])
            next_entry_coords = gff_entries[next_full_entry]
            current_entry_coords = gff_entries[entry]

            new_coords = [min(next_entry_coords + current_entry_coords), max(next_entry_coords + current_entry_coords)]
            new_entry = entry + "_" + next_full_entry

            gff_entries[new_entry] = new_coords

            # remove entries already tested
            gff_entries.pop(entry)
            gff_entries.pop(next_full_entry)

            # remove from entry_list, as this is what we are iterating over
            entry_list.remove(entry)
            entry_list.remove(next_full_entry)

            # print("next entry triggered")

if args.extend:
    for key in gff_entries:
        gff_entries[key][0] = gff_entries[key][0] - args.extend_len - 1 # -1 b/c bed format is 0-based
        gff_entries[key][1] = gff_entries[key][1] + args.extend_len - 1 # -1 b/c bed format is 0-based

# read genome (reference) fasta file
ref_seq = list(SeqIO.parse(args.reference, "fasta"))[0]

# generate fasta to write
target_seqs = []
for entry in gff_entries:
    target_seq = ">" + entry + "\n" + str(ref_seq.seq[gff_entries[entry][0]:gff_entries[entry][1]]).upper()
    target_seqs.append(target_seq)

# Write output to file
if os.path.isfile(args.output_fasta):
    with open(args.output_fasta, "a") as outfile:
        outfile.write("\n" + "\n".join(target_seqs))
else:
    with open(args.output_fasta, "w") as outfile:
        outfile.write("\n".join(target_seqs))















#####
