import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--rgi_txt", help="card rgi .txt output file")
parser.add_argument("--genome", help="genome identifier/accession")
parser.add_argument("--summary_file", help="output summary file")
parser.add_argument("--single_copy_dir", help="output dir for single copy fastas")
parser.add_argument("--multi_copy_dir", help="output dir for multi copy fastas")
parser.add_argument("--single_copy_aberrant_len_dir")

args = parser.parse_args()

# initialize summary file list:
summary_file_list = []
summary_file_list.append(args.genome) # add genome ID as first element

# read rgi output file:
with open(args.rgi_txt, "r") as infile1:

    blaz_entries = []
    for line in infile1:
        if "blaZ" in line:
            blaz_entries.append(line)
        else:
            continue

    """
    if # entries == 0:
        append summary file that no entries present
    elif # entries == 1:
        append summary file that 1 entry present
        get entry elements
        write aa and dna fasta in single-copy dir
    elif # entries > 1:
        append summary file that >1 entry present
        write aa and dna fastas in multi-copy dir
    """
    if len(blaz_entries) == 0:
        summary_file_list.append("0") # append that 0 entries present (gene absent)

    elif len(blaz_entries) == 1:
        summary_file_list.append("1") # append that 1 entry present

        entry_elements = blaz_entries[0].strip().split("\t")

        # if length is proper:
        if len(entry_elements[17].strip()) == 846:
            summary_file_list.append(str(len(entry_elements[17].strip())))
            summary_file_list.append(str(len(entry_elements[18].strip())))

            gene = "_".join(entry_elements[16].split(" "))

            dna_seq = entry_elements[17] # get dna seq
            dna_seq_multiline = []
            while len(dna_seq) > 0: # convert to multiline dna sequence
                dna_seq_multiline.append(dna_seq[:70])
                dna_seq = dna_seq[70:]
            dna_seq_multiline = "\n".join(dna_seq_multiline)

            aa_seq = entry_elements[18] # get aa seq
            aa_seq_multiline = []
            while len(aa_seq) > 0: # convert to multiline aa sequence
                aa_seq_multiline.append(aa_seq[:70])
                aa_seq = aa_seq[70:]
            aa_seq_multiline = "\n".join(aa_seq_multiline)

            out_dna_handle = args.single_copy_dir + args.genome + "-blaZ.fna"
            out_aa_handle = args.single_copy_dir + args.genome + "-blaZ.faa"

            to_write_dna = ">" + args.genome + "_" + gene + "\n" + dna_seq_multiline
            to_write_aa = ">" + args.genome + "_" + gene + "\n" + aa_seq_multiline

            with open(out_dna_handle, "w") as dna_fasta:
                dna_fasta.write(to_write_dna)
            with open(out_aa_handle, "w") as aa_fasta:
                aa_fasta.write(to_write_aa)

            ######### blaZ typing of AA seq below #############
            aa_1 = entry_elements[18][118:119]
            aa_2 = entry_elements[18][206:207]

            blaz_type = ""
            if aa_1 == "T" and aa_2 == "S":
                blaz_type = "A"
            elif aa_1 == "K" and aa_2 == "N":
                blaz_type = "B"
            elif aa_1 == "T" and aa_2 == "N":
                blaz_type = "C"
            elif aa_1 == "A" and aa_2 == "S":
                blaz_type = "D"
            else:
                blaz_type = "Unknown"

            summary_file_list.append("") # append blank for note field
            summary_file_list.append(aa_1)
            summary_file_list.append(aa_2)
            summary_file_list.append(blaz_type)

        # if length is improper:
        elif len(entry_elements[17].strip()) != 846:
            summary_file_list.append(str(len(entry_elements[17].strip())))
            summary_file_list.append(str(len(entry_elements[18].strip())))
            summary_file_list.append("aberrant length")

            gene = "_".join(entry_elements[16].split(" "))

            dna_seq = entry_elements[17] # get dna seq
            dna_seq_multiline = []
            while len(dna_seq) > 0: # convert to multiline dna sequence
                dna_seq_multiline.append(dna_seq[:70])
                dna_seq = dna_seq[70:]
            dna_seq_multiline = "\n".join(dna_seq_multiline)

            aa_seq = entry_elements[18] # get aa seq
            aa_seq_multiline = []
            while len(aa_seq) > 0: # convert to multiline aa sequence
                aa_seq_multiline.append(aa_seq[:70])
                aa_seq = aa_seq[70:]
            aa_seq_multiline = "\n".join(aa_seq_multiline)

            out_dna_handle = args.single_copy_aberrant_len_dir + args.genome + "-blaZ.fna"
            out_aa_handle = args.single_copy_aberrant_len_dir + args.genome + "-blaZ.faa"

            to_write_dna = ">" + args.genome + "_" + gene + "\n" + dna_seq_multiline
            to_write_aa = ">" + args.genome + "_" + gene + "\n" + aa_seq_multiline

            with open(out_dna_handle, "w") as dna_fasta:
                dna_fasta.write(to_write_dna)
            with open(out_aa_handle, "w") as aa_fasta:
                aa_fasta.write(to_write_aa)

    # if >1 entry, don't worry about lengths...
    elif len(blaz_entries) > 1:
        summary_file_list.append(str(len(blaz_entries)))

        count = 1 # count for numbering output files; increases with each subsequent entry
        tmp_dna_len_list = []
        tmp_aa_len_list = []
        for blaz_entry in blaz_entries:
            entry_elements = blaz_entry.strip().split("\t")

            gene = "_".join(entry_elements[16].split(" "))

            dna_seq = entry_elements[17] # get dna seq
            tmp_dna_len_list.append(str(len(dna_seq.strip())))

            ####### Biostars page for creating multiline fasta: https://www.biostars.org/p/205745/
            dna_seq_multiline = []
            while len(dna_seq) > 0: # convert to multiline dna sequence
                dna_seq_multiline.append(dna_seq[:70])
                dna_seq = dna_seq[70:]
            dna_seq_multiline = "\n".join(dna_seq_multiline)

            aa_seq = entry_elements[18] # get aa seq
            tmp_aa_len_list.append(str(len(aa_seq.strip())))

            aa_seq_multiline = []
            while len(aa_seq) > 0: # convert to multiline aa sequence
                aa_seq_multiline.append(aa_seq[:70])
                aa_seq = aa_seq[70:]
            aa_seq_multiline = "\n".join(aa_seq_multiline)

            out_dna_handle = args.multi_copy_dir + args.genome + "-blaZ" + "_" + str(count) + ".fna"
            out_aa_handle = args.multi_copy_dir + args.genome + "-blaZ" + "_" + str(count) + ".faa"

            to_write_dna = ">" + args.genome + "_" + gene + "\n" + dna_seq_multiline
            to_write_aa = ">" + args.genome + "_" + gene + "\n" + aa_seq_multiline

            with open(out_dna_handle, "w") as dna_fasta:
                dna_fasta.write(to_write_dna)
            with open(out_aa_handle, "w") as aa_fasta:
                aa_fasta.write(to_write_aa)

            count += 1

        summary_file_list.append("_".join(tmp_dna_len_list))
        summary_file_list.append("_".join(tmp_aa_len_list))

        summary_file_list.append("multicopy")

# Write summary file:
if not os.path.isfile(args.summary_file):
    with open(args.summary_file, "w") as summary_out:
        header = "Genome\tNum_blaZ_Genes\tDNA_Len\tAA_Len\tNotes\tAA119\tAA207\tBlaZ_Type"
        summary_out.write(header + "\n" + "\t".join(summary_file_list))
else:
    with open(args.summary_file, "a") as summary_out:
        summary_out.write("\n" + "\t".join(summary_file_list))
