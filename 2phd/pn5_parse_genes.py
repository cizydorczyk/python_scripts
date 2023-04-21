#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-03-27
Purpose: parse mutated genes to get PLES number
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="parse mutated genes to get PLES number",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_genes",
        help="genes file from bedtools intersect (see pn5)",
        metavar="str",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--output_genes",
        help="output genes file",
        metavar="str",
        type=str,
        required=True,
    )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Get PLES numbers; get genes without PLES"""

    args = get_args()

    output_lines = []
    with open(args.input_genes, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            if line_elements[-1] != "0":
                if line_elements[2] == "gene": # this skips entries without a corresponding gene entry, such as regulatory regions!
                    gene_entry = next(infile1).strip().split("\t")
                    gene_type = gene_entry[2]
                    gene_start = gene_entry[3]
                    gene_end = gene_entry[4]
                    gene_info = gene_entry[8]

                    gene_info_list = gene_info.split(";")
                    A180_annot = gene_info_list[0].split("=")[-1]

                    note_entry = [note for note in gene_info_list if "Note=" in note]

                    if len(note_entry) != 0:
                        note_entry_list = note_entry[0].split(",")
                        locus_tag = [tag for tag in note_entry_list if "PLES" in tag]

                        if len(locus_tag) == 1:
                            output_line = f'{A180_annot}\t{gene_type}\t{gene_start}\t{gene_end}\t{locus_tag[0].split(":")[1]}'
                            output_lines.append(output_line)
                        elif len(locus_tag) == 0:
                            output_line = f"{A180_annot}\t{gene_type}\t{gene_start}\t{gene_end}\tNO_PLES_TAG"
                            output_lines.append(output_line)
                    elif len(note_entry) == 0:
                        print("No Note entry in VCF record...")
                    else:
                        print("Something else is wrong...")

    to_write = "\n".join(output_lines)
    with open(args.output_genes, "w") as outfile1:
        outfile1.write(to_write)


# --------------------------------------------------
if __name__ == "__main__":
    main()
