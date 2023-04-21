#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk
Date   : 2023-04-11
Purpose: get new coords for plotting from jasmine VCF for a single isolate's vcf file
"""

import argparse
import statistics
import portion
import pandas as pd


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Rock the Casbah",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-o",
        "--output_coords",
        help="tsv with new coords for each indel in the isolates vcf",
        metavar="str",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-od",
        "--output_dels",
        help="tsv with new coords for each del in the isolates vcf",
        metavar="str",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-oi",
        "--output_ins",
        help="tsv with new coords for each ins in the isolates vcf",
        metavar="str",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-f",
        "--fofn",
        help="file of file names, one vcf for each isolate. tab sep, should be isolate<tab>/path/to/file.vcf",
        metavar="str",
        type=str,
        required=True,
    )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Get new coords"""

    args = get_args()

    insertions = {}
    insertions_start_positions = []
    insertion_lengths = []
    insertion_supports = []

    deletions = {}
    deletion_lengths = []
    deletion_supports = []

    inversions = {}
    inversion_lengths = []
    inversion_supports = []

    dup_tandem = {}
    dup_tandem_lengths = []
    dup_tandem_supports = []

    other = []

    # Parse individual isolate VCFs to get intervals:
    print("Parsing individual isolate files...")
    isolate_vcfs = {}
    with open(args.fofn, "r") as infile2:
        for line in infile2:
            isolate = line.strip().split("\t")[0]
            vcf_file = line.strip().split("\t")[1]
            isolate_vcfs[isolate] = vcf_file

    # Parse isolate VCF files to get master list of insertions:
    for isolate in isolate_vcfs:
        with open(isolate_vcfs[isolate], "r") as infile1:
            for line in infile1:
                if not line.startswith("#"):
                    line_elements = line.strip().split("\t")
                    if line_elements[2].split(".")[2] == "INS":
                        ins_start = int(line_elements[1])
                        ins_length = int(line_elements[7].split(";")[2].split("=")[-1])
                        ins = (ins_start, ins_length)

                        if not ins_start in insertions:
                            insertions[ins_start] = ins_length
                        elif ins_start in insertions:
                            if ins_length > insertions[ins_start]:
                                insertions[ins_start] = ins_length

    
    # Generate corrected reference positions dict:
    ref_positions = get_corrected_ref_positions(ref_length=6601757, insertions_dict=insertions)
    
    # Create corrected insertions dict:
    insertions_corrected = {}
    for insertion in insertions:
        corrected_ins_start = ref_positions[insertion]
        insertions_corrected[corrected_ins_start] = insertions[insertion]
    
    # Create output data list to hold all data for all isolates:
    output_data = []
    del_output = []
    ins_output = []

    # Create interval for pseudo reference (reference + insertions):
    full_interval = portion.open(1, sum([6601757, sum([insertions_corrected[i] for i in insertions_corrected])]))


    # Add reference to output data:
    ref_gaps = []
    for insertion in insertions_corrected:
        insertion_end = insertion + insertions_corrected[insertion]
        ref_gap_interval = portion.closed(insertion, insertion_end)
        ref_gaps.append(ref_gap_interval)
    ref_gaps = portion.Interval(*ref_gaps)
    ref_intervals = list(
        full_interval - ref_gaps
    )  # get difference between full interval & ref, represented by insertions
    for interval in ref_intervals:
        output_data.append(
            ["LESB58", interval.lower, interval.upper]
        )  # append intervals to final output

    # Parse individual isolate VCFs:
    isolate_vcfs = {}
    with open(args.fofn, "r") as infile2:
        for line in infile2:
            isolate = line.strip().split("\t")[0]
            vcf_file = line.strip().split("\t")[1]
            isolate_vcfs[isolate] = vcf_file
    for isolate in isolate_vcfs:
        isolate_deletions = []
        isolate_insertions = {}
        isolate_uncorrected_deletions = []
        with open(isolate_vcfs[isolate], "r") as infile3:
            
            for line in infile3:
                
                if not line.startswith("#"):
                    line_elements = line.strip().split("\t")
                    
                    if line_elements[2].split(".")[2] == "DEL": # if deletion
                        new_start_pos = ref_positions[int(line_elements[1])]
                        new_end_pos = new_start_pos + abs(int(line_elements[7].split(";")[2].split("=")[-1]))
                        isolate_deletions.append((new_start_pos, new_end_pos))

                        isolate_uncorrected_deletions.append((int(line_elements[1]), int(line_elements[1]) + abs(int(line_elements[7].split(";")[2].split("=")[-1]))))
                    
                    elif line_elements[2].split(".")[2] == "INS": # if insertion
                        ins_start_pos = ref_positions[int(line_elements[1])]
                        ins_length = int(line_elements[7].split(";")[2].split("=")[-1])
                        isolate_insertions[ins_start_pos] = ins_length
                        
                        if ins_start_pos not in insertions_corrected:
                            print("ERROR ERROR ERROR")
                            print(isolate, ins_start_pos)

        missing_insertions = []
        for insertion_start in sorted(insertions_corrected):
            if insertion_start not in isolate_insertions:
                missing_insertions.append((insertion_start, insertion_start + insertions_corrected[insertion_start]))
            
            elif insertion_start in isolate_insertions:
                
                if insertions_corrected[insertion_start] == isolate_insertions[insertion_start]:
                    continue
                
                elif insertions_corrected[insertion_start] > isolate_insertions[insertion_start]:
                    length_dif = insertions_corrected[insertion_start] - isolate_insertions[insertion_start]
                    missing_portion_end = insertion_start + insertions_corrected[insertion_start]
                    missing_portion_start = missing_portion_end - length_dif
                    missing_insertions.append((missing_portion_start, missing_portion_end))
        
        missing_regions = isolate_deletions + missing_insertions
        missing_regions = isolate_deletions + missing_regions
        missing_regions = [portion.closed(i[0], i[1]) for i in missing_regions]
        missing_regions = portion.Interval(*missing_regions)
        output_intervals = list(full_interval - missing_regions)

        for interval in output_intervals:
            output_data.append([isolate, interval.lower, interval.upper])

        # Deletions only:
        del_only = [portion.closed(i[0], i[1]) for i in isolate_uncorrected_deletions]
        del_only = portion.Interval(*del_only)
        ref_interval = portion.closed(1, 6601757)
        del_intervals = list(ref_interval - del_only)

        for interval in del_intervals:
            del_output.append([isolate, interval.lower, interval.upper])
        
        # Insertions only:
        for ins in isolate_insertions:
            ins_output.append([isolate, ins, ins+isolate_insertions[ins]])
        
    # Convert output to Pandas df:
    output_data_df = pd.DataFrame(output_data, columns=["Isolate", "Start", "End"])
    del_df = pd.DataFrame(del_output, columns=["Isolate", "Start", "End"])
    ins_df = pd.DataFrame(ins_output, columns=["Isolate", "Start", "End"])

    # Write output df to file:
    print("Writing final data to file...")
    output_data_df.to_csv(
        path_or_buf=args.output_coords, sep="\t", header=True, index=False
    )

    del_df.to_csv(path_or_buf=args.output_dels, sep="\t", header=True, index=False)
    ins_df.to_csv(path_or_buf=args.output_ins, sep="\t", header=True, index=False)
    
    print("Done.")

    print(ref_positions[665561], ref_positions[680385])
    print(ref_positions[863875], ref_positions[906018])
    print(ref_positions[1433825], ref_positions[1476616])
    print(ref_positions[1684114], ref_positions[1720919])
    print(ref_positions[2691759], ref_positions[2741659])
    print(ref_positions[4546499], ref_positions[4554097])
    print(ref_positions[2504700], ref_positions[2552409])
    print(ref_positions[2753109], ref_positions[2784809])
    print(ref_positions[2798145], ref_positions[2908715])
    print(ref_positions[3394109], ref_positions[3433537])
    print(ref_positions[4932837], ref_positions[4962250])
    print(ref_positions[280017], ref_positions[292442])
    print(ref_positions[2054752], ref_positions[2071349])
    print(ref_positions[2611035], ref_positions[2622606])
    print(ref_positions[2640878], ref_positions[2649200])
    print(ref_positions[3148400], ref_positions[3194166])
    print(ref_positions[3744468], ref_positions[3798235])
    print(ref_positions[3947016], ref_positions[3955409])
    print(ref_positions[4770083], ref_positions[4791709])
    print(ref_positions[4967101], ref_positions[4974424])
    print(ref_positions[5574296], ref_positions[5585225])
    print(ref_positions[5658978], ref_positions[5686227])
    print(ref_positions[6129698], ref_positions[6166692])


def get_corrected_ref_positions(ref_length, insertions_dict):
    ref_positions = {}
    
    for i in range(1,ref_length+1):
        ref_positions[i] = i

    for position in ref_positions:
        adjusted_position = position
        for insertion in sorted(insertions_dict):
            if not position > insertion:
                break
            else:
                adjusted_position += insertions_dict[insertion]
        ref_positions[position] = adjusted_position
    
    return(ref_positions)


# --------------------------------------------------
if __name__ == "__main__":
    main()
