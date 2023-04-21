#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-04-20
Purpose: get vcf stats; parse vcf variant types
"""

import argparse
import statistics
import portion


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='get vcf stats; parse vcf variant types',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f',
                        '--fofn',
                        help='File of file names: "isolate<tab>/path/to/file.vcf", one per line.',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-d',
                        '--deletions_vcf',
                        help='vcf with deletions only',
                        metavar='str',
                        type=str,
                        default='')
    
    parser.add_argument('-i',
                        '--insertions_vcf',
                        help='vcf with insertions only',
                        metavar='str',
                        type=str,
                        default='')
    
    parser.add_argument('-v',
                        '--inversions_vcf',
                        help='vcf with inversions only',
                        metavar='str',
                        type=str,
                        default='')
    
    parser.add_argument('-t',
                        '--tandemdup_vcf',
                        help='vcf with tandem dup only',
                        metavar='str',
                        type=str,
                        default='')
    
    parser.add_argument('-db',
                        '--deletions_bed',
                        help='bed with deletions only',
                        metavar='str',
                        type=str,
                        default='')
    
    parser.add_argument('-ib',
                        '--insertions_bed',
                        help='bed with insertions only',
                        metavar='str',
                        type=str,
                        default='')
    
    parser.add_argument('-vb',
                        '--inversions_bed',
                        help='bed with inversions only',
                        metavar='str',
                        type=str,
                        default='')
    
    parser.add_argument('-tdb',
                        '--tandemdup_bed',
                        help='bed with tandem dup only',
                        metavar='str',
                        type=str,
                        default='')

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    vcf_variants = get_variants_from_vcf(args.fofn)

    # Get summary stats for each variant type:
    insertion_length_stats = get_stats(vcf_variants[1])
    insertion_support_stats = get_stats(vcf_variants[2])

    deletion_length_stats = get_stats(vcf_variants[4])
    deletion_support_stats = get_stats(vcf_variants[5])

    inversion_length_stats = get_stats(vcf_variants[7])
    inversion_support_stats = get_stats(vcf_variants[8])

    dup_tandem_length_stats = get_stats(vcf_variants[10])
    dup_tandem_support_stats = get_stats(vcf_variants[11])

    # Print summary stats to screen:
    print("\nStats order: num_variant, mean, min, max, Q1, Q2 (median), Q3")
    print(
        f"Insertion stats:\n\tLength stats:\n\t\t{insertion_length_stats}\n\tSupport stats:\n\t\t{insertion_support_stats}\n"
        f"Deletion stats:\n\tLength stats:\n\t\t{deletion_length_stats}\n\tSupport stats:\n\t\t{deletion_support_stats}\n"
        f"Inversion stats:\n\tLength stats:\n\t\t{inversion_length_stats}\n\tSupport stats:\n\t\t{inversion_support_stats}\n"
        f"Tandem duplication stats:\n\tLength stats:\n\t\t{dup_tandem_length_stats}\n\tSupport stats:\n\t\t{dup_tandem_support_stats}"
    )

    # Print other variants that were not covered by ins, del, inv, or tandem_dup
    # These should be manually investigated
    print("Other unidentified variants:\n", vcf_variants[12])

    # Header for writing vcfs:
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"

    # Write deletions VCF if requested:
    print("Writing deletions VCF file...")
    if args.deletions_vcf != '':
        to_write_del = [vcf_variants[3][entry] for entry in vcf_variants[3]]
        to_write_del = "\n".join(to_write_del)
        to_write_del = header + to_write_del
        with open(args.deletions_vcf, "w") as outfile1:
            outfile1.write(to_write_del)

    # Write deletions BED if requested:
    if args.deletions_bed != '':
        base_path = "/".join(args.deletions_bed.strip().split("/")[0:-1])
        
        # If add a '/' if base_path is a path; otherwise, exclude & file will be written in current directory
        if base_path != "":
            base_path = base_path + "/"
        base_handle = ".".join(args.deletions_bed.strip().split("/")[-1].split(".")[0:-1])
        
        # Create handles for merged and unmerged bed files
        merged_handle = f"{base_path}{base_handle}.MERGED.bed"
        unmerged_handle = f"{base_path}{base_handle}.UNMERGED.bed"
        
        # Write merged BED:
        print("Writing ***MERGED*** deletions bed file...")
        del_coords = []
        for i in vcf_variants[4]:
            del_start = int(i.strip().split("_")[0])
            del_end = int(i.strip().split("_")[0]) + int(i.strip().split("_")[1])
            del_coords.append(portion.closed(del_start, del_end))
        del_coords = portion.Interval(*del_coords)
        del_coords = [f"{i.lower}\t{i.upper}" for i in del_coords]
        del_coords = "\n".join(del_coords)
        with open(merged_handle, "w") as infile1:
            infile1.write(del_coords)
        
        # Write unmerged BED:
        print("Writing ***UNMERGED*** deletions bed file...")
        del_coords = []
        for i in vcf_variants[4]:
            del_start = int(i.strip().split("_")[0])
            del_end = int(i.strip().split("_")[0]) + int(i.strip().split("_")[1])
            del_coords.append((del_start, del_end))
        del_coords = sorted(del_coords)
        del_coords = [f"{i[0]}\t{i[1]}" for i in del_coords]
        del_coords = "\n".join(del_coords)
        with open(unmerged_handle, "w") as infile1:
            infile1.write(del_coords)

    # Write insertions VCF if requested:
    print("Writing insertions VCF file...")
    if args.insertions_vcf != '':
        to_write_ins = [vcf_variants[0][entry] for entry in vcf_variants[0]]
        to_write_ins = "\n".join(to_write_ins)
        to_write_ins = header + to_write_ins
        with open(args.insertions_vcf, "w") as outfile1:
            outfile1.write(to_write_ins)
    
    # Write insertions BED if requested:
    print("Writing ***UN-MERGED*** insertions bed file...")
    if args.insertions_bed != '':
        ins_coords = []
        for i in vcf_variants[1]:
            ins_start = int(i.strip().split("_")[0])
            ins_end = int(i.strip().split("_")[0]) + int(i.strip().split("_")[1])
            ins_coords.append((ins_start, ins_end))
        ins_coords = sorted(ins_coords)
        ins_coords = [f"{i[0]}\t{i[1]}" for i in ins_coords]
        ins_coords = "\n".join(ins_coords)
        with open(args.insertions_bed, "w") as infile1:
            infile1.write(ins_coords)

    # Write inversions VCF if requested:
    print("Writing inversions VCF file...")
    if args.inversions_vcf != '':
        to_write_inv = [vcf_variants[6][entry] for entry in vcf_variants[6]]
        to_write_inv = "\n".join(to_write_inv)
        to_write_inv = header + to_write_inv
        with open(args.inversions_vcf, "w") as outfile1:
            outfile1.write(to_write_inv)
    
    # Write inversions BED if requested:
    print("Writing ***UN-MERGED*** inversions bed file...")
    if args.inversions_bed != '':
        inv_coords = []
        for i in vcf_variants[7]:
            inv_start = int(i.strip().split("_")[0])
            inv_end = int(i.strip().split("_")[0]) + int(i.strip().split("_")[1])
            inv_coords.append((inv_start, inv_end))
        inv_coords = sorted(inv_coords)
        inv_coords = [f"{i[0]}\t{i[1]}" for i in inv_coords]
        inv_coords = "\n".join(inv_coords)
        with open(args.inversions_bed, "w") as infile1:
            infile1.write(inv_coords)
    
    # Write tandem dup VCF if requested:
    print("Writing tandem duplications VCF file...")
    if args.tandemdup_vcf != '':
        to_write_td = [vcf_variants[9][entry] for entry in vcf_variants[9]]
        to_write_td = "\n".join(to_write_td)
        to_write_td = header + to_write_td
        with open(args.tandemdup_vcf, "w") as outfile1:
            outfile1.write(to_write_td)
    
    # Write tandem dup BED if requested:
    print("Writing ***UN-MERGED*** tandem duplications bed file...")
    if args.tandemdup_bed != '':
        td_coords = []
        for i in vcf_variants[10]:
            td_start = int(i.strip().split("_")[0])
            td_end = int(i.strip().split("_")[0]) + int(i.strip().split("_")[1])
            td_coords.append((td_start, td_end))
        td_coords = sorted(td_coords)
        td_coords = [f"{i[0]}\t{i[1]}" for i in td_coords]
        td_coords = "\n".join(td_coords)
        with open(args.tandemdup_bed, "w") as infile1:
            infile1.write(td_coords)
    
    print("Done.")


def get_variants_from_vcf(vcf_fofn):

    # Create Dict of isolate VCFs:
    isolate_vcfs = {}
    with open(vcf_fofn, "r") as infile1:
        for line in infile1:
            isolate = line.strip().split("\t")[0]
            vcf_file = line.strip().split("\t")[1]
            isolate_vcfs[isolate] = vcf_file
    
    # Parse variant types:
    insertions = {}
    insertions_lengths = {}
    insertions_supports = {}

    deletions = {}
    deletions_lengths = {}
    deletions_supports = {}

    inversions = {}
    inversions_lengths = {}
    inversions_supports = {}

    dup_tandem = {}
    dup_tandem_lengths = {}
    dup_tandem_supports = {}

    other = []

    for isolate in isolate_vcfs:
        
        with open(isolate_vcfs[isolate], "r") as infile3:
            
            for line in infile3:
                
                if not line.startswith("#"):
                    variant_type = line.strip().split("\t")[2].split(".")[2]

                    if variant_type == "INS":
                        
                        variant_stats = get_variant_stats(line)

                        insertion_key = f"{variant_stats[0]}_{variant_stats[1]}"

                        # Add insertion vcf entry:
                        if insertion_key not in insertions:
                            insertions[insertion_key] = line.strip() + "\t" + isolate
                        elif insertion_key in insertions:
                            insertions[insertion_key] = insertions[insertion_key] + "," + isolate
                        
                        # Store support for variant:
                        if insertion_key not in insertions_supports:
                            insertions_supports[insertion_key] = 1
                        
                        elif insertion_key in insertions_supports:
                            insertions_supports[insertion_key] += 1
                        
                        # Store length for variant:
                        if insertion_key not in insertions_lengths:
                            insertions_lengths[insertion_key] = variant_stats[1]
                    
                    elif variant_type == "DEL":

                        variant_stats = get_variant_stats(line)

                        deletion_key = f"{variant_stats[0]}_{variant_stats[1]}"

                        # Add deletion vcf entry:
                        if deletion_key not in deletions:
                            deletions[deletion_key] = line.strip() + "\t" + isolate
                        
                        elif deletion_key in deletions:
                            deletions[deletion_key] = deletions[deletion_key] + "," + isolate
                        
                        # Store support for variant:
                        if deletion_key not in deletions_supports:
                            deletions_supports[deletion_key] = 1
                        
                        elif deletion_key in insertions_supports:
                            insertions_supports[deletion_key] += 1
                        
                        # Store length for variant:
                        if deletion_key not in deletions_lengths:
                            deletions_lengths[deletion_key] = variant_stats[1]
                    
                    elif variant_type == "INV":

                        inv_start = int(line.strip().split("\t")[1])
                        inv_length = int(line.strip().split("\t")[7].split(";")[1].split("=")[-1]) - inv_start

                        inversion_key = f"{inv_start}_{inv_length}"

                        # Add inversion vcf entry:
                        if inversion_key not in inversions:
                            inversions[inversion_key] = line.strip() + "\t" + isolate
                        
                        elif inversion_key in inversions:
                            inversions[inversion_key] = inversions[inversion_key] + "," + isolate

                        # Store support for variant:
                        if inversion_key not in inversions_supports:
                            inversions_supports[inversion_key] = 1
                        
                        elif inversion_key in inversions_supports:
                            inversions_supports[inversion_key] += 1

                        # Store length for variant:
                        if inversion_key not in inversions_lengths:
                            inversions_lengths[inversion_key] = inv_length
                    
                    elif variant_type == "DUP_TANDEM":

                        variant_stats = get_variant_stats(line)

                        dup_tandem_key = f"{variant_stats[0]}_{variant_stats[1]}"

                        # Add tandem dup vcf entry:
                        if dup_tandem_key not in dup_tandem:
                            dup_tandem[dup_tandem_key] = line.strip() + "\t" + isolate
                        
                        elif dup_tandem_key in dup_tandem:
                            dup_tandem[dup_tandem_key] = dup_tandem[dup_tandem_key] + "," + isolate
                        
                        # Store support for variant:
                        if dup_tandem_key not in dup_tandem_supports:
                            dup_tandem_supports[dup_tandem_key] = 1
                        
                        elif dup_tandem_key in dup_tandem_supports:
                            dup_tandem_supports[dup_tandem_key] += 1
                        
                        # Store length for variant:
                        if dup_tandem_key not in dup_tandem_lengths:
                            dup_tandem_lengths[dup_tandem_key] = variant_stats[1]
                    
                    else:
                        if line.strip().split("\t")[2].split(".")[2] == "BND":
                            continue
                        else:
                            other.append(line.strip())

    return([insertions, insertions_lengths, insertions_supports, deletions, deletions_lengths, deletions_supports,
            inversions, inversions_lengths, inversions_supports, dup_tandem, dup_tandem_lengths, dup_tandem_supports, other])
    
def get_variant_stats(vcf_line):
    """Parse vcf line to get basic values for each SV."""

    variant_start = int(vcf_line.strip().split("\t")[1])
    variant_length = abs(int(vcf_line.strip().split("\t")[7].split(";")[2].split("=")[-1]))

    return [variant_start, variant_length]

def get_stats(stat_element):
    """Get basic stats for a vector of values."""

    output_string = ""

    # If function receives list (length stats):
    if type(stat_element) is list:

        if len(stat_element) > 1:

            vector_len = len(stat_element)
            vector_avg = round(statistics.mean(stat_element), 2)
            vector_min = min(stat_element)
            vector_max = max(stat_element)
            vector_quantiles = "\t".join(
                [str(i) for i in statistics.quantiles(stat_element, n=4)]
            )

            output_string = f"{vector_len}\t{vector_avg}\t{vector_min}\t{vector_max}\t{vector_quantiles}"

        elif len(stat_element) == 1:

            vector_len = 1
            vector_avg = round(statistics.mean(stat_element), 2)
            vector_min = min(stat_element)
            vector_max = max(stat_element)
            vector_median = statistics.median(stat_element)

            output_string = f"{vector_len}\t{vector_avg}\t{vector_min}\t{vector_max}\tNA\t{vector_median}\tNA"
    
    # If function receives dict (support stats):
    elif type(stat_element) is dict:

        stat_list = [stat_element[i] for i in stat_element]

        if len(stat_list) > 1:

            vector_len = len(stat_list)
            vector_avg = round(statistics.mean(stat_list), 2)
            vector_min = min(stat_list)
            vector_max = max(stat_list)
            vector_quantiles = "\t".join(
                [str(i) for i in statistics.quantiles(stat_list, n=4)]
            )

            output_string = f"{vector_len}\t{vector_avg}\t{vector_min}\t{vector_max}\t{vector_quantiles}"

        elif len(stat_list) == 1:

            vector_len = 1
            vector_avg = round(statistics.mean(stat_list), 2)
            vector_min = min(stat_list)
            vector_max = max(stat_list)
            vector_median = statistics.median(stat_list)

            output_string = f"{vector_len}\t{vector_avg}\t{vector_min}\t{vector_max}\tNA\t{vector_median}\tNA"

    return output_string

# --------------------------------------------------
if __name__ == '__main__':
    main()
