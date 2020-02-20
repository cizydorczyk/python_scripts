import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--input_no_density_snp_tsv", help="snvphyl snvTable.tsv file")
parser.add_argument("--density", help="# snps in window to filter")
parser.add_argument("--window", help="Window size in which to filter 'density' number of SNPs")
parser.add_argument("--ref", help="reference genome")

args = parser.parse_args()

vcf_dict = {}
vcf_pos_list = []
isolate_dict = {}
with open(args.input_no_density_snp_tsv, 'r') as infile1:
    for line in infile1:
        if not line.startswith("#"):
            # if "filtered-coverage" not in line:
            #     if "filtered-mpileup" not in line:
            #         if "filtered-invalid" not in line:
            line_elements = line.strip().split("\t")
            vcf_dict[line_elements[1]] = line_elements
            vcf_pos_list.append(int(line_elements[1]))

with open(args.input_no_density_snp_tsv, 'r') as infile2:
    for line in infile2:
        if line.startswith("#"):
            line_el = line.strip().split("\t")
            isolates = line_el[3:]
            for i in isolates:
                isolate_dict[i] = [isolates.index(i) + 3]
        elif not line.startswith("#"): #and "filtered" not in line:
            line_el2 = line.strip().split("\t")
            for i in isolates:
                isolate_dict[i].append(line_el2[isolate_dict[i][0]])

from Bio import SeqIO
ref = list(SeqIO.parse(args.ref, "fasta"))[0].seq

def merge_regions(regions):

    if len(regions) == 0:
        return regions
    regions = sorted(regions)
    merged_regions = list()
    merged_regions.append(regions[0])
    for region in regions[1:]:
        last_merged_region = merged_regions[-1]
        last_merged_region_start, last_merged_region_end = last_merged_region
        region_start, region_end = region
        if region_start >= last_merged_region_start and region_end <= last_merged_region_end:
            pass # discard region contained in the last region
        elif region_start <= (last_merged_region_end + 1) and region_end > last_merged_region_end:
            merged_regions[-1] = (last_merged_region_start, region_end) # extend last region by overlapping or adjacent region
        else:
            merged_regions.append(region) # add non-overlapping region to sorted list
    return merged_regions

def find_dense_regions(max_allowed_snps, window_size, snps):

    snp_count = len(snps)
    dense_region_list = []
    for idx, pos_start in enumerate(snps):
        # print(idx, pos_start)
        if (idx + max_allowed_snps) < snp_count: # if not out of range
            pos_end = snps[idx + max_allowed_snps] # get snp at index corresponding to max # allowed SNPs in window
            print(idx, pos_start, pos_end)
            if (pos_start + window_size - 1) >= pos_end:
                dense_region_list.append((pos_start, pos_end))
                print(pos_start, pos_end)
    dense_region_list = merge_regions(dense_region_list)
    return dense_region_list

f = find_dense_regions(2, 25, [1,24,25,26,100,101, 155])
print(f)

def in_region(pos, regions):
    """Find whether a position is included in a region.
    Parameters
    ----------
    pos : int
        DNA base position.
    regions : list of tuples
        List of (start, end) position integers.
    Returns
    -------
    bool
        True if the position is within an of the regions, False otherwise.
    Examples
    --------
    # Empty list
    >>> in_region(1, [])
    False
    # In list
    >>> in_region(10, [(3, 5), (9, 12)])
    True
    # Not in list
    >>> in_region(10, [(3, 5), (11, 12)])
    False
    """
    for region in regions:
        if (pos >= region[0]) and (pos <= region[1]):
            return True

    return False


# to_mask = []
# for i in range(0,len(ref)-25+1+1):
#     start = i+1
#     stop = i+1+25-1
#
#     num_snps = 0
#     for i in vcf_pos_list:
#         if start <= i <= stop:
#             num_snps += 1
#
#     if num_snps > 2:
#         to_mask.append((start, stop))
#
# for i in vcf_pos_list:
#     for j in to_mask:
#         if j[0] <= i <= j[1]:
#             continue
#         else:
#             print(i)















# print(len(vcf_pos_list))
# for i in isolate_dict:
#     print(len(isolate_dict[i][1:]))
#
# isolate_output_pos = {}


# for isolate in isolate_dict:
#     isolate_output_pos[isolate] = set()
#     if isolate.lower() != "reference":
#         isolate_seq = isolate_dict[isolate][1:]
#         ref_seq = isolate_dict["Reference"][1:]
#
#         ref_bases = []
#         alt_bases = []
#         removed = []
#
#         for i, j, k in zip(ref_seq, isolate_seq, vcf_pos_list):
#             if j in "N-" or i in "N-":
#                 removed.append(k)
#             else:
#                 if i == j:
#                     ref_bases.append(k)
#                 elif i != j:
#                     alt_bases.append(k)
#
#         print(isolate, len(ref_bases), len(alt_bases), len(removed))







        # ref_positions = [k for i, j, k in zip(ref_seq, isolate_seq, vcf_pos_list) if i == j]
        # for i in ref_positions:
        #     isolate_output_pos[isolate].add(i)
        #
        # non_ref_bases = [j for i, j in zip(ref_seq, isolate_seq) if i != j]
        # non_ref_pos = [k for i, j, k in zip(ref_seq, isolate_seq, vcf_pos_list) if i != j]
        #
        # for i, j in zip(non_ref_bases, non_ref_pos):
        #
        #     try:
        #         next_pos = non_ref_pos[non_ref_pos.index(j)+1]
        #     except IndexError:
        #         next_pos = j + 10000
        #
        #     try:
        #         prev_pos = non_ref_pos[non_ref_pos.index(j)-1]
        #     except IndexError:
        #         prev_pos = j - 10000
        #
        #     if (next_pos - j) > 25 and (j - prev_pos) > 25 and i not in "N-": # test threshold; easy in case of 2 in XX window
        #         isolate_output_pos[isolate].add("j")
        #     else:
        #         continue
# final_set = set()
# for i in isolate_output_pos:
#     if len(final_set) == 0:
#         final_set = isolate_output_pos[i]
#     else:
#         final_set = final_set & isolate_output_pos[i]
#
# for i in sorted(list(final_set)):
#     print(i)
#
# print(len(final_set))










# for j in isolate_dict:
    # print(j) # print isolates in dict, including reference!
    # print(j + "\t" + "".join(isolate_dict[i][1:])) # print isolate SNP sequences (non-filtered)
    # print(len(isolate_dict[i][1:])) # print isolate SNP sequence length (non-filtered)
#
# for key in isolate_dict:
#     if key.lower() != "reference":
#         ref_seq = isolate_dict["Reference"][1:]
#         isolate_seq = isolate_dict[key][1:]
#
#         output_seq = []
#
#         for i in range(0,len(ref_seq)):
#             ref_base = ref_seq[i]
#             isolate_base = isolate_seq[i]
#             pos = vcf_pos_list[i]
#
#             try:
#                 next_pos = vcf_pos_list[i+1]
#                 nref_base = ref_seq[i+1]
#                 nisolate_base = isolate_seq[i+1]
#             except IndexError:
#                 next_pos = pos
#                 nref_base = ref_base
#                 nisolate_base = isolate_base
#
#             if ref_base in "N-":
#                 continue
#             elif isolate_base in "N-":
#                 continue
#             elif ref_base == isolate_base:
#                 output_seq.append(isolate_base)
#             elif ref_base != isolate_base:
#                 if (next_pos - pos) <= 25 and nref_base != nisolate_base:
#                     continue
#                 else:
#                     # print(vcf_dict[str(pos)])
#                     output_seq.append(isolate_base)
#         print(len(output_seq))
