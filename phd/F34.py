import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--input_vcf", help="isolate-specific vcf file")
parser.add_argument("--density", help="# snps in window to filter", type=int)
parser.add_argument("--window", help="Window size in which to filter 'density' number of SNPs", type=int)
#parser.add_argument("--ref", help="reference genome")
parser.add_argument("--fasta", help="fasta from which vcf was generated")
parser.add_argument("--sample", help="sample that is being analyzed")
parser.add_argument("--mask_char", help="character with which to mask sequence, default = 'N'", default="N")
parser.add_argument("--output_fasta", help="output fasta file -- can be appended to")

args = parser.parse_args()

# for each entry in vcf:
#     set isolate base to proper alt/ref base
#     if base is not same as ref:
#         keep entry
#     else:
#         skip entry


def parse_vcf(vcf_file):
    """
    Parse a VCF file & return a dictionary with VCF entry positions as keys
    and the entire VCF line (i.e. the VCF entry) as the value.

    The function parses the VCF, keeping only entries that are NOT identical
    to the reference and do not contain ambiguous sites.

    -------------------
    Parameters
    -------------------
        vcf_file = single isolate VCF file with genotype indicated by 0/1/2/etc.
        where 0 = reference allele, 1 = first alt allele, 2 = second alt allele,
        etc.
    """
    vcf_pos_dict = {}
    with open(vcf_file, 'r') as infile1:
        for line in infile1:
            if not line.startswith("#"):

                line_elements = line.strip().split("\t") # get vcf line elements
                pos = line_elements[1]

                ref_base = [line_elements[3]] # Get reference base
                alt_bases = line_elements[4].strip().split(",") # get all alt bases
                possible_bases = ref_base + alt_bases # create list of ref + alt bases
                genotype = int(line_elements[9]) # get genotype
                alt_base = possible_bases[genotype] # assign alt base based on the index of all bases genotype corresponds to

                # print(line_elements[1], possible_bases, genotype, alt_base)

                if alt_base != "*":
                    if alt_base == ref_base[0]:
                        continue
                    elif alt_base != ref_base[0]:
                        vcf_pos_dict[int(pos)] = line_elements
                else:
                    # print(line_elements[1], possible_bases, genotype, alt_base)
                    continue
    return vcf_pos_dict

def find_dense_regions(max_allowed_snps, window_size, snps): # Function taken from CF-SAN
    """Scan a list of snp positions to find regions where the snp density exceeds the
    allowed thershold.
    Parameters
    ----------
    max_allowed_snps : int
        Maximum allowed number of snps in a given rolling window.
    window_size : int
        Size of rolling window along the length of a genome which
        is scanned for excessive snps.
    snps : list of int
        Sorted list of snp positions
    Returns
    -------
    dense_region_list : list of tuples
        List of (start_position, end_position) tuples identifying the list of
        dense snp regions.
    Examples
    --------
    # Empty list
    >>> find_dense_regions(3, 1000, [])
    []
    # Not dense window
    >>> find_dense_regions(3, 1000, [1, 2, 3, 1001])
    []
    # One more than max_allowed_snps at window boundaries
    >>> find_dense_regions(3, 1000, [1, 20, 30, 1000])
    [(1, 1000)]
    # Two more than max_allowed_snps at window boundaries
    >>> find_dense_regions(3, 1000, [1, 20, 30, 40, 1000])
    [(1, 1000)]
    # Overlapping dense regions with combined size greater than window_size
    >>> find_dense_regions(3, 1000, [1, 20, 30, 40, 501, 600, 1000, 1500])
    [(1, 1500)]
    # Multiple dense regions
    >>> find_dense_regions(3, 1000, [1, 2, 3, 1000, 1500, 3001, 3002, 3003, 4000])
    [(1, 1000), (3001, 4000)]
    """
    snp_count = len(snps)
    dense_region_list = []
    for idx, pos_start in enumerate(snps):
        if (idx + max_allowed_snps) < snp_count:
            pos_end = snps[idx + max_allowed_snps]
            if (pos_start + window_size - 1) >= pos_end:
                dense_region_list.append((pos_start, pos_end))
    dense_region_list = merge_regions(dense_region_list)
    return dense_region_list

def merge_regions(regions): # Function taken from CF-SAN
    """Coalesce regions.
    Scans a sorted list of region starting and ending positions looking
    for the outer-most start and end positions to coalesce overlapping
    and contained regions into a smaller list of larger regions.
    Parameters
    ----------
    regions : list of tuples
        List of (start, end) position integers.
    Returns
    -------
    regions : list of tuples
        List of merged (start, end) position integers.
    Examples
    --------
    >>> # Empty list
    >>> merge_regions([])
    []
    >>> # Only one region
    >>> merge_regions([(10,20)])
    [(10, 20)]
    >>> # Discard contained region at left
    >>> merge_regions([(10,20), (10,15)])
    [(10, 20)]
    >>> # Discard contained region at right
    >>> merge_regions([(10,20), (15,20)])
    [(10, 20)]
    >>> # Discard contained region exact match
    >>> merge_regions([(10,20), (10,20)])
    [(10, 20)]
    >>> # Discard contained region fully contained
    >>> merge_regions([(10,20), (11,19)])
    [(10, 20)]
    >>> # Extend region by overlap right
    >>> merge_regions([(10,20), (15,25)])
    [(10, 25)]
    >>> # Extend region by overlap left
    >>> merge_regions([(10,20), (5,15)])
    [(5, 20)]
    >>> # Extend immediately adjacent region by extension
    >>> merge_regions([(10,20), (21,30)])
    [(10, 30)]
    >>> # No overlap
    >>> merge_regions([(40,50), (25,30)])
    [(25, 30), (40, 50)]
    >>> # Single position region : discard contained region
    >>> merge_regions([(40,50), (40,40)])
    [(40, 50)]
    >>> # Single position region : discard contained region
    >>> merge_regions([(40,50), (50,50)])
    [(40, 50)]
    >>> # Single position region : discard contained region
    >>> merge_regions([(40,50), (41,41)])
    [(40, 50)]
    >>> # Single position region : discard contained region
    >>> merge_regions([(40,50), (49,49)])
    [(40, 50)]
    >>> # Single position region : extend immediately adjacent region by extension
    >>> merge_regions([(10,10), (11,21)])
    [(10, 21)]
    >>> # Single position region : extend immediately adjacent region by extension
    >>> merge_regions([(10,20), (21,21)])
    [(10, 21)]
    >>> # Single position region : merge two immediately adjacent single-position regions
    >>> merge_regions([(20,20), (21,21)])
    [(20, 21)]
    >>> # Single position region : no overlap
    >>> merge_regions([(40,50), (60,60)])
    [(40, 50), (60, 60)]
    >>> # Single position region : no overlap
    >>> merge_regions([(40,40), (50,60)])
    [(40, 40), (50, 60)]
    >>> # Single position region : no overlap
    >>> merge_regions([(40,40), (50,50)])
    [(40, 40), (50, 50)]
    """
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

def in_region(pos, regions): # Function taken from CF-SAN
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
        True if the position is within any of the regions, False otherwise.
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

def mask_fasta(fasta_sequence, invalid_vcf_entries):
    """
    Function to mask dense snp regions.
    -------------------
    Parameters

        fasta_sequence = fasta sequence object from BioPython SeqIO that will
                         be masked.

        invalid_vcf_entries = dict of vcf entries (generated above); sites in
                              this dict will be masked with "N" in the fasta
                              sequence.
    """



    sequence_list = [char for char in str(fasta_sequence.seq)] # convert seq to list
    sequence_header = fasta_sequence.id # get seq header
    pos_list = [key for key in invalid_vcf_entries] # get list of invalid positions

    for pos in pos_list:
        sequence_list[pos-1] = args.mask_char # set base to N for invalid positions; -1 b/c of 0-based indexing (0 = 1st position)

    output_seq_record = SeqRecord(Seq("".join(sequence_list)), id=sequence_header) # create new seq record object with masked sequence

    return(output_seq_record)

def check_masking(unmasked_seq_record, masked_seq_record, invalid_vcf_entries):
    if not ((masked_seq_record.seq.count(args.mask_char)) - (unmasked_seq_record.seq.count(args.mask_char))) == len(invalid_vcf_entries):

        masked_dif = ((masked_seq_record.seq.count(args.mask_char)) - (unmasked_seq_record.seq.count(args.mask_char)))
        len_to_mask = len(invalid_vcf_entries)
        print("Error: number of masked positions does not equal number of invalid positions. Something is wrong with the script...")
        print("Count of masked positions is: ", masked_dif)
        print("Count should be: ", len_to_mask)
        sys.exit()

# Main

print("Working on sample " + args.sample + "...")

# Get fasta sequence of isolate being analyzed
print("Parsing fasta file...")
fasta_records = list(SeqIO.parse(args.fasta, "fasta"))

sample_record = [record for record in fasta_records if record.id == args.sample][0]

print("Parsing vcf file...")
vcf_entries = parse_vcf(args.input_vcf)

vcf_positions_list = [key for key in vcf_entries]

print(len(vcf_positions_list))

print("Getting dense positions...")
dense_regions = find_dense_regions(args.density, args.window, vcf_positions_list)

valid_entries = {} # contains entries passing density filter
invalid_entries = {} # contains entries failing density filter (to be masked)
for entry in vcf_entries:
    if in_region(entry, dense_regions):
        invalid_entries[entry] = vcf_entries[entry]
    else:
        valid_entries[entry] = vcf_entries[entry]

print("Masking dense positions...")
masked_sample_record = mask_fasta(sample_record, invalid_entries)

## Sanity check: difference should = number of invalid positions

# print(len(invalid_entries))
# print(sample_record.seq.count("N"))
# print(masked_sample_record.seq.count("N"))
# print((masked_sample_record.seq.count("N")) - (sample_record.seq.count("N"))) # should equal len(invalid_entries)
check_masking(sample_record, masked_sample_record, invalid_entries)

print("Writing masked sequence to file...")

output_record_id = ">" + masked_sample_record.id
output_record_seq = str(masked_sample_record.seq)

if os.path.isfile(args.output_fasta):
    with open(args.output_fasta, 'a') as outfile1:
        outfile1.write("\n" + output_record_id + "\n" + output_record_seq)

elif not os.path.isfile(args.output_fasta):
    with open(args.output_fasta, 'w') as infile1:
        infile1.write(output_record_id + "\n" + output_record_seq)

# Uncomment lines below to print removed positions to stdout; useful to compare to CF-SAN "removed" VCF file
# for i in invalid_entries:
#     print(i)
