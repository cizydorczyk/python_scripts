import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--vcf", help="vcf file to get segregating sites from")
parser.add_argument("--num_isolates", help="number of isolates, not including reference, in vcf/fasta file")
parser.add_argument("--output_vcf", help="vcf of segregating sites only; final output of script")
parser.add_argument("--summary_table", help="summary table of vcf results; appropriate for paper supplementary table")

args = parser.parse_args()

# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format, pa, record):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = format
        self.pa = pa #presence absence of isolates; variable in length so all gets captured together
        self.record = record

# Define function to parse VCF files and create VcfVariant class objects for each variant:
def parse_vcf(vcf_file, header_list):
    variant_objects_list_1 = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                line_elements = line.strip().split("\t")
                elements_length = len(line_elements)
                # For example, 28 total vcf elements/fields, 18 isolates not including reference
                # pa therefore is line_elements[28-18:]

                variant_objects_list_1.append(VcfVariant(line_elements[0], int(line_elements[1]), line_elements[2], line_elements[3], line_elements[4],
                line_elements[5], line_elements[6], line_elements[7], line_elements[8], line_elements[elements_length-int(args.num_isolates):],
                "\t".join(line_elements)))
            elif line.startswith("#"): # append header lines to a list; will be used for writing header later
                header_list.append(line.strip())

    return variant_objects_list_1

# Read VCF file
header_list = []
vcf_elements = parse_vcf(args.vcf, header_list)

# Function to check if all elements in list equal:
def checkEqual(list):
    return len(set(list)) <= 1

# Function to determine whether site is segregating or simply absent in some isolates but different from ref:
def checkSegregating(alt_bases, pa):
    alt_list = alt_bases.split(",")

    if "*" not in alt_list:
        return_value = True
    elif "*" in alt_list:
        if len(alt_list) == 1: # shouldn't see this from snp-sites anyway...but here as a sanity check
            return_value = False
        elif len(alt_list) > 1:
            if alt_list[0] == "*":
                pa = [y for y in pa if y != "1"]
                if checkEqual(pa):
                    return_value = False
                elif checkEqual(pa) == False:
                     return_value = True
            elif alt_list[1] == "*":
                pa = [y for y in pa if y != "2"]
                if checkEqual(pa):
                    return_value = False
                elif checkEqual(pa) == False:
                    return_value = True

    return return_value

output_vcf_elements = []

# Remove vcf records that are all identical among my isolates but differ from reference:
for i in vcf_elements:
    if checkEqual(i.pa):
        continue
    else:
        output_vcf_elements.append(i)

# Remove vcf records where the only difference in 1/more isolates is a missing base (*):
final_vcf_elements = []
for i in output_vcf_elements:
    if checkSegregating(i.alt, i.pa):
        final_vcf_elements.append(i)
    else:
        continue

# Create output strings:
output_list = [i.record for i in final_vcf_elements]
output_string = "\n".join(output_list)
header_string = "\n".join(header_list)

# Create summary table:
summary_table = []

for i in final_vcf_elements:
    ann_elements = i.info.split("|")
    gene = ann_elements[4]
    gene_name = ann_elements[3]
    variant_effect = ann_elements[1]

    summary_string = "\t".join([i.chrom, str(i.pos), i.ref, i.alt, variant_effect, gene, gene_name, "\t".join(i.pa)])
    summary_table.append(summary_string)

summary_table_string = "\n".join(summary_table)

# Write output vcf to file w/ header:
with open(args.output_vcf, 'w') as outfile1:
    outfile1.write("\n".join([header_string, output_string]))

# Write output summary table to file:
with open(args.summary_table, 'w') as outfile2:
    outfile2.write(summary_table_string)
