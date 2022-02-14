import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--vcf")
parser.add_argument("--num_isolates", help="number of isolates, not including reference, in vcf/fasta file")


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

for entry in vcf_elements:
    print(entry.pa)
