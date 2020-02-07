import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--ref", help="reference fasta sequence")
parser.add_argument("--seg_vcf", help="vcf with only segregating sites (all sites, not just core!)")

args = parser.parse_args()

###############################################
#
# Note: this script works to parse a vcf file and insert the proper  ref bases based on a reference fasta file
# when the vcf contains ref positions from the first isolate; i.e. generated from a fasta file using snp-sites
# without the reference included in the fasta file.
#
# HOWEVER: IT DOES NOT CURRENTLY FIX THE PRESENCE/ABSENCE FIELD OF THE VCF FILE. As such, it's output is not
# yet fully functional.
#
# Use with caution. But the code that is here is good at least.
#
###############################################

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
def parse_vcf(vcf_file):
    variant_objects_list_1 = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                line_elements = line.strip().split("\t")

                variant_objects_list_1.append(VcfVariant(line_elements[0], int(line_elements[1]), line_elements[2], line_elements[3], line_elements[4], line_elements[5], line_elements[6], line_elements[7], line_elements[8], "\t".join(line_elements[9:]), "\t".join(line_elements)))

    return variant_objects_list_1

# Read VCF file
vcf_elements = parse_vcf(args.seg_vcf)

# Read reference
ref = SeqIO.read(args.ref, "fasta")

# Heavy lifting:
output_vcf_elements = []

for record in vcf_elements:
    ref_base = ref.seq[record.pos - 1]
    # print(ref_base + "\t" + record.ref + "\t" + record.alt)
    if record.ref == ref_base:
        complete_record = "\t".join([record.chrom, str(record.pos), record.id, record.ref, record.alt, record.qual, record.filter, ".", record.format, record.pa])
        new_record = VcfVariant(record.chrom, str(record.pos), record.id, record.ref, record.alt, record.qual, record.filter, ".", record.format, record.pa, complete_record)
        output_vcf_elements.append(new_record)

    elif record.ref != ref_base:
        original_alt_bases = record.alt.split(",")
        new_alt_bases = original_alt_bases.copy()
        new_alt_bases.remove(ref_base)
        new_alt_bases.append(record.ref)
        # print(str(record.pos) + "\t" + ref_base + "\t" + record.ref + "\t" + record.alt + "\t" + ",".join(new_alt_bases))

        complete_record = "\t".join([record.chrom, str(record.pos), record.id, ref_base, ",".join(new_alt_bases), record.qual, record.filter, ".", record.format, record.pa])
        new_record = VcfVariant(record.chrom, str(record.pos), record.id, ref_base, ",".join(new_alt_bases), record.qual, record.filter, ".", record.format, record.pa, complete_record)
        output_vcf_elements.append(new_record)

for i in output_vcf_elements:
    print(i.record)
