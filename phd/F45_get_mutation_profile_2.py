import argparse
import os.path
from Bio import SeqIO

### NOTE
# The script ignores complex variants currently. It also reports the GC content for the reference genomes and NOT the isolate genomes.


parser = argparse.ArgumentParser()

parser.add_argument("--input_vcf", help="input vcf file")
parser.add_argument("--ref", help="reference genome")
parser.add_argument("--out", help="output file")
parser.add_argument("--isolate", help="isolate number")

args = parser.parse_args()

print(args.isolate)

snps = []
insertions = []
deletions = []

# Read reference fasta & get sequence:
reference_sequence = list(SeqIO.parse(args.ref, "fasta"))[0].seq

# Get length of reference sequence:
ref_length = len(reference_sequence)

# Get %GC content of reference sequence:
c_count = reference_sequence.count("C")
g_count = reference_sequence.count("G")

gc_percent = round(((c_count + g_count) / ref_length) * 100, 2)

print("Percent GC\t%s" % gc_percent)

# Define VcfVariants class;
class VcfVariant(object):
    def __init__(self, isolateid, chrom, pos, ref, alt, filename):
        self.isolateid = isolateid
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.filename = filename

# Classify each variant position in VCF file by type of variant (SNP, insertion, or deletion):
with open(args.input_vcf, 'r') as infile1:
    #isolate_id = os.path.basename(args.input_vcf)[0:-4]
    for line in infile1:
        if not line.startswith("#"):
            line_list = line.strip().split('\t')
            if len(line_list[3]) == 1 and len(line_list[4]) == 1:
                vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], args.input_vcf)
                snps.append(vcf_variant_obj)
            elif len(line_list[3]) > 1:
                vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], args.input_vcf)
                deletions.append(vcf_variant_obj)
            elif len(line_list[4]) > 1:
                if "," in line_list[4]:
                    bases = line_list[4].split(",")
                    bases.remove("*")

                    for base in bases:
                        vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], base, args.input_vcf)
                        snps.append(vcf_variant_obj)
                        # print(vcf_variant_obj.ref, vcf_variant_obj.alt)

                else:
                    # print("indel")
                    pass

# Define numbers of each type of variant (SNP, insertion, or deletion):
num_snps = len(snps)
num_insertions = len(insertions)
num_deletions = len(deletions)

# Get types of mutations
nuc_mutations = {"A->C":0, "A->G":0, "A->T":0, "C->G":0, "C->A":0, "C->T":0, "G->C":0, "G->A":0, "G->T":0, "T->A":0, "T->C":0, "T->G":0}

for variant in snps:
    if len(variant.alt) > 1:
        print(variant.alt)
    if variant.ref == "A" and variant.alt == "C":
        nuc_mutations["A->C"] += 1
    elif variant.ref == "A" and variant.alt == "G":
        nuc_mutations["A->G"] += 1
    elif variant.ref == "A" and variant.alt == "T":
        nuc_mutations["A->T"] += 1

    elif variant.ref == "C" and variant.alt == "G":
        nuc_mutations["C->G"] += 1
    elif variant.ref == "C" and variant.alt == "A":
        nuc_mutations["C->A"] += 1
    elif variant.ref == "C" and variant.alt == "T":
        nuc_mutations["C->T"] += 1

    elif variant.ref == "G" and variant.alt == "C":
        nuc_mutations["G->C"] += 1
    elif variant.ref == "G" and variant.alt == "A":
        nuc_mutations["G->A"] += 1
    elif variant.ref == "G" and variant.alt == "T":
        nuc_mutations["G->T"] += 1

    elif variant.ref == "T" and variant.alt == "A":
        nuc_mutations["T->A"] += 1
    elif variant.ref == "T" and variant.alt == "C":
        nuc_mutations["T->C"] += 1
    elif variant.ref == "T" and variant.alt == "G":
        nuc_mutations["T->G"] += 1


# Get raw mutation rates:
print("Raw mutations profile:")
for i in nuc_mutations:
    print("\tMutation:\t" + i + "\tRaw number of:\t" + str(nuc_mutations[i]))

raw_nuc_mutation_rates = {}

raw_nuc_mutation_rates["A->C"] = nuc_mutations["A->C"] / num_snps
raw_nuc_mutation_rates["A->G"] = nuc_mutations["A->G"] / num_snps
raw_nuc_mutation_rates["A->T"] = nuc_mutations["A->T"] / num_snps

raw_nuc_mutation_rates["C->A"] = nuc_mutations["C->A"] / num_snps
raw_nuc_mutation_rates["C->G"] = nuc_mutations["C->G"] / num_snps
raw_nuc_mutation_rates["C->T"] = nuc_mutations["C->T"] / num_snps

raw_nuc_mutation_rates["G->A"] = nuc_mutations["G->A"] / num_snps
raw_nuc_mutation_rates["G->C"] = nuc_mutations["G->C"] / num_snps
raw_nuc_mutation_rates["G->T"] = nuc_mutations["G->T"] / num_snps

raw_nuc_mutation_rates["T->A"] = nuc_mutations["T->A"] / num_snps
raw_nuc_mutation_rates["T->C"] = nuc_mutations["T->C"] / num_snps
raw_nuc_mutation_rates["T->G"] = nuc_mutations["T->G"] / num_snps

print("Raw mutation rates:")
for i in raw_nuc_mutation_rates:
    print("\tMutation:\t" + i + "\tRaw number of:\t" + str(raw_nuc_mutation_rates[i]))

# Get corrected mutation rates as per Transition bias influences the evolution of antibiotic resistance in Mycobacterium tuberculosis:
gc_correction = gc_percent/(100 - gc_percent)

corrected_mutation_rates = {}

corrected_num_snps = (nuc_mutations["A->C"] * gc_correction) + (nuc_mutations["A->G"] * gc_correction) + (nuc_mutations["A->T"] * gc_correction) + (nuc_mutations["T->A"] * gc_correction) + (nuc_mutations["T->C"] * gc_correction) + (nuc_mutations["T->G"] * gc_correction) + nuc_mutations["C->A"] + nuc_mutations["C->G"] + nuc_mutations["C->T"] + nuc_mutations["G->A"] + nuc_mutations["G->C"] + nuc_mutations["G->T"]
# print(corrected_num_snps)

corrected_mutation_rates["A->C"] = (nuc_mutations["A->C"] * gc_correction) / corrected_num_snps
corrected_mutation_rates["A->G"] = (nuc_mutations["A->G"] * gc_correction) / corrected_num_snps
corrected_mutation_rates["A->T"] = (nuc_mutations["A->T"] * gc_correction) / corrected_num_snps

corrected_mutation_rates["T->A"] = (nuc_mutations["T->A"] * gc_correction) / corrected_num_snps
corrected_mutation_rates["T->C"] = (nuc_mutations["T->C"] * gc_correction) / corrected_num_snps
corrected_mutation_rates["T->G"] = (nuc_mutations["T->G"] * gc_correction) / corrected_num_snps

corrected_mutation_rates["C->A"] = nuc_mutations["C->A"] / corrected_num_snps
corrected_mutation_rates["C->G"] = nuc_mutations["C->G"] / corrected_num_snps
corrected_mutation_rates["C->T"] = nuc_mutations["C->T"] / corrected_num_snps

corrected_mutation_rates["G->A"] = nuc_mutations["G->A"] / corrected_num_snps
corrected_mutation_rates["G->C"] = nuc_mutations["G->C"] / corrected_num_snps
corrected_mutation_rates["G->T"] = nuc_mutations["G->T"] / corrected_num_snps

print("Corrected mutation rates:")
for i in corrected_mutation_rates:
    print("\tMutation:\t" + i + "\tRaw number of:\t" + str(corrected_mutation_rates[i]))

# Get & print corrected (relative) transition & transversion rates:
transitions_rate = corrected_mutation_rates["A->G"] + corrected_mutation_rates["G->A"] + corrected_mutation_rates["C->T"] + corrected_mutation_rates["T->C"]

transversions_rate = corrected_mutation_rates["A->C"] + corrected_mutation_rates["A->T"] + corrected_mutation_rates["G->C"] + corrected_mutation_rates["G->T"] + corrected_mutation_rates["C->A"] + corrected_mutation_rates["C->G"] + corrected_mutation_rates["T->A"] + corrected_mutation_rates["T->G"]

print("Transitions rate:\t%s" % transitions_rate)
print("Transversions rate:\t%s" % transversions_rate)

transitions_transversions_rate_ratio = transitions_rate / transversions_rate

print("Transition:Transversion rate ratio:\t%s" % transitions_transversions_rate_ratio)

# Write data to file:
header = "isolate" + '\t' + "%GC" + '\t' + "num_snps" + '\t' + "num_insertions" + '\t' + "num_deletions" + '\t' + \
    "corrected_num_snps" + '\t' + "ts_tv_rel_rate_ratio" + '\t' + \
    "corrected_AC_rate" + '\t' + "corrected_AG_rate" + '\t' + "corrected_AT_rate" + '\t' + \
    "corrected_CA_rate" + '\t' + "corrected_CG_rate" + '\t' + "corrected_CT_rate" + '\t' + \
    "corrected_GA_rate" + '\t' + "corrected_GC_rate" + '\t' + "corrected_GT_rate" + '\t' + \
    "corrected_TA_rate" + '\t' + "corrected_TC_rate" + '\t' + "corrected_TC_rate" + '\t' + \
    "num_AC_snps" + '\t' + "num_AG_snps" + '\t' + "num_AT_snps" + '\t' +\
    "num_CA_snps" + '\t' + "num_CG_snps" + '\t' + "num_CT_snps" + '\t' +\
    "num_GA_snps" + '\t' + "num_GC_snps" + '\t' + "num_GT_snps" + '\t' +\
    "num_TA_snps" + '\t' + "num_TC_snps" + '\t' + "num_TG_snps" + '\t' +\
    "raw_AC_rate" + '\t' + "raw_AG_rate" + '\t' + "raw_AT_rate" + '\t' +\
    "raw_CA_rate" + '\t' + "raw_CG_rate" + '\t' + "raw_CT_rate" + '\t' +\
    "raw_GA_rate" + '\t' + "raw_GC_rate" + '\t' + "raw_GT_rate" + '\t' +\
    "raw_TA_rate" + '\t' + "raw_TC_rate" + '\t' + "raw_TG_rate"

to_write = args.isolate + '\t' + str(gc_percent) + '\t' + str(num_snps) + '\t' + str(num_insertions) + '\t' + str(num_deletions) + '\t' + \
    str(corrected_num_snps) + '\t' + str(transitions_transversions_rate_ratio) + '\t' +\
    str(corrected_mutation_rates["A->C"]) + '\t' + str(corrected_mutation_rates["A->G"]) + '\t' + str(corrected_mutation_rates["A->T"]) + '\t' + \
    str(corrected_mutation_rates["C->A"]) + '\t' + str(corrected_mutation_rates["C->G"]) + '\t' + str(corrected_mutation_rates["C->T"]) + '\t' + \
    str(corrected_mutation_rates["G->A"]) + '\t' + str(corrected_mutation_rates["G->C"]) + '\t' + str(corrected_mutation_rates["G->T"]) + '\t' + \
    str(corrected_mutation_rates["T->A"]) + '\t' + str(corrected_mutation_rates["T->C"]) + '\t' + str(corrected_mutation_rates["T->G"]) + '\t' + \
    str(nuc_mutations["A->C"]) + '\t' + str(nuc_mutations["A->G"]) + '\t' + str(nuc_mutations["A->T"]) + '\t' + str(nuc_mutations["C->A"]) + '\t' + \
    str(nuc_mutations["C->G"]) + '\t' + str(nuc_mutations["C->T"]) + '\t' + str(nuc_mutations["G->A"]) + '\t' + str(nuc_mutations["G->C"]) + '\t' + \
    str(nuc_mutations["G->T"]) + '\t' + str(nuc_mutations["T->A"]) + '\t' + str(nuc_mutations["T->C"]) + '\t' + str(nuc_mutations["T->G"]) + '\t' + \
    str(raw_nuc_mutation_rates["A->C"]) + '\t' + str(raw_nuc_mutation_rates["A->G"]) + '\t' + str(raw_nuc_mutation_rates["A->T"]) + '\t' + \
    str(raw_nuc_mutation_rates["C->A"]) + '\t' + str(raw_nuc_mutation_rates["C->G"]) + '\t' + str(raw_nuc_mutation_rates["C->T"]) + '\t' + \
    str(raw_nuc_mutation_rates["G->A"]) + '\t' + str(raw_nuc_mutation_rates["G->C"]) + '\t' + str(raw_nuc_mutation_rates["G->T"]) + '\t' + \
    str(raw_nuc_mutation_rates["T->A"]) + '\t' + str(raw_nuc_mutation_rates["T->C"]) + '\t' + str(raw_nuc_mutation_rates["T->G"])

# write results to file:
if os.path.isfile(args.out):
    with open(args.out, 'a') as outfile:
        outfile.write('\n' + to_write)
else:
    with open(args.out, 'w') as outfile:
        outfile.write(header + '\n' + to_write)
