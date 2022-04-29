#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 14:41:37 2022

@author: conrad
"""

###############################################################################

# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format, genotypes, record):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = format
        self.genotypes = genotypes #presence absence of isolates; variable in length so all gets captured together
        self.record = record
        
# Define gff3 format class to allow easy access to entries of gff3 formatted file:
class GffEntry(object):
    def __init__(self, seq, source, feature_type, feature_start, feature_end, score, strand, phase, attributes, record):
        self.seq = seq
        self.source = source
        self.feature_type = feature_type
        self.feature_start = feature_start
        self.feature_end = feature_end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
        self.record = record

# Define function to parse VCF files and create VcfVariant class objects for each variant:
def parse_vcf(vcf_file, header_list, num_isolates):
    variant_objects_list_1 = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                line_elements = line.strip().split("\t")
                elements_length = len(line_elements)
                # For example, 28 total vcf elements/fields, 18 isolates not including reference
                # pa therefore is line_elements[28-18:]
                
                

                variant_objects_list_1.append(VcfVariant(line_elements[0], int(line_elements[1]), line_elements[2], line_elements[3], line_elements[4],
                line_elements[5], line_elements[6], line_elements[7], line_elements[8], line_elements[elements_length-int(num_isolates):],
                "\t".join(line_elements)))
            elif line.startswith("#"): # append header lines to a list; will be used for writing header later
                header_list.append(line.strip())

    return variant_objects_list_1

# Function to get SNPs occurring DURING infection; requires coords of isolates' column in vcf file:
def filter_vcf(vcf_elements_list, genotypes_coords_list):
    ambiguous_index_dict = {0:"1", 1:"2", 2:"3", 3:"4", 4:"5"}
    
    records_to_keep = []
    
    for entry in vcf_elements_list:
        
        genotypes_alleles = []
        for coord in genotypes_coords_list:
            genotypes_alleles.append(entry.genotypes[int(coord)])
        
        if len(set(genotypes_alleles)) != 1:
            if "*" in entry.alt:
                ambiguous_index = entry.alt.split(",").index("*")
                ambiguous_allele = ambiguous_index_dict[ambiguous_index]
                
                while ambiguous_allele in genotypes_alleles:
                    genotypes_alleles.remove(ambiguous_allele)
                    
                if len(set(genotypes_alleles)) > 1:
                    records_to_keep.append(entry)
            else:
                records_to_keep.append(entry)
    
    return(records_to_keep)

# Function to parse gff3 formated file:
def parse_gff3(gff3_file):
    gff_records = []
    with open(gff3_file, "r") as infile1:
        for line in infile1:
            if not line.startswith("#"):
                line_elements = line.strip().split("\t")
                gff_records.append(GffEntry(line_elements[0], line_elements[1], line_elements[2], line_elements[3],
                                            line_elements[4], line_elements[5], line_elements[6], line_elements[7],
                                            line_elements[8], "\t".join(line_elements)))
    return(gff_records)

# Function to get peg name from gff dict using peg from vcf variant dict:
def get_gff_names(st_variant_dict, gff_dict):
    for record in st_variant_dict:
        peg = st_variant_dict[record][0].info.split("|")[4]
        
        if "RAST2" not in peg:
            peg_name = gff_dict[peg]
            st_variant_dict[record].append(peg)
            st_variant_dict[record].append(peg_name)
        elif "RAST2" in peg:
            peg = st_variant_dict[record][0].info.split("|")[4:6]
            peg1 = peg[0].split("-RAST2")[0]
            peg2 = peg[1]
            
            peg1_name = gff_dict[peg1]
            peg2_name = gff_dict[peg2]
            
            # st_variant_dict[record].append(peg1_name)
            # st_variant_dict[record].append(peg2_name)
            
            st_variant_dict[record].append(peg1 + "-" + peg2)
            st_variant_dict[record].append(peg1_name + "-" + peg2_name)
            
    return(st_variant_dict)

###############################################################################

#### ST-5 ####

###############################################################################

# Set VCF file handle & number of isolates (including reference):
st5_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st5/mmgene-testing/5.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st5_num_isolates = 8

# Read VCF file
st5_header_list = []
vcf_elements = parse_vcf(st5_vcf, st5_header_list, st5_num_isolates)

# A057 isolate pair:
a057_pa_coords_list = [1,2] # 0-based, ref should always have coord 0
a057_records = filter_vcf(vcf_elements, a057_pa_coords_list)

# A376 isolate pair:
a376_pa_coords_list = [6,7]
a376_records = filter_vcf(vcf_elements, a376_pa_coords_list)

# Create combined positions dict:
st5_dict = {}
for record in a057_records:
    if record.pos not in st5_dict:
        st5_dict[record.pos] = [record]

for record in a376_records:
    if record.pos not in st5_dict:
        st5_dict[record.pos] = [record]

# Parse gff for obtaining gene names based on peg:
st5_gff = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st5/mmgene-testing/5.concatenated-reference.gff"
st5_gff_records = parse_gff3(st5_gff)

st5_gff_dict = {}
for record in st5_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg
    name = record.attributes.split(";")[1] # name
    st5_gff_dict[peg] = name

# Append names to st5 vcf variant dict:
st5_dict = get_gff_names(st5_dict, st5_gff_dict)

# Write output:
## Uncomment lines below to write to file:
# to_write = []
# for entry in st5_dict:
#     to_write.append(st5_dict[entry][0].record + "\t" + st5_dict[entry][1] + "\t" + st5_dict[entry][2])

# with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st5/mmgene-testing/5.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     st5_header = "\n".join(st5_header_list)
#     st5_records = "\n".join(to_write)
#     outfile1.write(st5_header + "\n" + st5_records)

###############################################################################

#### ST-23 ####

###############################################################################

st23_vcf = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st23/mmgene-testing/23.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st23_num_isolates = 5

st23_header_list = []
st23_vcf_elements = parse_vcf(st23_vcf, st23_header_list, st23_num_isolates)

# A130
a130_pa_coords_list = [1,2,3,4]
a130_records = filter_vcf(st23_vcf_elements, a130_pa_coords_list)

# Create positions dict:
st23_dict = {}
for record in a130_records:
    if record.pos not in st23_dict:
        st23_dict[record.pos] = [record]

# Parse gff:
st23_gff = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st23/mmgene-testing/23.concatenated-reference.gff"
st23_gff_records = parse_gff3(st23_gff)

st23_gff_dict = {}
for record in st23_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg]
    name = record.attributes.split(";")[1] # name
    st23_gff_dict[peg] = name

# Append names:
st23_dict = get_gff_names(st23_dict, st23_gff_dict)

# Write output:
## Uncomment lines below to write to file:
# to_write = []
# for entry in st23_dict:
#     to_write.append(st23_dict[entry][0].record + "\t" + st23_dict[entry][1] + "\t" + st23_dict[entry][2])

# with open("/home/conrad/m2-sm-notebook/m2-2-snp-calling/st23/mmgene-testing/23.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     st23_header = "\n".join(st23_header_list)
#     st23_records = "\n".join(to_write)
#     outfile1.write(st23_header + "\n" + st23_records)

###############################################################################

#### ST-39_n1 ####

###############################################################################



# ST39_n1 has no pairs of isolates from the same patient, so nothing to parse.



###############################################################################

#### ST-91 ####

###############################################################################

st91_vcf = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st91/mmgene-testing/91.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st91_num_isolates = 6

st91_header_list = []
st91_vcf_elements = parse_vcf(st91_vcf, st91_header_list, st91_num_isolates)

# A374
a374_pa_coords_list = [1,2,3,4,5]
a374_records = filter_vcf(st91_vcf_elements, a374_pa_coords_list)

# Create positions dict:
st91_dict = {}
for record in a374_records:
    if record.pos not in st91_dict:
        st91_dict[record.pos] = [record]
        
# Parse gff:
st91_gff = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st91/mmgene-testing/91.concatenated-reference.gff"
st91_gff_records = parse_gff3(st91_gff)

st91_gff_dict = {}
for record in st91_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg]
    name = record.attributes.split(";")[1] # name
    st91_gff_dict[peg] = name

# Append names:
st91_dict = get_gff_names(st91_dict, st91_gff_dict)

# Write output:
## Uncomment lines below to write to file:
# to_write = []
# for entry in st91_dict:
#     to_write.append(st91_dict[entry][0].record + "\t" + st91_dict[entry][1] + "\t" + st91_dict[entry][2])

# with open("/home/conrad/m2-sm-notebook/m2-2-snp-calling/st91/mmgene-testing/91.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     st91_header = "\n".join(st91_header_list)
#     st91_records = "\n".join(to_write)
#     outfile1.write(st91_header + "\n" + st91_records)

###############################################################################

#### ST-199 ####

###############################################################################

st199_vcf = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st199/mmgene-testing/199.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st199_num_isolates = 20

st199_header_list = []
st199_vcf_elements = parse_vcf(st199_vcf, st199_header_list, st199_num_isolates)

# A055:
a055_pa_coords_list = [2,3,4,5,6,7,8,9,10,11,12,13,14,15] # omit SM003 (index 1) due to distance to rest of isolates
a055_records = filter_vcf(st199_vcf_elements, a055_pa_coords_list)

# A090:
a090_pa_coords_list = [17,18,19]
a090_records = filter_vcf(st199_vcf_elements, a090_pa_coords_list)
    
# Create positions dict:
st199_dict = {}
for record in a055_records:
    if record.pos not in st199_dict:
        st199_dict[record.pos] = [record]

for record in a090_records:
    if record.pos not in st199_dict:
        st199_dict[record.pos] = [record]

# Parse gff:
st199_gff = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st199/mmgene-testing/199.concatenated-reference.gff"
st199_gff_records = parse_gff3(st199_gff)

st199_gff_dict = {}
for record in st199_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg]
    name = record.attributes.split(";")[1] # name
    st199_gff_dict[peg] = name

# Append names:
st199_dict = get_gff_names(st199_dict, st199_gff_dict)

# Write output:
# ## Uncomment lines below to write to file:
# to_write = []
# for entry in st199_dict:
#     to_write.append(st199_dict[entry][0].record + "\t" + st199_dict[entry][1] + "\t" + st199_dict[entry][2])

# with open("/home/conrad/m2-sm-notebook/m2-2-snp-calling/st199/mmgene-testing/199.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     st199_header = "\n".join(st199_header_list)
#     st199_records = "\n".join(to_write)
#     outfile1.write(st199_header + "\n" + st199_records)

###############################################################################

#### ST-220 ####

###############################################################################

st220_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st220/mmgene-testing/220.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st220_num_isolates = 4

st220_header_list = []
st220_vcf_elements = parse_vcf(st220_vcf, st220_header_list, st220_num_isolates)

# A371:
a371_pa_coords_list = [2,3]
a371_records = filter_vcf(st220_vcf_elements, a371_pa_coords_list)

# Create positions dict:
st220_dict = {}
for record in a371_records:
    if record.pos not in st220_dict:
        st220_dict[record.pos] = [record]

# Parse gff:
st220_gff = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st220/mmgene-testing/220.concatenated-reference.gff"
st220_gff_records = parse_gff3(st220_gff)

st220_gff_dict = {}
for record in st220_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg]
    name = record.attributes.split(";")[1] # name
    st220_gff_dict[peg] = name

# Append names:
st220_dict = get_gff_names(st220_dict, st220_gff_dict)

# Write output:
## Uncomment lines below to write to file:
# to_write = []
# for entry in st220_dict:
#     to_write.append(st220_dict[entry][0].record + "\t" + st220_dict[entry][1] + "\t" + st220_dict[entry][2])

# with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st220/mmgene-testing/220.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     st220_header = "\n".join(st220_header_list)
#     st220_records = "\n".join(to_write)
#     outfile1.write(st220_header + "\n" + st220_records)

###############################################################################

#### ST-224 ####

###############################################################################

st224_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st224/mmgene-testing/224.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st224_num_isolates = 4

st224_header_list = []
st224_vcf_elements = parse_vcf(st224_vcf, st224_header_list, st224_num_isolates)

# A013:
a013_pa_coords_list = [1,2]
a013_records = filter_vcf(st224_vcf_elements, a013_pa_coords_list)

# Create positions dict:
st224_dict = {}
for record in a013_records:
    if record.pos not in st224_dict:
        st224_dict[record.pos] = [record]

# Parse gff:
st224_gff = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st224/mmgene-testing/224.concatenated-reference.gff"
st224_gff_records = parse_gff3(st224_gff)

st224_gff_dict = {}
for record in st224_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg]
    name = record.attributes.split(";")[1] # name
    st224_gff_dict[peg] = name

# Append names:
st224_dict = get_gff_names(st224_dict, st224_gff_dict)

# Write output:
## Uncomment lines below to write to file:
# to_write = []
# for entry in st224_dict:
#     to_write.append(st224_dict[entry][0].record + "\t" + st224_dict[entry][1] + "\t" + st224_dict[entry][2])

# with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st224/mmgene-testing/224.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     st224_header = "\n".join(st224_header_list)
#     st224_records = "\n".join(to_write)
#     outfile1.write(st224_header + "\n" + st224_records)

###############################################################################

#### ST-246 ####

###############################################################################

st246_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st246/mmgene-testing/246.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st246_num_isolates = 3

st246_header_list = []
st246_vcf_elements = parse_vcf(st246_vcf, st246_header_list, st246_num_isolates)

# A057:
a057_pa_coords_list = [1,2]
a057_records = filter_vcf(st246_vcf_elements, a057_pa_coords_list)

# Create positions dict:
st246_dict = {}
for record in a057_records:
    if record.pos not in st246_dict:
        st246_dict[record.pos] = [record]

# Parse gff:
st246_gff = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st246/mmgene-testing/246.concatenated-reference.gff"
st246_gff_records = parse_gff3(st246_gff)

st246_gff_dict = {}
for record in st246_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg]
    name = record.attributes.split(";")[1] # name
    st246_gff_dict[peg] = name

# Append names:
st246_dict = get_gff_names(st246_dict, st246_gff_dict)

# Write output:
## Uncomment lines below to write to file:
# to_write = []
# for entry in st246_dict:
#     to_write.append(st246_dict[entry][0].record + "\t" + st246_dict[entry][1] + "\t" + st246_dict[entry][2])

# with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st246/mmgene-testing/246.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     st246_header = "\n".join(st246_header_list)
#     st246_records = "\n".join(to_write)
#     outfile1.write(st246_header + "\n" + st246_records)

###############################################################################

#### ST-365 ####

###############################################################################

st365_vcf = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st365/mmgene-testing/365.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st365_num_isolates = 4

st365_header_list = []
st365_vcf_elements = parse_vcf(st365_vcf, st365_header_list, st365_num_isolates)

# A057:
a344_pa_coords_list = [1,2,3]
a344_records = filter_vcf(st365_vcf_elements, a344_pa_coords_list)

# Create positions dict:
st365_dict = {}
for record in a344_records:
    if record.pos not in st365_dict:
        st365_dict[record.pos] = [record]

# Parse gff:
st365_gff = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st365/mmgene-testing/365.concatenated-reference.gff"
st365_gff_records = parse_gff3(st365_gff)

st365_gff_dict = {}
for record in st365_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg]
    name = record.attributes.split(";")[1] # name
    st365_gff_dict[peg] = name

# Append names:
st365_dict = get_gff_names(st365_dict, st365_gff_dict)

# Write output:
## Uncomment lines below to write to file:
# to_write = []
# for entry in st365_dict:
#     to_write.append(st365_dict[entry][0].record + "\t" + st365_dict[entry][1] + "\t" + st365_dict[entry][2])

# with open("/home/conrad/m2-sm-notebook/m2-2-snp-calling/st365/mmgene-testing/365.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     st365_header = "\n".join(st365_header_list)
#     st365_records = "\n".join(to_write)
#     outfile1.write(st365_header + "\n" + st365_records)

###############################################################################

#### ST-N2 ####

###############################################################################

stn2_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/stn2/mmgene-testing/n2.annotated.chrom-renamed.rc-masked.with-ref.vcf"
stn2_num_isolates = 4

stn2_header_list = []
stn2_vcf_elements = parse_vcf(stn2_vcf, stn2_header_list, stn2_num_isolates)

# A057:
a090_pa_coords_list = [1,2,3]
a090_records = filter_vcf(stn2_vcf_elements, a090_pa_coords_list)

# Create positions dict:
stn2_dict = {}
for record in a090_records:
    if record.pos not in stn2_dict:
        stn2_dict[record.pos] = [record]

# Parse gff:
stn2_gff = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/stn2/mmgene-testing/n2.concatenated-reference.gff"
stn2_gff_records = parse_gff3(stn2_gff)

stn2_gff_dict = {}
for record in stn2_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg]
    name = record.attributes.split(";")[1] # name
    stn2_gff_dict[peg] = name

# Append names:
stn2_dict = get_gff_names(stn2_dict, stn2_gff_dict)

# Write output:
## Uncomment lines below to write to file:
# to_write = []
# for entry in stn2_dict:
#     to_write.append(stn2_dict[entry][0].record + "\t" + stn2_dict[entry][1] + "\t" + stn2_dict[entry][2])

# with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/stn2/mmgene-testing/n2.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf", "w") as outfile1:
#     stn2_header = "\n".join(stn2_header_list)
#     stn2_records = "\n".join(to_write)
#     outfile1.write(stn2_header + "\n" + stn2_records)

###############################################################################

# Added February 22, 2022:

###############################################################################

# GET NON-CF MUTATIONS SEGREGATING BETWEEN PATIENTS

###############################################################################

# Function to filter mutations, considering ALL genotypes (except reference):
def filter_vcf2(vcf_elements_list):
    ambiguous_index_dict = {0:"1", 1:"2", 2:"3", 3:"4", 4:"5"}
    
    records_to_keep = []
    
    for entry in vcf_elements_list:
        
        genotypes_alleles = entry.genotypes[1:] # ignore reference '0' genotype
        
        if len(set(genotypes_alleles)) != 1:
            if "*" in entry.alt:
                ambiguous_index = entry.alt.split(",").index("*")
                ambiguous_allele = ambiguous_index_dict[ambiguous_index]
                
                while ambiguous_allele in genotypes_alleles:
                    genotypes_alleles.remove(ambiguous_allele)
                    
                if len(set(genotypes_alleles)) > 1:
                    records_to_keep.append(entry)
            else:
                records_to_keep.append(entry)
    
    return(records_to_keep)

# Define function to append variant type to vcf file:
def get_variant_type(variants_list):
    return_dict = {}
    for record in variants_list:
        mutation_type = record.info.split("|")[1]
        
        if mutation_type == "synonymous_variant":
            mutation_record = "1" + "\t" + "0" + "\t" + "0" + "\t" + "0"
            return_dict[record.pos] = [record, mutation_record]
        elif mutation_type == "missense_variant":
            mutation_record = "0" + "\t" + "1" + "\t" + "0" + "\t" + "0"
            return_dict[record.pos] = [record, mutation_record]
        elif mutation_type == "stop_gained":
            mutation_record = "0" + "\t" + "0" + "\t" + "1" + "\t" + "0"
            return_dict[record.pos] = [record, mutation_record]
        elif mutation_type == "intergenic_region":
            mutation_record = "0" + "\t" + "0" + "\t" + "0" + "\t" + "1"
            return_dict[record.pos] = [record, mutation_record]
        else:
            mutation_record = "0" + "\t" + "0" + "\t" + "0" + "\t" + "1"
            return_dict[record.pos] = [record, mutation_record]
            
    return(return_dict)

################### ST-5 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
st5_filtered_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st5/mmgene-testing/5.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st5_filtered_vcf_header_list = []
st5_filtered_vcf_records = parse_vcf(st5_filtered_vcf, st5_filtered_vcf_header_list, st5_num_isolates)

st5_filtered_positions = [i.pos for i in st5_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-5 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
st5_non_cf_mutations = []
for record in vcf_elements: # st5 records are in vcf_elements; other STs are in st##_vcf_elements
    if record.pos not in st5_filtered_positions:
        st5_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
st5_non_cf_filtered_records = filter_vcf2(st5_non_cf_mutations)

# Append type of mutation for easy calculations:    
st5_non_cf_filtered_records_dict = get_variant_type(st5_non_cf_filtered_records)

# Append names to st5 vcf variant dict:
st5_non_cf_filtered_records_dict = get_gff_names(st5_non_cf_filtered_records_dict, st5_gff_dict)

# Write output to file:
to_write = []
for record in st5_non_cf_filtered_records_dict:
    to_write.append(st5_non_cf_filtered_records_dict[record][0].record + "\t" + st5_non_cf_filtered_records_dict[record][1] + "\t" + st5_non_cf_filtered_records_dict[record][2] + "\t" + st5_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st5/mmgene-testing/st5.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(st5_filtered_vcf_header_list) + "\n" + "\n".join(to_write))

################### ST-23 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
st23_filtered_vcf = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st23/mmgene-testing/23.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st23_filtered_vcf_header_list = []
st23_filtered_vcf_records = parse_vcf(st23_filtered_vcf, st23_filtered_vcf_header_list, st23_num_isolates)

st23_filtered_positions = [i.pos for i in st23_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-23 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
st23_non_cf_mutations = []
for record in st23_vcf_elements:
    if record.pos not in st23_filtered_positions:
        st23_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
st23_non_cf_filtered_records = filter_vcf2(st23_non_cf_mutations)

# Append type of mutation for easy calculations:    
st23_non_cf_filtered_records_dict = get_variant_type(st23_non_cf_filtered_records)

# Append names to st23 vcf variant dict:
st23_non_cf_filtered_records_dict = get_gff_names(st23_non_cf_filtered_records_dict, st23_gff_dict)

# Write output to file:
to_write = []
for record in st23_non_cf_filtered_records_dict:
    to_write.append(st23_non_cf_filtered_records_dict[record][0].record + "\t" + st23_non_cf_filtered_records_dict[record][1] + "\t" + st23_non_cf_filtered_records_dict[record][2] + "\t" + st23_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m2-sm-notebook/m2-2-snp-calling/st23/mmgene-testing/st23.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(st23_filtered_vcf_header_list) + "\n" + "\n".join(to_write))

################### ST-91 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
st91_filtered_vcf = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st91/mmgene-testing/91.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st91_filtered_vcf_header_list = []
st91_filtered_vcf_records = parse_vcf(st91_filtered_vcf, st91_filtered_vcf_header_list, st91_num_isolates)

st91_filtered_positions = [i.pos for i in st91_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-91 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
st91_non_cf_mutations = []
for record in st91_vcf_elements:
    if record.pos not in st91_filtered_positions:
        st91_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
st91_non_cf_filtered_records = filter_vcf2(st91_non_cf_mutations)

# Append type of mutation for easy calculations:    
st91_non_cf_filtered_records_dict = get_variant_type(st91_non_cf_filtered_records)

# Append names to st91 vcf variant dict:
st91_non_cf_filtered_records_dict = get_gff_names(st91_non_cf_filtered_records_dict, st91_gff_dict)

# Write output to file:
to_write = []
for record in st91_non_cf_filtered_records_dict:
    to_write.append(st91_non_cf_filtered_records_dict[record][0].record + "\t" + st91_non_cf_filtered_records_dict[record][1] + "\t" + st91_non_cf_filtered_records_dict[record][2] + "\t" + st91_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m2-sm-notebook/m2-2-snp-calling/st91/mmgene-testing/st91.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(st91_filtered_vcf_header_list) + "\n" + "\n".join(to_write))

################### ST-199 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
st199_filtered_vcf = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st199/mmgene-testing/199.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st199_filtered_vcf_header_list = []
st199_filtered_vcf_records = parse_vcf(st199_filtered_vcf, st199_filtered_vcf_header_list, st199_num_isolates)

st199_filtered_positions = [i.pos for i in st199_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-199 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
st199_non_cf_mutations = []
for record in st199_vcf_elements:
    if record.pos not in st199_filtered_positions:
        st199_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
st199_non_cf_filtered_records = filter_vcf2(st199_non_cf_mutations)

# Append type of mutation for easy calculations:    
st199_non_cf_filtered_records_dict = get_variant_type(st199_non_cf_filtered_records)

# Append names to st91 vcf variant dict:
st199_non_cf_filtered_records_dict = get_gff_names(st199_non_cf_filtered_records_dict, st199_gff_dict)

# Write output to file:
to_write = []
for record in st199_non_cf_filtered_records_dict:
    to_write.append(st199_non_cf_filtered_records_dict[record][0].record + "\t" + st199_non_cf_filtered_records_dict[record][1] + "\t" + st199_non_cf_filtered_records_dict[record][2] + "\t" + st199_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m2-sm-notebook/m2-2-snp-calling/st199/mmgene-testing/st199.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(st199_filtered_vcf_header_list) + "\n" + "\n".join(to_write))

################### ST-220 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
st220_filtered_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st220/mmgene-testing/220.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st220_filtered_vcf_header_list = []
st220_filtered_vcf_records = parse_vcf(st220_filtered_vcf, st220_filtered_vcf_header_list, st220_num_isolates)

st220_filtered_positions = [i.pos for i in st220_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-220 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
st220_non_cf_mutations = []
for record in st220_vcf_elements:
    if record.pos not in st220_filtered_positions:
        st220_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
st220_non_cf_filtered_records = filter_vcf2(st220_non_cf_mutations)

# Append type of mutation for easy calculations:    
st220_non_cf_filtered_records_dict = get_variant_type(st220_non_cf_filtered_records)

# Append names to st220 vcf variant dict:
st220_non_cf_filtered_records_dict = get_gff_names(st220_non_cf_filtered_records_dict, st220_gff_dict)

# Write output to file:
to_write = []
for record in st220_non_cf_filtered_records_dict:
    to_write.append(st220_non_cf_filtered_records_dict[record][0].record + "\t" + st220_non_cf_filtered_records_dict[record][1] + "\t" + st220_non_cf_filtered_records_dict[record][2] + "\t" + st220_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st220/mmgene-testing/st220.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(st220_filtered_vcf_header_list) + "\n" + "\n".join(to_write))

################### ST-224 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
st224_filtered_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st224/mmgene-testing/224.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st224_filtered_vcf_header_list = []
st224_filtered_vcf_records = parse_vcf(st224_filtered_vcf, st224_filtered_vcf_header_list, st224_num_isolates)

st224_filtered_positions = [i.pos for i in st224_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-224 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
st224_non_cf_mutations = []
for record in st224_vcf_elements:
    if record.pos not in st224_filtered_positions:
        st224_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
st224_non_cf_filtered_records = filter_vcf2(st224_non_cf_mutations)

# Append type of mutation for easy calculations:    
st224_non_cf_filtered_records_dict = get_variant_type(st224_non_cf_filtered_records)

# Append names to st224 vcf variant dict:
st224_non_cf_filtered_records_dict = get_gff_names(st224_non_cf_filtered_records_dict, st224_gff_dict)

# Write output to file:
to_write = []
for record in st224_non_cf_filtered_records_dict:
    to_write.append(st224_non_cf_filtered_records_dict[record][0].record + "\t" + st224_non_cf_filtered_records_dict[record][1] + "\t" + st224_non_cf_filtered_records_dict[record][2] + "\t" + st224_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st224/mmgene-testing/st224.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(st224_filtered_vcf_header_list) + "\n" + "\n".join(to_write))

################### ST-246 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
st246_filtered_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st246/mmgene-testing/246.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st246_filtered_vcf_header_list = []
st246_filtered_vcf_records = parse_vcf(st246_filtered_vcf, st246_filtered_vcf_header_list, st246_num_isolates)

st246_filtered_positions = [i.pos for i in st246_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-246 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
st246_non_cf_mutations = []
for record in st246_vcf_elements:
    if record.pos not in st246_filtered_positions:
        st246_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
st246_non_cf_filtered_records = filter_vcf2(st246_non_cf_mutations)

# Append type of mutation for easy calculations:    
st246_non_cf_filtered_records_dict = get_variant_type(st246_non_cf_filtered_records)

# Append names to st246 vcf variant dict:
st246_non_cf_filtered_records_dict = get_gff_names(st246_non_cf_filtered_records_dict, st246_gff_dict)

# Write output to file:
to_write = []
for record in st246_non_cf_filtered_records_dict:
    to_write.append(st246_non_cf_filtered_records_dict[record][0].record + "\t" + st246_non_cf_filtered_records_dict[record][1] + "\t" + st246_non_cf_filtered_records_dict[record][2] + "\t" + st246_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st246/mmgene-testing/st246.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(st246_filtered_vcf_header_list) + "\n" + "\n".join(to_write))

################### ST-365 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
st365_filtered_vcf = "/home/conrad/m2-sm-notebook/m2-2-snp-calling/st365/mmgene-testing/365.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
st365_filtered_vcf_header_list = []
st365_filtered_vcf_records = parse_vcf(st365_filtered_vcf, st365_filtered_vcf_header_list, st365_num_isolates)

st365_filtered_positions = [i.pos for i in st365_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-365 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
st365_non_cf_mutations = []
for record in st365_vcf_elements:
    if record.pos not in st365_filtered_positions:
        st365_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
st365_non_cf_filtered_records = filter_vcf2(st365_non_cf_mutations)

# Append type of mutation for easy calculations:    
st365_non_cf_filtered_records_dict = get_variant_type(st365_non_cf_filtered_records)

# Append names to st365 vcf variant dict:
st365_non_cf_filtered_records_dict = get_gff_names(st365_non_cf_filtered_records_dict, st365_gff_dict)

# Write output to file:
to_write = []
for record in st365_non_cf_filtered_records_dict:
    to_write.append(st365_non_cf_filtered_records_dict[record][0].record + "\t" + st365_non_cf_filtered_records_dict[record][1] + "\t" + st365_non_cf_filtered_records_dict[record][2] + "\t" + st365_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m2-sm-notebook/m2-2-snp-calling/st365/mmgene-testing/st365.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(st365_filtered_vcf_header_list) + "\n" + "\n".join(to_write))

################### ST-n2 #########################

# Read output vcf from above:
# NOTE: parse_vcf() does not parse these VCF files properly, because it includes the additional annotation info in the presence/absence
# since it takes ALL columns after the start of presence/absence data. This does not matter here, however.
stn2_filtered_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/stn2/mmgene-testing/n2.filtered.annotated.chrom-renamed.rc-masked.with-ref.vcf"
stn2_filtered_vcf_header_list = []
stn2_filtered_vcf_records = parse_vcf(stn2_filtered_vcf, stn2_filtered_vcf_header_list, stn2_num_isolates)

stn2_filtered_positions = [i.pos for i in stn2_filtered_vcf_records]


# Next, make sure the original, unfiltered ST-n2 vcf has been read in above

# Identify vcf entries in original file not included in filtered file:
stn2_non_cf_mutations = []
for record in stn2_vcf_elements:
    if record.pos not in stn2_filtered_positions:
        stn2_non_cf_mutations.append(record)
    else:
        continue

# Filter non-CF mutations so that they are at least segregating within the ST
# Ignoring differences from only the ref (which is not RC masked):
stn2_non_cf_filtered_records = filter_vcf2(stn2_non_cf_mutations)

# Append type of mutation for easy calculations:    
stn2_non_cf_filtered_records_dict = get_variant_type(stn2_non_cf_filtered_records)

# Append names to stn2vcf variant dict:
stn2_non_cf_filtered_records_dict = get_gff_names(stn2_non_cf_filtered_records_dict, stn2_gff_dict)

# Write output to file:
to_write = []
for record in stn2_non_cf_filtered_records_dict:
    to_write.append(stn2_non_cf_filtered_records_dict[record][0].record + "\t" + stn2_non_cf_filtered_records_dict[record][1] + "\t" + stn2_non_cf_filtered_records_dict[record][2] + "\t" + stn2_non_cf_filtered_records_dict[record][3])

with open("/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/stn2/mmgene-testing/stn2.non-cf-filtered-mutations.vcf", "w") as outfile1:
    outfile1.write("\n".join(stn2_filtered_vcf_header_list) + "\n" + "\n".join(to_write))




