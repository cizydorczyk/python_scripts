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
    
    for entry in vcf_elements:
        
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
            
            st_variant_dict[record].append(peg1_name)
            st_variant_dict[record].append(peg2_name)
            
    return(st_variant_dict)

###############################################################################

#### ST-5 ####

###############################################################################

# Set VCF file handle & number of isolates:
st5_vcf = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st5/mmgene-testing/ann.chrom-renamed.st5.rc-masked.with-ref.vcf"
st5_num_isolates = 8

# Read VCF file
st5_header_list = []
vcf_elements = parse_vcf(st5_vcf, st5_header_list, st5_num_isolates)

# A057 isolate pair:
a057_pa_coords_list = [1,2]
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
st5_gff = "/home/conrad/m-sm-notebook/snp-calling/m5-snp-calling/snp-pipeline/st5/mmgene-testing/st5.concatenated-reference.gff"
st5_gff_records = parse_gff3(st5_gff)

st5_gff_dict = {}
for record in st5_gff_records:
    peg = record.attributes.split(";")[0].split("|")[1] # peg
    name = record.attributes.split(";")[1] # name
    st5_gff_dict[peg] = name

# Append names to st5 vcf variant dict:
st5_dict = get_gff_names(st5_dict, st5_gff_dict)

###############################################################################

#### ST-23 ####

###############################################################################

st23_vcf = 



        
        
        
        
        
        
        
        
        