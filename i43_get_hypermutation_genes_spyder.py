#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 10:36:12 2021

@author: conrad
"""
            
def parse_gff(gff_file):
    gff_records = {}
    with open(gff_file, "r") as infile1:
        for line in infile1:
            if not line.startswith("#"):
                line_elements = line.strip().split("\t")
                chromosome = line_elements[0]
                start = line_elements[3]
                stop = line_elements[4]
                annotation_element = ".".join(line_elements[8].split(";")[0].split(".")[-2:])
                gff_records[annotation_element] = [chromosome, start, stop]
    return(gff_records)

def parse_blast_output(input_file):
    with open(input_file) as infile1:
        for line in infile1:
            if not line.startswith("#"):
                line_elements = line.strip().split("\t")
                annotation_element = ".".join(line_elements[1].split(".")[-2:])
    return(annotation_element)

def read_vcf(input_file):
    vcf_records = {}
    with open(input_file, "r") as infile1:
        for line in infile1:
            if not line.startswith("#"):
                line_elements = line.strip().split("\t")
                entry_id = line_elements[0] + "_" + line_elements[1]
                vcf_records[entry_id] = line.strip()
    return(vcf_records)                
   
def read_snps_tab(input_file):
    tab_records = {}
    with open(input_file, "r") as infile1:
        for line in infile1:
            if not line.startswith("CHROM"):
                line_elements = line.strip().split("\t")
                entry_id = line_elements[0] + "_" + line_elements[1]
                tab_records[entry_id] = line.strip()
    return(tab_records)

hypermutator_genes_list = ["dam", "dnaQ", "mutH",  
                           "mutL", "mutM", "mutS", 
                           "mutT", "mutY", "polA", 
                           "uvrD"]

# =============================================================================
# ST-11
# =============================================================================

st11_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A323-H173-21-03-2011.gff")

st11_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st11/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st11_pegs[i] = annotation_peg


st11_isolate_list = []
with open("/home/conrad/hinfluenzae/st11_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st11_isolate_list.append(line.strip())

for isolate in st11_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st11_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st11_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {}
    
    for key in st11_pegs:
        gene_loc = st11_reference_gff[st11_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st11/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))
        
# =============================================================================
# ST-12
# =============================================================================

st12_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A274-H194-06-03-2006.gff")

st12_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st12/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st12_pegs[i] = annotation_peg


st12_isolate_list = []
with open("/home/conrad/hinfluenzae/st12_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st12_isolate_list.append(line.strip())

for isolate in st12_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st12_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st12_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st12_pegs:
        gene_loc = st12_reference_gff[st12_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st12/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))
        
# =============================================================================
# ST-103
# =============================================================================

st103_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A276-H42-02-09-2009.gff")

st103_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st103/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st103_pegs[i] = annotation_peg


st103_isolate_list = []
with open("/home/conrad/hinfluenzae/st103_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st103_isolate_list.append(line.strip())

for isolate in st103_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st103_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st103_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st103_pegs:
        gene_loc = st103_reference_gff[st103_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st103/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-14
# =============================================================================

st14_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A367-H205-09-12-2010.gff")

st14_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st14/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st14_pegs[i] = annotation_peg


st14_isolate_list = []
with open("/home/conrad/hinfluenzae/st14_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st14_isolate_list.append(line.strip())

for isolate in st14_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st14_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st14_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st14_pegs:
        gene_loc = st14_reference_gff[st14_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st14/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-105
# =============================================================================

st105_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A367-H211-26-10-2011.gff")

st105_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st105/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st105_pegs[i] = annotation_peg


st105_isolate_list = []
with open("/home/conrad/hinfluenzae/st105_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st105_isolate_list.append(line.strip())

for isolate in st105_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st105_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st105_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st105_pegs:
        gene_loc = st105_reference_gff[st105_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st105/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-124
# =============================================================================

st124_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A274-H191-02-05-2005.gff")

st124_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st124/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st124_pegs[i] = annotation_peg


st124_isolate_list = []
with open("/home/conrad/hinfluenzae/st124_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st124_isolate_list.append(line.strip())

for isolate in st124_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st124_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st124_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st124_pegs:
        gene_loc = st124_reference_gff[st124_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st124/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-145
# =============================================================================

st145_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A406-H224-11-05-2015.gff")

st145_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st145/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st145_pegs[i] = annotation_peg


st145_isolate_list = []
with open("/home/conrad/hinfluenzae/st145_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st145_isolate_list.append(line.strip())

for isolate in st145_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st145_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st145_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st145_pegs:
        gene_loc = st145_reference_gff[st145_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st145/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-203
# =============================================================================

st203_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A367-H204-07-12-2010.gff")

st203_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st203/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st203_pegs[i] = annotation_peg


st203_isolate_list = []
with open("/home/conrad/hinfluenzae/st203_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st203_isolate_list.append(line.strip())

for isolate in st203_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st203_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st203_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st203_pegs:
        gene_loc = st203_reference_gff[st203_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st203/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-321
# =============================================================================

st321_reference_gff = parse_gff("/home/conrad/hinfluenzae/pangenome_analyses/i41_public_private_pangenome/st321/rasttk_annotations/raw_annotated_gff/SAMN03450976.gff")

st321_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st321/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st321_pegs[i] = annotation_peg


st321_isolate_list = []
with open("/home/conrad/hinfluenzae/st321_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st321_isolate_list.append(line.strip())

for isolate in st321_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i39_st321_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i39_st321_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st321_pegs:
        gene_loc = st321_reference_gff[st321_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st321/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-393
# =============================================================================

st393_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A360-H04-26-03-2012.gff")

st393_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st393/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st393_pegs[i] = annotation_peg


st393_isolate_list = []
with open("/home/conrad/hinfluenzae/st393_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st393_isolate_list.append(line.strip())

for isolate in st393_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st393_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st393_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st393_pegs:
        gene_loc = st393_reference_gff[st393_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st393/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-1034
# =============================================================================

st1034_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A058-H216-28-09-2012.gff")

st1034_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st1034/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st1034_pegs[i] = annotation_peg


st1034_isolate_list = []
with open("/home/conrad/hinfluenzae/st1034_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st1034_isolate_list.append(line.strip())

for isolate in st1034_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st1034_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/st1034_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st1034_pegs:
        gene_loc = st1034_reference_gff[st1034_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/st1034/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))

# =============================================================================
# ST-n1
# =============================================================================

stn1_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A370-H186-20-04-2016.gff")

stn1_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/stn1/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    stn1_pegs[i] = annotation_peg


stn1_isolate_list = []
with open("/home/conrad/hinfluenzae/stn1_isolate_list.txt", "r") as infile2:
    for line in infile2:
        stn1_isolate_list.append(line.strip())

for isolate in stn1_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/stn1_snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i22_snp_calling/stn1_snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in stn1_pegs:
        gene_loc = stn1_reference_gff[stn1_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/blast_output/stn1/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))
        
# =============================================================================
# ST-n1 ALTERNATE REFERENCE (A370-H14)
# =============================================================================

stn1_alt_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A370-H14-21-11-2012.gff")

stn1_alt_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/stn1_alt_ref/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    stn1_alt_pegs[i] = annotation_peg


stn1_alt_isolate_list = []
with open("/home/conrad/hinfluenzae/stn1_isolate_list.txt", "r") as infile2:
    for line in infile2:
        stn1_alt_isolate_list.append(line.strip())

for isolate in stn1_alt_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i27_alt_ref_snp_calling/stn1/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i27_alt_ref_snp_calling/stn1/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in stn1_alt_pegs:
        gene_loc = stn1_alt_reference_gff[stn1_alt_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/stn1_alt_ref/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))
        
# =============================================================================
# ST-393 ALTERNATE REFERENCE (A037-H11)
# =============================================================================

st393_alt_reference_gff = parse_gff("/home/conrad/hinfluenzae/genome_annotations/i33_rast_test/exported_gff/A037-H11-04-08-2004.gff")

st393_alt_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/st393_alt_ref/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st393_alt_pegs[i] = annotation_peg


st393_alt_isolate_list = []
with open("/home/conrad/hinfluenzae/st393_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st393_alt_isolate_list.append(line.strip())

for isolate in st393_alt_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/snp_calling/i27_alt_ref_snp_calling/st393/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/snp_calling/i27_alt_ref_snp_calling/st393/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st393_alt_pegs:
        gene_loc = st393_alt_reference_gff[st393_alt_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/st393_alt_ref/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))
        
# =============================================================================
# ST-393 PUBLIC REFERENCE (GCF_003203035.1_ASM320303v1_genomic.fna)
# =============================================================================

st393_pub_reference_gff = parse_gff("/home/conrad/hinfluenzae/i43_hypermutator_testing/st393_public_ref/ncbi-genomes-2021-09-17/rasttk-annotation/GCF_003203035.1.gff")

st393_pub_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i43_hypermutator_testing/st393_public_ref/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st393_pub_pegs[i] = annotation_peg


st393_pub_isolate_list = []
with open("/home/conrad/hinfluenzae/st393_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st393_pub_isolate_list.append(line.strip())

for isolate in st393_pub_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/i43_hypermutator_testing/st393_public_ref/snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/i43_hypermutator_testing/st393_public_ref/snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st393_pub_pegs:
        gene_loc = st393_pub_reference_gff[st393_pub_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = "_".join(record.split("_")[0:2])
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i43_hypermutator_testing/st393_public_ref/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))
        
# =============================================================================
# i44 ST-203
# =============================================================================

st203_reference_gff = parse_gff("/home/conrad/hinfluenzae/i44/st203/rasttk_annotation/SRR9847108.gff")

st203_pegs = {}
for i in hypermutator_genes_list:
    blast_file = "/home/conrad/hinfluenzae/i44/st203/hypermutator_genes/" + i + ".blast_output.txt"
    annotation_peg = parse_blast_output(blast_file)
    st203_pegs[i] = annotation_peg


st203_isolate_list = []
with open("/home/conrad/hinfluenzae/st203_isolate_list.txt", "r") as infile2:
    for line in infile2:
        st203_isolate_list.append(line.strip())


for isolate in st203_isolate_list:
    isolate_vcf = "/home/conrad/hinfluenzae/i44/st203/snp_calling/raw_snippy_output/" + isolate + "/snps.vcf"
    isolate_tab = "/home/conrad/hinfluenzae/i44/st203/snp_calling/raw_snippy_output/" + isolate + "/snps.tab"
    
    isolate_vcf_records = read_vcf(isolate_vcf)
    isolate_tab_records = read_snps_tab(isolate_tab)
    record_keys = [key for key in isolate_vcf_records]

    tab_records_to_keep = {"dam":[], "dnaQ":[], "mutH":[], "mutL":[], "mutM":[], "mutS":[], "mutT":[], "mutY":[], "polA":[], "uvrD":[]}
    
    for key in st203_pegs:
        gene_loc = st203_reference_gff[st203_pegs[key]]
        chrom = gene_loc[0]
        start = int(gene_loc[1])
        stop = int(gene_loc[2])
    
        count = 0
        for record in record_keys:
            record_chrom = record.split("_")[0] # record chromosome format needed channging compared to above since using public reference
            record_position = int(record.split("_")[-1])
        
            if record_chrom == chrom and record_position >= start and record_position <= stop:
                count += 1
                tab_record = isolate_tab_records[record].split("\t")
                vcf_record = isolate_vcf_records[record]
                vcf_impact = vcf_record.split("\t")[7].split(";")[7].split("|")[2]
                tab_record.append(vcf_impact)
                if "LOF" in vcf_record:
                    vcf_lof = "Predicted loss of function"
                else:
                    vcf_lof = ""
                tab_record.append(vcf_lof)
                        
                tab_records_to_keep[key].append("\t".join(tab_record))
    
    
    
    output_handle = "/home/conrad/hinfluenzae/i44/st203/hypermutator_genes/" + isolate + ".hypermutator_genes.txt"
    with open(output_handle, "w") as outfile1:
        to_write_list = []
        for key in tab_records_to_keep:
            if len(tab_records_to_keep[key]) > 0:
                record = "\n".join(tab_records_to_keep[key])
                to_write_list.append(record)
        outfile1.write("\n".join(to_write_list))