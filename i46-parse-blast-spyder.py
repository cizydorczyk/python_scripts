#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 20:53:43 2021

@author: conrad
"""

from Bio.Blast import NCBIXML
from Bio import SeqIO

# Uses fna blast results to get pegs for relevant gene(s); then uses id of that/those pegs to get corresponding faa seqs
# Pegs are same for fna/faa
# blastp find multiple results & makes parsing difficult...blastn finds less variable seq & makes it easy
def getSeqs(xml_file, fna_file, faa_file):
    xml_handle = open(xml_file)
    xml_records = NCBIXML.parse(xml_handle)
    
    headers = []
    for record in xml_records:
        for alignment in record.alignments:
            header = alignment.title.split(" ")[1]
            if header not in headers:
                headers.append(header)
    xml_handle.close()
    
    fna_records = list(SeqIO.parse(fna_file, "fasta"))
    faa_records = list(SeqIO.parse(faa_file, "fasta"))
    
    output_fna_records = []
    output_faa_records = []
    
    for fnarecord in fna_records:
        if fnarecord.id in headers:
            output_fna_records.append(fnarecord)
    
    output_fna_records_ids = [fnarec.id for fnarec in output_fna_records]
    for faarecord in faa_records:
        if faarecord.id in output_fna_records_ids:
            output_faa_records.append(faarecord)
    
    return([output_fna_records, output_faa_records])
            
### Test function ###
# Should produce 2 seq records each for fna/faa, ending in peg.2186 and peg.2185
# test_output = getSeqs("/home/conrad/hinfluenzae/i46-tbpab/blast-output/tbp1/fna/A013-H12-09-07-2003-tbp1-raw-blast-output.xml",
#                       "/home/conrad/hinfluenzae/i46-tbpab/isolate-fna/A013-H12-09-07-2003.fna",
#                       "/home/conrad/hinfluenzae/i46-tbpab/isolate-faa/A013-H12-09-07-2003.faa")
# print(test_output)
### ###

# Read in isolate list:
isolate_list = []
with open("/home/conrad/hinfluenzae/hi_isolate_list.txt", "r") as infile1:
    for line in infile1:
        isolate_list.append(line.strip())

# Get seqs for each isolate, getting fna & faa seqs one isolate at a time, 
# doing both tbp1 and tbp2 before moving on to next isolate...
for isolate in isolate_list:
    tbp1_fna_blast_xml = "/home/conrad/hinfluenzae/i46-tbpab/blast-output/tbp1/fna/" + isolate + "-tbp1-raw-blast-output.xml"
    tbp2_fna_blast_xml = "/home/conrad/hinfluenzae/i46-tbpab/blast-output/tbp2/fna/" + isolate + "-tbp2-raw-blast-output.xml"
    
    isolate_fna = "/home/conrad/hinfluenzae/i46-tbpab/isolate-fna/" + isolate + ".fna"
    isolate_faa = "/home/conrad/hinfluenzae/i46-tbpab/isolate-faa/" + isolate + ".faa"
    
    # For tbp1:
    tbp1_seqs = getSeqs(tbp1_fna_blast_xml, isolate_fna, isolate_faa)
    
    for nuc_record in tbp1_seqs[0]:
        nuc_record.id = isolate + " tbpA"
    for aa_record in tbp1_seqs[1]:
        aa_record.id = isolate + " tbpA"
    
    # For tbp2:
    tbp2_seqs = getSeqs(tbp2_fna_blast_xml, isolate_fna, isolate_faa)
    
    for nuc_record in tbp2_seqs[0]:
        nuc_record.id = isolate + " tbpB"
    for aa_record in tbp2_seqs[1]:
        aa_record.id = isolate + " tbpB"
    
    # Write output:
    tbp1_fna_output_file = "/home/conrad/hinfluenzae/i46-tbpab/output-sequences/tbp1/fna/" + isolate + "-tbp1.fna"
    tbp1_faa_output_file = "/home/conrad/hinfluenzae/i46-tbpab/output-sequences/tbp1/faa/" + isolate + "-tbp1.faa"
    
    tbp2_fna_output_file = "/home/conrad/hinfluenzae/i46-tbpab/output-sequences/tbp2/fna/" + isolate + "-tbp2.fna"
    tbp2_faa_output_file = "/home/conrad/hinfluenzae/i46-tbpab/output-sequences/tbp2/faa/" + isolate + "-tbp2.faa"
    
    with open(tbp1_fna_output_file, "w") as outfile1:
        SeqIO.write(tbp1_seqs[0], outfile1, "fasta")
    
    with open(tbp1_faa_output_file, "w") as outfile2:
        SeqIO.write(tbp1_seqs[1], outfile2, "fasta")
        
    with open(tbp2_fna_output_file, "w") as outfile3:
        SeqIO.write(tbp2_seqs[0], outfile3, "fasta")
    
    with open(tbp2_faa_output_file, "w") as outfile4:
        SeqIO.write(tbp2_seqs[1], outfile4, "fasta")























