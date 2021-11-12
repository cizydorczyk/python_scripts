import argparse
from Bio import SeqIO
import sys
import os.path
import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

parser = argparse.ArgumentParser()

parser.add_argument("--input_annotations", help="input gff file from RAST")
parser.add_argument("--fasta", help="fasta file of contigs used in RAST annotation")
parser.add_argument("--output_file", help="Formatted GFF file ready to be used with PIRATE/Panaroo")
parser.add_argument("--contigs_to_remove", help="comma sep list of contigs to remove from gff", default="")
parser.add_argument("--filetype", help="gff or gbk to parse, default=gff", default="gff")

args = parser.parse_args()

contigs_to_remove = args.contigs_to_remove.strip().split(",")

if args.filetype not in ["gff", "gbk"]:
    print("\nPlease specify filetype as either 'gff' or 'gbk'. Exiting.")
    sys.exit()
elif args.filetype == "":
    print("\nNo filetype specified. Parsing input file as 'gff'.\n")

print(f"\nRunning parsing script with the following options:\n\n"
    f"\tFiletype to parse: {args.filetype}\n"
    f"\tInput file: {os.path.basename(args.input_annotations)}\n"
    f"\tInput contigs: {os.path.basename(args.fasta)}\n"
    f"\tNumber of contigs to remove: {len(contigs_to_remove)}\n"
    f"\tOutput file: {os.path.basename(args.output_file)}\n")

######################
### Main script 
######################

if args.filetype == "gff":
    contigs = list(SeqIO.parse(args.fasta, "fasta"))
    print(f"\nTotal number of contigs in input contigs fasta: {len(contigs)}\n")
    
    sequence_regions = []
    fasta_seqs = []
    for i in contigs:
        if i.id not in contigs_to_remove:
            sequence_region = "##sequence-region " + str(i.id) + "\t1\t" + str(len(i.seq))# Create sequence-region strings
            sequence_regions.append(sequence_region)

            fasta_seq = ">" + str(i.id) + "\n" + str(i.seq).upper() # Create fasta sequences to add to end of file; uppercase all letters
            fasta_seqs.append(fasta_seq)
        else:
            print(f"Skipping contig {i.id} because it was requested to be removed.")
            continue

    sequence_region_to_write = "\n".join(sequence_regions)
    fasta_seqs_to_write = "\n".join(fasta_seqs)
    fasta_line = "##FASTA"

    # Read original gff file:
    # orig_gff = "\n".join([line.rstrip("\n") for line in open(args.input_annotations)][1:])
    gff_header = "##gff-version 3"

    filt_gff_list = [] # filtered gff entries list
    
    orig_gff_annot_by_contig = {} # keep track of # of annotations/contig
    for contig in contigs:
        orig_gff_annot_by_contig[contig.id] = 0 # initialize at 0 before reading gff
    
    with open(args.input_annotations, 'r') as infile1:
        # orig_num_annotations = 0 # initialize count of ALL gff entries
        for line in infile1:
            if not line.startswith("#"):
                # orig_num_annotations += 1
                contig = line.strip().split("\t")[0]
                
                orig_gff_annot_by_contig[contig] += 1 # increase annotation count for contig

                if contig not in contigs_to_remove:
                    filt_gff_list.append(line.strip())
                else:
                    # print(f"Omitting annotations on contig {i.id} as it was requested to be removed.")
                    continue
        # print(f"\nTotal number of annotations in gff file: {orig_num_annotations}")
        # print(f"Total number of annotations after filtering: {len(filt_gff_list)}")
    
    ####################################
    ### Print some useful information...
    ####################################
    print(f"\nBreakdown of number of annotations per contig in original gff file:\n")
    for contig in orig_gff_annot_by_contig:
        sum_annot = sum([orig_gff_annot_by_contig[contig] for contig in orig_gff_annot_by_contig])        
        print(f"\t{contig}\t{orig_gff_annot_by_contig[contig]}")
    print(f"\nTotal number of annotations in original gff file: {sum_annot}")

    print(f"Total of number of annotations per contig in filtered gff file: {len(filt_gff_list)}\n")
    #####################################

    filt_gff = "\n".join(filt_gff_list) # join output gff entries by newlines to facilitate writing to file
    
    # Organize output file:
    to_write = gff_header + "\n" + sequence_region_to_write + "\n" + filt_gff + "\n" + fasta_line + "\n" + fasta_seqs_to_write

    # Write output file:
    with open(args.output_file, 'w') as outfile:
        outfile.write(to_write)

if args.filetype == "gbk":
    input_gbk_records = SeqIO.parse(args.input_annotations, "genbank")
    
    output_gbk_records = (record for record in input_gbk_records if record.id not in contigs_to_remove)

    SeqIO.write(output_gbk_records, args.output_file, "genbank")    