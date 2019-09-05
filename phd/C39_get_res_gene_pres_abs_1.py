import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--card_res_genes", help="file with all genes for which to get pres/abs; output of C39_get_all_res_genes_1.py script")
parser.add_argument("--input_file_list", help="input file listing all card output files, one per line")
parser.add_argument("--output_prefix", help="output file prefix")
parser.add_argument("--wd", help="working directory for output files")

args = parser.parse_args()

# Set working directory:
os.chdir(args.wd)

# Get list of output files:
file_list =[]
with open(args.input_file_list, "r") as infile1:
    for line in infile1:
        file_list.append(line.strip())

# Get resistance genes:
res_genes = {}
with open(args.card_res_genes, 'r') as infile2:
    for line in infile2:
        line_elements = line.strip().split("\t")
        res_genes[line_elements[0]] = line_elements[1:3]


# Initialize output dictionaries:

## presence absence dictionary with % identities recorded per gene:
pres_abs_ident_dict = {}
for i in file_list:
    pres_abs_ident_dict[i.split('/')[-1].split('_')[0]] = []

# Create class for card tab output format entries:
class ResGene(object):
    def __init__(self, ORF_ID, Contig, Start, Stop, Orientation, Cut_Off, Pass_Bitscore, Best_Hit_Bitscore, Best_Hit_ARO, Best_Identities, ARO, Model_type, SNPs_in_Best_Hit_ARO, Other_SNPs, Drug_Class, Resistance_Mechanism, AMR_Gene_Family, Predicted_DNA, Predicted_Protein, CARD_Protein_Sequence, Percentage_Length_of_Reference_Sequence, ID, Model_ID, Nudged, Note):
        self.ORF_ID = ORF_ID
        self.Contig = Contig
        self.Start = Start
        self.Stop = Stop
        self.Orientation = Orientation
        self.Cut_Off = Cut_Off
        self.Pass_Bitscore = Pass_Bitscore
        self.Best_Hit_Bitscore = Best_Hit_Bitscore
        self.Best_Hit_ARO = Best_Hit_ARO
        self.Best_Identities = Best_Identities
        self.ARO = ARO
        self.Model_type = Model_type
        self.SNPs_in_Best_Hit_ARO = SNPs_in_Best_Hit_ARO
        self.Other_SNPs = Other_SNPs
        self.Drug_Class = Drug_Class
        self.Resistance_Mechanism = Resistance_Mechanism
        self.AMR_Gene_Family = AMR_Gene_Family
        self.Predicted_DNA = Predicted_DNA
        self.Predicted_Protein = Predicted_Protein
        self.CARD_Protein_Sequence = CARD_Protein_Sequence
        self.Percentage_Length_of_Reference_Sequence = Percentage_Length_of_Reference_Sequence
        self.ID = ID
        self.Model_ID = Model_ID
        self.Nudged = Nudged
        self.Note = Note

# Main algorithm:

# Create list for recording order of genes in output files:
output_gene_order = []
for i in sorted(res_genes):
    output_gene_order.append(i)

for isolate in file_list:
    isolate_id = isolate.strip().split('/')[-1].split('_')[0]
    isolate_res_genes = {}
    with open(isolate, 'r') as infile3:
        for line in infile3:
            if not line.startswith("ORF"):
                line_elements = line.strip().split('\t')
                if len(line_elements) == 23:
                    resgene_obj = ResGene(line_elements[0], line_elements[1], line_elements[2], line_elements[3], line_elements[4], line_elements[5], line_elements[6], line_elements[7], line_elements[8], line_elements[9], line_elements[10], line_elements[11], line_elements[12], line_elements[13], line_elements[14], line_elements[15], line_elements[16], line_elements[17], line_elements[18], line_elements[19], line_elements[20], line_elements[21], line_elements[22], "na", "na")
                elif len(line_elements) == 25:
                    resgene_obj = ResGene(line_elements[0], line_elements[1], line_elements[2], line_elements[3], line_elements[4], line_elements[5], line_elements[6], line_elements[7], line_elements[8], line_elements[9], line_elements[10], line_elements[11], line_elements[12], line_elements[13], line_elements[14], line_elements[15], line_elements[16], line_elements[17], line_elements[18], line_elements[19], line_elements[20], line_elements[21], line_elements[22], line_elements[23], line_elements[24])

                isolate_res_genes[line_elements[8]] = resgene_obj

    for gene in sorted(res_genes):
        if gene in isolate_res_genes:

            if isolate_res_genes[gene].Cut_Off == "Perfect":

                if isolate_res_genes[gene].Model_type == "protein homolog model":
                    pm = "PHM"
                elif isolate_res_genes[gene].Model_type == "protein variant model":
                    pm = "PVM"
                elif isolate_res_genes[gene].Model_type == "protein overexpression model":
                    pm = "POM"
                elif isolate_res_genes[gene].Model_type == "toxin system meta model":
                    pm = "TSMM"
                elif isolate_res_genes[gene].Model_type == "rRNA gene variant model":
                    pm = "RGVM"
                elif isolate_res_genes[gene].Model_type == "pilus meta-model":
                    pm = "PMM"
                elif isolate_res_genes[gene].Model_type == "protein knockout model":
                    pm = "PKM"
                elif isolate_res_genes[gene].Model_type == "nonfunctional insertion model":
                    pm = "NFIM"
                elif isolate_res_genes[gene].Model_type == "gene cluster meta-model":
                    pm = "GCMM"
                elif isolate_res_genes[gene].Model_type == "efflux pump system meta-model":
                    pm = "EPSMM"
                elif isolate_res_genes[gene].Model_type == "protein domain meta-model":
                    pm = "PDMM"

                pres_abs_ident_dict[isolate_id].append(isolate_res_genes[gene].Best_Identities + "_" + isolate_res_genes[gene].Percentage_Length_of_Reference_Sequence + "_" + pm + "*")

            elif isolate_res_genes[gene].Cut_Off == "Strict" and isolate_res_genes[gene].Nudged == "na":

                if isolate_res_genes[gene].Model_type == "protein homolog model":
                    pm = "PHM"
                elif isolate_res_genes[gene].Model_type == "protein variant model":
                    pm = "PVM"
                elif isolate_res_genes[gene].Model_type == "protein overexpression model":
                    pm = "POM"
                elif isolate_res_genes[gene].Model_type == "toxin system meta model":
                    pm = "TSMM"
                elif isolate_res_genes[gene].Model_type == "rRNA gene variant model":
                    pm = "RGVM"
                elif isolate_res_genes[gene].Model_type == "pilus meta-model":
                    pm = "PMM"
                elif isolate_res_genes[gene].Model_type == "protein knockout model":
                    pm = "PKM"
                elif isolate_res_genes[gene].Model_type == "nonfunctional insertion model":
                    pm = "NFIM"
                elif isolate_res_genes[gene].Model_type == "gene cluster meta-model":
                    pm = "GCMM"
                elif isolate_res_genes[gene].Model_type == "efflux pump system meta-model":
                    pm = "EPSMM"
                elif isolate_res_genes[gene].Model_type == "protein domain meta-model":
                    pm = "PDMM"

                pres_abs_ident_dict[isolate_id].append(isolate_res_genes[gene].Best_Identities + "_" + isolate_res_genes[gene].Percentage_Length_of_Reference_Sequence + "_" + pm)

            elif isolate_res_genes[gene].Cut_Off == "Strict" and isolate_res_genes[gene].Nudged == "True":

                if isolate_res_genes[gene].Model_type == "protein homolog model":
                    pm = "PHM"
                elif isolate_res_genes[gene].Model_type == "protein variant model":
                    pm = "PVM"
                elif isolate_res_genes[gene].Model_type == "protein overexpression model":
                    pm = "POM"
                elif isolate_res_genes[gene].Model_type == "toxin system meta model":
                    pm = "TSMM"
                elif isolate_res_genes[gene].Model_type == "rRNA gene variant model":
                    pm = "RGVM"
                elif isolate_res_genes[gene].Model_type == "pilus meta-model":
                    pm = "PMM"
                elif isolate_res_genes[gene].Model_type == "protein knockout model":
                    pm = "PKM"
                elif isolate_res_genes[gene].Model_type == "nonfunctional insertion model":
                    pm = "NFIM"
                elif isolate_res_genes[gene].Model_type == "gene cluster meta-model":
                    pm = "GCMM"
                elif isolate_res_genes[gene].Model_type == "efflux pump system meta-model":
                    pm = "EPSMM"
                elif isolate_res_genes[gene].Model_type == "protein domain meta-model":
                    pm = "PDMM"

                pres_abs_ident_dict[isolate_id].append(isolate_res_genes[gene].Best_Identities + "_" + isolate_res_genes[gene].Percentage_Length_of_Reference_Sequence + "_" + pm + "**")

        elif gene not in isolate_res_genes:
            pres_abs_ident_dict[isolate_id].append("0") # record '0' % identity if no hit

# Write output to file:
header = ['\t'.join(["sample"] + output_gene_order)]

to_write_identities = []
for i in pres_abs_ident_dict:
    output_line = '\t'.join([i] + pres_abs_ident_dict[i])
    to_write_identities.append(output_line)
to_write_identities_ = '\n'.join(header + to_write_identities)

# to_write_lengths = []
# for i in pres_abs_len_dict:
#     output_line = '\t'.join([i] + pres_abs_len_dict[i])
#     to_write_lengths.append(output_line)
# to_write_lengths_ = '\n'.join(header + to_write_lengths)

output_identities_file = args.output_prefix + "_pres_abs_by_identity.txt"
# output_lengths_file = args.output_prefix + "_pres_abs_by_lengths.txt"

with open(output_identities_file, 'w') as outfile1:
    outfile1.write(to_write_identities_)

# with open(output_lengths_file, 'w') as outfile2:
#     outfile2.write(to_write_lengths_)
