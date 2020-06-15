import argparse
import os
from operator import itemgetter
import pandas

parser = argparse.ArgumentParser()

parser.add_argument("--rgi_output_dir", help="directory with all .txt files from rgi")
parser.add_argument("--output_file", help="parsed output file with presence absence data")
parser.add_argument("--input_source", help="protein or nucleotide (default); type of input used for rgi", default=
                    "nucleotide")

args = parser.parse_args()

output_file_list = [os.path.join(args.rgi_output_dir, i) for i in os.listdir(args.rgi_output_dir) if i.endswith(".txt")]

resistance_elements = {}
isolate_records = {}
isolates_list = []

for ifile in output_file_list:
    isolate = ifile.split("/")[-1].split(".")[0].split("_")[0]
    isolates_list.append(isolate)
    isolate_records[isolate] = []

    with open(ifile, 'r') as infile1:
        for line in infile1:
            if not line.startswith("ORF_ID"):
                line_elements = line.strip().split("\t")
                if line_elements[8] not in resistance_elements:
                    resistance_elements[line_elements[8]] = []

                isolate_records[isolate].append(line)

            elif line.startswith("ORF_ID"):
                continue

for element in resistance_elements:
    for isolate in isolates_list:
        iso_record = ""
        for record in isolate_records[isolate]:
            record_elements = record.split("\t")
            if element == record_elements[8]:
                element_presence = "yes"
                if args.input_source == "nucleotide":
                    iso_record = itemgetter(1, 5, 9, 11, 12, 14, 15, 16, 20, 23)(record_elements)
                elif args.input_source == "protein":
                    iso_record = itemgetter(0, 5, 9, 11, 12, 14, 15, 16, 20, 23)(record_elements)
        if iso_record != "":
            resistance_elements[element].append(";".join(iso_record))
        elif iso_record == "":
            resistance_elements[element].append("")

res_element_pres_abs = pandas.DataFrame.from_dict(resistance_elements, orient="index", columns=isolates_list)

# Comment for output df:
comment_string = "# Presence format: contig_rgiID;cutoff;%_identity;model_type;snps;drug_class;resistance_mechanism;amr_gene_\
family;%_length_of_ref;nudged\n"

output_file = open(args.output_file, "a")
output_file.write(comment_string)
res_element_pres_abs.to_csv(path_or_buf=output_file, sep="\t")
output_file.close()
