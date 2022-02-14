import argparse
from Bio import SeqIO
from collections import OrderedDict

parser = argparse.ArgumentParser()

parser.add_argument("--vf_file", help="text file with virulence factors; 3 columns with 3rd containing locus tags")
parser.add_argument("--lesb58_gbk", help="lesb58 genbank file, NC_011770 accession")
parser.add_argument("--output_coords", help="output coords only, one set per line, dash separated")
parser.add_argument("--output_records", help="output coords + locus tags")

args = parser.parse_args()


# Read and parse gbk
new_locustag_dict = {}
old_locustag_dict = {}

gbk_record = SeqIO.read(args.lesb58_gbk, "genbank")

for feature in gbk_record.features:

    if feature.type == "gene":
        coord_start = str(feature.location).split("(")[0].strip("[]").split(":")[0]
        coord_end = str(feature.location).split("(")[0].strip("[]").split(":")[1]

        new_locus_tag = feature.qualifiers.get('locus_tag')[0]
        new_locustag_dict[new_locus_tag] = coord_start + "-" + coord_end

        old_locus_tag = feature.qualifiers.get('old_locus_tag')

        if old_locus_tag:
            old_locustag_dict[old_locus_tag[0]] = coord_start + "-" + coord_end
        else:
            continue

# Read and parse virulence file:
output_dict = OrderedDict()
feature_list = []
with open(args.vf_file, "r") as infile1:
    for line in infile1:
        line_elements = line.strip().split("\t")
        locus = line_elements[3]
        if locus in feature_list:
            print(locus)
        else:
            feature_list.append(locus)

        coords = old_locustag_dict.get(locus)

        if coords:
            # print(locus, coords)
            output_dict[locus] = coords
        else:
            coords = new_locustag_dict[locus]
            # print(locus, coords)
            output_dict[locus] = coords

print(len(output_dict))


# Write coords to file with and without locus tags (without to facilitate samtools coverage run):

with open(args.output_records, "w") as outfile1:
    to_write = []
    for record in output_dict:
        to_write.append(record + "\t" + output_dict[record])
    outfile1.write("\n".join(to_write))

with open(args.output_coords, "w") as outfile2:
    to_write = []
    for record in output_dict:
        to_write.append(output_dict[record])
    outfile2.write("\n".join(to_write))
