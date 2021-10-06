import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_blast_file", help="outfmt 7 required")

args = parser.parse_args()

coordinates = []
with open(args.input_blast_file, "r") as infile1:
    for line in infile1:
        if not line.startswith("#"):
            line_elements = line.strip().split("\t")
            annotation_element = ".".join(line_elements[1].split(".")[-2:])
            gene = args.input_blast_file.split("/")[-1].split(".")[0]
            print(gene, annotation_element)
