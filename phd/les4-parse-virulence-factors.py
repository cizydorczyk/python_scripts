import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_vf", help="input samtools coverage file for VFs")
parser.add_argument("--output_vf", help="file with output")

args = parser.parse_args()

vf_list = []

with open(args.input_vf, "r") as infile1:
    for line in infile1:
        if not line.startswith("#"):
            line_elements = line.strip().split("\t")
            cov = round(float(line_elements[5]), 2)
            mean_depth = round(float(line_elements[6]), 2)

            vf_list.append(str(cov) + "_" + str(mean_depth))

to_write = "\n".join(vf_list)
with open(args.output_vf, "w") as outfile1:
    outfile1.write(to_write)
