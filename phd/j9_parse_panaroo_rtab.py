import argparse
import pandas as pd
import scipy.spatial

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", help="anvio summary file")
parser.add_argument("--output_file", help="output fasta")

args = parser.parse_args()

input_data = pd.read_csv(filepath_or_buffer = args.input_file, sep="\t", header=0, index_col=0)

input_df = input_data.to_dict(orient="list")

output_list = []
for key in input_df:
    header = ">" + key
    seq = "".join(str(i) for i in input_df[key])
    seq = seq.replace("1", "A")
    seq = seq.replace("0", "C")
    output_list.append(header + "\n" + seq)

with open(args.output_file, 'w') as outfile:
    outfile.write("\n".join(output_list))
