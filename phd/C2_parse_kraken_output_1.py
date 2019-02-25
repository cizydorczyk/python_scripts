import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input_kraken_output", help="file to parse")
parser.add_argument("parsed_kraken_file", help="parsed file ready to be input to krona")

args = parser.parse_args()

output_lines = []
with open(args.input_kraken_output, 'r') as infile:
    for line in infile:
        line = line.strip().split('\t')
        query_id = line[1]
        tax_id = line[2].split('(')[-1][:-1].split(' ')[1]
        output_lines.append(query_id + '\t' + tax_id)

with open(args.parsed_kraken_file, 'w') as outfile:
    to_write = '\n'.join(output_lines)
    outfile.write(to_write)
