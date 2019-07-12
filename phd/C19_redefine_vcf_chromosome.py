import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_vcf", help="file to parse")
parser.add_argument("--output_vcf", help="output vcf file with chromosome redefined")
parser.add_argument("--new_chrom", help='new chromosome name to add')

args = parser.parse_args()

to_write = []
with open(args.input_vcf, 'r') as infile1:
    for line in infile1:
        
        if line.startswith("#"):
            to_write.append(line.strip())

        elif not line.startswith("#"):
            line_elements = line.strip().split('\t')
            line_elements[0] = args.new_chrom

            to_write.append('\t'.join(line_elements))

with open(args.output_vcf, 'w') as outfile1:
    outfile1.write('\n'.join(to_write))
