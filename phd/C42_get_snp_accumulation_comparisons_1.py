import argparse
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--input_vcf", help="list of vcf files")
parser.add_argument("--isolates", nargs='*', help="isolates to analyze")

args = parser.parse_args()
