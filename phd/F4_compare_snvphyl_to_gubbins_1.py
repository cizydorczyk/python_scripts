import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--gubbins_vcf", help="gubbins vcf file")
parser.add_argument("--snvtable", help="snvphyl snv table file")
parser.add_argument("--no_density_filter_snvtable", help="snv table from snvphyl\
                    calculated with no density filter in place (used to get true\
                    negatives)")
parser.add_argument("--window_size", help="window size used in comparison")
parser.add_argument("--density_threshold", help="density threshold used in snvphyl")
parser.add_argument("--output_file", help="output file to contain calculations\
                    can be appended to)")

args = parser.parse_args()

# Get gubbins vcf positions:
gubbins_positions = []

with open(args.gubbins_vcf, 'r') as infile1:
    for line in infile1:
        if line.startswith("#"):
            continue
        else:
            line_elements = line.strip().split('\t')
            ref = line_elements[3]
            alt = line_elements[4]
            position = line_elements[1]

            if 'N' not in ref and 'N' not in alt:
                gubbins_positions.append(position)
            else:
                continue

# Convert gubbins positions to sets for easy comparisons
gubbins_pos_set = set(gubbins_positions)

# Get SNVPhyl positions
snvphyl_positions = []
with open(args.snvtable, 'r') as infile2:
    for line in infile2:
        if line.startswith("#"):
            continue
        else:
            line_elements = line.strip().split('\t')
            position = line_elements[1]

            if 'filtered' in line_elements[2]:
                continue
            else:
                snvphyl_positions.append(position)
# Convert snvphyl positions to set for easy comparisons
snvphyl_pos_set = set(snvphyl_positions)

# Get true negative positions:
true_negatives = []

with open(args.no_density_filter_snvtable, 'r') as infile3:
    for line in infile3:
        if line.startswith("#"):
            continue
        else:
            line_elements = line.strip().split("\t")
            position = line_elements[1]

            if 'filtered' in line_elements[2]:
                continue
            else: # true negatives are positions in the unfiltered snp dataset
                  # that are absent from the gubbins SNP positions and snvphyl
                  # SNP positions

                if position in gubbins_pos_set:
                    continue
                elif position in snvphyl_pos_set:
                    continue
                else:
                    true_negatives.append(position)

# Get True Positives:
tp = gubbins_pos_set & snvphyl_pos_set
fp = snvphyl_pos_set - gubbins_pos_set
fn = gubbins_pos_set - snvphyl_pos_set
tn = set(true_negatives)

# Sensitivity = TP/(TP+FN)
# = how many SNPs that Gubbins recovered were recovered by snvphyl
sensitivity = len(tp)/(len(tp) + len(fn))
print("Sensitivity (tp, fp, fn, tn): ", sensitivity, str(len(tp)), str(len(fp)),\
    str(len(fn)), str(len(tn)))

# Specificity = TN/(TN+FP)
# = how many SNPs that Gubbins excluded were excluded by snvphyl
specificity = len(tn)/(len(tn)+len(fp))
print("Specificity: ", specificity)

header = "Comparison\tWindow_Size\tDensity_Threshold\tTP\tFP\tFN\tTN\tSensitivity\tSpecificity"
comparison_name = args.density_threshold + " SNPs in " + args.window_size + " bp window"
output_line = comparison_name + "\t" + args.window_size + "\t" + args.density_threshold + "\t" +\
            str(len(tp)) + "\t" + str(len(fp)) + "\t" + str(len(fn)) + "\t" + str(len(tn)) +\
            "\t" + str(sensitivity) + "\t" + str(specificity)

if not os.path.isfile(args.output_file):
    with open(args.output_file, 'w') as outfile1:
        output = header + "\n" + output_line
        outfile1.write(output)
elif os.path.isfile(args.output_file):
    with open(args.output_file, 'a') as outfile1:
        outfile1.write("\n" + output_line)
