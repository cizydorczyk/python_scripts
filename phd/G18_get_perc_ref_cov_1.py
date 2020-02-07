import argparse
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--flagstat_dir", help="dir with flagstat files")
parser.add_argument("--depth_dir", help="dir with depth files")
parser.add_argument("--out_prefix", help="output file prefix; should be a full path")
parser.add_argument("--isolate_list", help="isolate list")
parser.add_argument("--pairs", help="suspected pairs, tab separated")

args = parser.parse_args()

with open(args.isolate_list, 'r') as infile0:
    isolates = [line.strip() for line in infile0]

output_dict = {}
called_positions_dict = {}


for isolate in isolates:

    output_dict[isolate] = []
    called_positions_dict[isolate] = []

    flagstat_file = isolate + "_flagstat.txt"
    flagstat_file_path = os.path.join(args.flagstat_dir, flagstat_file)

    with open(flagstat_file_path, 'r') as infile1:
        for line in infile1:
            if "mapped (" in line:
                perc_mapped = line.strip().split(" ")[-3].strip("(").strip("%")

                output_dict[isolate].append(perc_mapped)

    depth_file = isolate + "_depths.txt"
    depth_file_path = os.path.join(args.depth_dir, depth_file)

    with open(depth_file_path, 'r') as infile2:
        total = 0
        base_calls = 0
        for line in infile2:
            total += 1
            line_elements = line.strip().split("\t")
            if int(line_elements[2]) > 0:
                base_calls += 1
                called_positions_dict[isolate].append(line_elements[1])

        proportion_covered = round(base_calls / total * 100, 2)

        output_dict[isolate].append(str(base_calls))
        output_dict[isolate].append(str(total))
        output_dict[isolate].append(str(proportion_covered))

# write per-isolate statistics to file:
to_write = []
header = "Isolate\tPerc_reads_mapped\tNum_pos_covered\tTotal_pos_in_ref\tPerc_ref_covered"
to_write.append(header)

for key in output_dict:
    to_write.append(key + "\t" + "\t".join(output_dict[key]))

output_individual_stats_file = os.path.join(args.out_prefix, "individual_mapping_stats.txt")
with open(output_individual_stats_file, 'w') as outfile1:
    outfile1.write("\n".join(to_write))

# Get pairwise stats:
pairs_list = []
with open(args.pairs, 'r') as infile3:
    for line in infile3:
        line_elements = line.strip().split("\t")
        pairs_list.append((line_elements[0], line_elements[1]))

output_dict2 = {}

for pair in pairs_list:
    a = pair[0]
    b = pair[1]

    dict_key = a + "\t" + b
    output_dict2[dict_key] = []

    a_positions = set(called_positions_dict[a]) # positions in a with base calls
    b_positions = set(called_positions_dict[b]) # positions in b with base calls

    a_length = len(a_positions) # number positions a covers in ref
    b_length = len(b_positions) # number positions b covers in ref
    a_perc = round(a_length/total*100, 2) # percentage of ref covered by positions in a
    b_perc = round(b_length/total*100, 2) # percentage of ref covered by positions in b

    output_dict2[dict_key].append(str(a_length))
    output_dict2[dict_key].append(str(a_perc))
    output_dict2[dict_key].append(str(b_length))
    output_dict2[dict_key].append(str(b_perc))

    a_and_b_of_ref = len(a_positions & b_positions)
    perc_a_and_b_of_ref = round(len(a_positions & b_positions) / total * 100, 2) # percent of reference covered by a and b

    output_dict2[dict_key].append(str(a_and_b_of_ref))
    output_dict2[dict_key].append(str(perc_a_and_b_of_ref))

    a_unique = len(a_positions - b_positions) # number of a's positions absent in b
    b_unique = len(b_positions - a_positions) # number of b's positions absent in a
    perc_a_unique = round(len(a_positions - b_positions) / a_length * 100, 2) # percent of a unique to a
    perc_b_unique = round(len(b_positions - a_positions) / b_length * 100, 2) # percent of b unique to b

    output_dict2[dict_key].append(str(a_unique))
    output_dict2[dict_key].append(str(perc_a_unique))
    output_dict2[dict_key].append(str(b_unique))
    output_dict2[dict_key].append(str(perc_b_unique))

    #print(a_length, b_length, a_perc, b_perc, a_and_b_of_ref, perc_a_and_b_of_ref, a_unique, b_unique, perc_a_unique, perc_b_unique)

# Write pairwise stats to file:
to_write = []
header = "Isolate_1\tIsolate_2\tI1_num_pos\tI1_perc_ref\tI2_num_pos\tI2_perc_ref\tShared_num_pos\tShared_perc_ref\tI1_unique_pos\tI1_unique_perc\tI2_unique_pos\tI2_unique_perc"
to_write.append(header)

for key in output_dict2:
    to_write.append(key + "\t" + "\t".join(output_dict2[key]))

output_pairwise_stats_file = os.path.join(args.out_prefix, "pairwise_mapping_stats.txt")

with open(output_pairwise_stats_file, 'w') as outfile2:
    outfile2.write("\n".join(to_write))
