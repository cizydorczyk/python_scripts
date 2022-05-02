import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("--input_file")
parser.add_argument("--output_file")

parser.add_argument("--list_mutation_types", action='store_true', default=False)
parser.add_argument("--analyse", action='store_true', default=False)

args = parser.parse_args()

# Define parser function:
def ParseMutationType(mutation_type, mutation_types_list):
    output_list = []
    num_mutation_types = len(mutation_types_list)

    for i in range(0, num_mutation_types):
        if mutation_type == mutation_types_list[i]:
            output_list.append("1")
        else:
            output_list.append("0")

    return(output_list)

# Control flow:
if not args.list_mutation_types and not args.analyse:
    print("Please specify either --list_mutation_types or --analyse. Exiting.")
    sys.exit()

if args.list_mutation_types and args.analyse:
    print("Please specify EITHER --list_mutation_types or --analyse, not both. Exiting.")
    sys.exit()

# If want to get only mutation types without analysis:
if args.list_mutation_types and not args.analyse:
    mutation_types = []
    with open(args.input_file, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            mutation_type = "\t" + line_elements[5].split("|")[1]
            if mutation_type not in mutation_types:
                mutation_types.append(mutation_type)
    print("\nMutation types:\n\n" + "\n".join(mutation_types) + "\n\nPlease edit line 49 in the script to include all of these mutation types and re-run the script with the --analyse option.")

# If analyzing data:
if not args.list_mutation_types and args.analyse:
    #mutation_types_order = 'synonymous_variant', 'stop_gained', 'missense_variant', 'intergenic_region', 'stop_lost&splice_region_variant'
    mutation_types_list = ["synonymous_variant", "stop_gained", "missense_variant", "intergenic_region", "stop_lost&splice_region_variant", "start_lost", "splice_region_variant&stop_retained_variant"]
    print("\nOrder of mutation types in output file is: " + ", ".join(mutation_types_list) + "\n")
    output_lines = []
    with open(args.input_file, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            mutation_type = line_elements[5].split("|")[1]
            parsed_mutation_types_list = ParseMutationType(mutation_type, mutation_types_list)
            output_line = line_elements + parsed_mutation_types_list
            output_lines.append("\t".join(output_line))
        
    with open(args.output_file, "w") as outfile1:
        outfile1.write("\n".join(output_lines))
