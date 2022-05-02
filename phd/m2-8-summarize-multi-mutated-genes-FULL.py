import argparse
import sys

#####################################################################
# Arguments
#####################################################################

parser = argparse.ArgumentParser()

# Blast parsing options:

## Control option:
parser.add_argument("--blast_parse", help="Run the blast parsing part of the script.", action='store_true', default=False)

## Input files:
parser.add_argument("--blast_raw_input", help="Blast output file. Custom format: -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore' ")

## Output files:
parser.add_argument("--blast_unique_pegs", help="Output file with non-multi-mutated pegs.")
parser.add_argument("--blast_groups_output", help="Output file with blast groups. These genes are multi-mutated ACROSS STs.")

## Group filtering options:
parser.add_argument("--pid", help="pid threshold; default = 70.0", default=70.0, type=float)
parser.add_argument("--perclen", help="perclen threshold; default = 70.0", default=70.0, type=float)



# Summarize VCFs options:

## Control option:
parser.add_argument("--summarize_vcfs", help="Run the summarize VCFs analysis part of the script.", action='store_true', default=False)

## Input files:
parser.add_argument("--summarize_vcfs_input_files", help="comma-separated list (no spaces) of all vcf files to be analyzed. Required.")

## Output files:
parser.add_argument("--summarize_vcfs_output_file", help="Summary file with all VCFs summarized together.")



# Get mutation types options:

## Control option:
parser.add_argument("--muttypes", help="Run the mutation types of the analysis. This includes assigning a mutation type (S, NS, etc.) to each mutation.",
action='store_true', default=False)

## Input files:
parser.add_argument("--muttypes_input_summary_file", help="Summary file of VCF output from all STs being analyzed produced by this script. Should include header line starting with '#'. Required.")
parser.add_argument("--muttypes_blast_groups", help="Blast groups file produced by this script. Required.")

## Output files:
parser.add_argument("--muttypes_output_summary_file", help="Output summary file with mutation types and groups appended.")


args = parser.parse_args()

# Control flow:
if sum([args.blast_parse, args.summarize_vcfs, args.muttypes]) > 1:
    print("Please pick only one analysis to run. Exiting.")
    sys.exit()
elif sum([args.blast_parse, args.summarize_vcfs, args.muttypes]) == 0:
    print("Please select one analysis to run from: --blast_parse, --summarize_vcfs, --muttypes. Exiting.")
    sys.exit()

## Blast parse analysis:
if args.blast_parse:
    if args.blast_raw_input is None:
        print("Please supply raw blast input file using --blast_raw_input. Exiting.")
        sys.exit()
    if args.blast_unique_pegs is None:
        print("Please supply file for unique (not belonging to a group) pegs using --blast_unique_pegs. Exiting.")
        sys.exit()
    if args.blast_groups_output is None:
        print("Please supply file for blast groups using --blast_groups_output. Exiting.")
        sys.exit()

elif args.summarize_vcfs:
    if args.summarize_vcfs_input_files is None:
        print("Please supply input VCFs to summarize using --summarize_vcfs_input_files. Exiting.")
        sys.exit()
    if args.summarize_vcfs_output_file is None:
        print("Please supply file for summarized VCFs output using --summarize_vcfs_output_file. Exiting.")
        sys.exit()

elif args.muttypes:
    if args.muttypes_input_summary_file is None:
        print("Please supply input summarized VCFs file using --muttypes_input_summary_file. This is created using the --summarize_vcfs analysis option in this script. Exiting.")
        sys.exit()
    if args.muttypes_blast_groups is None:
        print("Please supply input blast groups file using --muttypes_blast_groups. This is created using the --blast_parse analysis option in this script. Exiting.")
        sys.exit()
    if args.muttypes_output_summary_file is None:
        print("Please supply output summary file using --muttypes_output_summary_file. Exiting.")
        sys.exit()



#####################################################################
# Classes & functions:
#####################################################################
 
#  For parsing blast output:

## Blast hit class:
class BlastHit(object):
    def __init__(self, qaccver, saccver, pident, length, mismatches, gapopens, qstart, qend, qlen, sstart, send, slen, evalue, bitscore, perclen, record):
        self.qaccver = qaccver
        self.saccver = saccver
        self.pident = pident
        self.length = length
        self.mismatches = mismatches
        self.gapopens = gapopens
        self.qstart = qstart
        self.qend = qend
        self.qlen = qlen
        self.sstart = sstart
        self.send = send
        self.slen = slen
        self.evalue = evalue
        self.bitscore = bitscore
        self.perclen = perclen
        self.record = record

## Blast helper function:
def find_key(input_dict, value):
    return {k for k, v in input_dict.items() if value in v}

# For obtaining mutation types:

## Define parser function:
def ParseMutationType(mutation_type, mutation_types_list):
    output_list = []
    num_mutation_types = len(mutation_types_list)

    for i in range(0, num_mutation_types):
        if mutation_type == mutation_types_list[i]:
            output_list.append("1")
        else:
            output_list.append("0")

    return(output_list)

#####################################################################
# Main analysis functions:
#####################################################################

# Main parse blast analysis function:
def ParseBlast():
    groups_dict = {}
    unique_pegs_list = []
    groups_counter = 1
    with open(args.blast_raw_input, "r") as infile1:
        for line in infile1:

            le = line.strip().split("\t")
            perclen = round((int(le[3])/max(int(le[8]), int(le[11]))) * 100, 2)
            blast_hit = BlastHit(le[0], le[1], le[2], le[3], le[4], le[5], le[6], le[7], le[8], le[9], le[10], le[11], le[12], le[13], perclen, line.strip())

            if blast_hit.qaccver != blast_hit.saccver:
                if float(blast_hit.pident) > args.pid and blast_hit.perclen > args.perclen:
                    key = find_key(groups_dict, blast_hit.qaccver)
                    if len(key) == 0: # if fail to find query in any keys
                        # print("1")
                        key = find_key(groups_dict, blast_hit.saccver) # try finding subject
                        if len(key) == 0: # if subject not found, start new group
                            # print("2")
                            groups_dict[groups_counter] = [blast_hit.qaccver, blast_hit.saccver]
                            groups_counter += 1
                        elif len(key) == 1: # if subject found, append query to dict entry
                            # print("3")
                            groups_dict[list(key)[0]].append(blast_hit.qaccver)
                        elif len(key) > 1: # if multiple keys with subject found, print error and break
                            print("4")
                            print(key, "Something is weird...multiple groups with value found...")
                            break
                    elif len(key) == 1: # if query found
                        # print("5")
                        key2 = find_key(groups_dict, blast_hit.saccver) # test if subject found
                        if len(key2) == 0: # if subject not found
                            # print("6")
                            groups_dict[list(key)[0]].append(blast_hit.saccver) # append subject to query dict entry
                        elif len(key2) == 1: # if subject found
                            # print("7")
                            if list(key)[0] == list(key2)[0]: # if query key == subject key
                                # print("8")
                                continue # do nothing; repeat entry for pair
                            elif key != key2: # if query != subject key; may happen if order of hits is weird...
                            # aka. 1-2, 3-4, 3-1; 1/2 grouped, 3/4 grouped, then 3/1 throw this error...not sure what sort order blast uses but this can technically occur...
                                print("9" + "\t" + "Multiple keys for pair...combining groups...\n\tCHECK THIS ENTRY MANUALLY\n\t" + blast_hit.qaccver + "\t" + blast_hit.saccver + "\n\tGroups " + str(key) + "\t" + str(key2))
                                groups_dict[list(key)[0]] = groups_dict[list(key)[0]] + groups_dict[list(key2)[0]]
                                del groups_dict[list(key2)[0]]
                        elif len(key2) > 1:
                            print("10")
                            print("Something is weird...multiple groups with subject found...")
                            break


            elif blast_hit.qaccver == blast_hit.saccver:
                unique_pegs_list.append(blast_hit.qaccver.split("|")[1])

    ## Write groups file:
    to_write = []
    for group in groups_dict:
        output_string = str(group) + "\t" + "\t".join(groups_dict[group])
        to_write.append(output_string)
    with open(args.blast_groups_output, "w") as outfile1:
        outfile1.write("\n".join(to_write))

    ## Add one representative peg from each group to unique pegs list:
    for group in groups_dict:
        unique_pegs_list.append(groups_dict[group][0].strip().split("|")[1])

    ## Write unique pegs to file:
    with open(args.blast_unique_pegs, "w") as outfile2:
        outfile2.write("\n".join(unique_pegs_list))

# Main summarize VCFs analysis:

def SummarizeVcfs():
    vcfs_list = args.summarize_vcfs_input_files.split(",")
    vcf_lines = ["ST\tCHROM\tPOS\tREF\tALT\tINFO\tPEG\tNAME"]
    for vcf_file in vcfs_list:
        st = vcf_file.strip().split("/")[-1].split(".")[0]
        with open(vcf_file, "r") as infile2:
            for line in infile2:
                if not line.startswith("#"):
                    line_elements = line.strip().split("\t")
                    new_line = st + "\t" + line_elements[0] + "\t" + line_elements[1] + "\t" + line_elements[3] + "\t" + line_elements[4] + "\t" + line_elements[7] + "\t" + line_elements[-2] + "\t" + line_elements[-1]
                    vcf_lines.append(new_line)

    with open(args.summarize_vcfs_output_file, "w") as outfile1:
        outfile1.write("\n".join(vcf_lines))

# Main get mutation types analysis:

def GetMutationTypes():
    ## Read in blast groups:
    blast_groups = {}
    with open(args.muttypes_blast_groups, "r") as infile3:
        for line in infile3:
            line_elements = line.strip().split("\t")
            pegs = []
            for i in line_elements[1:]:
                pegs.append(i.strip().split("|")[1])
            blast_groups[line_elements[0]] = pegs

    ## Mutation types/blast groups analysis:
    output_lines = []
    mutation_types_list = []
    with open(args.muttypes_input_summary_file, "r") as infile4:
        for line in infile4:
            if not line.startswith("ST"):
                line_elements = line.strip().split("\t")
                mutation_type = line_elements[5].split("|")[1]
                if mutation_type not in mutation_types_list:
                    mutation_types_list.append(mutation_type)
            elif line.startswith("ST"):
                header = line.strip()
                #output_lines.append(header)
        header = header + "\t" + "\t".join(mutation_types_list) + "\tBlast_Group"
        output_lines.append(header)

        infile4.seek(0)
        for line in infile4:
            if not line.startswith("ST"):
                line_elements = line.strip().split("\t")
                mutation_type = line_elements[5].split("|")[1]

                # Get mutation type
                parsed_mutation_types_list = ParseMutationType(mutation_type, mutation_types_list)
                
                # Get blast group (if any)
                peg = line_elements[6]
                blast_group = ""
                for key, value in blast_groups.items():
                    if peg in value:
                        blast_group = key
                    # else:
                    #     blast_group = ""
                
                # Concatenate output
                output_line = line_elements + parsed_mutation_types_list + [blast_group]
                output_lines.append("\t".join(output_line))

    # Write output  
    with open(args.muttypes_output_summary_file, "w") as outfile1:
        outfile1.write("\n".join(output_lines))

#####################################################################
# Run main analyses:
#####################################################################

if args.blast_parse:
    print("Parsing blast results.\n\tOutput blast groups file is: " + args.blast_groups_output + "\n\tOutput unique pegs file is: " + args.blast_unique_pegs)
    ParseBlast()

if args.summarize_vcfs:
    print("Summarizing VCFs.\n\tOutput summary file is: " + args.summarize_vcfs_output_file)
    SummarizeVcfs()

if args.muttypes:
    print("Parsing mutation types and blast groups.\n\tOutput summary file is: " + args.muttypes_output_summary_file + "\nPLEASE NOTE that the script does not assign groups to intergenic SNPs, even if the up/downstream gene is in a group. This must be checked manually. Intergenic-only groups are not recorded. Groups that consist of only one SNP are those for which the remaining SNPs are intergenic. Groups with multiple SNPs may or may not have additional intergenic SNPs.")
    GetMutationTypes()