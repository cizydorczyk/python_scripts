import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--accessions_list", help="original names of fastq")
parser.add_argument("--strains_list", help="new names for isolates; MUST be sorted in corresponding (same) order as accessions list")
parser.add_argument("--arc_fastq_dir", help="fastq dir on arc (or wherever files reside")
parser.add_argument("--rename_commands_file", help="output text file with one mv command per line")

args = parser.parse_args()

# Read isolate list:
with open(args.accessions_list, "r") as infile1:
    accessions_list = [line.strip() for line in infile1]

# Read new isolate names list:
with open(args.strains_list, "r") as infile2:
    strains_list = [line.strip() for line in infile2]

# Create fastq rename commands:
mv_cmds = []
for old_isolate, new_isolate in zip(accessions_list, strains_list):
    mv_cmd_1 = "mv " + args.arc_fastq_dir + old_isolate + "_1.fastq.gz " + args.arc_fastq_dir + new_isolate + "_1.fastq.gz"

    mv_cmd_2 = "mv " + args.arc_fastq_dir + old_isolate + "_2.fastq.gz " + args.arc_fastq_dir + new_isolate + "_2.fastq.gz"

    mv_cmds.append(mv_cmd_1)
    mv_cmds.append(mv_cmd_2)

# Write to file:
with open(args.rename_commands_file, "w") as outfile1:
    outfile1.write("\n".join(mv_cmds))
