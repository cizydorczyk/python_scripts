import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--isolate_list")
parser.add_argument("--sample_list", help="this one contains dates")
parser.add_argument("--output_commands")

args = parser.parse_args()

with open(args.isolate_list, "r") as infile1:
    isolate_list = [line.strip() for line in infile1]

with open(args.sample_list, "r") as infile2:
    sample_list = [line.strip() for line in infile2]

to_write = []
for isolate, sample in zip(isolate_list, sample_list):
    mv_r1_cmd = "mv /home/conrad/saureus/fastq-files/sa1-trimmed-fastq/" + isolate + "_1.fastq.gz /home/conrad/saureus/fastq-files/sa1-trimmed-fastq/" + sample + "_1.fastq.gz"
    mv_r2_cmd = "mv /home/conrad/saureus/fastq-files/sa1-trimmed-fastq/" + isolate + "_2.fastq.gz /home/conrad/saureus/fastq-files/sa1-trimmed-fastq/" + sample + "_2.fastq.gz"

    to_write.append(mv_r1_cmd + "\n" + mv_r2_cmd)

with open(args.output_commands, "w") as outfile1:
    outfile1.write("\n".join(to_write))