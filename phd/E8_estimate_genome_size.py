import argparse
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser()

parser.add_argument("--R1", help="input R1 fastq file")
parser.add_argument("--R2", help="input R2 fastq file")
parser.add_argument("--fastq_dir", help='dir where subsampled fastq files will be output')
parser.add_argument("--isolate", help='isolate number/id')

args = parser.parse_args()

print("\nWorking on isolate " + args.isolate + "...")

R1_file = args.R1.strip().split("/")[-1]
R2_file = args.R2.strip().split("/")[-1]

def GetBasesFromFq(inputfq):
    base_count = 0
    with open(inputfq, 'r') as infile1:
        for title, seq, qual in FastqGeneralIterator(infile1):
            base_count += len(seq)
    return base_count

print("\nGetting total bases from fastq files...")
R1_bases = GetBasesFromFq(args.R1)
R2_bases = GetBasesFromFq(args.R2)
total_bases = R1_bases + R2_bases

print("\nEstimating genome size with 'mash'...")

mash_cmd = ["mash", "sketch", "-o", "temp.R1.mash.sketch", "-k", "32", "-m", "3", "-s", "1000", "-r", args.R1]
mash_cmd_raw_output = subprocess.run(mash_cmd, stderr=subprocess.PIPE)
subprocess.run(["rm", "temp.R1.mash.sketch.msh"])

mash_cmd_output = mash_cmd_raw_output.stderr.decode('utf-8')
mash_gsize_est = float(mash_cmd_output.strip().split('\n')[0].split(' ')[-1])
print("Estimated genome size is ~" + str(mash_gsize_est))
# mash_cov_est = mash_cmd_output.strip().split('\n')[1].split(' ')[-1]

print('\nEstimating sequencing depth...')
orig_depth = total_bases/mash_gsize_est

print("Original sequencing depth is ~" + str(round(orig_depth,3)))

if orig_depth > 110:
    depth_factor = round(100/orig_depth, 3)
    print("\nSubsampling factor to go from ~" + str(round(orig_depth,3)) + "x depth to ~100x depth is ~" + str(depth_factor))

    print("Subsampling reads (R1 & R2) by factor of " + str(depth_factor) + "..")

    p1 = subprocess.Popen("seqtk sample -s 42 " + args.R1 + " " + str(depth_factor) + " > " + args.fastq_dir + args.isolate + "_1.sub.fastq", shell=True)
    p2 = subprocess.Popen("seqtk sample -s 42 " + args.R2 + " " + str(depth_factor) + " > " + args.fastq_dir + args.isolate + "_2.sub.fastq", shell=True)

    print("Done.")

else:
    print("\nNo subsampling necessary for isolate " + args.isolate)
    print("Copying reads to subsampled directory and adding '.sub.' between _1/_2 and fastq...")

    subprocess.run(["cp", args.R1, args.fastq_dir])
    subprocess.run(["mv", args.fastq_dir + R1_file, args.fastq_dir + args.isolate + "_1.sub.fastq"])

    subprocess.run(["cp", args.R2, args.fastq_dir])
    subprocess.run(["mv", args.fastq_dir + R2_file, args.fastq_dir + args.isolate + "_2.sub.fastq"])

    print("Done.")
