import argparse
import itertools
import subprocess
import os

# import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--ref", help="reference genome")
parser.add_argument("--snippy_output_dir", help="directory with snippy output directories")
parser.add_argument("--isolate_list_file", help="file with isolates - one per line")
parser.add_argument("--wd", help="working directory (where snippy will output files)")

args = parser.parse_args()

# set wd:
os.chdir(args.wd)

# read isolates file:
isolates = []
with open(args.isolate_list_file, 'r') as infile1:
    for line in infile1:
        isolates.append(line.strip())

# create a list of all possible unordered pairs of isolates:
isolate_pairs = list(itertools.combinations(isolates, 2))

#snippy_core_cmd = 'cd ' + os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling") + ' && snippy-core --ref ' + args.reference + ' ' + snp_calling_directory_list
#subprocess.run(snippy_core_cmd, shell=True)

for pair in isolate_pairs:

    isolate_1 = pair[0]
    isolate_2 = pair[1]

    snippy_core_cmd = "snippy-core --prefix " + isolate_1 + "_" + isolate_2 + " --ref " + args.ref + " " + args.snippy_output_dir + isolate_1 + "/ " + args.snippy_output_dir + isolate_2 + "/"
    subprocess.run(snippy_core_cmd, shell = True)

    remove_unnec_files_cmd = "rm " + isolate_1  + "_" + isolate_2 + ".vcf " + isolate_1  + "_" + isolate_2 + ".txt " + isolate_1  + "_" + isolate_2 + ".tab " + isolate_1  + "_" + isolate_2 + ".ref.fa " + isolate_1 + "_" + isolate_2 + ".aln"
    subprocess.run(remove_unnec_files_cmd, shell = True)
