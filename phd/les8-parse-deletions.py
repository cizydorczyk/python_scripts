import argparse
from operator import itemgetter
from itertools import groupby


parser = argparse.ArgumentParser()

parser.add_argument("--depth_file", help="samtools depth output file")
parser.add_argument("--deletions_file", help="output file")

args = parser.parse_args()

deletions = []

with open(args.depth_file, "r") as infile1:
    prev_line_depth = 1
    for line in infile1:
        if not line.startswith("#"):
            pos = int(line.strip().split("\t")[1])
            depth = int(line.strip().split("\t")[2])

            if depth > 0:
                continue
            elif depth == 0:
                deletions.append(pos)

ranges = []
for k,g in groupby(enumerate(deletions),lambda x:x[0]-x[1]):
    group = (map(itemgetter(1),g))
    group = list(map(int,group))
    ranges.append((group[0],group[-1]))

to_write = []
for range in ranges:
    to_write.append(str(range[0]) + "\t" + str(range[1]) + "\t" + str(range[1]-range[0]))

with open(args.deletions_file, "w") as outfile1:
    outfile1.write("\n".join(to_write))
