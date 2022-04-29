import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_spa_typing_file", help="file with stdout from spaTyper (isolate must be echoed to file before spaTyper output!)")
parser.add_argument("--parsed_spa_typing_output", help="output file with parsed output")

args = parser.parse_args()

isolate_list = []
type_list = []
with open(args.input_spa_typing_file, "r") as infile1:
    for line in infile1:
        if line.startswith("A") or line.startswith("SA"):
            isolate = line.strip()
            isolate_list.append(isolate)
        if "length" in line:
            type_line = line.strip()
            type_list.append(type_line)

to_write = []
for i, j in zip(isolate_list, type_list):
    to_write_line = i + "\t" + j
    to_write.append(to_write_line)

with open(args.parsed_spa_typing_output, "w") as outfile1:
    outfile1.write("\n".join(to_write))