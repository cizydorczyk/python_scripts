from sys import argv

script, input_file, output_file = argv

input_file_lines = []
with open(input_file, 'r') as infile1:
	for line in infile1:
		input_file_lines.append(line.strip())

lines_with_na = []
for line in input_file_lines:
	if "NA" in line:
		lines_with_na.append(line)

with open(output_file, 'w') as outfile1:
	to_write = '\n'.join(lines_with_na)
	header = input_file_lines[0]
	outfile1.write(header + '\n')
	outfile1.write(to_write)