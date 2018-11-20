from sys import argv

script, inputfile, outputfile = argv

isolate = inputfile.strip().split('/')[-2]

with open(inputfile, 'r') as infile1:
	for line in infile1:
		if not line.startswith("ST"):
			data = line.strip()

with open(outputfile, 'a+') as outfile1:
	outfile1.write(isolate + '\t' + data + '\n')