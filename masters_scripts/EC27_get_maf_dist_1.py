from sys import argv
import numpy as np
import itertools

script, inputallelicdepth, outputfile = argv
print "Working on file: " + inputallelicdepth.split('/')[-1]

with open(inputallelicdepth, 'r') as infile1:
    lines = infile1.read().splitlines()
    del lines[0]

all_maf = []
major_maf = []

for i in lines:
	line = i.strip().split('\t')
	ad = [float(j) for j in line[-1].split(',')]
	adsum = sum(ad)
	ad.remove(max(ad))

	numbases = len(ad[0:-1])

	if adsum > 0.0:
		for k in ad[0:-1]:
				all_maf.append(round((k/adsum),2))

		if max(ad) > 0.0:
			major_maf.append(round((max(ad)/adsum),2))

# Count total proportions:
all_maf_dict = {}
major_maf_dict = {}

proportions = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50]

for i in proportions:
	all_maf_dict[str(i)] = all_maf.count(i)
	major_maf_dict[str(i)] = major_maf.count(i)

with open(outputfile, 'w') as outfile1:
	outfile1.write('maf\ttotal_count\tmajor_count\n')
	for keyt, key2 in itertools.izip(sorted(all_maf_dict.keys()), sorted(major_maf_dict.keys())):
		outfile1.write(str(keyt) + '\t' + str(all_maf_dict[keyt]) + '\t' + str(major_maf_dict[key2]) + '\n')


