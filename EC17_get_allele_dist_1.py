from sys import argv
import numpy as np
import itertools

script, inputallelicdepth, outputfile = argv
print "Working on file: " + inputallelicdepth.split('/')[-1]


with open(inputallelicdepth, 'r') as infile1:
    lines = infile1.read().splitlines()
    del lines[0]

proportions_breakdown = {1:[], 2:[], 3:[], 4:[]}
proportions = []

for i in lines:
	line = i.strip().split('\t')
	ad = [float(j) for j in line[-1].split(',')]
	adsum = sum(ad)

	numbases = len(ad[0:-1])

	if adsum != 0.0:
		for k in ad[0:-1]:
			proportions_breakdown[numbases].append(round((k/adsum),2))
			proportions.append(round((k/adsum),2))
	
	elif adsum == 0.0:
	# 	proportions[numbases].append(0.00)
		continue

# Count total proportions:
proportions_dict = {}

for i in np.arange(0,1.01, 0.01):
	proportions_dict[str(i)] = proportions.count(i)

# Count proportions with 2, 3, and 4 bases separately:
proportions_2_dict = {}
proportions_3_dict = {}
proportions_4_dict = {}

for i in np.arange(0,1.01, 0.01):
	proportions_2_dict[str(i)] = proportions_breakdown[2].count(i)

for i in np.arange(0,1.01, 0.01):
	proportions_3_dict[str(i)] = proportions_breakdown[3].count(i)

for i in np.arange(0,1.01, 0.01):
	proportions_4_dict[str(i)] = proportions_breakdown[4].count(i)

with open(outputfile, 'w') as outfile1:
	outfile1.write('proportion\ttotal_count\tcount_2\tcount_3\tcount_4\n')
	for keyt, key2, key3, key4 in itertools.izip(sorted(proportions_dict.keys()), sorted(proportions_2_dict.keys()), sorted(proportions_3_dict.keys()), sorted(proportions_4_dict.keys())):
		outfile1.write(str(keyt) + '\t' + str(proportions_dict[keyt]) + '\t' + str(proportions_2_dict[key2]) + '\t' + str(proportions_3_dict[key3]) + '\t' + str(proportions_4_dict[key4]) + '\n')

	# for key, value in sorted(proportions_dict.iteritems()):
	# 	outfile1.write(str(key) + '\t' + str(value) + '\n')