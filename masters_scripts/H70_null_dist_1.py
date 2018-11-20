from sys import argv
import itertools

script, maxtest, mintest, phenoorder, nulldistrout = argv

print "Reading file: ", maxtest.strip().split('.')[0]

maxlist = []
minlist = []
phenoorderlist = []

with open(maxtest, 'r') as infile1:
	for line in infile1:
		maxlist.append(line.strip())

with open(mintest, 'r') as infile2:
	for line in infile2:
		minlist.append(line.strip())

with open(phenoorder, 'r') as infile3:
	for line in infile3:
		phenoorderlist.append(line.strip().split(".")[1])

for pheno in phenoorderlist:
	if pheno == 'PHENO1':
		del maxlist[phenoorderlist.index(pheno)]
		del minlist[phenoorderlist.index(pheno)]

absmaxlist = []
for max_, min_ in itertools.izip(maxlist, minlist):
	abslist = [abs(float(min_)), float(max_)]
	absmaxval = max(abslist)
	absmaxlist.append(absmaxval)

absmaxlist = map(str, absmaxlist)

with open(nulldistrout, 'a') as outfile1:

	outfile1.write('\n'.join(absmaxlist))
	outfile1.write('\n')

