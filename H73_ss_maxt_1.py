from sys import argv
import itertools

script, nulldist, pval, output_file = argv

nulldist_list = []
with open(nulldist, 'r') as infile1:
	for line in infile1:
		nulldist_list.append(float(line.strip()))
nulldist_list.sort(reverse=True)

pval_list = []
snp_list = []
with open(pval, 'r') as infile2:
	for line in infile2:
		if not line.startswith('#'):
			snp_list.append(line.strip().split('\t')[1])
			if line.strip().split('\t')[-1] != 'NA':
				pval_list.append(float(line.strip().split('\t')[-1]))
			else:
				pval_list.append('NA')

ssmaxt_pval = []

for i in pval_list:
	if i != 'NA':
		count = 1
		for j in nulldist_list:
			if j > i:
				count += 1
			else:
				break
		denom = len(nulldist_list) + 1
		cor_pval = count / denom
		ssmaxt_pval.append(cor_pval)

with open(output_file, 'w') as outfile1:
	to_write = []
	for snp, corpval in itertools.izip(snp_list, ssmaxt_pval):
		to_write.append(str(snp) + '\t' + str(corpval))
	outfile1.write('\n'.join(to_write))

# for i in pval_list:
# 	if i != 'NA':
# 		numerator = sum(j > i for j in nulldist_list) + 1
# 		denominator = len(nulldist_list) + 1
# 		cor_pval = float(numerator) / float(denominator)
# 		ssmaxt_pval.append(cor_pval)

# print len(ssmaxt_pval)

# for i in pval_list:
# 	if i != 'NA':
# 		count = 1
# 		for j in nulldist_list:
# 			if j > i:
# 				count += 1
# 			else:
# 				continue
# 		denom = len(nulldist_list) + 1
# 		cor_pval = count / denom
# 		ssmaxt_pval.append(cor_pval)

# print len(ssmaxt_pval)