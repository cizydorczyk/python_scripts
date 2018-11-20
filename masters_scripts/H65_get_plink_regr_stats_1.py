from sys import argv
import os.path
import itertools

script, input_log_file, output_file = argv

header = 'num_pcs' + '\t' + 'lambda' + '\t' + 'proportion_NA' + '\t' + 'warnings'
total_num_var = 33886.0

inlines = []
with open(input_log_file, 'r') as infile1:
	for line in infile1:
		inlines.append(line.strip())
		if '(based on median chisq)' in line:
			if len(line.strip().split(' ')[-1]) == 2:
				lambda_est = '1'
				
			else:
				lambda_est = line.strip().split(' ')[-1][:-2]
				
		elif '--adjust values (' in line:
			num_var = line.strip().split(' ')[2][1:]
			
			proportions = str((total_num_var - float(num_var))/total_num_var*100)
		elif 'covariates loaded from' in line:
			pc_num = line.strip().split(' ')[0]
		elif 'covariate loaded from' in line:
			pc_num = line.strip().split(' ')[0]

for i in inlines:
	if "Warning" in i:
		warning = i + '...'
		break
	else:
		warning = 'NA'

print pc_num, lambda_est, num_var, proportions, warning

if os.path.isfile(output_file) == True:
	with open(output_file, 'a') as outfile1:
		outfile1.write('\n' + pc_num + '\t' + lambda_est + '\t' + num_var + '\t' + proportions + '\t' + warning)

if os.path.isfile(output_file) != True:
	with open(output_file, 'w') as outfile1:
		outfile1.write(header + '\n')
		outfile1.write(pc_num + '\t' + lambda_est + '\t' + num_var + '\t' + proportions + '\t' + warning)

