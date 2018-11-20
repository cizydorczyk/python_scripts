from sys import argv
import itertools

script, input_perm_files_list, output_jobs_dir = argv

scinetdir = '/scratch/d/dguttman/cizydor/plink_perms/'

input_files = []
perm_nums = []
with open(input_perm_files_list, 'r') as infile1:
	for line in infile1:
		input_files.append(line.strip().split('/')[-1])
		perm_nums.append(line.strip().split('/')[-1].split('.')[0])

for permfile, num in itertools.izip(input_files, perm_nums):
	outfile = output_jobs_dir + str(num) + ".perm1.job"
	with open(outfile, 'w') as outfile1:

		outfile1.write('#!/bin/bash' + '\n' + '#PBS -l nodes=1:ppn=8' + '\n' + '#PBS -l walltime=1:00:00' + '\n' + '#PBS -N ' + str(num) + '_1kperm_job' + '\n' + '#PBS -A zke-503-ab' + '\n' + '\n')
		outfile1.write('module load plink2/plink20' + '\n' + '\n')
		
		outfile1.write("cd /scratch/d/dguttman/cizydor/plink_perms/raw_output_files/" + '\n')
		outfile1.write('mkdir perm' + str(num) + '\n')
		outfile1.write('cd perm' + str(num) + '\n')

		outfile1.write('\n' + 'plink2 --bfile ' + scinetdir + 'input_files/snps.unique_seg.maf005.pheno --pheno ' + scinetdir + 'H68_plink_permutations/' + permfile + ' --chr-set -1 --glm hide-covar firth-fallback --covar ' + scinetdir + 'input_files/snps.unique_seg.maf005.pheno.pca.eigenvec --covar-name PC1 PC2 PC4 PC5 PC6 PC7' + '\n')
		outfile1.write('\n' + "for i in *.hybrid; do awk 'BEGIN{a=0}{if ($12>0+a) a=$12} END{print a}' $i >> " + scinetdir + "null_distr_files/" + str(num) + ".maxtest.txt; awk 'BEGIN{a=1000}{if ($12<0+a) a=$12} END{print a}' $i >> " + scinetdir + "null_distr_files/" + str(num) + ".mintest.txt; echo $i >> " + scinetdir + "null_distr_files/" + str(num) + ".pheno_order.txt; done" + '\n')

		outfile1.write('\n' + 'cd ..' + '\n')
		outfile1.write('\n' + 'rm -r perm' + str(num) + '\n')


#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -N 36_pae
#PBS -A zke-503-ab