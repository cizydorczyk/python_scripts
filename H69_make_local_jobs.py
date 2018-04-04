from sys import argv
import itertools

script, input_perm_files_list, output_jobs_dir = argv

inputdir = '/home/conrad/Data/primary_project_3/gwas/March_8_2018_up_to_date_gwas_input_files'

input_files = []
perm_nums = []
with open(input_perm_files_list, 'r') as infile1:
	for line in infile1:
		input_files.append(line.strip().split('/')[-1])
		perm_nums.append(line.strip().split('/')[-1].split('.')[0])

for permfile, num in itertools.izip(input_files, perm_nums):
	outfile = output_jobs_dir + str(num) + ".perm1.job"
	with open(outfile, 'w') as outfile1:

		outfile1.write('#!/bin/bash' + '\n' + '\n')
		
		outfile1.write("cd /home/conrad/Data/primary_project_3/H69_local_output/" + '\n' + '\n')
		outfile1.write('mkdir perm' + str(num) + '\n')
		outfile1.write('cd perm' + str(num) + '\n')

		outfile1.write('\n' + '~/Software/plink2 --bfile ' + inputdir + '/bed_bim_fam_files/snps.unique_seg.maf005.nopheno --pheno /home/conrad/Data/primary_project_3/gwas/H68_plink_permutations/' + permfile + ' --chr-set -1 --glm hide-covar firth-fallback --covar /home/conrad/Data/primary_project_3/gwas/H66_snps_firth_regression_w_PCs/plink_pca/snps.unique_seg.maf005.pheno.pca.eigenvec --covar-name PC1 PC2 PC4 PC5 PC6 PC7' + '\n')
		outfile1.write('\n' + "for i in *.hybrid; do awk 'BEGIN{a=0}{if ($12>0+a) a=$12} END{print a}' $i >> /home/conrad/Data/primary_project_3/gwas/H69_null_distr_files/" + str(num) + ".maxtest.txt; awk 'BEGIN{a=1000}{if ($12<0+a) a=$12} END{print a}' $i >> /home/conrad/Data/primary_project_3/gwas/H69_null_distr_files/" + str(num) + ".mintest.txt; done" + '\n')

		outfile1.write('\n' + 'cd ..' + '\n')
		# outfile1.write('\n' + 'rm -r perm' + str(num) + '\n')


#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -N 36_pae
#PBS -A zke-503-ab