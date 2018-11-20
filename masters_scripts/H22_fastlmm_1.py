# Load FaST-LMM basic association test:
from fastlmm.association import single_snp
from pysnptools.snpreader import Ped
from pysnptools.snpreader import Pheno
from pysnptools.snpreader import wrap_plink_parser
import numpy as np
from sys import argv

script, inped_file, inpheno_file, results_dataframe = argv

# Load snp data:
print "Loading variant data..."
ped_file = Ped(inped_file)
print "Loading phenotype data..."
pheno_fn = Pheno(inpheno_file)

# Run basic association test:
print "Running FaST-LMM single_snp test..."
results_df = single_snp(test_snps=ped_file, pheno=pheno_fn, leave_out_one_chrom=0, output_file_name=results_dataframe) 