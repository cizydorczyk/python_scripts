# Load FaST-LMM basic association test:
from fastlmm.association import single_snp
from pysnptools.snpreader import Ped
from pysnptools.snpreader import Pheno
from pysnptools.snpreader import wrap_plink_parser
import numpy as np
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fastlmm.util.util as flutil

script, inped_file, inpheno_file, results_dataframe, output_manhattan = argv

# Load snp data:
print "Loading variant data..."
ped_file = Ped(inped_file)
print "Loading phenotype data..."
pheno_fn = Pheno(inpheno_file)

# Run basic association test:
print "Running FaST-LMM single_snp test..."
results_df = single_snp(test_snps=ped_file, pheno=pheno_fn, leave_out_one_chrom=0, output_file_name=results_dataframe)

chromosome_starts = flutil.manhattan_plot(results_df.as_matrix(["Chr", "ChrPos", "PValue"]), pvalue_line=4.4e-7, xaxis_unit_bp=True)
plt.show()
# fig = plt.figure()
# fig.savefig(output_manhattan)