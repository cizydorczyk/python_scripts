import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()

parser.add_argument("--chewbacca_table", help="table exported from chewbbaca")
parser.add_argument("--st_list", help="isolate list for specific ST")
parser.add_argument("--snp_aln", help="snp alignment")
parser.add_argument("--output_fasta", help="combined snp/pres-abs fasta")
parser.add_argument("--output_nexus", help="combined snp/pres-abs nexus")

args = parser.parse_args()

st_list = []
with open(args.st_list, "r") as infile1:
    for line in infile1:
        st_list.append(line.strip())

presData = pd.read_csv(filepath_or_buffer="/home/conrad/m-sm-notebook/mlst/m3-wgmlst/extract-cgmlst-alleles-threshold-0/cgMLST-edited.tsv",
    sep="\t",
    index_col=0,
    header=0)

presData.astype("int") # convert dataframe entries to integers

presData = presData.loc[st_list] # keep only rows corresponding to specified ST

presData[presData > 0] = 1 # replace all values >1 with 1 (presence)

presData = presData[presData.columns[presData.sum()>0]] # remove columns that sum to 0 (all absent)

presDict = presData.transpose().to_dict(orient="list")

# convert list value to string sequence:
for i in presDict:
    presDict[i] = "".join([str(int) for int in presDict[i]])

# Read in SNP alignment:
snpDict = SeqIO.to_dict(SeqIO.parse(args.snp_aln, "fasta"))

# Convert snpDict to match format of presDict and combine with presDict data:
to_write = []
for i in snpDict:
    len_snps = len(snpDict[i])
    len_pres = len(presDict[i])

    seq = str(str(snpDict[i].seq) + presDict[i])
    total_len = len(seq)

    to_write.append(">" + i.replace("-", "."))
    to_write.append(seq)

# Write fasta to file:
with open(args.output_fasta, "w") as outfile1:
    outfile1.write("\n".join(to_write))

##########
# Creates non-interleaved nexus...which MrBayes doesn't like.
# So use seqmagick convert to get interleaved nexus and manually adjust for partitioned data.
# Like so: seqmagick convert --alphabet dna-ambiguous test_combined_fasta.fasta test_combined_nexus.nex

##########

# Create nexus...
header = "#NEXUS\n\nBegin data;\n\tDimensions ntax=" + str(len(st_list)) + " nchar=" + str(total_len) +\
    ";\n\tFormat datatype=mixed(DNA:1-" + str(len_snps) + ",Restriction:" + str(1+len_snps) + "-" + str(total_len) +\
    ") interleave=no gap=- missing=?;\n\tMatrix\n"
print(header) # print for reference for manually generating partitioned nexus from seqmagick





#seq_body = []
#for i in snpDict:
#    seq = str(str(snpDict[i].seq) + presDict[i])
#    seq_body.append(i)
#    seq_body.append(seq)

#seq_body = "\n".join(seq_body) + "\n\t;\nend;"

# Write nexus:
#with open(args.output_nexus, "w") as outfile2:
#    outfile2.write(header + seq_body)
