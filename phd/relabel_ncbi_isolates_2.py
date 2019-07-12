isolates = []

with open("/home/conrad/ncbi_ecoli_genomes/ncbi-genomes-2019-06-07/in_silico_st131_typing/ncbi_st131_c1_isolates.txt", 'r') as infile1:
    for line in infile1:
        accession = line.strip().split("_")[0] + "_" + line.strip().split("_")[1]
        isolates.append(accession)

for i in isolates:
    print(i)

# to_write = '\n'.join(isolates)

# with open("/home/conrad/ncbi_ecoli_genomes/mlst/ncbi_st1193_isolates.txt", 'w') as outfile1:
#     outfile1.write(to_write)

#/home/conrad/ncbi_ecoli_genomes/ncbi-genomes-2019-06-07/fasta_files/GCA_000010485.1_ASM1048v1_genomic.fna.gz
