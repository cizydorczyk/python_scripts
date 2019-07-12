file1 = "/home/conrad/ncbi_ecoli_genomes/clade_assemblies/all_clades_strains.txt"

isolate_files = []

with open(file1, 'r') as infile1:
    for line in infile1:
        isolate_files.append('_'.join(line.strip().split('_')[0:2]))

isolates_set = sorted(list(set(isolate_files)))

with open("/home/conrad/ncbi_ecoli_genomes/clade_assemblies/all_clades.txt", 'w') as outfile1:
    to_write = '\n'.join(isolates_set)
    outfile1.write(to_write)
