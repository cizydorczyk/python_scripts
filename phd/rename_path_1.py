
lines1 = []
with open("/home/conrad/ncbi_ecoli_genomes/mlst/ncbi_st1193_genomes.txt", 'r') as infile1:
    for line in infile1:
        lines1.append(line.strip())

new_lines = []
for i in lines1:
    line_list = i.split('/')
    line_list.insert(5, "fasta_files")
    new_lines.append('/'.join(line_list))

with open("/home/conrad/ncbi_ecoli_genomes/mlst/ncbi_st1193_genomes_2.txt", 'w') as outfile1:
    to_write = '\n'.join(new_lines)
    outfile1.write(to_write)
