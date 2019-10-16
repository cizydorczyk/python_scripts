isolates = []
with open("/home/conrad/smaltophilia/ncbi_genomes/partial_and_complete_genomes/sm_partial_and_complete_genomes_list.txt", 'r') as infile1:
    for line in infile1:
        if line.startswith("GCF"):
            line_elements = line.strip().split("_")
            isolate = line_elements[0] + "_" + line_elements[1] + "_" + line_elements[2]
            isolates.append(isolate)

with open("/home/conrad/smaltophilia/ncbi_genomes/partial_and_complete_genomes/sm_partial_and_complete_genomes_list2.txt", 'w') as outfile1:
    to_write = "\n".join(isolates)
    outfile1.write(to_write)
