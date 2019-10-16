import os.path

file_list = []
contigs_list = []
with open("/home/conrad/hinfluenzae/unicycler_testing/read_mapping_test/uncorrected/uncor.sorted.depths.txt", "r") as infile1:
    for line in infile1:
        file_list.append(line.strip().split("\t"))
        contigs_list.append(line.strip().split("\t")[0])

contigs_list_set = set(contigs_list)

for i in contigs_list_set:
    sum = 0
    length = int(i.split("_")[3])
    for j in file_list:
        if i == j[0]:
            print(j)
            sum += int(j[2])

    average_cov = sum/length
    print(i, average_cov)
