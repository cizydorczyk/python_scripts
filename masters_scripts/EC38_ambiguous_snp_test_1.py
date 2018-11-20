from sys import argv
import os.path

script, input_snp_alignment, input_snp_positions, input_checked_snps, output_file = argv

# Read snp alignment:
alignment_dict = {}

with open(input_snp_alignment, 'r') as infile1:
    for line in infile1:
        if line.startswith(">"):
            alignment_dict[line.strip()[1:]] = next(infile1)

# Read snp positions:
snp_positions = []

with open(input_snp_positions, 'r') as infile2:
    for line in infile2:
        snp_positions.append(line.strip().split('-')[1])

# Read checked snps:
checked_snps_dict = {}

with open(input_checked_snps, 'r') as infile3:
    for line in infile3:
        line_isolate_number = line.strip().split('\t')[0].split('-')[1].split('/')[0]
        line_position = line.strip().split('\t')[2]

        if line_isolate_number not in set([key for key in checked_snps_dict]):

            checked_snps_dict[line_isolate_number] = {}
            checked_snps_dict[line_isolate_number][line_position] = line.strip().split('\t')[3:]

        elif line_isolate_number in set([key for key in checked_snps_dict]):
            checked_snps_dict[line_isolate_number][line_position] = line.strip().split('\t')[3:]

for isolate in alignment_dict:
    
    records = {"too_few_alt_f_r_reads":0, "low_ref_alt_ratio":0, "other":0, "low_quality":0}
    
    seq_list = list(alignment_dict[isolate])
    ncount = seq_list.count('N')
    indices = [i for i , x in enumerate(seq_list) if x == 'N']
    
    for index in indices:
        position = snp_positions[index]
        print checked_snps_dict[isolate][position]
        ref_dp4_f = int(checked_snps_dict[isolate][position][6])
        ref_dp4_r = int(checked_snps_dict[isolate][position][7])
        alt_dp4_f = int(checked_snps_dict[isolate][position][9])
        alt_dp4_r = int(checked_snps_dict[isolate][position][10])

        if alt_dp4_f == 0 or alt_dp4_r == 0:
            records["too_few_alt_f_r_reads"] += 1
        
        elif (float(alt_dp4_f) + float(alt_dp4_r))/(float(ref_dp4_f) + float(ref_dp4_r) + float(alt_dp4_f) + float(alt_dp4_r)) < 0.80:
            records["low_ref_alt_ratio"] += 1
        
        elif float(checked_snps_dict[isolate][position][2]) < 25.0:
            records["low_quality"] += 1
        else:
            records["other"] += 1

    if not os.path.isfile(output_file):
        with open(output_file, 'w') as outfile1:
            record = isolate + "\t" + str(records["too_few_alt_f_r_reads"]) + '\t' + str(records["low_ref_alt_ratio"]) + '\t' + str(records["low_quality"]) + '\t' + str(records["other"])
            outfile1.write("Isolate\ttoo_few_alt_reads\tlow_ref_alt_ratio\tlow_quality\tother\n" + record)

    elif os.path.isfile(output_file):
        with open(output_file, 'a') as outfile1:
            record = isolate + "\t" + str(records["too_few_alt_f_r_reads"]) + '\t' + str(records["low_ref_alt_ratio"]) + '\t' + str(records["low_quality"]) + '\t' + str(records["other"])
            outfile1.write("\n" + record)

