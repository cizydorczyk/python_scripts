from sys import argv
import os.path
import itertools

script, input_snp_alignment, input_snp_positions, input_checked_snps, output_dir = argv

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

n_total_depths = {}
n_hq_depths = {}
snp_total_depths = {}
snp_hq_depths = {}
n_unbal_total_depths = {}
n_unbal_hq_depths = {}
n_polym_total_depths = {}
n_polym_hq_depths = {}
n_low_qual_total_depths = {}
n_low_qual_hq_depths = {}
n_other_total_depths = {}
n_other_hq_depths = {}

for isolate in alignment_dict:
    
    records = {"too_few_alt_f_r_reads":0, "low_ref_alt_ratio":0, "other":0, "low_quality":0}
    
    seq_list = list(alignment_dict[isolate])
    ncount = seq_list.count('N')
    
    # Get positions called ambiguous (N):
    nindices = [i for i , x in enumerate(seq_list) if x == 'N']
    npositions = [snp_positions[i] for i in nindices]
        
    # Get SNP positions that are not called ambiguous (N):
    all_isolate_positions = [key for key in checked_snps_dict[isolate]]
    npositions_set = set(npositions)
    all_isolate_positions_set = set(all_isolate_positions)
    non_n_positions = list(all_isolate_positions_set.difference(npositions_set))

    # Get total and hq depth distributions:
    
    ## Get overall N depth distributions:
    for nposition in npositions:
        total_depth = int(checked_snps_dict[isolate][nposition][3])
        hq_depth = int(checked_snps_dict[isolate][nposition][4])

        if total_depth not in n_total_depths:
            n_total_depths[total_depth] = 1
        elif total_depth in n_total_depths:
            n_total_depths[total_depth] += 1
        
        if hq_depth not in n_hq_depths:
            n_hq_depths[hq_depth] = 1
        elif hq_depth in n_hq_depths:
            n_hq_depths[hq_depth] += 1

    ## Get distributions among categories of N sites:

        ref_dp4_f = int(checked_snps_dict[isolate][nposition][6])
        ref_dp4_r = int(checked_snps_dict[isolate][nposition][7])
        alt_dp4_f = int(checked_snps_dict[isolate][nposition][9])
        alt_dp4_r = int(checked_snps_dict[isolate][nposition][10])

        ### Among unbalanced forward/reverse read sites:
        if alt_dp4_f == 0 or alt_dp4_r == 0:
            
            if total_depth not in n_unbal_total_depths:
                n_unbal_total_depths[total_depth] = 1
            elif total_depth in n_unbal_total_depths:
                n_unbal_total_depths[total_depth] += 1
            
            if hq_depth not in n_unbal_hq_depths:
                n_unbal_hq_depths[hq_depth] = 1
            elif hq_depth in n_unbal_hq_depths:
                n_unbal_hq_depths[hq_depth] += 1
        
        # Among polymorphic sites:
        elif (float(alt_dp4_f) + float(alt_dp4_r))/(float(ref_dp4_f) + float(ref_dp4_r) + float(alt_dp4_f) + float(alt_dp4_r)) < 0.80:
            
            if total_depth not in n_polym_total_depths:
                n_polym_total_depths[total_depth] = 1
            elif total_depth in n_polym_total_depths:
                n_polym_total_depths[total_depth] += 1

            if hq_depth not in n_polym_hq_depths:
                n_polym_hq_depths[hq_depth] = 1
            elif hq_depth in n_polym_hq_depths:
                n_polym_hq_depths[hq_depth] += 1
        
        # Among low quality sites:
        elif float(checked_snps_dict[isolate][nposition][2]) < 25.0:

            if total_depth not in n_low_qual_total_depths:
                n_low_qual_total_depths[total_depth] = 1
            elif total_depth in n_low_qual_total_depths:
                n_low_qual_total_depths[total_depth] += 1
            
            if hq_depth not in n_low_qual_hq_depths:
                n_low_qual_hq_depths[hq_depth] = 1
            elif hq_depth in n_low_qual_hq_depths:
                n_low_qual_hq_depths[hq_depth] += 1
        
        # Among "other" sites:
        else:

            if total_depth not in n_other_total_depths:
                n_other_total_depths[total_depth] = 1
            elif total_depth in n_other_total_depths:
                n_other_total_depths[total_depth] += 1
            
            if hq_depth not in n_other_hq_depths:
                n_other_hq_depths[hq_depth] = 1
            elif hq_depth in n_other_hq_depths:
                n_other_hq_depths[hq_depth] += 1

    ## Get overall SNP depth distributions:
    for snpposition in non_n_positions:
        total_depth = int(checked_snps_dict[isolate][snpposition][3])
        hq_depth = int(checked_snps_dict[isolate][snpposition][4])

        if total_depth not in snp_total_depths:
            snp_total_depths[total_depth] = 1
        elif total_depth in snp_total_depths:
            snp_total_depths[total_depth] += 1
        
        if hq_depth not in snp_hq_depths:
            snp_hq_depths[hq_depth] = 1
        elif hq_depth in snp_hq_depths:
            snp_hq_depths[hq_depth] += 1

with open(output_dir + "/snp_total_depths_dist.txt", 'w') as outfile1:
    output_list = []
    for key, value in snp_total_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile1.write(to_write)

with open(output_dir + "/snp_hq_depths_dist.txt", 'w') as outfile2:
    output_list = []
    for key, value in snp_hq_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile2.write(to_write)

with open(output_dir + "/n_total_depths_dist.txt", 'w') as outfile3:
    output_list = []
    for key, value in n_total_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile3.write(to_write)

with open(output_dir + "/n_hq_depths_dist.txt", 'w') as outfile4:
    output_list = []
    for key, value in n_hq_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile4.write(to_write)

with open(output_dir + "/n_unbal_total_depths_dist.txt", 'w') as outfile5:
    output_list = []
    for key, value in n_unbal_total_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile5.write(to_write)

with open(output_dir + "/n_unbal_hq_depths_dist.txt", 'w') as outfile6:
    output_list = []
    for key, value in n_unbal_hq_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile6.write(to_write)

with open(output_dir + "/n_polym_total_depths_dist.txt", 'w') as outfile7:
    output_list = []
    for key, value in n_polym_total_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile7.write(to_write)

with open(output_dir + "/n_polym_hq_depths_dist.txt", 'w') as outfile8:
    output_list = []
    for key, value in n_polym_hq_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile8.write(to_write)

with open(output_dir + "/n_low_qual_total_depths_dist.txt", 'w') as outfile9:
    output_list = []
    for key, value in n_low_qual_total_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile9.write(to_write)

with open(output_dir + "/n_low_qual_hq_depths_dist.txt", 'w') as outfile10:
    output_list = []
    for key, value in n_low_qual_hq_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile10.write(to_write)

with open(output_dir + "/n_other_total_depths_dist.txt", 'w') as outfile11:
    output_list = []
    for key, value in n_other_total_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile11.write(to_write)

with open(output_dir + "/n_other_hq_depths_dist.txt", 'w') as outfile12:
    output_list = []
    for key, value in n_other_hq_depths.iteritems():
        output_list.append(str(key) + '\t' + str(value))
    to_write = '\n'.join(output_list)
    outfile12.write(to_write)
