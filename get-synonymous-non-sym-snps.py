from sys import argv
# script, input_vcf_list = argv

# input_vcf_files = []
# with open(input_vcf_list, 'r') as infile1:
#     for line in infile1:
#         input_vcf_files.append(line.strip())

script, input_vcf = argv

vcf_snp_count = 0

sequence_ontology_list = {"coding_sequence_variant":0, "chromosome":0, "duplication":0, "inversion":0, "inframe_insertion":0, "disruptive_inframe_insertion":0, "inframe_deletion":0, "disruptive_inframe_deletion":0, "exon_variant":0, "exon_loss_variant":0, "frameshift_variant":0, "feature_ablation":0, "gene_fusion":0, "bidirectional_gene_fusion":0, "rearranged_at_DNA_level":0,  "intergenic_region":0, "conserved_intergenic_variant":0, "intragenic_variant":0, "miRNA":0, "missense_variant":0, "initiator_codon_variant":0, "stop_retained_variant":0, "protein_protein_contact":0, "structural_interaction_variant":0, "rare_amino_acid_variant":0, "splice_acceptor_variant":0, "splice_donor_variant":0, "splice_region_variant":0, "stop_lost":0, "start_lost":0, "stop_gained":0, "synonymous_variant":0, "start_retained":0, "stop_retained_variant":0, "transcript_variant":0, "regulatory_region_variant":0}


with open(input_vcf, 'r') as infile2:
    for line in infile2:
        if not line.startswith("#"):
            ann_line = line.strip().split('\t')[7]
            if len(ann_line) > 1:
                vcf_snp_count += 1
            for term in sequence_ontology_list:
                if term in ann_line:
                    sequence_ontology_list[term] += 1

print sequence_ontology_list
print "total snps: ", vcf_snp_count
snps_sum = 0
for term in sequence_ontology_list:
    snps_sum += sequence_ontology_list[term]
print "snps annotated: ", snps_sum




# start_lost = 0
# stop_gained = 0
# synonymous_variant = 0
# start_retained = 0
# stop_retained_variant = 0
# transcript_variant = 0
# regulatory_region_variant = 0
# upstream_gene_variant = 0
# coding_sequence_variant = 0
# chromosome = 0
# duplication = 0
# inversion = 0
# inframe_insertion = 0
# disruptive_inframe_insertion = 0
# inframe_deletion = 0
# disruptive_inframe_deletion = 0
# downstream_gene_variant = 0
# exon_variant = 0
# exon_loss_variant = 0
# frameshift_variant = 0
# gene_variant = 0
# feature_ablation = 0
# gene_fusion = 0
# bidirectional_gene_fusion = 0
# rearranged_at_DNA_level = 0
# intergenic_region = 0
# conserved_intergenic_variant = 0
# intragenic_variant = 0
# miRNA = 0
# missense_variant = 0
# initiator_codon_variant = 0
# stop_retained_variant = 0
# protein_protein_contact = 0
# structural_interaction_variant = 0
# rare_amino_acid_variant = 0
# splice_acceptor_variant = 0
# splice_donor_variant = 0
# splice_region_variant = 0
# stop_lost = 0



# with open(input_vcf, 'r') as infile2:
#     for line in infile2:
#         if not line.startswith("#"):
#             ann_line = line.strip().split('\t')[7]
#             if len(ann_line) > 1:
#                 vcf_snp_count += 1
#             if 'intergenic' in ann_line:
#                 intergenic_snps += 1
#             if 'synonymous' in ann_line:
#                 synonymous_snps += 1

# print "total snps ",vcf_snp_count
# print "intergenic snps ", intergenic_snps
# print "synonymous snps ", synonymous_snps