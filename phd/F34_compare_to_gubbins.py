import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--gubbins_fasta", help="rc masked fasta file produced by gubbins & masked with maskrc-svg")
parser.add_argument("--filtered_fasta", help="fasta with sequences obtained by density filtering using F34.py scriptt")
parser.add_argument("--ref_fasta", help="reference fasta sequence")
parser.add_argument("--unmasked_fasta", help="unmasked fasta alignment")

args = parser.parse_args()

gubbins_seqs = list(SeqIO.parse(args.gubbins_fasta, "fasta"))
gubbins_dict = {}
for sequence in gubbins_seqs:
    gubbins_dict[sequence.id] = sequence

filtered_seqs = list(SeqIO.parse(args.filtered_fasta, "fasta"))
filtered_dict = {}
for sequence in filtered_seqs:
    filtered_dict[sequence.id] = sequence

ref_seq = SeqIO.read(args.ref_fasta, "fasta")

useq_seq = list(SeqIO.parse(args.unmasked_fasta, "fasta"))
useq_dict = {}
for sequence in useq_seq:
    useq_dict[sequence.id] = sequence

tp = 0
tn = 0
fp = 0
fn = 0

for record in filtered_seqs:
    gubbins_seq_record = [gubbins_dict[key] for key in gubbins_dict if key == record.id][0]
    useq_record = [useq_dict[key] for key in useq_dict if key == record.id][0]

    filtered_sequence = [i for i in record.seq]
    gubbins_sequence = [i for i in gubbins_seq_record.seq]
    useq_sequence = [i for i in useq_record.seq]
    ref_sequence = [i for i in ref_seq.seq]

    for fb, gb, ub, rb in zip(filtered_sequence, gubbins_sequence, useq_sequence, ref_sequence):
        if ub != rb and ub not in "N-": # SNP site
            if gb not in "N-" and fb not in "N-": # we don't care about masking at ambiguous positions; redundant to saying 'ub not in "N-"' in line above...
                if gb == fb and gb != "X" and fb != "X":
                    # True Negative
                    # Not masked by either gubbins nor filtering
                    tn += 1
                elif gb == fb and gb == "X" and fb == "X":
                    # True Positive
                    # Masked by both Gubbins and density filtering
                    tp += 1
                elif gb != fb: # if either gb or fb is masked:
                    if gb == "X" and fb != "X":
                        # False Negative
                        # Masked by Gubbins but not density filtering
                        fn += 1
                    elif gb != "X" and fb == "X":
                        # False Positive
                        # Masked by density filtering but not Gubbins
                        fp += 1
                    else:
                        print(gb, fb, rb)


print(tp, tn, fp, fn)




# for gseq in gubbins_seqs:
#     for fseq in filtered_seqs:
#
#
#         tp = 0
#         tn = 0
#         fp = 0
#         fn = 0
#         non_snp = 0
#
#         if gseq.id == fseq.id:
#             gseq_sequence = [i for i in gseq.seq]
#             fseq_sequence = [i for i in fseq.seq]
#             rseq_sequence = [i for i in ref_seq.seq]
#             useq_sequence = [i for i in useq_seq.seq]
#             print(len(fseq))
#
#
#             ref_count = 0
#             for rb, ub, gb, fb in zip(rseq_sequence, gseq_sequence, fseq_sequence):
#                 if gb == fb and gb != rb and fb != rb: # This is a SNP/masked position
#                     if gb not in "XN-" and fb not in "XN-": # This ensures it is NOT a masked position
#                         print(rb, gb, fb) # Thus this is a SNP



                # if gb == fb and gb == rb:
                #     # Invariant site
                #     ref_count += 1
                # elif gb == fb and gb != rb:
                #     # True Negative
                #     tn += 1
                # elif





            #     if gb != fb:
            #         if gb == "X" and fb != "X":
            #             fn += 1 # FN
            #         elif gb != "X" and fb == "X":
            #             fp += 1 # FP
            #         elif gb != "X" and fb != "X":
            #             tn += 1 # TN
            #     elif gb == "X" and fb == "X":
            #         tp += 1 # TP
            #     else:
            #         non_snp += 1 # non-SNPs or one sequence has ambiguous base (shouldn't happen; both should have N or -)
            #
            # sum_stats = tp+tn+fp+fn+non_snp
            # print(tp, tn, fp, fn, non_snp, sum_stats)
