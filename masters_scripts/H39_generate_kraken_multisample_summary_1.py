from sys import argv
import pandas as pd
import os.path

script, input_kraken_report, output_tsv = argv

isolate_number = input_kraken_report.strip().split('_')[0]
print isolate_number

class kraken_object(object):
    def __init__(self, perc_reads, root_reads, taxon_reads, rank, taxid, sciname, summary):
        self.perc_reads = perc_reads
        self.root_reads = root_reads
        self.taxon_reads = taxon_reads
        self.rank = rank
        self.taxid = taxid
        self.sciname = sciname
        self.summary = summary

kraken_objects = []
with open(input_kraken_report, 'r') as infile1:
    for line in infile1:
        line_ = line.strip().split('\t')
        line_object = kraken_object(float(line_[0]), int(line_[1]), int(line_[2]), line_[3], int(line_[4]), line_[5].strip(), ", ".join(line_))
        kraken_objects.append(line_object)

sum_unclassified = 0
sum_aeruginosa = 0
sum_pseudomonads = 0

to_plot_species = []
for i in kraken_objects:
    if i.rank == "U":
        if i.perc_reads > 0.00:
            print i.rank, i.perc_reads, i.sciname
            sum_unclassified += i.perc_reads
    elif i.rank == "S":
        if "aeruginosa" in i.sciname:
            if i.perc_reads > 0.00:
                print i.rank, i.perc_reads, i.sciname
                sum_aeruginosa += i.perc_reads
        elif "Pseudomonas" in i.sciname and "aeruginosa" not in i.sciname:
            if i.perc_reads > 0.00:
                print i.rank, i.perc_reads, i.sciname
                sum_pseudomonads += i.perc_reads
    else:
        pass

sum_other = (100-sum_unclassified-sum_aeruginosa-sum_pseudomonads)

to_write = '\t'.join([str(isolate_number), str(sum_aeruginosa), str(sum_pseudomonads), str(sum_other), str(sum_unclassified)])
print isolate_number, sum_aeruginosa, sum_pseudomonads, sum_other, sum_unclassified

if os.path.isfile(output_tsv) == True:
    with open(output_tsv, 'a+') as outfile1:
        outfile1.write(to_write + '\n')
else:
    with open(output_tsv, 'w') as outfile1:
        outfile1.write("Isolate" + '\t' + "%_p_aeruginosa" + '\t' + "%_pseudomonads" + '\t' + "%_other" + '\t' "%_unclassified" + '\n')
        outfile1.write(to_write + '\n')