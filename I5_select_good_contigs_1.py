from sys import argv
import pandas
from Bio import SeqIO

script, inputcontigs, input_parsedblastoutput, output_stats, output_contigs = argv

print "Working on sample: " + inputcontigs.split('/')[-1]

class stats_object(object):
    def __init__(self, name, length, coverage, color, outline, record):
        self.name = name
        self.length = length
        self.coverage = coverage
        self.color = color
        self.outline = outline
        self.record = record

stats_objects = []
with open(input_parsedblastoutput, 'r') as infile1:
    for line in infile1:
        line_ = line.strip().split('\t')
        stats_objects.append(stats_object(line_[0], line_[1], line_[2], line_[3], line_[4], '\t'.join(line_)))
stats_objects.pop(0)

records = list(SeqIO.parse(inputcontigs, "fasta"))

print "Total number of contigs: ", len(stats_objects)

################
# Create list of good stats objects (> 10x coverage, > 1000 bp, color=deepskyblue):

good_stats_objects = []
for obj in stats_objects:
    # print obj.record
    if obj.color == 'deepskyblue':
        if float(obj.coverage) > 10.0:
            if int(obj.length) > 1000:
                good_stats_objects.append(obj)

print "Number of good contigs: ",len(good_stats_objects)

# Create list of good contig fasta objects:

good_contig_names = [i.name for i in good_stats_objects]

good_contigs = []
for i in records:
    if i.name in good_contig_names:
        good_contigs.append(i)

# Write good contigs and good stats objects to file:

with open(output_contigs, 'w') as outfile1:
    SeqIO.write(good_contigs, outfile1, "fasta")

with open(output_stats, 'w') as outfile2:
    outfile2.write("Name\tLength\tCoverage\tColor\tOutline\n")
    for i in good_stats_objects:
        outfile2.write(i.record + '\n')

