from sys import argv

script, input_vcf = argv

# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, idd, ref, alt, quality, ffilter, iinfo):
        self.chrom = chrom
        self.pos = pos
        self.idd = idd
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.ffilter = ffilter
        self.iinfo = iinfo

header_lines = []
vcf_object_list = []
with open(input_vcf, 'r') as infile1:
    for line in infile1:
        if not line.startswith("#"):
            split_line = line.strip().split('\t')
            # print split_line
            vcf_object = VcfVariant(split_line[0], split_line[1], split_line[2], split_line[3], split_line[4], split_line[5], split_line[6], split_line[7])
            vcf_object_list.append(vcf_object)
        elif line.startswith("#"):
            header_lines.append(line.strip())

print "number of unique indels: ", len(vcf_object_list)

indel_length_sum = 0
min_indel = 100
max_indel = 0
for i in vcf_object_list:
    max_ = max(len(i.ref), len(i.alt))
    min_ = min(len(i.ref), len(i.alt))
    indel_length = max_ - min_

    if indel_length > max_indel:
        max_indel = indel_length
    
    if indel_length < min_indel:
        min_indel = indel_length
    # Add to total indel length
    indel_length_sum += indel_length


print "average indel size: ", indel_length_sum/len(vcf_object_list)
print "max indel size: ", max_indel 
print "min indel size: ", min_indel