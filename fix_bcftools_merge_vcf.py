from sys import argv

script, input_vcf, output_vcf = argv

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

corrected_vcf_objects = []
for item in vcf_object_list:
    ref_base = item.ref
    alt_list = item.alt.split(',')
    # print ref_base, alt_list
    if len(alt_list) > 1:
        if ref_base in alt_list:
            ref_index = alt_list.index(ref_base)
            alt_list.pop(ref_index)
    
    new_vcf_object = VcfVariant(item.chrom, item.pos, item.idd, item.ref, ','.join(alt_list), item.quality, item.ffilter, item.iinfo)
    corrected_vcf_objects.append(new_vcf_object)

# for i in corrected_vcf_objects:
#     print i.chrom, i.pos, i.idd, i.ref, i.alt, i.quality, i.ffilter, i.iinfo

# print '\n'.join(header_lines)

to_write_vcf_records = []
for i in corrected_vcf_objects:
    record = i.chrom + '\t' + i.pos + '\t' +  i.idd + '\t' +  i.ref + '\t' +  i.alt + '\t' +  i.quality + '\t' +  i.ffilter + '\t' +  i.iinfo
    to_write_vcf_records.append(record)

with open(output_vcf, 'w') as outfile1:
    outfile1.write('\n'.join(header_lines) + '\n' + '\n'.join(to_write_vcf_records))