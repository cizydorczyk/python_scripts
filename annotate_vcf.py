from sys import argv

script, tab_annotation, GTF_intron_annotation, input_vcf, output_annotations_file = argv

# Create table of codons:
codons = {'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'TGC':'C', 'TGT':'C', 'GAC':'D', 'GAT':'D',\
'GAA':'E', 'GAG':'E', 'TTC':'F', 'TTT':'F', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'CAC':'H', \
'CAT':'H', 'ATA':'I', 'ATC':'I', 'ATT':'I', 'AAA':'K', 'AAG':'K', 'TTA':'L', 'TTG':'L', 'CTA':'L', \
'CTC':'L', 'CTG':'L', 'CTT':'L', 'ATG':'M', 'AAC':'N', 'AAT':'N', 'CCA':'P', 'CCC':'P', 'CCG':'P', \
'CCT':'P', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R', 'AGG':'R', \
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'AGC':'S', 'AGT':'S', 'ACA':'T', 'ACC':'T', 'ACG':'T', \
'ACT':'T', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'TGG':'W', 'TAC':'Y', 'TAT':'Y', 'TAA':'St1', \
'TAG':'St2', 'TGA':'St3'}


# Intergenic annotation format:
# type (intergenic), start, end, length of intergenic region, left gene, right gene
# Parse intergenic annotation file:
GTF_intergenic_dict = {}
GTF_CDS_dict = {}
GTF_pseudo_dict = {}
GTF_rRNA_dict = {}
GTF_tRNA_dict = {}
GTF_tmRNA_dict = {}
GTF_ncRNA_dict = {}

with open(GTF_intron_annotation, 'r') as infile1:
    for line in infile1:
        splitline = line.strip().split('\t')
        line_type = splitline[2]
        key = (int(splitline[3]), int(splitline[4]))
        if line_type == 'CDS':
            GTF_CDS_dict[key] = splitline[8].split(';')[0].split(' ')[1].strip('"')
        elif line_type == 'intergenic':
            GTF_intergenic_dict[key] = splitline[8].split(';')[3] + "_" + splitline[8].split(';')[4]
        elif line_type == 'pseudo':
            GTF_pseudo_dict[key] = splitline[8].split(';')[0].split(' ')[1].strip('"')
        elif line_type == 'rRNA':
            GTF_rRNA_dict[key] = splitline[8].split(';')[0].split(' ')[1].strip('"')
        elif line_type == 'tRNA':
            GTF_tRNA_dict[key] = splitline[8].split(';')[0].split(' ')[1].strip('"')
        elif line_type == 'tmRNA':
            GTF_tmRNA_dict[key] = splitline[8].split(';')[0].split(' ')[1].strip('"')
        elif line_type == 'ncRNA':
            GTF_ncRNA_dict[key] = splitline[8].split(';')[0].split(' ')[1].strip('"')

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

snps_dict = {}
for i in vcf_object_list:
    snps_dict[i.pos+'_'+i.alt] = [i.ref, i.alt]

# Create function for testing if snp is in range of a gene/intergenic region:
def lookup_bounds(bounds, value):
    for min_, max_ in bounds:
        # Don't have to test value < min because this is being run on a sorted
        # list of tuples, so a value < min would have already been accepted in the
        # previous tuple, and if it's not in any tuples, the elif statement
        # covers that and returns something in the same step. It's broken down
        # into two statements to speed things up,
        # so that if the value is larger than the max of the tuple, it doesn't
        # even bother testing the minimum (making two comparisons):
        if value > max_:
            continue
        elif min_ <= value:
            return min_, max_
        else:
            return

# Create sorted tuple lists of keys for each annotation dictionary:

sorted_GTF_intergenic_dict = sorted(GTF_intergenic_dict)
sorted_GTF_CDS_dict = sorted(GTF_CDS_dict)
sorted_GTF_pseudo_dict = sorted(GTF_pseudo_dict)
sorted_GTF_rRNA_dict = sorted(GTF_rRNA_dict)
sorted_GTF_tRNA_dict = sorted(GTF_tRNA_dict)
sorted_GTF_tmRNA_dict = sorted(GTF_tmRNA_dict)
sorted_GTF_ncRNA_dict = sorted(GTF_ncRNA_dict)

# Create tab annotation dict:
tab_ann_dict = {}
with open(tab_annotation, 'r') as infile2:
    for line in infile2:
        if not line.startswith("#") and not line.startswith("Sequence"):
            splitline = line.strip().split('\t')
            gene = splitline[1]
            start = splitline[2]
            end = splitline[3]
            seq = splitline[4]
            tab_ann_dict[gene] = [gene, start, end, seq]
        

# print tab_ann_dict      
#~~~~~~~~~~~~~~~~~~~~Actual heavy lifting script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
intergeniccount = 0
pseudocount = 0
rrnacount = 0
trnacount = 0
tmrnacount = 0
ncrnacount = 0
errorcount = 0
synonymouscount = 0
nonsynonymouscount = 0
ambiguouscount = 0


# Output list:
to_file = []
# For every snp in snps dictionary:
for i in sorted(snps_dict):
    # Identify genic region snp falls in:
    snppos = int(i.split("_")[0])
    bound = lookup_bounds(sorted_GTF_CDS_dict, snppos)

    # If falls in genic region:
    if bound is not None:
        gene = GTF_CDS_dict[bound]
        alt_base = i.split("_")[1]
        sample_sequence = list(tab_ann_dict[gene][3])
        # Zero-based index of the SNP in the gene:
        snp_index = snppos - int(tab_ann_dict[gene][2])
        codon = ''
         # Identify reference codon, depending on position in codon alternate
        # base falls in, but recording the reference codon:

        if (snp_index + 1) % 3 == 0:
            try:
                # print 't1'
                codon = ''.join([sample_sequence[snp_index-2], sample_sequence[snp_index-1], sample_sequence[snp_index]])
            except IndexError:
                print 'e1'
                codon = 'NN' + sample_sequence[snp_index]
            # print codon

        elif (snp_index + 2) % 3 == 0:
            try:
                # print 't2'
                codon = ''.join([sample_sequence[snp_index-1], sample_sequence[snp_index], sample_sequence[snp_index+1]])
            except IndexError:
                print 'e2'
                codon = 'N' + sample_sequence[snp_index] + 'N'
            # print codon

        elif (snp_index + 3) % 3 == 0:
            try:
                # print 't3'
                codon = ''.join([sample_sequence[snp_index], sample_sequence[snp_index+1], sample_sequence[snp_index+2]])
            except IndexError:
                print 'e3'
                codon = sample_sequence[snp_index] + 'NN'
            # print codon

    # Insert alternate base into reference sequence:
        sample_sequence[snp_index] = alt_base

        # Alt codon:
        codon_alt = ''
        # Position in codon:
        position_in_codon = ''
         # Identify alternate codon based on position in codon of alternate base,
         # this time including alternate base:
        if (snp_index + 1) % 3 == 0:
            try:
                # print 't1'
                codon_alt = ''.join([sample_sequence[snp_index-2], sample_sequence[snp_index-1], sample_sequence[snp_index]])
                position_in_codon = 3
            except IndexError:
                print 'e1'
                codon_alt = 'NN' + sample_sequence[snp_index]
                position_in_codon = 3
            # print codon_alt

        elif (snp_index + 2) % 3 == 0:
            try:
                # print 't2'
                codon_alt = ''.join([sample_sequence[snp_index-1], sample_sequence[snp_index], sample_sequence[snp_index+1]])
                position_in_codon = 2
            except IndexError:
                print 'e2'
                codon_alt = 'N' + sample_sequence[snp_index] + 'N'
                position_in_codon = 2
            # print codon_alt

        elif (snp_index + 3) % 3 == 0:
            try:
                # print 't3'
                codon_alt = ''.join([sample_sequence[snp_index], sample_sequence[snp_index+1], sample_sequence[snp_index+2]])
                position_in_codon = 1
            except IndexError:
                print 'e3'
                codon_alt = sample_sequence[snp_index] + 'NN'
                position_in_codon = 1
            # print codon_alt


        # print codon, codon_alt, i

        # Identify ref and alt aa:
        ref_aa = ''
        if 'N' not in codon:
            ref_aa = codons[codon]
        elif 'N' in codon:
            ref_aa = 'X'
        alt_aa = ''
        if 'N' not in codon_alt:
            alt_aa = codons[codon_alt]
        elif 'N' in codon_alt:
            alt_aa = 'X'

        # Type of codon/aa change:
        codon_change = ''

        if ref_aa == alt_aa:
            if ref_aa == 'X' and alt_aa == 'X':
                codon_change = "ambiguous"
                ambiguouscount += 1
            else:
                codon_change = "synonymous"
                synonymouscount += 1

        elif ref_aa != alt_aa:
            if alt_aa == '*': #doesn't do anything anymore
                codon_change = "stop" #doesn't do anything anymore
            elif alt_aa == 'X':
                codon_change = "ambiguous"
                ambiguouscount += 1
            else:
                codon_change = "non-synonymous"
                nonsynonymouscount += 1
        
        outputline = str(snppos) + '\t' + alt_base + '\t' + str(position_in_codon) + '\t' + codon + '\t' + codon_alt + '\t' + codon_change
        # print outputline
        to_file.append(outputline)
    
    elif bound is None:
        bound = lookup_bounds(sorted_GTF_intergenic_dict, snppos)
        altbase = i.split("_")[1]

        if bound is not None:
            outputline = str(snppos) + '\t' + altbase + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + 'intergenic'
            to_file.append(outputline)
            intergeniccount += 1
            # print outputline
        
        elif bound is None:
            bound = lookup_bounds(sorted_GTF_pseudo_dict, snppos)

            if bound is not None:
                outputline = str(snppos) + '\t' + altbase + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + 'pseudogene'
                to_file.append(outputline)
                pseudocount += 1

            elif bound is None:
                bound = lookup_bounds(sorted_GTF_rRNA_dict, snppos)

                if bound is not None:
                    outputline = str(snppos) + '\t' + altbase + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + 'rRNA'
                    to_file.append(outputline)
                    rrnacount += 1

                elif bound is None:
                    bound = lookup_bounds(sorted_GTF_tRNA_dict, snppos)

                    if bound is not None:
                        outputline = str(snppos) + '\t' + altbase + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + 'tRNA'
                        to_file.append(outputline)
                        trnacount += 1
                    
                    elif bound is None:
                        bound = lookup_bounds(sorted_GTF_tmRNA_dict, snppos)

                        if bound is not None:
                            outputline = str(snppos) + '\t' + altbase + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + 'tmRNA'
                            to_file.append(outputline)
                            tmrnacount += 1
                        
                        elif bound is None:
                            bound = lookup_bounds(sorted_GTF_ncRNA_dict, snppos)

                            if bound is not None:
                                outputline = str(snppos) + '\t' + altbase + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + 'ncRNA'
                                to_file.append(outputline)
                                ncrnacount += 1
                            
                            elif bound is None:
                                outputline = str(snppos) + '\t' + altbase + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + 'error'
                                to_file.append(outputline)
                                print outputline
                                errorcount += 1

header = "Position" + '\t' + "Alt" + '\t' + "Pos_in_Codon" + '\t' + "Ref_Codon" + '\t' + "Alt_Codon" + '\t' + "Type"
with open(output_annotations_file, 'w') as outfile1:
    outfile1.write(header + '\n' + '\n'.join(to_file))

print "total # of SNPs: ", len(snps_dict)
print "synonymous snps: ", synonymouscount
print "non-synonymous snps: ", nonsynonymouscount
print "intergenic snps: ", intergeniccount
print "rRNA snps: ", rrnacount
print "tRNA snps: ", trnacount
print "tmRNA snps: ", tmrnacount
print "ncRNA snps: ", ncrnacount
print "pseudogene snps: ", pseudocount
print "ambiguous snps: ", ambiguouscount
print "error snps (intergenic): ", errorcount




