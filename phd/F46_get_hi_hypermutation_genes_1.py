import argparse
import os.path
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_vcf", help="input vcf file")
parser.add_argument("--out", help="output file")
parser.add_argument("--isolate", help="isolate number")
parser.add_argument("--st", help="sequence type (just the number; 131, 73, or 1193)")

args = parser.parse_args()

# Create empty list to hold VCF variants:
variants = []

# Define VcfVariants class;
class VcfVariant(object):
    def __init__(self, hypermutationgene, isolateid, chrom, pos, ref, alt, info, filename):
        self.hypermutationgene = hypermutationgene
        self.isolateid = isolateid
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info
        self.filename = filename

# Classify each variant position in VCF file by type of variant (SNP, insertion, or deletion):
with open(args.input_vcf, 'r') as infile1:
    #isolate_id = os.path.basename(args.input_vcf)[0:-4]
    for line in infile1:
        if not line.startswith("#"):
            line_list = line.strip().split('\t')
            vcf_variant_obj = VcfVariant("NA", args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], line_list[7], args.input_vcf)
            variants.append(vcf_variant_obj)

# Create dictionaries with gene corodinates for each ST:
st12_hi_genes = {"mutL":["NZ_CP031688",1015068,1016957], "mutS":["NZ_CP031688",326916,329501],
                "mutH":["NZ_CP031688",615497,616168], "uvrD":["NZ_CP031688",1665611,1667791],
                "dam":["NZ_CP031688",858356,859216], "dnaQ":["NZ_CP031688",929746,930513],
                "polA":["NZ_CP031688",143359,146166], "recA":["NZ_CP031688",420694,421758],
                "mfd":["NZ_CP031688",1493257,1496697], "rep":["NZ_CP031688",382403,384418],
                "lexA":["NZ_CP031688",237202,237825], "mutT":["NZ_CP031688",89960,90370],
                "mutM":["NZ_CP031688",48103,48918], "mutY":["NZ_CP031688",225199,226335],
                "ung":["NZ_CP031688",1069387,1070046], "sodA":["NZ_CP031688",1762801,1763448],
                "oxyR":["NZ_CP031688",449062,449967]}
st12_all_hi_genes = {"mutL":["NZ_CP031688.1",1015068,1016957], "mutS":["NZ_CP031688.1",326916,329501],
                "mutH":["NZ_CP031688.1",615497,616168], "uvrD":["NZ_CP031688.1",1665611,1667791],
                "dam":["NZ_CP031688.1",858356,859216], "dnaQ":["NZ_CP031688.1",929746,930513],
                "polA":["NZ_CP031688.1",143359,146166], "recA":["NZ_CP031688.1",420694,421758],
                "mfd":["NZ_CP031688.1",1493257,1496697], "rep":["NZ_CP031688.1",382403,384418],
                "lexA":["NZ_CP031688.1",237202,237825], "mutT":["NZ_CP031688.1",89960,90370],
                "mutM":["NZ_CP031688.1",48103,48918], "mutY":["NZ_CP031688.1",225199,226335],
                "ung":["NZ_CP031688.1",1069387,1070046], "sodA":["NZ_CP031688.1",1762801,1763448],
                "oxyR":["NZ_CP031688.1",449062,449967]}
st14_hi_genes = {"mutL":["NZ_CP007472",1656263,1658152], "mutS":["NZ_CP007472",1004243,1006828],
                "mutH":["NZ_CP007472",1292805,1293476], "uvrD":["NZ_CP007472",464452,466632],
                "dam":["NZ_CP007472",1499078,1499938], "dnaQ":["NZ_CP007472",1570691,1571458],
                "polA":["NZ_CP007472",783624,786416], "recA":["NZ_CP007472",1103973,1105037],
                "mfd":["NZ_CP007472",292329,295769], "rep":["NZ_CP007472",1064622,1066640],
                "lexA":["NZ_CP007472",914402,915025], "mutT":["NZ_CP007472",730347,730757],
                "mutM":["NZ_CP007472",687683,688498], "mutY":["NZ_CP007472",865401,866537],
                "ung":["NZ_CP007472",1710438,1711097], "sodA":["NZ_CP007472",561571,562218],
                "oxyR":["NZ_CP007472",1127141,1128046]}
st103_hi_genes = {"mutL":["NZ_CP020008",70192,72081], "mutS":["NZ_CP020008",770635,773220],
                "mutH":["NZ_CP020008",475718,476389], "uvrD":["NZ_CP020008",1290431,1292614],
                "dam":["NZ_CP020008",227682,228542], "dnaQ":["NZ_CP020008",156410,157177],
                "polA":["NZ_CP020008",954069,956876], "recA":["NZ_CP020008",674210,675274],
                "mfd (BV121_RS07240)":["NZ_CP020008",1461871,1465311], "rep":["NZ_CP020008",711536,713551],
                "lexA":["NZ_CP020008",861365,861988], "mutT":["NZ_CP020008",1017263,1017673],
                "mutM":["NZ_CP020008",968475,969287], "mutY":["NZ_CP020008",873881,875017],
                "ung (BV121_RS00095)":["NZ_CP020008",18688,19347], "sodA":["NZ_CP020008",1194799,1195437],
                "oxyR":["NZ_CP020008",646582,647487]}
st124_hi_genes = {"mutL":["NC_022356",1598577,1600478], "mutS":["NC_022356",385957,388542],
                "mutH":["NC_022356",1792145,1792816], "uvrD":["NC_022356",1468917,1471100],
                "dam":["NC_022356",1702272,1703132], "dnaQ":["NC_022356",1748825,1749595],
                "polA":["NC_022356",1685446,1688238], "recA":["NC_022356",241995,243059],
                "mfd (HIFGL_RS04735)":["NC_022356",972003,975443], "rep":["NC_022356",286327,288342],
                "lexA":["NC_022356",191450,192073], "mutT":["NC_022356",537115,537525],
                "mutM":["NC_022356",577844,578659], "mutY":["NC_022356",435062,436198],
                "ung (HifGL_001516)":["NC_022356",1543908,1544567], "sodA":["NC_022356",757674,758321],
                "oxyR":["NC_022356",219618,220523]}
st145_hi_genes = {"mutL":["NZ_CP031686",1041118,1043007], "mutS":["NZ_CP031686",334976,337561],
                "mutH":["NZ_CP031686",631273,631944], "uvrD":["NZ_CP031686",1713073,1715253],
                "dam":["NZ_CP031686",884796,885656], "dnaQ":["NZ_CP031686",955793,956560],
                "polA":["NZ_CP031686",151097,153889], "recA":["NZ_CP031686",432344,433408],
                "mfd":["NZ_CP031686",1274960,1278400], "rep":["NZ_CP031686",394183,396198],
                "lexA":["NZ_CP031686",245940,246563], "mutT":["NZ_CP031686",90197,90607],
                "mutM":["NZ_CP031686",48584,49399], "mutY":["NZ_CP031686",232911,234047],
                "ung":["NZ_CP031686",1093852,1094511], "sodA":["NZ_CP031686",1810292,1810939],
                "oxyR":["NZ_CP031686",460829,461734]}
st321_hi_genes = {"mutL":["NZ_CP008740",73185,75074], "mutS":["NZ_CP008740",828324,830909],
                "mutH":["NZ_CP008740",549248,549919], "uvrD":["NZ_CP008740",1376152,1378335],
                "dam":["NZ_CP008740",271347,272207], "dnaQ":["NZ_CP008740",193350,194120],
                "polA":["NZ_CP008740",1008846,1011653], "recA":["NZ_CP008740",744974,746038],
                "mfd (C645_RS09185)":["NZ_CP008740",1803268,1806708], "rep":["NZ_CP008740",784024,786036],
                "lexA":["NZ_CP008740",907796,908419], "mutT":["NZ_CP008740",1074662,1075072],
                "mutM":["NZ_CP008740",1118891,1119706], "mutY":["NZ_CP008740",919356,920492],
                "ung (C645_RS00095)":["NZ_CP008740",18681,19340], "sodA":["NZ_CP008740",1246616,1247251],
                "oxyR":["NZ_CP008740",716226,717131]}
st393_A037_H11_hi_genes = {"mutL":["11",48503,50383], "mutS":["20",14137,16722],
                "mutH":["1",113712,114383], "uvrD":["4",106609,108792],
                "dam":["7",17599,18459], "dnaQ":["15",8660,9430],
                "polA":["3",1274,4081], "recA":["23",7985,9049],
                "mfd":["2",173012,176494], "rep":["9",22552,24567],
                "lexA":["5",32281,32904], "mutT":["3",54214,54624],
                "mutM":["3",98449,99264], "mutY":["5",44789,45925],
                "ung":["10",53557,54216], "sodA":["14",162,797],
                "oxyR":["18",4177,5082]}
st393_A360_H04_hi_genes = {"mutL":["5",51847,53727], "mutS":["16",13922,16507],
                "mutH":["4",179372,179656], "uvrD":["3",43520,45703],
                "dam":["6",17640,18500], "dnaQ":["6",91188,91958],
                "polA":["1",120436,123243], "recA":["25",5991,7055],
                "mfd":["2",172860,176342], "rep":["13",27205,29220],
                "lexA":["1",32281,32904], "mutT":["1",173378,173788],
                "mutM":["1",217613,218428], "mutY":["1",44789,45925],
                "ung":["5",107600,108259], "sodA":["3",144479,145114],
                "oxyR":["4",4177,5082]}


# Check if any SNPs/insertions/deletions occur in hypermutation genes:

# Create empty lists to hold results:
hypermutation_variants = []

if args.st == "12":
    for variant in variants:
        for gene in st12_hi_genes:
            if st12_hi_genes[gene][0] == variant.chrom:
                if int(variant.pos) >= st12_hi_genes[gene][1] and int(variant.pos) <= st12_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

if args.st == "12all": # For running this script vs. F31 all isolates vs. ST-12
    for variant in variants:
        for gene in st12_all_hi_genes:
            if st12_all_hi_genes[gene][0] == variant.chrom:
                if int(variant.pos) >= st12_all_hi_genes[gene][1] and int(variant.pos) <= st12_all_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

if args.st == "14":
    for variant in variants:
        for gene in st14_hi_genes:
            if st14_hi_genes[gene][0] == variant.chrom:
                if int(variant.pos) >= st14_hi_genes[gene][1] and int(variant.pos) <= st14_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

if args.st == "103":
    for variant in variants:
        for gene in st103_hi_genes:
            if st103_hi_genes[gene][0] == variant.chrom:
                if int(variant.pos) >= st103_hi_genes[gene][1] and int(variant.pos) <= st103_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

if args.st == "124":
    for variant in variants:
        for gene in st124_hi_genes:
            if st124_hi_genes[gene][0] == variant.chrom:
                if int(variant.pos) >= st124_hi_genes[gene][1] and int(variant.pos) <= st124_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

if args.st == "145":
    for variant in variants:
        for gene in st145_hi_genes:
            if st145_hi_genes[gene][0] == variant.chrom:
                if int(variant.pos) >= st145_hi_genes[gene][1] and int(variant.pos) <= st145_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

if args.st == "321":
    for variant in variants:
        for gene in st321_hi_genes:
            if st321_hi_genes[gene][0] == variant.chrom:
                if int(variant.pos) >= st321_hi_genes[gene][1] and int(variant.pos) <= st321_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

if args.st == "393_A037":
    for variant in variants:
        for gene in st393_A037_H11_hi_genes:
            # print(st393_A037_H11_hi_genes[gene][0], variant.chrom)
            if st393_A037_H11_hi_genes[gene][0] == variant.chrom:
                print(st393_A037_H11_hi_genes[gene][0], variant.chrom)
                if int(variant.pos) >= st393_A037_H11_hi_genes[gene][1] and int(variant.pos) <= st393_A037_H11_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

if args.st == "393_A360":
    for variant in variants:
        for gene in st393_A360_H04_hi_genes:
            if st393_A360_H04_hi_genes[gene][0] == variant.chrom:
                print(st393_A360_H04_hi_genes[gene][0], variant.chrom)
                if int(variant.pos) >= st393_A360_H04_hi_genes[gene][1] and int(variant.pos) <= st393_A360_H04_hi_genes[gene][2]:
                    variant.hypermutationgene = gene
                    hypermutation_variants.append(variant)
                else:
                    continue

# Write results to file:

header = "hypermutation_gene\tisolate\tchrom\tpos\tref\talt\tinfo"

if len(hypermutation_variants) == 0:
    to_write = 'NA' + '\t' + args.isolate + '\tNA\tNA\tNA\tNA\tNA'

    if os.path.isfile(args.out):
        with open(args.out, 'a') as outfile:
            outfile.write('\n' + to_write)

    else:
        with open(args.out, 'w') as outfile:
            outfile.write(header + '\n' + to_write)

elif len(hypermutation_variants) != 0:
    to_write = []
    for i in hypermutation_variants:
        write_str = i.hypermutationgene + '\t' + i.isolateid + '\t' + i.chrom + '\t' + i.pos + '\t' + i.ref + '\t' + i.alt + '\t' + i.info
        to_write.append(write_str)

    if os.path.isfile(args.out):
        with open(args.out, 'a') as outfile:
            outfile.write('\n' + '\n'.join(to_write))
    else:
        with open(args.out, 'w') as outfile:
            outfile.write(header + '\n' + '\n'.join(to_write))
