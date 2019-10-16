import pandas
import collections
import argparse
import statistics

parser = argparse.ArgumentParser()

# General required arguments:
parser.add_argument("--raw_blast_input_file", help="raw blast output to parse")
parser.add_argument("--parsed_blast_output_file", help="parsed blast output file")

args = parser.parse_args()

def ParseBlastOutput(raw_blast_input_file, parsed_blast_output_file):
    """Function to parse raw contig blast output on Synergy."""

    contig_list = []
    no_result_contig_list = []
    query_list = []
    fill_color = []
    out_color = []
    contig_lengths = []
    contig_depths = []

    with open(args.raw_blast_input_file, 'r') as infile:
        contents = list(infile)

    for lnum, line in enumerate(contents):
        if "# Query" in line:
            if (lnum+4) < len(contents): # ensures we don't read past end of file

                line_elements = line.strip().split(" ")
                contig_num = line_elements[2]
                query_list.append(contig_num)

                contig_lengths.append(line_elements[3].split("=")[1])
                contig_depths.append(float(line_elements[4].split("=")[1][0:-1]))

                for line in contents:
                    if line.startswith(contig_num + '\t') and ("Haemophilus" in line):
                        fill_color.append("deepskyblue")
                        out_color.append("blue")
                        break

                else:
                    fill_color.append("red")
                    out_color.append("red4")

            else:
                query_list.append(line.strip().split(' ')[-1])
                fill_color.append("red")
                out_color.append("red4")

    data = collections.OrderedDict([('Contig', query_list), ('Length', contig_lengths), ('Coverage', contig_depths),
                                    ('Color', fill_color), ('Outline', out_color)])

    df = pandas.DataFrame(data)

    Q1 = df['Coverage'].quantile(0.25)
    Q3 = df['Coverage'].quantile(0.75)
    IQR = round(Q3 - Q1,2)
    MEAN = round(statistics.mean(df['Coverage']),2)
    MEDIAN = statistics.median(df['Coverage'])
    MIN = min(df['Coverage'])
    MAX = max(df['Coverage'])
    OUTLIER_MIN = round(Q1-1.5*IQR,2)
    OUTLIER_MAX = round(Q3+1.5*IQR,2)

    header1 = "#\tMEAN\tMEDIAN\tMIN\tMAX\tQ1\tQ3\tIQR\tQ1-1.5*IQR\tQ3+1.5*IQR"
    header2 = "\t".join(["#", str(MEAN), str(MEDIAN), str(MIN), str(MAX), str(Q1), str(Q3), str(IQR), str(OUTLIER_MIN), str(OUTLIER_MAX)])

    # Uncomment lines below to mark contigs as outliers based on traditional values (i.e. Q1-1.5*IQR, Q3+1.5*IQR).
    # Contig lengths almost certainly do not follow a normal distribution, however, so it is questionable how useful this would be...
    # df.loc[(df["Coverage"] > OUTLIER_MAX) & (df["Color"] == "deepskyblue"), "Color"] = "chartreuse3"
    # df.loc[df["Color"] == "chartreuse3", "Outline"] = "forestgreen"
    # df.loc[(df["Coverage"] < OUTLIER_MIN) & (df["Color"] == "deepskyblue"), "Color"] = "chartreuse3"
    # df.loc[df["Color"] == "chartreuse3", "Outline"] = "forestgreen"
    #
    # df.loc[(df["Coverage"] > OUTLIER_MAX) & (df["Color"] == "red"), "Color"] = "darkorange2"
    # df.loc[df["Color"] == "darkorange2", "Outline"] = "darkorange4"
    # df.loc[(df["Coverage"] < OUTLIER_MIN) & (df["Color"] == "red"), "Color"] = "darkorange2"
    # df.loc[df["Color"] == "darkorange2", "Outline"] = "darkorange4"

    outputfile = open(args.parsed_blast_output_file, 'a')
    outputfile.write(header1 + "\n")
    outputfile.write(header2 + "\n")
    df.to_csv(outputfile, sep="\t", header=True, index=False)
    outputfile.close()

ParseBlastOutput(args.raw_blast_input_file, args.parsed_blast_output_file)

##### Testing code outside function (for ease of testing) below #####
# contig_list = []
# no_result_contig_list = []
#
# with open("/home/conrad/hinfluenzae/assembly_pipeline_testing/blast/raw_blast_output/A058-H09-28-09-2012_raw_blast_output.txt", 'r') as infile:
#     contents = list(infile)
#
# query_list = []
# fill_color = []
# out_color = []
# contig_lengths = []
# contig_depths = []
#
# for lnum, line in enumerate(contents):
#     if "# Query" in line:
#         if (lnum+4) < len(contents):
#             contig_num = line.strip().split(' ')[2]
#             query_list.append(contig_num)
#
#             contig_lengths.append(line.strip().split(" ")[3].split("=")[1])
#             contig_depths.append(line.strip().split(" ")[4].split("=")[1][0:-1])
#
#             for line in contents:
#                 if line.startswith(contig_num + '\t') and ("Haemophilus" in line):
#                     fill_color.append("deepskyblue")
#                     out_color.append("blue")
#                     break
#
#             else:
#                 fill_color.append("red")
#                 out_color.append("red4")
#
#         else:
#             query_list.append(line.strip().split(' ')[-1])
#             fill_color.append("red")
#             out_color.append("red4")
#
# data = collections.OrderedDict([('Contig', query_list), ('Length', contig_lengths), ('Coverage', contig_depths),
#                                 ('Color', fill_color), ('Outline', out_color)])
#
# df = pandas.DataFrame(data)
#
# df.to_csv(outputfile, sep='\t', header=True, index=False)
