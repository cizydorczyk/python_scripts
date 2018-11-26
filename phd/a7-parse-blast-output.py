from sys import argv
import pandas
import collections

script, inputblastfile, outputfile = argv

contig_list = []
no_result_contig_list = []

with open(inputblastfile, 'r') as infile:
    contents = list(infile)

print("Working on file: " + str(inputblastfile.split('/')[-1]))
# data = pandas.read_table(inputcontigstats, sep="\t")
# length_list = list(data["Consensus length"])
# coverage_list = list(data["Average coverage"])

query_list = []
fill_color = []
out_color = []
contig_lengths = []
contig_depths = []

for lnum, line in enumerate(contents):
    if "# Query" in line:
        if (lnum+4) < len(contents):
            contig_num = line.strip().split(' ')[2]
            query_list.append(contig_num)

            contig_lengths.append(line.strip().split(" ")[3].split("=")[1])
            contig_depths.append(line.strip().split(" ")[4].split("=")[1])

            for line in contents:
                if line.startswith(contig_num + '\t') and ("Escherichia" in line or "Shigella" in line):
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

data = collections.OrderedDict([('Filename', query_list), ('Length', contig_lengths), ('Coverage', contig_depths),
                                ('Color', fill_color), ('Outline', out_color)])

df = pandas.DataFrame(data)

df.to_csv(outputfile, sep='\t', header=True, index=False)
