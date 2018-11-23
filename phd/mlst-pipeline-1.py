import argparse
import subprocess
import os
from Bio.Seq import Seq

parser = argparse.ArgumentParser()

##################################################################################
#  .-----------------. .----------------.  .----------------.  .----------------.
# | .--------------. || .--------------. || .--------------. || .--------------. |
# | | ____  _____  | || |     ____     | || |  _________   | || |  _________   | |
# | ||_   \|_   _| | || |   .'    `.   | || | |  _   _  |  | || | |_   ___  |  | |
# | |  |   \ | |   | || |  /  .--.  \  | || | |_/ | | \_|  | || |   | |_  \_|  | |
# | |  | |\ \| |   | || |  | |    | |  | || |     | |      | || |   |  _|  _   | |
# | | _| |_\   |_  | || |  \  `--'  /  | || |    _| |_     | || |  _| |___/ |  | |
# | ||_____|\____| | || |   `.____.'   | || |   |_____|    | || | |_________|  | |
# | |              | || |              | || |              | || |              | |
# | '--------------' || '--------------' || '--------------' || '--------------' |
#  '----------------'  '----------------'  '----------------'  '----------------'

# All directories must be specified as full paths!

# Caveats:
#   The pipeline may fail if there are multiple hits per gene in an assembly. Haven't tested it out.
#   The pipeline may fail if none of the existing genes are similar enough to novel genes; can add the '--type blastn'
#       option to the blast command below to use more relaxed blasting criteria. Usually this shows up as truncations
#       in an MLST gene...
#   I should add a message in the help menu about requiring absolute paths...

###################################################################################

# Required arguments:
parser.add_argument("input_dir", help="Directory with input assemblies")
parser.add_argument("project_dir", help="Directory for all output files. Will be created if doesn't exist.")

# Optional arguments:
parser.add_argument("--mlsa", help="Create MLSA/MLST in addition to TS mlst?", action="store_true")
parser.add_argument("--blastdbdir", default=None, help="Directory with blast databases (1 db/MLST gene) for \
your species. DBs should be named with simply the gene name and nothing else, e.g. 'adk'. These MUST be made prior to \
 running this script.", action="store")
parser.add_argument("--mlst_genes", default=None, help="Comma-separated list of MLST genes for your species. \
E.g. adk,fumc,gyrb,icd,mdh,pura,reca .")

args = parser.parse_args()

##############
# Check if Torsten Seemann's mlst program is in path:
try:
    subprocess.run(['mlst', '--check'], check=True)
except subprocess.CalledProcessError as err:
    print('ERROR:', err)

# Check if BLAST is in path:
try:
    subprocess.run(['blastn', '-version'], check=True)
except subprocess.CalledProcessError as err:
    print('ERROR:', err)

# Check if input dir exists:
inputdirexists = os.path.isdir(args.input_dir)
if not inputdirexists:
    raise OSError("input directory does not exist: %s" % args.input_dir)
else:
    print('Input directory:', args.input_dir)

# Check if ouput dir exists and create if not:
outputdirexists = os.path.isdir(args.project_dir)
if outputdirexists:
    print("Project directory: %s" % args.project_dir)
else:
    os.mkdir(args.project_dir)
    print("Project directory created: %s" % args.project_dir)

# Get list of files in input directory and remove any non-fasta files:
input_file_list = os.listdir(args.input_dir)
# print(input_file_list)
for file_ in input_file_list:
    if not file_.lower().endswith(('.fa', '.fasta', '.fna')):
        input_file_list.remove(file_)

input_fasta_files = []
for fa in input_file_list:
    fasta_ = os.path.join(args.input_dir, fa)
    input_fasta_files.append(fasta_)

# Create output file suffixes (isolate numbers):
isolate_numbers = [i.split('.')[0] for i in input_file_list]


##############
# Get STs using Torsten Seemann's mlst program:
print("Getting STs from Torsten Seemann's mlst...")

output_ts_mlst = []

for fasta_file in sorted(input_fasta_files):
    mlst_completed = subprocess.run(['mlst', '--quiet', '--nopath', fasta_file], stdout=subprocess.PIPE)
    output_ts_mlst.append(mlst_completed.stdout.decode('utf-8'))

# print(''.join(output_ts_mlst))
print("mlst done.")

##############
# Perform MLSA/MLST:

if args.mlsa is True:

    print("Performing MLSA/MLST...")

    # Check if blastdb directory and databases exist:
    # Blast will complain if dbs don't exist, so not doing it here...

    blastdb_dir = os.path.isdir(args.blastdbdir)
    if not blastdb_dir:
        raise OSError("Blastdb directory does not exist: %s" % args.blastdbdir)
    else:
        print('Blastdb directory:', args.blastdbdir)

    # Get paths to Blast gene dbs:

    genes_list = args.mlst_genes.split(',')
    gene_dbs = [os.path.join(args.blastdbdir, gene) for gene in genes_list]

    # Blast isolate assemblies agaisnt mlst gene dbs:

    for fasta in input_fasta_files:

        isolate = fasta.split('/')[-1].split('.')[0]
        print("Blasting isolate " + isolate + "...")

        os.mkdir(os.path.join(args.project_dir, isolate))

        for db in gene_dbs:
            gene = db.split('/')[-1]
            subprocess.run(['blastn', '-db', db, '-query', fasta, '-out',
                            os.path.join(args.project_dir, isolate, (isolate + "." + gene + ".txt")), '-outfmt',
                            '6 sseqid sstrand qseqid qseq', '-evalue', '0.00001', '-max_target_seqs', '1',
                            '-num_threads', '8'])

    isolate_out_dirs = os.listdir(args.project_dir)

    print("Generating MLSA fasta file...")

    for isolate in isolate_out_dirs:
        isolate_ = isolate.split('_')[0]

        # Combined genes list contains lines from each isolate's blast results for each gene:
        combined_genes_list = []
        gene_files = os.listdir(os.path.join(args.project_dir, isolate))

        for file_ in sorted(gene_files):
            with open(os.path.join(args.project_dir, isolate, file_), 'r') as infile1:
                for line in infile1:
                    combined_genes_list.append(line)

        header = []
        sequence = []
        # combined_genes_list.sort(key=lambda x: x[1])

        for line_ in sorted(combined_genes_list):
            split_line = line_.strip().split('\t')

            header.append((split_line[0] + "," + split_line[2]))
            if split_line[1] == 'minus':
                sequence_ = Seq(split_line[3])
                sequence.append(str(sequence_.reverse_complement()))
                # print(split_line[0], sequence_.reverse_complement()[1:10])
            else:
                sequence.append(split_line[3])
                # print('plus strand', split_line[0], split_line[3][1:10])

        with open(os.path.join(args.project_dir, isolate, (isolate + ".concat_mlst_genes.fasta")), 'w') as outfile2:
            outfile2.write(">" + '|'.join(header) + '\n' + ''.join(sequence))

        with open(os.path.join(args.project_dir, isolate, (isolate + ".combined_mlst_genes.fasta")), 'w') as outfile1:
            to_write = []
            for line in combined_genes_list:
                split_line = line.strip().split('\t')
                if split_line[1] == 'minus':
                    sequence_ = Seq(split_line[3]).reverse_complement()
                else:
                    sequence_ = split_line[3]
                to_write.append(">" + split_line[0] + "|" + split_line[2] + '\n' + str(sequence_) + '\n')
            outfile1.write(''.join(to_write))

        with open(os.path.join(args.project_dir, "mlsa_alignment.fasta"), 'a+') as outfile3:
            outfile3.write(">" + str(isolate_) + '|' + '|'.join(header) + '\n' + ''.join(sequence) + '\n')

# Write TS's mlst results to file here so as to avoid problems navigating directories when creating 'gene_files' list...
with open(os.path.join(args.project_dir, "TS_mlst_results.txt"), 'w') as outfile5:
    outfile5.write(''.join(output_ts_mlst))

print("MLST pipeline finished.")

####################################################################################
#  .----------------.  .----------------.  .----------------.  .-----------------. #
# | .--------------. || .--------------. || .--------------. || .--------------. | #
# | |     ______   | || |     _____    | || |     ____     | || | ____  _____  | | #
# | |   .' ___  |  | || |    |_   _|   | || |   .'    `.   | || ||_   \|_   _| | | #
# | |  / .'   \_|  | || |      | |     | || |  /  .--.  \  | || |  |   \ | |   | | #
# | |  | |         | || |      | |     | || |  | |    | |  | || |  | |\ \| |   | | #
# | |  \ `.___.'\  | || |     _| |_    | || |  \  `--'  /  | || | _| |_\   |_  | | #
# | |   `._____.'  | || |    |_____|   | || |   `.____.'   | || ||_____|\____| | | #
# | |              | || |              | || |              | || |              | | #
# | '--------------' || '--------------' || '--------------' || '--------------' | #
#  '----------------'  '----------------'  '----------------'  '----------------'  #
####################################################################################
