import argparse
import subprocess
import os
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from progress.bar import Bar


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
# blastn and mlst must be in path!!

# This pipeline uses Torsten Seemann's mlst program to get STs. It can also create a multiple sequence alignment of MLST
#   loci, if given a fasta file for each gene in the MLST scheme containing all (or as many as you want) alleles.
#   In short, this pipeline blasts alleles against an isolate database, keeping one hit per allele. It then takes the
#   best hit based on the Blast SCORE and adds that into a multiple sequence alignment.
#   Obviously there may be situations where this fails..but in simple cases, it should work.
#
# The pipeline now (should) handle isolates with missing/truncated alleles. It also has progress bars...

# Caveats:
#   The pipeline may fail if there are multiple hits per gene in an assembly. Haven't tested it out.
#   The pipeline may fail if none of the existing genes are similar enough to novel genes; can add the '--type blastn'
#       option to the blast command below to use more relaxed blasting criteria. Usually this shows up as truncations
#       in an MLST gene...
#   I should add a message in the help menu about requiring absolute paths...
#   A gene split over 2 contigs should be okay now...haven't tested it though, and don't know why this would happen in
#       the first place...

###################################################################################

# Required arguments:
parser.add_argument("input_dir", help="Directory with input assemblies")
parser.add_argument("project_dir", help="Directory for all output files. Will be created if doesn't exist.")

# Optional arguments:
parser.add_argument("--mlsa", help="Create MLSA/MLST in addition to TS mlst?", action="store_true")
# parser.add_argument("--blastdbdir", default=None, help="Directory with blast databases (1 db/MLST gene) for \
# your species. DBs should be named with simply the gene name and nothing else, e.g. 'adk'. These MUST be made prior to \
#  running this script.", action="store")
parser.add_argument("--mlst_genes", default=None, help="Comma-separated list of MLST genes for your species. \
E.g. adk,fumc,gyrb,icd,mdh,pura,reca .")
parser.add_argument("--mlst_alleles_directory", help="Directory containing fasta files with all alleles per gene. 1 file"
                                                     " per gene, named simply 'gene.fasta'.")
parser.add_argument("--scheme", default='', help="TS's mlst scheme to use (optional - otherwise it uses autodetect).")

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
# for file_ in input_file_list:
#     if not file_.lower().endswith(('.fa', '.fasta', '.fna')):
#         input_file_list.remove(file_)
# print(input_file_list)

# Get full paths to fasta files:
input_fasta_files = []
for fa in input_file_list:
    fasta_ = os.path.join(args.input_dir, fa)
    input_fasta_files.append(fasta_)

# Create output file suffixes (isolate numbers):
isolate_numbers = [i.split('.')[0] for i in input_file_list]


##############
# Get STs using Torsten Seemann's mlst program:
print("Getting STs from Torsten Seemann's mlst...")

bar = Bar('Progress:', max=len(input_fasta_files), width=100)

output_ts_mlst = []

# for i in range(mlst_iter):
    # Do some work
for fasta_file in sorted(input_fasta_files):
    if args.scheme == '':
        mlst_completed = subprocess.run(['mlst', '--threads', '8', '--quiet', '--nopath', fasta_file],
                                        stdout=subprocess.PIPE)
        output_ts_mlst.append(mlst_completed.stdout.decode('utf-8'))
    else:
        mlst_completed = subprocess.run(
            ['mlst', '--threads', '8', '--scheme', args.scheme, '--quiet', '--nopath', fasta_file],
            stdout=subprocess.PIPE)
        output_ts_mlst.append(mlst_completed.stdout.decode('utf-8'))
    bar.next()
bar.finish()

print("mlst done.")

##############
# Perform MLSA/MLST:

if args.mlsa is True:

    print("Performing MLSA...")

    print("Creating isolate Blast DBs...")

    for isolate_fasta in input_fasta_files:
        subprocess.run(['makeblastdb', '-in', isolate_fasta, '-dbtype', 'nucl', '-out',
                        os.path.join(args.input_dir, isolate_fasta.split('.')[0])])

    print("Performing MLSA Blast Analysis...")

    bar2 = Bar("Progress:", max=len(isolate_numbers), width=100)

    mlst_genes = args.mlst_genes.split(',')
    for isolate in isolate_numbers:

        # Make directory for each isolate for intermediate output files:
        os.mkdir(os.path.join(args.project_dir, isolate))

        for gene in mlst_genes:
            subprocess.run(['blastn', '-db', os.path.join(args.input_dir, isolate), '-query',
                            os.path.join(args.mlst_alleles_directory, (gene + ".fasta")), '-out',
                            os.path.join(args.project_dir, isolate, (isolate + "." + gene + ".xml")), '-outfmt',
                            '5', '-num_threads', '8', '-max_hsps', '1'])

        # print("Parsing blast results for isolate %s" % isolate)

        isolate_blast_files = os.listdir(os.path.join(args.project_dir, isolate))

        combined_genes_list = []

        for file_ in sorted(isolate_blast_files):
            # print("\t" + file_)
            result_handle = open(os.path.join(args.project_dir, isolate, file_), 'r')
            blast_records = NCBIXML.parse(result_handle)
            best_score = 0
            top_hit = ''
            top_hit_seq = ''
            contig = ''
            for blast_record in blast_records:
                if blast_record.alignments:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.score > best_score:
                                best_score = hsp.score
                                top_hit = blast_record.query
                                top_hit_seq = hsp.sbjct

                                # Test for truncations/indels:
                                gene_length = blast_record.query_length
                                if len(str(top_hit_seq)) < gene_length:
                                    top_hit = top_hit + "_TRUNCATED"

                                contig = alignment.hit_def.split(' ')[0]
                elif not blast_record.alignments and best_score is 0:
                    top_hit = blast_record.query + "_MISSING"
                    contig = "NA"

            combined_genes_list.append(str(top_hit + "\t" + contig + "\t" + top_hit_seq))

        header = []
        sequence = []

        for line_ in sorted(combined_genes_list):
            split_line = line_.strip().split('\t')
            if len(split_line) == 3:
                header.append((split_line[0] + "," + split_line[1]))
                sequence.append(split_line[2])
            elif len(split_line) == 2:
                header.append((split_line[0] + "," + split_line[1]))

        with open(os.path.join(args.project_dir, isolate, (isolate + ".concat_mlst_genes.fasta")), 'w') as outfile2:
            outfile2.write(">" + '|'.join(header) + '\n' + ''.join(sequence))

        with open(os.path.join(args.project_dir, isolate, (isolate + ".combined_mlst_genes.fasta")), 'w') as outfile1:

            to_write = []
            for line in combined_genes_list:
                split_line = line.strip().split('\t')
                if len(split_line) == 3:
                    to_write.append(">" + split_line[0] + "|" + split_line[1] + '\n' + split_line[2] + '\n')
                elif len(split_line) == 2:
                    to_write.append(">" + split_line[0] + "|" + split_line[1] + '\n')

            outfile1.write(''.join(to_write))

        with open(os.path.join(args.project_dir, "mlsa_alignment.fasta"), 'a+') as outfile3:
            outfile3.write(">" + isolate + "\t" + '|'.join(header) + '\n' + ''.join(sequence) + '\n')

        bar2.next()
bar2.finish()

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
