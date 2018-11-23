import argparse
import subprocess
import os
from Bio.Blast import NCBIXML

parser = argparse.ArgumentParser()

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

# Get list of files in input directory and remove any non-fasta files (not working for some reason...):
input_file_list = os.listdir(args.input_dir)
# print(input_file_list)
# for file_ in input_file_list:
#     if not (file_.lower().endswith(('.fa', '.fasta', '.fna'))):
#         input_file_list.remove(file_)
# print(input_file_list)

# Get full paths to input fasta files:
input_fasta_files = []
for fa in input_file_list:
    fasta_ = os.path.join(args.input_dir, fa)
    input_fasta_files.append(fasta_)

# Create output file suffixes (isolate numbers):
isolate_numbers = [i.split('.')[0] for i in input_file_list]


##############
# Get STs from Torsten Seemann's mlst program:
print("Getting STs from Torsten Seemann's mlst...")

# output_ts_mlst = []
#
# for fasta_file in sorted(input_fasta_files):
#     mlst_completed = subprocess.run(['mlst', '--quiet', '--nopath', fasta_file], stdout=subprocess.PIPE)
#     output_ts_mlst.append(mlst_completed.stdout.decode('utf-8'))
#
# print(''.join(output_ts_mlst))

##############
# Perform MLSA/MLST:

if args.mlsa is True:

    print("Performing MLSA/MLST...")

    # Check if blastdb directory and databases exist:

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

        os.mkdir(os.path.join(args.project_dir, isolate))

        for db in gene_dbs:
            gene = db.split('/')[-1]
            subprocess.run(['blastn', '-db', db, '-query', fasta, '-out',
                            os.path.join(args.project_dir, isolate, (isolate + "." + gene + ".xml")), '-outfmt',
                            '5', '-evalue', '0.00001', '-max_target_seqs', '1', '-num_threads', '8'])

    for isolate in isolate_numbers:
        gene_files = os.listdir(os.path.join(args.project_dir, isolate))

        for file_ in sorted(gene_files):
            result_handle = open(os.path.join(args.project_dir, isolate, file_), 'r')
            blast_records_generator = NCBIXML.parse(result_handle)
            for blast_record in blast_records_generator:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        print(alignment, hsp.strand)

        #     with open(os.path.join(args.project_dir, isolate, file_), 'r') as infile1:
        #         for line in infile1:
        #             pass

    # result_handle = open(blastxml, 'r')
    #
    # # blast_record = NCBIXML.read(result_handle)
    # # for i in blast_record.alignments:
    # #     for hsp in i.hsps:
    # #         print len(hsp.sbjct)
    #
    # blast_records = NCBIXML.parse(result_handle)
    # list_records = list(blast_records)
    #
    # sequences_list = []
    # header_list = []
    # for blast_record in list_records:
    #     for record in blast_record.alignments:
    #         header = ">" + record.title.split(" ")[0]
    #         for hsp in record.hsps:
    #             sequences_list.append(hsp.sbjct)




    # isolate_out_dirs = os.listdir(args.project_dir)
    #
    # for isolate in isolate_out_dirs:
    #     isolate_ = isolate.split('_')[0]
    #
    #     # Combined genes list contains lines from each isolate's blast results for each gene:
    #     combined_genes_list = []
    #     gene_files = os.listdir(os.path.join(args.project_dir, isolate))

        # for file_ in sorted(gene_files):
        #     with open(os.path.join(args.project_dir, isolate, file_), 'r') as infile1:
        #         for line in infile1:
        #             combined_genes_list.append(line)
        # with open(os.path.join(args.project_dir, isolate, (isolate + ".combined_mlst_genes.fasta")), 'w') as outfile1:
        #     to_write = []
        #     for line_ in combined_genes_list:
        #         to_write.append(">" + line_.strip().split('\t')[0] + "|" + line_.strip().split('\t')[1] + '\n' +
        #                         line_.strip().split('\t')[2] + '\n')
        #     outfile1.write(''.join(to_write))
        #
        # header = []
        # sequence = []
        # # combined_genes_list.sort(key=lambda x: x[1])
        #
        # for line_ in sorted(combined_genes_list):
        #     header.append(line_.strip().split('\t')[0] + "," + line_.strip().split('\t')[1])
        #     sequence.append(line_.strip().split('\t')[2])
        #
        # with open(os.path.join(args.project_dir, isolate, (isolate + ".concat_mlst_genes.fasta")), 'w') as outfile2:
        #     outfile2.write(">" + '|'.join(header) + '\n' + ''.join(sequence))
        #
        # with open(os.path.join(args.project_dir, "mlsa_alignment.fasta"), 'a+') as outfile3:
        #     outfile3.write(">" + str(isolate_) + '|' + '|'.join(header) + '\n' + ''.join(sequence) + '\n')
