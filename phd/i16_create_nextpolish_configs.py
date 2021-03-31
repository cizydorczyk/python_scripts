import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--job_type", help='set to "local"')
parser.add_argument("--job_prefix", help='job prefix to use for outputs')
parser.add_argument("--task", help='algorithms to run; consult manual. default \
= "best"', default='best')
parser.add_argument("--rewrite", help='rewrite existing dir? yes or no; default\
= no', default='no')
parser.add_argument("--rerun", help='# reruns to run unfinished jobs until \
finished or run for this many loops; int; default=3', default='3')
parser.add_argument("--parallel_jobs", help='# tasks to run in parallel; req.')
parser.add_argument("--multithread_jobs", help='# threads per task; set to 1')
parser.add_argument("--genome", help='full path to assembly')
parser.add_argument("--polish_options", help='options; see manual; default = \
-p {multithread_jobs}', default='-p {multithread_jobs}')
parser.add_argument("--sgs_fofn", help='sgs file of file names')
parser.add_argument("--sgs_options", help='sgs options; see manual; default='\
'-max_depth 100 -bwa', default='-max_depth 100 -bwa')

parser.add_argument("--config_file", help='output config file')

args = parser.parse_args()


output_str = "[General]\njob_type = " + args.job_type + "\ntask = " + args.task\
+ "\nrewrite = " + args.rewrite + "\nrerun = " + args.rerun + "\nparallel_jobs\
 = " + args.parallel_jobs + "\nmultithread_jobs = " + args.multithread_jobs + \
 "\ngenome = " + args.genome + "\ngenome_size = auto\nworkdir = ./01_rundir" \
 + "\npolish_options = " + args.polish_options + "\n\n[sgs_option]\nsgs_fofn ="\
 + " " + args.sgs_fofn + "\nsgs_options = " + args.sgs_options

with open(args.config_file, 'w') as outfile:
    outfile.write(output_str)
