#!/usr/bin/env python3
"""
Author : conrad <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-09-27
Purpose: Concatenate results from CARD, AMRFinderPlus, and ResFinder.
"""

import argparse
import pandas as pd


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Concatenate results from CARD, AMRFinderPlus, and ResFinder.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-iso',
                        '--isolate_list',
                        help='Isolate list containing all isolates for which AMR finding was run. One isolate per line.',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-os',
                        '--output_summary',
                        help='Output tsv file summarizing results with presence/absence of genes recorded. Short format output file.',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-od',
                        '--output_detailed',
                        help='Output tsv file containing detailed summary results; more detailed version of --output_summary.',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-afp',
                        '--amrfinderplus_path',
                        help='Path to directory containing AMRFinderPlus output files.',
                        metavar='str',
                        type=str,
                        required=True)

    parser.add_argument('-afs',
                        '--amrfinderplus_suffix',
                        help='Name of AMRFinderPlus output files not including the isolate name found in the isolate list. E.g. JMB00935-amrfinderplus-output.txt ==> the suffix is "-amrfinderplus-output.txt".',
                        metavar='str',
                        type=str,
                        required=True)

    parser.add_argument('-cp',
                        '--card_path',
                        help='Path to directory containing CARD output files.',
                        metavar='str',
                        type=str,
                        required=True)

    parser.add_argument('-cs',
                        '--card_suffix',
                        help='Name of CARD output files not including the isolate name found in the isolate list. E.g. JMB00935-card-output.txt ==> the suffix is "-card-output.txt".',
                        metavar='str',
                        type=str,
                        required=True)

    parser.add_argument('-rfp',
                        '--resfinder_path',
                        help='Path to directory containing ResFinder output files.',
                        metavar='str',
                        type=str,
                        required=True)

    parser.add_argument('-rfs',
                        '--resfinder_suffix',
                        help='Name of ResFinder output files not including the isolate name found in the isolate list. E.g. JMB00935-resfinder-output.txt ==> the suffix is "-resfinder-output.txt".',
                        metavar='str',
                        type=str,
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    isolate_list = args.isolate_list
    output_summary = args.output_summary
    output_detailed = args.output_detailed
    amrfinderplus_path = args.amrfinderplus_path
    amrfinderplus_suffix = args.amrfinderplus_suffix
    card_path = args.card_path
    card_suffix = args.card_suffix
    resfinder_path = args.resfinder_path
    resfinder_suffix = args.resfinder_suffix

    print(f'\nRunning script with the following arguments:\n')
    print(f'\tisolate list = "{isolate_list}"')
    print(f'\toutput summary file = "{output_summary}"')
    print(f'\toutput detailed file = "{output_detailed}"')
    print(f'\tamrfinderplus output directory path = "{amrfinderplus_path}"')
    print(f'\tamrfinderplus suffix = "{amrfinderplus_suffix}"')
    print(f'\tcard output directory path = "{card_path}"')
    print(f'\tcard suffix = "{card_suffix}"')
    print(f'\tresfinder output directory path = "{resfinder_path}"')
    print(f'\tresfinder suffix = "{resfinder_suffix}"\n')

    # # Read isolate list:
    # with open(isolate_list, "r") as infile1:
    #     isolates = [line.strip() for line in infile1]
    
    parseCard("/home/conrad/sa-ncfb/sa11-card/JMB00950-card-output.txt")

# --------------------------------------------------
def parseCard(card_file):
    """ Parse CARD output file & get relevant info."""

    with open(card_file, "r") as infile2:
        for line in infile2:
            if not line.startswith("ORF_ID"):
                line_elements = line.strip().split("\t")
                gene = line_elements[8].lower()

                print(gene, perc_coverage, gene_length, nudge, entry_note)


# --------------------------------------------------
if __name__ == '__main__':
    main()
