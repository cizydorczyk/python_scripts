#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conradizydorczk@outlook.com>
Date   : 2024-01-30
Purpose: get fastq file lengths
"""

import argparse
import json


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='get fastq file lengths',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--multiqc_input',
                        help='multiqc json file',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output_file',
                        help='output tsv with read lengths',
                        metavar='str',
                        type=str,
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    read_lengths = {}
    
    with open(args.multiqc_input, "r") as infile1:
        data = json.load(infile1)
    
    data_sub = data["report_saved_raw_data"]["multiqc_fastqc"]

    for isolate in data_sub:
        isolate_ = isolate.split("_")[0]
        lengths = data_sub[isolate]['Sequence length']

        if isolate_ not in read_lengths:
            read_lengths[isolate_] = [lengths]
        else:
            read_lengths[isolate_].append(lengths)
    
    to_write = []
    for key, value in sorted(read_lengths.items()):
        output_string = f"{key}\t{value[0]}\t{value[1]}"
        to_write.append(output_string)

    with open(args.output_file, "w") as outfile1:
        outfile1.write("\n".join(to_write))


        # print(isolate, data_sub[isolate]['Sequence length'])
        

    # for sub_dict in data_sub:
    #     for key in sub_dict:
    #         key_sub = sub_dict[key]
    #         # output_string = f"{key}\t{key_sub['median_sequence_length']}"
    #         output_string = (key, key_sub['median_sequence_length'])
    #         read_lengths.append(output_string)
    
    # read_lengths_2 = {}
    # for i in sorted(read_lengths):
    #     isolate = i[0].split("_")[0]
    #     length = i[1]

    #     if isolate not in read_lengths_2: 
    #         read_lengths_2[isolate] = [length]
    #     else:
    #         read_lengths_2[isolate].append(length)
    
    # for key in data:
    #     print(key)


            


# --------------------------------------------------
if __name__ == '__main__':
    main()
