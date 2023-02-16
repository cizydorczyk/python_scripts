#!/usr/bin/env python3
"""
Author : conrad <conrad@localhost>
Date   : 2022-12-14
Purpose: Get PA isolates to sequence to fill in gaps. Jan 2023
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Get isolates to sequence. Dec 2022',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input_file',
                        help='input file',
                        metavar='str',
                        type=str,
                        )

    parser.add_argument('-o',
                        '--output_file',
                        help='output file',
                        metavar='str',
                        type=str,
                        )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Get isolates to sequence. Dec 2022"""

    args = get_args()
    
    patient_years_dict_even = {}
    patient_years_2017_2020 = {}

    even_year_isolates_to_seq = []
    isolates_to_seq_2017_2020 = []

    with open(args.input_file, "r") as infile1:
        for line in infile1:
            if not line.startswith("PID"):
                line_elements = line.strip().split("\t")
                patient_anumber = line_elements[1].split("_")[0]
                patient_year = line_elements[2].split("-")[0]
                
                ## we want every other year isolates from 2002-2016
                # if year is even and 2016 or earlier (min 2002):
                if int(patient_year) % 2 == 0 and int(patient_year) <= 2016:
                    # if patient not already in dict:
                    if patient_anumber not in patient_years_dict_even:
                        # add patient to dict, add current year to patient-specific year dict, and current line (biobank entry) to year dict
                        patient_years_dict_even[patient_anumber] = {patient_year:[line]}
                    # else if patient already in dict:
                    else:
                        # if current year not in dict:
                        if patient_year not in patient_years_dict_even[patient_anumber]:
                            # add current year to patient-specific year dict and current line (biobank entry to year dict
                            patient_years_dict_even[patient_anumber][patient_year] = [line]
                        # if current line not in dict but patient and current year are:
                        if line not in patient_years_dict_even[patient_anumber][patient_year]:
                            # append current line to current year dict in current patient dict:
                            patient_years_dict_even[patient_anumber][patient_year].append(line)
                
                ## we want yearly isolates for 2017-2020
                # else if year is greater than 2016 (i.e. 2017-2020):
                elif int(patient_year) > 2016:
                    # if patient not in dict:
                    if patient_anumber not in patient_years_2017_2020:
                        # add patient to dict, add current year to patient-specific year dict, and current line (biobank entry) to year dict
                        patient_years_2017_2020[patient_anumber] = {patient_year:[line]}
                    # else if patient already in dict:
                    else:
                        # if current year not in dict:
                        if patient_year not in patient_years_2017_2020[patient_anumber]:
                            # add current year to patient-specific year dict and current line (biobank entry to year dict
                            patient_years_2017_2020[patient_anumber][patient_year] = [line]
                        # if current line not in dict but patient and current year are:
                        if line not in patient_years_2017_2020[patient_anumber][patient_year]:
                            # append current line to current year dict in current patient dict:
                            patient_years_2017_2020[patient_anumber][patient_year].append(line)
    
    # go through 2002-2016 patients/years and select isolates from even years where none exist:
    for patient in patient_years_dict_even:
        for year in patient_years_dict_even[patient]:
            print(f"Checking patient {patient} and year {year}...")
            wgs_codes = []
            pfge_codes = []
            isolate_dates = []
            isolate_lines = []
            for isolate in patient_years_dict_even[patient][year]:
                isolate_elements = isolate.strip().split("\t")
                wgs_codes.append(int(isolate_elements[5]))
                pfge_codes.append(int(isolate_elements[4]))
                isolate_dates.append(isolate_elements[2])
                isolate_lines.append(isolate)
            
            if any(wgs_codes):
                print(f"\tWGS present, moving on to next year/patient.")
                continue
            else:
                print(f"\tWGS missing, checking PFGE...\n\tPFGE presence/absence codes for {year} are: {pfge_codes}")
                if any(pfge_codes):
                    pfge_index = pfge_codes.index(1)
                    isolate_to_seq = isolate_dates[pfge_index]
                    isolate_line_to_seq = isolate_lines[pfge_index]
                    print(f"\tPFGE present, selecting first isolate with PFGE to sequence...\n\tIsolate to sequence is ---> {isolate_to_seq} <--- out of possible isolates {isolate_dates}.")
                    even_year_isolates_to_seq.append(isolate_line_to_seq)
                else:
                    print(f"\tPFGE missing...\n\tSelecting 1st available isolate for sequencing ---> {isolate_dates[0]} <--- out of available isolates {isolate_dates}.")
                    even_year_isolates_to_seq.append(isolate_lines[0])

    # go through 2017-2020 patients/years and select yearly isolates where none exist:
    for patient in patient_years_2017_2020:
        for year in patient_years_2017_2020[patient]:
            print(f"Checking patient {patient} and year {year}...")
            wgs_codes = []
            pfge_codes = []
            isolate_dates = []
            isolate_lines = []
            for isolate in patient_years_2017_2020[patient][year]:
                isolate_elements = isolate.strip().split("\t")
                wgs_codes.append(int(isolate_elements[5]))
                pfge_codes.append(int(isolate_elements[4]))
                isolate_dates.append(isolate_elements[2])
                isolate_lines.append(isolate)
            
            if any(wgs_codes):
                print(f"\tWGS present, moving on to next year/patient.")
                continue
            else:
                print(f"\tWGS missing, checking PFGE...\n\tPFGE presence/absence codes for {year} are: {pfge_codes}")
                if any(pfge_codes):
                    pfge_index = pfge_codes.index(1) # index of first pfge positive isolate
                    isolate_to_seq = isolate_dates[pfge_index]
                    isolate_line_to_seq = isolate_lines[pfge_index]
                    print(f"\tPFGE present, selecting first isolate with PFGE to sequence...\n\tIsolate to sequence is ---> {isolate_to_seq} <--- out of possible isolates {isolate_dates}.")
                    isolates_to_seq_2017_2020.append(isolate_line_to_seq)
                else:
                    print(f"\tPFGE missing...\n\tSelecting 1st available isolate for sequencing ---> {isolate_dates[0]} <--- out of available isolates {isolate_dates}.")
                    isolates_to_seq_2017_2020.append(isolate_lines[0])


    header = f"PID\tADATE\tCULTDAT\tvalues\tPFGE_Isolates\tSequenced_Isolates\n"      
    to_write = sorted(even_year_isolates_to_seq + isolates_to_seq_2017_2020)
    to_write = "".join(to_write)
    print("Writing output file...")
    with open(args.output_file, "w") as outfile1:
        outfile1.write(header + to_write)


    

# --------------------------------------------------
if __name__ == '__main__':
    main()
