import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("--isolate_list", help="isolate list, one isolate per line")
parser.add_argument("--pirate_matrix", help="PIRATE binary presence/absence matrix")
parser.add_argument("--output_matrix", help="output per patient binary presence absence matrix", default="")

args = parser.parse_args()

pres_abs_df = pd.read_table(args.pirate_matrix, sep="\t", header=0, index_col=0)

patients_list = []
with open(args.isolate_list, 'r') as infile1:
    for line in infile1:
        patient = line.strip().split("-")[0]
        if patient not in patients_list:
            patients_list.append(patient)


df_columns = list(pres_abs_df.columns)
df_rows = list(pres_abs_df.index)

patient_pres_abs_dict = {}

for patient in patients_list:
    patient_columns = [i for i in df_columns if i.startswith(patient)]
    patient_df = pres_abs_df[patient_columns]

    patient_seq = []
    for index, row in patient_df.iterrows():
        row_list = list(row)
        if 1 in row_list:
            patient_seq.append(1)
        else:
            patient_seq.append(0)

    patient_pres_abs_dict[patient] = patient_seq

patient_pres_abs_df = pd.DataFrame.from_dict(data=patient_pres_abs_dict, orient='columns', dtype=int)
patient_pres_abs_df.index = df_rows

# Write output dict to file:
if args.output_matrix != "":
    print("Writing matrix to file...")
    patient_pres_abs_df.to_csv(path_or_buf=args.output_matrix, sep="\t")
else:
    print("No output file provided! Provide one using '--output_matrix <output matrix file path>'.")
