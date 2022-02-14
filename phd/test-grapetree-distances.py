import argparse

parser = argparse.ArgumentParser()

# General script options:
parser.add_argument("--profile_file")

args = parser.parse_args()

profiles_dict = {}
with open(args.profile_file, 'r') as infile1:
    for line in infile1:
        if not line.startswith("FILE"):
            line_elements = line.strip().split("\t")
            isolate = line_elements[0]
            alleles = line_elements[1:]

            profiles_dict[isolate] = alleles

A371_SM076_list = profiles_dict["A371-SM076-30-11-2015"]
A371_SM078_list = profiles_dict["A371-SM078-18-01-2016"]
A090_SM080_list = profiles_dict["A090-SM080-8-2-2016"]

A090_SM043_list = profiles_dict["A090-SM043-26-09-2012"]
A090_SM159_list = profiles_dict["A090-SM159-29-04-2013"]

dif_count2 = 0
sim_count2 = 0
for i, j in zip(A371_SM076_list, A090_SM080_list):
    if i == j:
        sim_count2 += 1
    elif i != j:
        dif_count2 += 1
print("--missing 2: ", sim_count2, dif_count2)
# --missing 2 calculates the p-distance as (number differences)/(number differences + number similarities)
# where number differences includes comparisons of allele vs. missing
# and number of similarities does as well

dif_count3 = 0
sim_count3 = 0
for i, j in zip(A371_SM076_list, A090_SM080_list):
    if i != "0" and j != "0":
        if i == j:
            sim_count3 += 1
        elif i != j:
            dif_count3 += 1
print("--missing 3: ", sim_count3, dif_count3)
# --missing 0 calculates p-distance as (number of differences)/(number of differences + number similarities)
# where number differences does NOT include comparisons of allele vs. missing
# and number of similarities does NOT include such comparisons as well
# --missing 3 simply outputs the number of differences that --missing 0 uses in its calculation

#############
# Take a look at ST-N2 isolates, which appear much more closely related using both distance calculations

A090_sim_count2 = 0
A090_dif_count2 = 0
for i, j in zip(A090_SM043_list, A090_SM159_list):
    if i == j:
        A090_sim_count2 += 1
    elif i != j:
        A090_dif_count2 += 1
print("--missing 2: ", A090_sim_count2, A090_dif_count2)

A090_sim_count3 = 0
A090_dif_count3 = 0
for i, j in zip(A090_SM043_list, A090_SM159_list):
    if i != "0" and j != "0":
        if i == j:
            A090_sim_count3 += 1
        elif i != j:
            A090_dif_count3 += 1
print("--missing 3: ", A090_sim_count3, A090_dif_count3)
