import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--strain")
parser.add_argument("--pgap_yaml")
parser.add_argument("--submol_yaml")

parser.add_argument("--assembly")
parser.add_argument("--organism", help="'Genus species' in single quotes")

parser.add_argument("--topology", help="topology for ALL chromosomes. Default = linear", default="linear")
parser.add_argument("--location", help="chromosome or plasmid. default = chromosome", default="chromosome")

args = parser.parse_args()

# Create input yaml:
input_yaml_text = "fasta:\n    class: File\n    location: " + args.assembly + "\nsubmol:\n    class: File\n    location: " + args.submol_yaml

# Creat submol yaml:
submol_yaml_text = "topology: '" + args.topology + "'\nlocation: '" + args.location + "'\norganism:\n    genus_species: '" + args.organism + "'\n    strain: '" + args.strain + "'\ncontact_info:\n    last_name: 'Izydorczyk'\n    first_name: 'Conrad'\n    email: 'conrad.izydorczyk@ucalgary.ca'\n    organization: 'University of Calgary'\n    department: 'Microbiology Immunology and Infectious Diseases'\n    street: '3330 Hospital Drive NW'\n    city: 'Calgary'\n    state: 'AB'\n    postal_code: 'T2N 4N1'\n    country: 'Canada'\nauthors:\n    - author:\n        last_name: 'Izydorczyk'\n        first_name: 'Conrad'\n    - author:\n        last_name: 'Parkins'\n        first_name: 'Michael'"

# Write input yaml to file:
with open(args.pgap_yaml, "w") as outfile1:
    outfile1.write(input_yaml_text)

with open(args.submol_yaml, "w") as outfile2:
    outfile2.write(submol_yaml_text)



