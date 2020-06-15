import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--genus", help="genus, capitalized!")
parser.add_argument("--species", help="species, lowercase!")
parser.add_argument("--fname")
parser.add_argument("--lname")
parser.add_argument("--email", help="can be dummy")
parser.add_argument("--organization")
parser.add_argument("--department")
parser.add_argument("--street")
parser.add_argument("--city")
parser.add_argument("--state")
parser.add_argument("--pcode")
parser.add_argument("--country")
parser.add_argument("--isolate")
parser.add_argument("--output_yaml")

args = parser.parse_args()

organism = args.genus + "_" + args.species

to_write = "organism:\n\tgenus_species: '" + organism + "'\n\tstrain: '" + args.isolate + "'\ncontact_info:\n\tlast_name: '" + \
            args.lname + "'\n\tfirst_name: '" + args.fname + "'\n\temail: '" + args.email + "'\n\torganization: '" + \
            args.organization + "'\n\tdepartment: '" + args.department + "'\n\tstreet: '" + args.street + "'\n\tcity: '" + \
            args.city + "'\n\tstate: '" + args.state + "'\n\tpostal_code: '" + args.pcode + "'\n\tcountry: '" + args.country + \
            "'\nauthors:\n\t- author:\n\t\tlast_name: '" + args.lname + "'\n\t\tfirst_name: '" + args.fname + "'"

with open(args.output_yaml, 'w') as outfile1:
    outfile1.write(to_write)
