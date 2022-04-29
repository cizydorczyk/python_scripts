import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_blast_file", help="(outfmt 6 output blast file, columns qaccver saccver pident length mismatches gapopens qstart qend qlen sstart send slen evalue bitscore)")
parser.add_argument("--blast_groups", help="blast groups output file, listing pegs & groups")
parser.add_argument("--unique_pegs", help="output file with unique pegs, all non-group pegs + 1/group")

parser.add_argument("--pid", help="pid threshold; default = 70.0", default=70.0, type=float)
parser.add_argument("--perclen", help="perclen threshold; default = 70.0", default=70.0, type=float)

args = parser.parse_args()

class BlastHit(object):
    def __init__(self, qaccver, saccver, pident, length, mismatches, gapopens, qstart, qend, qlen, sstart, send, slen, evalue, bitscore, perclen, record):
        self.qaccver = qaccver
        self.saccver = saccver
        self.pident = pident
        self.length = length
        self.mismatches = mismatches
        self.gapopens = gapopens
        self.qstart = qstart
        self.qend = qend
        self.qlen = qlen
        self.sstart = sstart
        self.send = send
        self.slen = slen
        self.evalue = evalue
        self.bitscore = bitscore
        self.perclen = perclen
        self.record = record

def find_key(input_dict, value):
    return {k for k, v in input_dict.items() if value in v}

# test_dict = {1:"!", 2:"@", 3:"#", 4:"$", 5:"%"}
# key = find_key(test_dict, "@")
# if len(key) == 0:
#     print("No key found")
# else:
#     print(key)

groups_dict = {}
unique_pegs_list = []
groups_counter = 1
with open(args.input_blast_file, "r") as infile1:
    for line in infile1:

        le = line.strip().split("\t")
        perclen = round((int(le[3])/max(int(le[8]), int(le[11]))) * 100, 2)
        blast_hit = BlastHit(le[0], le[1], le[2], le[3], le[4], le[5], le[6], le[7], le[8], le[9], le[10], le[11], le[12], le[13], perclen, line.strip())

        if blast_hit.qaccver != blast_hit.saccver:
            if float(blast_hit.pident) > args.pid and blast_hit.perclen > args.perclen:
                key = find_key(groups_dict, blast_hit.qaccver)
                if len(key) == 0: # if fail to find query in any keys
                    # print("1")
                    key = find_key(groups_dict, blast_hit.saccver) # try finding subject
                    if len(key) == 0: # if subject not found, start new group
                        # print("2")
                        groups_dict[groups_counter] = [blast_hit.qaccver, blast_hit.saccver]
                        groups_counter += 1
                    elif len(key) == 1: # if subject found, append query to dict entry
                        # print("3")
                        groups_dict[list(key)[0]].append(blast_hit.qaccver)
                    elif len(key) > 1: # if multiple keys with subject found, print error and break
                        print("4")
                        print(key, "Something is weird...multiple groups with value found...")
                        break
                elif len(key) == 1: # if query found
                    # print("5")
                    key2 = find_key(groups_dict, blast_hit.saccver) # test if subject found
                    if len(key2) == 0: # if subject not found
                        # print("6")
                        groups_dict[list(key)[0]].append(blast_hit.saccver) # append subject to query dict entry
                    elif len(key2) == 1: # if subject found
                        # print("7")
                        if list(key)[0] == list(key2)[0]: # if query key == subject key
                            # print("8")
                            continue # do nothing; repeat entry for pair
                        elif key != key2: # if query != subject key; may happen if order of hits is weird...
                        # aka. 1-2, 3-4, 3-1; 1/2 grouped, 3/4 grouped, then 3/1 throw this error...not sure what sort order blast uses but this can technically occur...
                            print("9" + "\t" + "Multiple keys for pair...combining groups...\n\tCHECK THIS ENTRY MANUALLY\n\t" + blast_hit.qaccver + "\t" + blast_hit.saccver + "\n\tGroups " + str(key) + "\t" + str(key2))
                            groups_dict[list(key)[0]] = groups_dict[list(key)[0]] + groups_dict[list(key2)[0]]
                            del groups_dict[list(key2)[0]]
                    elif len(key2) > 1:
                        print("10")
                        print("Something is weird...multiple groups with subject found...")
                        break


        elif blast_hit.qaccver == blast_hit.saccver:
            unique_pegs_list.append(blast_hit.qaccver.split("|")[1])

# Write groups file:
to_write = []
for group in groups_dict:
    output_string = str(group) + "\t" + "\t".join(groups_dict[group])
    to_write.append(output_string)
with open(args.blast_groups, "w") as outfile1:
    outfile1.write("\n".join(to_write))

# Add one representative peg from each group to unique pegs list:
for group in groups_dict:
    unique_pegs_list.append(groups_dict[group][0].strip().split("|")[1])

# Write unique pegs to file:
with open(args.unique_pegs, "w") as outfile2:
    outfile2.write("\n".join(unique_pegs_list))
