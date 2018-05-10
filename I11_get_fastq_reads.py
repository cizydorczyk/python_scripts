from sys import argv

script,sam_names = argv

names_list = []
with open(sam_names, 'r') as infile1:
    for line in infile1:
        names_list.append(line.strip())

print len(names_list)
print len(set(names_list))