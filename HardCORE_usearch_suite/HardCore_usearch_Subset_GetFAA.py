#! /usr/bin/python


from sys import argv

script, subset_file, AllProteinList = argv

subset = open(subset_file, "r")
AllProteins = open(AllProteinList, "r")

# First read the AllProteins file into a dictionary where the key is >Strain_#####[Strain]
SEQS = {}
first = 0
SEQS_key = ""
SEQS_value = ""
for line in AllProteins:
    if line.startswith(">"):
        if first != 0:
            SEQS[SEQS_key]= SEQS_value
            SEQS_key = ""
            SEQS_value = ""
        SEQS_key = ">" + line.split(" ")[1] + line.split(" ")[2]
        SEQS_key = SEQS_key.rstrip("\n")
        first = 1
    else:
        SEQS_value += line
    SEQS[SEQS_key] = SEQS_value
AllProteins.close()

# Read the subset file into a dict as well
SUBSET = {}
for line in subset:
    line = line.rstrip("\n")
    line = line.rstrip("\t")
    line_list = line.split('\t')
    for x in line_list:
        SUBSET[x] = "NOT_DONE"
subset.close()

out_file = open(subset_file+".faa", "w")
# Now use the subset file and cross ref with the SEQS hash to get the sequences
for key, value in SUBSET.items():
    if SEQS[key]:
        to_write = key + "\n" + SEQS[key]
        out_file.write(to_write)
out_file.close()