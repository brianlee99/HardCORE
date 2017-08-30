#! /usr/bin/python

from sys import argv
import os
import glob
import subprocess

# This script will take in a FAA file with the following format for a header: >SM01_00816[SM01]
# The user will need to supply the directory where the FFNs are located.

# This will go to the FFN dir, concat all the FFNs and move the AllFFN.all file to the working dir
script, subset_file, path_to_ffns = argv
subset = open(subset_file, "r")
home = os.getcwd()
os.chdir(path_to_ffns)
process = subprocess.Popen("cat *.ffn > AllFFN.all", shell=True)
process.wait()
process = subprocess.Popen("mv AllFFN.all " + home, shell=True)
process.wait()

# Read the AllFFN.all into a hash
os.chdir(home)
ffn_file = open("AllFFN.all", "r")
SEQS = {}
first = 0
SEQS_key = ""
SEQS_value = ""
for line in ffn_file:
    if line.startswith(">"):
        if first != 0:
            SEQS[SEQS_key]= SEQS_value
            SEQS_key = ""
            SEQS_value = ""
        SEQS_key = line.split(" ")[0] + "[" + line.split("_")[0].lstrip(">") + "]"
        first = 1
    else:
        SEQS_value += line
    SEQS[SEQS_key] = SEQS_value
ffn_file.close()

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