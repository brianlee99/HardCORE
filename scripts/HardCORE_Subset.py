#!/usr/bin/python
import subprocess
import os
import shutil
import glob

import sys
import random
import re

# This script should be run within the folder used for the main HardCORE analysis.
# The results will be written to CoreGeneCount.tsv, with each row being a distinct size
# ranging from step to the largest multiple of step at or below the total # of genomes.

# Note that the genomes should have a tab at the end, after which there is the genus name in square brackets
# e.g. [Citrobacter]

# Example usage:
# HardCORE_Subset.py --reps 10 --step 4
# or 
# HardCORE_Subset.py --step 5 --reps 8

# Incorrect usage:
# HardCORE_Subset.py --reps 6 --steps
# HardCORE_Subset.py --steps 32

# want to choose k genomes out of STRAINS
def CoreFinder(k):
    genomes = []
    subset = []
    genera = []
    genera_taken = {}
    taken = {}

    with open("STRAINS_LABELLED") as f:
        # count the number of lines
        for i, l in enumerate(f):
            l = l.rstrip()
            genomes.append(l)
    N = i + 1

    # 0 if not taken, 1 if taken
    for i in range(N):
        taken[i] = 0

    # find all the genera and put them into the genera list
    p = re.compile('.*(\[.+\]).*')
    for genome in genomes:
        m = p.match(genome)
        # the 1:-1 strips away the square brackets
        genus = m.groups(0)[0][1:-1]

        if genus not in genera:
            genera.append(genus)

    for genus in genera:
        genera_taken[genus] = 0

    # Here's the tricky part:
    # There are two cases to cover:
    # 1) The size is too small to cover all the possible categories, so we just select what we can
    if k < len(genera):
        upper_limit = k
    # 2) The size is big enough to cover at least one of every category
    # we find at least one genome in each category; rest can be random
    else:
        upper_limit = len(genera)

    for i in range(upper_limit):
        r = random.randint(0, N-1)
        m = p.match(genomes[r])
        genus = m.groups(0)[0][1:-1] # randomly chosen genus

        while genera_taken[genus] != 0:
            r = random.randint(0, N-1)
            m = p.match(genomes[r])
            genus = m.groups(0)[0][1:-1]

        # found one!
        genera_taken[genus] = 1
        taken[r] = 1
        subset.append(genomes[r])   

    if k >= len(genera):
        for i in range(k - len(genera)):
            r = random.randint(0, N-1)
            
            while taken[r] != 0:
                r = random.randint(0, N-1)
            
            taken[r] = 1
            subset.append(genomes[r])

    # strip the subset of the \t and the genus tag
    subset_stripped = []
    #q = re.compile('(.*)\\t\[.+\].*')
    q = re.compile('(\S*)(\s+)(\[.+\])')
    for i in range(len(subset)):
        m = q.match(subset[i])
        genome_stripped = m.groups(0)[0]
        subset_stripped.append(genome_stripped)

    return subset_stripped


def write_strains(subset):
    with open("STRAINS_SUBSET", 'w') as out:
        for genome in subset:
            out.write(genome)
            out.write('\n')


# if __name__ == "__main__":
#     if len(sys.argv) == 1:
#         print "You need to supply a k value!"
#         sys.exit(1)
#     k = int(sys.argv[1])
#     subset = CoreFinder(k)
#     write_strains(subset)


def strains_subset(k):
    subset = CoreFinder(k)
    write_strains(subset)

# This is going to run HardCORE on the subset, not the full thing
cwd = os.getcwd()

for i, arg in enumerate(sys.argv):
    if arg == '--reps':
        reps = int(sys.argv[i+1])
    elif arg == '--step':
        step = int(sys.argv[i+1])

print "Performing HardCORE Subset analysis with:"
print "Step size:", step
print "Repetitions:", reps

# Total number of genomes
numGenomes = 0
with open('STRAINS_LABELLED') as f:
    for line in f:
        numGenomes += 1

# if CoreGeneCount.tsv already exists, delete it
if 'CoreGeneCount.tsv' in os.listdir('.'):
    os.remove('CoreGeneCount.tsv')

if step == 1:
    size = 2  # cannot start with size = 1
else:
    size = step  # number of genomes to select from STRAINS; start with step/'k'

while size <= numGenomes:
    print "Analyzing a subset of size {}".format(size)
    with open('CoreGeneCount.tsv', 'a') as out:
        out.write(str(size))
        out.write('\t')
        
        for i in range(reps):
            # Step 1) Create STRAINS_SUBSET
            strains_subset(size)

            process = subprocess.Popen("StrainsSubset.py %d" % size, shell=True)
            process.wait()

            # Step 2) Create a tmp folder and move STRAINS_SUBSET
            # Also copy over relevant faa/ffn files
            os.makedirs('tmp/ffn')
            shutil.move("STRAINS_SUBSET", "tmp/STRAINS_SUBSET")

            strains = []
            with open('tmp/STRAINS_SUBSET') as file:
                for line in file:
                    strains.append(line.rstrip())

            # this is when we copy stuff over
            for strain in strains:
                shutil.copyfile("%s.faa" % strain, "tmp/%s.faa" % strain)
                shutil.copyfile("ffn/%s.ffn" % strain, "tmp/ffn/%s.ffn" % strain)

            print "--"
            print os.getcwd()

            # Step 3) The fun step! We get to run HardCORE_usearch.pl now.
            os.chdir('tmp')
            process = subprocess.Popen(("HardCORE_usearch.pl -dir ./ -ffn '%s/tmp/ffn' -strains STRAINS_SUBSET " 
            "-plen 0.9 -id 0.9 -threads 6") % (cwd), shell=True)  # 6 is just a placeholder for now
            process.wait()
            os.chdir('../')

            # Step 4) Once the analysis is complete, we will have to write to CoreGeneCount.tsv
            count = 0
            with open("tmp/Core_Genome/CoreGenome.usearch") as file:
                for line in file:
                    count += 1
            count = str(count)

            out.write(count)
            out.write('\t')

            #Step 5) Delete the tmp folder and repeat
            shutil.rmtree("tmp", ignore_errors=True)
        out.write('\n')  # new line for new set of entries
        size += step  # increment size

print "HardCore Subset analysis finished!"


