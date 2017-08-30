#!/usr/bin/python
import subprocess, sys, os
from ete3 import PhyloTree, Tree, NodeStyle, TreeStyle, TextFace, SeqMotifFace, CircleFace, random_color #, add_face_to_node, SeqMotifFace
import argparse
import pprint  # for debugging

# Tree_Visualizer.py 

# This is an extremely simple script that firs calls snp-sites
# and generates a SNP alignment file. This file can then be read
# and passed onto the tree-visualization section of this script,
# which will render the tree with the SNP alignment
# on the right-hand side.

# This script assumes you've already built the tree.

# example usage:
# Tree_Visualizer.py -t TREE_FILE -L LEGEND_FILE -out view

import random

def render_alignment(t, seqs):
    for genome, seq in seqs.items():
        seqFace = SeqMotifFace(seq, seq_format='seq')
        (t & genome[1:]).add_face(seqFace, 0, "aligned")  # get rid of the '>' with [1:]

def random_color():
    r1 = random.randint(0,255)
    r2 = random.randint(0,255)
    r3 = random.randint(0,255)  
    return '#%02X%02X%02X' % (r1,r2,r3)

# Draw the grouped genes as a large red circle (where clade collapsing occurred)
def label_groups(t, ngroups):
    nstyles = {}
    for i in range(1, ngroups+1):
        nstyle = NodeStyle()
        nstyle["size"] = 10
        color = random_color()
        nstyle["fgcolor"] = color
        # assign a node style for each group from 1...N
        nstyles["Group{}".format(i)] = nstyle  

    for n in t.traverse():
        if n.name.startswith('Group'):
            n.set_style(nstyles[n.name])

    return nstyles

def populate_legend(ts, nstyles):
    ts.legend_position = 4

    all_files = os.listdir('.')
    for file in all_files:
        # find the legend.txt file
        # this is a hack, should implement a real solution eventually
        if 'legend.txt' in file:
            #j = 0
            with open(file) as f:
                for line in f:
                    line = line.strip()
                    line_arr = line.split('\t')

                    x = CircleFace(5, nstyles[line_arr[0]]['fgcolor'])

                    ts.legend.add_face(x, column=0)
                    ts.legend.add_face(TextFace(line_arr[0]), column=1)
                    genome_list = ''
                    for i in range(1, len(line_arr)):
                        if i > 5:
                            genome_list += '...'  # too many genomes
                            break
                        if i == len(line_arr) - 1:  # last element shouldn't end with comma
                            genome_list += line_arr[i][1:]
                        else:
                            genome_list += line_arr[i][1:] + ', '
                    ts.legend.add_face(TextFace(genome_list), column=2)
                    #j += 1


# Main routine
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Specify a data type.')
    parser.add_argument('-i', help="input FASTA file",
    					metavar="filename")
    parser.add_argument('-t', help="input tree file (.phy/.phylip)", required=True,
    					metavar="filename")
    parser.add_argument('-b', help="show_branch_length/support", action='store_const', const=True)
    parser.add_argument('-L', help="supply legend for grouped clades", metavar="filename")

    parser.add_argument('-out', help="output format (choices: pdf or view); defaults to pdf")
    args = parser.parse_args()

    if args.i:
        aln_path = os.path.join(os.getcwd(), args.i)
    else:
        aln_path = None

    if args.L:
        legend_path = os.path.join(os.getcwd(), args.L)
    else:
        legend_path = None

    if args.out != 'pdf' and args.out != 'view':
        sys.exit("-out must be either pdf or view")

    tree = os.path.join(os.getcwd(), args.t)
    print "Calling Tree_Visualizer with -i {}, -t {}, -b {}, -L {}, -out {}".format(
        args.i, args.t, args.b, args.L, args.out)
    t = Tree(tree)

    # ------------------------SNP Alignment------------------------------
    #print "Generating snp-sites file..."

    # only render the alignment if the -i parameter was specified
    if aln_path:
        process = subprocess.Popen("snp-sites -m {}".format(aln_path), shell=True)
        process.wait()

        seqs = {}  # key: genome, value: sequence
        with open("{}.snp_sites.aln".format(args.i)) as f:
            aln = f.read()
            aln = aln.split('\n')
            # if the ith element is the genome, the i+1th element is its sequence
            for i in range(0, len(aln), 2):  
                if aln[i].startswith('>'):
                    seqs[aln[i]] = aln[i+1]
        # Give each sequence a motif (visual alignment)
        render_alignment(t, seqs)

    #print "Rendering tree..."
    ts = TreeStyle()
    ts.title.add_face(TextFace("Tree", fsize=20), column=0)

    # args.b : if the user wants branch length/support values shown on the tree
    if args.b:
        ts.show_branch_length = True
        ts.show_branch_support = True

    # look for a 'legend.txt' file
    if legend_path:
        legend_file = ''
        for file in os.listdir('.'):
            if 'legend.txt' in file:
                legend_file = file
                break
        
        if legend_file:
            ngroups = 0
            with open(legend_file) as f:
                for line in f:
                    ngroups += 1
            nstyles = label_groups(t, ngroups)  # returns the styles dictionary (used for constructing legend)
            populate_legend(ts, nstyles)
        else:
            raise IOError("Cannot locate legend. Please label with the extension legend.txt")

    ts.legend.fgcolor = "#0030c1"

    if args.out == 'view':
        t.show(tree_style=ts)
    elif args.out == 'pdf':
        t.render("tree.pdf", w=600, units="mm", tree_style=ts)
