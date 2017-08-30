#! /usr/bin/python

from sys import argv
import os
import glob
import subprocess
from subprocess import Popen
import multiprocessing
from itertools import islice


# This script will take in a FAA file with the following format for a header: >SM01_00816[SM01]
# The user will need to supply the directory where the FFNs are located.

# Requires full path to the CoreGenome.usearch


def AllFAA_descrip_dict(dir):
    # Read the AllFFN.all into a hash
    home = dir
    os.chdir(home)
    faa_file = open("AllFAA.all", "r")
    SEQS = {}
    for line in faa_file:
        if line.startswith(">"):
            SEQS_key = line.split(" ")[0].lstrip(">")
            SEQS_value = "\n\t\t\t\t\t\t\t\t\t\t\t\t\t" + line
            SEQS[SEQS_key] = SEQS_value
    faa_file.close()
    return SEQS




def AllFAA_seq_dict(dir):
    # Read the AllFFN.all into a hash
    home = dir
    os.chdir(home)
    faa_file = open("AllFAA.all", "r")

    SEQS = {}
    first = 0
    SEQS_key = ""
    SEQS_value = ""
    for line in faa_file:
        if line.startswith(">"):
            if first != 0:
                SEQS[SEQS_key] = SEQS_value
                SEQS_key = ""
                SEQS_value = ""
            SEQS_key = line.split(" ")[0] + "[" + line.split("_")[0].lstrip(">") + "]"
            first = 1
        else:
            SEQS_value += line
        SEQS[SEQS_key] = SEQS_value
    faa_file.close()
    return SEQS


def Create_ALL(core, dict_seqs):
    core_file = open(core, "r")
    home = os.getcwd()
    SEQS = dict_seqs

    # Read the CoreGenome.usearch file line by line, grab the sequence from the AllProteins dictionary
    # and output that to a single file
    os.chdir(home + "/" + "usearch_SNPanalysis_faa_files")
    for line in core_file:
        line = line.rstrip("\n")
        line = line.rstrip("\t")
        line_list = line.split("\t")
        file_rep_name = line.split("\t")[0].split("[")[0]
        file_rep_name = file_rep_name.lstrip(">") + ".fas"
        out_file = open(file_rep_name, "w")
        for x in line_list:
            new_seq_name = x.split("_")[0]
            new_seq_name = new_seq_name.lstrip(">")
            out_file.write(">" + new_seq_name + "\n" + SEQS[x])
        out_file.close()
    core_file.close()


def Muscle(dir):
    home = dir
    os.chdir(home + "/" + "usearch_SNPanalysis_faa_files")
    file_list = glob.glob(home + "/" + "usearch_SNPanalysis_faa_files/" + "*.fas")
    commands = []
    for genome in file_list:
        name = os.path.basename(genome)
        new_command = "muscle -in " + name + " -out " + name.rstrip(".fas") + ".aln"
        commands.append(new_command)

    max_workers = multiprocessing.cpu_count()
    processes = (Popen(cmd, shell=True) for cmd in commands)
    running_processes = list(islice(processes, max_workers))  # start new processes
    while running_processes:
        for i, process in enumerate(running_processes):
            if process.poll() is not None:  # the process has finished
                running_processes[i] = next(processes, None)  # start new process
                if running_processes[i] is None:  # no new processes
                    del running_processes[i]
                    break

    for genome in file_list:
        name = os.path.basename(genome)
        process = subprocess.Popen("mv " + name.rstrip(".fas") + ".aln " + name, shell=True)
        process.wait()


def BreakDown(dir, seq_dict_descrip):
    home = dir
    os.chdir(home + "/" + "usearch_SNPanalysis_faa_files")
    file_list = glob.glob(home + "/" + "usearch_SNPanalysis_faa_files/" + "*.fas")
    SEQS = seq_dict_descrip
    for genome in file_list:
        name = os.path.basename(genome)
        name_only = name.rstrip(".fas")
        process = subprocess.Popen("snp-sites -v -o " + name_only + ".vcf" + " " + name, shell=True)
        process.wait()
        vcf_file = name.rstrip(".fas") + ".vcf"
        vcf_file2 = name.rstrip(".fas") + ".final.report"
        vcf = open(vcf_file, "r")
        count = 0
        vcf_file_text = ""
        for line in vcf:
            count = count + 1
            if not line.startswith("##"):
                vcf_file_text = vcf_file_text + line
        vcf.close()
        if count > 3:
            new_vcf = open(vcf_file2, "w")
            new_vcf.write(SEQS[name_only] + vcf_file_text)
            new_vcf.close()
        process = subprocess.Popen("rm " + vcf_file, shell=True)
        process.wait()


###################################################################################################################

script, CORE, path_to_faa = argv

home = os.getcwd()
os.mkdir("usearch_SNPanalysis_faa_files")
os.chdir(path_to_faa)
process = subprocess.Popen("cat *.faa > AllFAA.all", shell=True)
process.wait()
process = subprocess.Popen("mv AllFAA.all " + home, shell=True)
process.wait()

Create_ALL(CORE, AllFAA_seq_dict(home))
Muscle(home)
BreakDown(home, AllFAA_descrip_dict(home))

os.chdir(home + "/" + "usearch_SNPanalysis_faa_files")
process = subprocess.Popen("cat *.final.report > usearch_CORE_faa.All_Genes.vcf", shell=True)
process.wait()
process = subprocess.Popen("rm *.final.report", shell=True)
process.wait()
process = subprocess.Popen("mv usearch_CORE_faa.All_Genes.vcf " + home, shell=True)
process.wait()
process = subprocess.Popen("FASconCAT_v1.0.pl -s -p", shell=True)
process.wait()
process = subprocess.Popen("mv -t " + home + " FcC_smatrix.phy FcC_smatrix.fas FcC_info.xls", shell=True)
process.wait()
os.chdir(home)
process = subprocess.Popen("mv FcC_smatrix.fas usearch_CORE_faa.fas", shell=True)
process.wait()
process = subprocess.Popen("mv FcC_smatrix.phy usearch_CORE_faa.phy", shell=True)
process.wait()
process = subprocess.Popen("snp-sites -v -o usearch_CORE_faa_SNPs.vcf usearch_CORE_faa.fas", shell=True)
process.wait()
process = subprocess.Popen("zip -r usearch_SNPanalysis_faa_files.zip usearch_SNPanalysis_faa_files", shell=True)
process.wait()
process = subprocess.Popen("rm -r usearch_SNPanalysis_faa_files", shell=True)
process.wait()
process = subprocess.Popen("rm FcC_info.xls", shell=True)
process.wait()
