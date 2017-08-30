#! /usr/bin/python
# Does most of the work in rendering the web pages for HardCORE.

from flask import render_template, flash, redirect, request, url_for, send_from_directory, abort
from app import app
from forms import *

from werkzeug.utils import secure_filename
import os
import shutil
import subprocess
from blessings import Terminal
import glob
import datetime
import pandas as pd
import json
import time 
import re
import collections
import xlsxwriter
import pprint
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

t = Terminal()

# app.config variables
DATA_FOLDER = app.config['DATA_FOLDER']
APP_FOLDER = app.config['APP_FOLDER']
OUTPUT_FOLDER = app.config['OUTPUT_FOLDER']
GRAPHS_FOLDER = app.config['GRAPHS_FOLDER']
RF_FOLDER = app.config['RF_FOLDER']
TREE_FOLDER = app.config['TREE_FOLDER']
VISUALS_FOLDER = app.config['VISUALS_FOLDER']
HC_FOLDER = app.config['HC_FOLDER']
SCRIPTS_FOLDER = app.config['SCRIPTS_FOLDER']
ACCEPTED_AA_ALGORITHMS = app.config['ACCEPTED_AA_ALGORITHMS']
ACCEPTED_NT_ALGORITHMS = app.config['ACCEPTED_NT_ALGORITHMS']
DEFAULT_AA_ALGORITHM = app.config['DEFAULT_AA_ALGORITHM']
DEFAULT_NT_ALGORITHM = app.config['DEFAULT_NT_ALGORITHM']


#########################################################################################################################
###                                                                                                                   ###
###                                                    ROUTES                                                         ###
###                                                                                                                   ###
#########################################################################################################################
# The index page; primarily focused on running the main HardCORE analysis
@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = HardCoreForm()

    if request.method == 'POST':
        if form.validate_on_submit():
            # timestamp uniquely identifies a run, even if the run name doesn't
            folder_name = create_timestamp()
            # if run name is supplied, use it, otherwise use the date
            if form.run_name.data:
                run_name = form.run_name.data
            else:
                run_name = folder_name
            
            # parameters for hardcore.pl
            ident = form.ident.data
            plen = form.plen.data
            threads = form.threads.data

            # num threads default to max number of cores that user has
            cpu_count = multiprocessing.cpu_count()
            if not threads or threads > cpu_count:
                threads = cpu_count

            try:
                os.chdir(DATA_FOLDER)
                os.mkdir(folder_name)
                upload_fastas(form, folder_name)
                check_fastas(folder_name)
                # create STRAINS file for HardCore script
                create_strains(folder_name)
                # measure how long it takes to complete the analysis
                start_time = time.time()
                print 'ok'
                run_hardcore(folder_name, ident, plen, threads)
                run_faa_snp_analysis(folder_name)
                run_ffn_snp_analysis(folder_name)
                end_time = time.time()
                time_elapsed = round((end_time - start_time) / 60, 2)
            except ValueError as e:
                # go back to the index page without doing any of the analysis
                flash(str(e))
                delete_incomplete_run(folder_name)
                return render_template('index.html', form=form)

            # update runs.txt and runs.json
            update_runs_txt()
            update_runs_json(folder_name, run_name, ident, plen, time_elapsed)

            # do a quick check on the "integrity" of the coregenome and pangenome files
            errors = check_pan_core_usearch(folder_name)
            for error_type, error in errors.items():
                if error:
                    flash(error_type)

            flash('Your file has been successfully uploaded and processed.')
            flash('The analysis took {} minutes.'.format(time_elapsed))
            return redirect(url_for('summary_table'))
        else:
            # Provide meaningful error messages for user
            if not form.threads.data:
                flash("threads must be an integer")
            elif form.threads.data <= 0:
                flash("threads must be positive (or blank)")
            if form.ident.data < 0.5 or form.ident.data > 1:
                flash("ident must be between 0.5 and 1")
            if form.plen.data < 0.5 or form.plen.data > 1:
                flash("plen must be between 0.5 and 1")
    
    return render_template('index.html', form=form)


# Create a unique timestamp, which will be used as the folder name for the run
# format: YYYYMMDD_HHMMSS
def create_timestamp():
    dt = datetime.datetime.now()
    dt_str = dt.strftime("%Y%m%d_%H%M%S")
    return dt_str


# Uploads and unzips the faa/ffn archives to the server side
def upload_fastas(form, folder_name):
    faa = form.faa.data
    ffn = form.ffn.data
    faa_filename = secure_filename(faa.filename)
    ffn_filename = secure_filename(ffn.filename)
    
    # check if the filenames are distinct
    if faa_filename == ffn_filename:
        raise ValueError("The .ffa and .ffn archives have the same name")

    # upload faa/ffn archives
    print t.blue("Uploading faa/ffn archives")
    wd = os.path.join(DATA_FOLDER, folder_name)
    os.chdir(wd)
    faa.save(os.path.join(wd, faa_filename))
    ffn.save(os.path.join(wd, ffn_filename))

    print t.blue("Upload complete.")
    print t.blue("Unzipping archives")

    # unzip archives
    extension = faa_filename.split('.')[-1]
    # e.g. .zip, .gz, .tar
    if extension == "zip":
        process = subprocess.Popen("unzip " + faa_filename, shell=True)
        process.wait()
        process = subprocess.Popen("unzip -d ffn " + ffn_filename, shell=True)
        process.wait()
    elif extension == "tar":
        os.mkdir("ffn")
        process = subprocess.Popen("tar -xf %s -C ffn" % ffn_filename, shell=True)
        process.wait()
        process = subprocess.Popen("tar -xf %s" % faa_filename, shell=True)
        process.wait()
    elif extension == "gz":
        os.mkdir("ffn")
        process = subprocess.Popen("tar -xzf %s -C ffn" % ffn_filename, shell=True)
        process.wait()
        process = subprocess.Popen("tar -xzf %s" % faa_filename, shell=True)
        process.wait()

    # delete zip files when finished
    os.remove(os.path.join(wd, faa_filename))
    os.remove(os.path.join(wd, ffn_filename))
    print t.blue("Unzip complete.")


# Checks if the faa and ffn archive contents are valid
def check_fastas(folder_name):
    goto_folder(folder_name)
    print t.green("Checking fastas:")

    # check for empty archives
    faa_list = glob.glob('*')
    if len(faa_list) == 0:
        print t.red("Error: .faa archive is empty.\nExiting...")
        raise ValueError(".faa archive is empty")
    if 'ffn' not in faa_list:
        print t.red("Error: .ffn archive is empty.\nExiting...")
        raise ValueError(".ffn archive is empty")

    faa_list.remove('ffn')  # remove 'ffn' folder from the glob
    
    # check if the files in the faa archive all terminate in .faa
    for file in faa_list:
        if not file.endswith(".faa"):
            print t.red("Error: .faa archive contains non-faa files.\nExiting...")
            raise ValueError(".faa archive contains non-faa files")
        # check just the first header of an .faa file to confirm
        # it was annotated properly by Prokka
        with open(file) as f:
            for line in file:
                if line.startswith('>'):
                    line = line.strip()[1:]
                    check_fasta_header(line)
            line = f.readline().strip()[1:]  # get rid of >
            check_fasta_header(line)

    # check if the files in the ffn archive all terminate in .ffn
    os.chdir('ffn')
    ffn_list = glob.glob('*')
    for file in ffn_list:
        if not file.endswith(".ffn"):
            print t.red("Error: .ffn archive contains non-ffn files.\nExiting...")
            raise ValueError(".ffn archive contains non-ffn files")
        # do the same with ffn files
        with open(file) as f:
            for line in file:
                if line.startswith('>'):
                    line = line.strip()[1:]
                    check_fasta_header(line)
            line = f.readline().strip()[1:]  # get rid of >
            check_fasta_header(line)

    os.chdir('../')

    # check if there are an equal number of faa's and ffn's
    if len(faa_list) != len(ffn_list):
        print t.red("Error: Lengths of faa and ffn archives differ.\nExiting...")
        raise ValueError("Lengths of faa and ffn archives differ")

    # sort the faa's/ffn's so it's easier to compare their names
    faa_list.sort()
    ffn_list.sort()
    accepted_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.-'

    for i in range(len(faa_list)):
        faa_string = faa_list[i]
        ffn_string = ffn_list[i]
        # taking out .faa/.ffn from the file names
        faa_prefix = faa_string[:-4]
        ffn_prefix = ffn_string[:-4]

        # double check that the filenames only contain accepted letters
        for letter in faa_prefix:
            if letter not in accepted_letters:
                print t.red("Error: Invalid faa filename")
                raise ValueError("Invalid faa filename")
        for letter in ffn_prefix:
            if letter not in accepted_letters:
                print t.red("Error: Invalid ffn filename")
                raise ValueError("Invalid ffn filename")
        # compare names against each other to confirm they are the same
        if faa_prefix != ffn_prefix:
            print t.red("The faa and ffn files do not correspond")
            raise ValueError("The faa and ffn files do not correspond")

    os.chdir('ffn')
    # Go inside each FFN file and ensure no sequence starts with ">NODE"
    for ffn in ffn_list:
        with open(ffn) as f:
            for line in f:
                if line.startswith('>NODE'):
                    print t.red("Error: Invalid PROKKA annotation")
                    raise ValueError("Invalid PROKKA annotation")

    print t.green("FFN/FAA check completed! Files are valid")


def check_fasta_header(line):
    accepted_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.-'
    split_underline = line.split('_')

    # split[0] is the genome, split[1] is the gene number
    if len(split_underline) != 2:
        print t.red("Error in file header labelling")
        raise ValueError("Error in file header labelling")
    
    # check if split_underline[0] only contains accepted letters
    for letter in split_underline[0]:
        if letter not in accepted_letters:
            print t.red("Error in file header labelling")
            raise ValueError("Error in file header labelling")
    
    # check if the prokka gene number is 5 digits, followed by a space
    for i, letter in enumerate(split_underline[1][:6]):
        if i == 5:
            if letter != ' ':
                print t.red("Error in file header labelling")
                raise ValueError("Error in file header labelling")
        else:
            if letter not in '0123456789':
                print t.red("Error in file header labelling")
                raise ValueError("Error in file header labelling")


# Creates a STRAINS file by reading the list of faa files
def create_strains(folder_name):
    os.chdir(os.path.join(DATA_FOLDER, folder_name))
    faa_list = sorted(glob.glob('*.faa'))
    with open('STRAINS', 'w') as strains:
        for faa in faa_list:
            faa = faa[:-4]  # get rid of prefix
            strains.write(faa + '\n')


# Calls the main hardcore script
def run_hardcore(folder_name, ident, plen, threads):
    os.chdir(os.path.join(DATA_FOLDER, folder_name))
    # need the full path for the ffn folder
    wd = os.path.join(DATA_FOLDER, folder_name)
    ffn_path = os.path.join(wd, 'ffn')

    script_dir = os.path.join(HC_FOLDER, "HardCORE_usearch.pl")

    process = subprocess.Popen("perl {} -dir ./ -strains STRAINS -ffn '{}' -id {} -plen {} -threads {}"
        .format(script_dir, ffn_path, ident, plen, threads), shell=True)
    process.wait()


# Runs the FAA "SNP" analysis script
def run_faa_snp_analysis(folder_name):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, "Core_Genome"))
    wd = os.path.join(DATA_FOLDER, folder_name)
    core_path = os.path.join(wd, 'Core_Genome/CoreGenome.usearch')
    prokka_path = os.path.join(wd, 'PROKKA_backups')  # Have to use original prokka-style files

    script_dir = os.path.join(HC_FOLDER, "HardCore_usearch_SNP_FAA_analysis.py")

    process = subprocess.Popen("python {} {} {}".format(script_dir, core_path, prokka_path), shell=True)
    process.wait()


# Runs the FFN SNP analysis script
def run_ffn_snp_analysis(folder_name):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, "Core_Genome"))
    wd = os.path.join(DATA_FOLDER, folder_name)
    core_path = os.path.join(wd, 'Core_Genome/CoreGenome.usearch')
    ffn_path = os.path.join(wd, 'ffn')

    script_dir = os.path.join(HC_FOLDER, "HardCore_usearch_SNP_FFN_analysis.py")

    process = subprocess.Popen("python {} {} {}".format(script_dir, core_path, ffn_path), shell=True)
    process.wait()


# Updates runs.txt after a successful run -> this is done by reading the
# list of folders (excluding Tree) and writing the directories as a text file
def update_runs_txt():
    os.chdir(DATA_FOLDER)
    folders = next(os.walk('.'))[1]
    if 'Tree' in folders:
        folders.remove('Tree')

    with open('runs.txt', 'w') as out:
        for folder in sorted(folders):  # folder is not sorted by default
            out.write(folder + '\n')

    print t.blue('runs.txt has been updated to reflect the new list of runs')


# Updates runs.json which contains a brief summary of each run
def update_runs_json(folder_name, run_name, ident, plen, time_elapsed):
    counts = calculate_counts(folder_name)
    os.chdir(DATA_FOLDER)
    data = {}
    # if json file exists already, load it in
    if 'runs.json' in os.listdir('.'):
        with open('runs.json') as json_file:
            data = json.load(json_file)

    # new entry for runs.json
    dt_obj = datetime.datetime.strptime(folder_name, '%Y%m%d_%H%M%S')
    pretty_date = dt_obj.strftime("%b %d %Y, %H:%M:%S")
    run_data = counts

    # add more information to the run_data entry
    run_data.update({'date': pretty_date,
                     'run_name': run_name,
                     'ident': float(ident*100),
                     'plen': float(plen*100),
                     'time_elapsed': time_elapsed
                     })
    data.update({folder_name: run_data})
    
    with open('runs.json', 'w') as out:
        json.dump(data, out, indent=4)


# Calculate statistics based on the run 'folder_name'
def calculate_counts(folder_name):
    core_count = 0
    pan_count = 0
    singleton_count = 0
    snp_count = 0
    mutation_count = 0
    genome_count = 0
    # =reads all .faa files to calculate # genes / genome (avg)
    genes_per_genome = average_genes_per_genome(folder_name)  

    goto_folder(folder_name)
    with open('Core_Genome/CoreGenome.usearch') as file:
        for line in file:
            core_count += 1
    with open('Pan_Genome/PanGenome.usearch') as file:
        for line in file:
            pan_count += 1
    with open('Singletons/Singletons.usearch') as file:
        for line in file:
            singleton_count += 1
    with open('Core_Genome/usearch_CORE_faa_SNPs.vcf') as file:
        for line in file:
            mutation_count += 1
        mutation_count -= 3  # First 3 lines are comments
    with open('Core_Genome/usearch_CORE_ffn_SNPs.vcf') as file:
        for line in file:
            snp_count += 1
        snp_count -= 3  # First 3 lines are comments
    with open('STRAINS') as file:
        for line in file:
            genome_count += 1
    
    return {'core_count': core_count, 'pan_count': pan_count, 'singleton_count': singleton_count,
            'snp_count': snp_count, 'mutation_count': mutation_count,
            'genome_count': genome_count, 'genes_per_genome': genes_per_genome}


# in calcuating genes/genome, # genes is meant as the total number of protein sequences
# (so not counting things like tRNA) over the number of genomes
def average_genes_per_genome(folder_name):
    num_genomes = 0
    gene_count = 0
    goto_folder(folder_name)
    faa_list = glob.glob('*.faa')
    for faa in faa_list:
        with open(faa) as file:
            num_genomes += 1
            for line in file:
                if line.startswith('>'):
                    gene_count += 1
    return round(float(gene_count) / num_genomes, 2)


# If the analysis was unsuccessful (usually due to invalid fails),
# deletes the folder that was created that would've contained the run
def delete_incomplete_run(folder_name):
    os.chdir(DATA_FOLDER)
    shutil.rmtree(folder_name)

# ================================= SUMMARY_TABLE =========================================
@app.route('/summary_table', methods=['GET', 'POST'])
def summary_table():
    form = SummaryTableForm()

    # load the runs into an ordered dictionary, ordered_runs
    with open(os.path.join(DATA_FOLDER, 'runs.json')) as f:
        runs = json.load(f)

    if request.method == 'POST':
        if form.validate_on_submit():
            to_process = request.form.getlist("check")

            # Download files (for a particular set of runs)
            if form.download.data: 
                if not form.singles.data and not form.core_genome.data and not form.pan_genome.data:
                    flash("You must select at least one download parameter!")
                elif len(to_process) == 0:
                    flash("You must select at least one run to download!")
                else:
                    prepare_summary_zip(form, to_process)
                    return redirect(url_for('download', filename="output/download.zip"))
            
            # Delete run(s)
            elif form.delete.data:
                delete_runs(to_process)
                deleted_runs = ""
                for run in to_process:
                    deleted_runs += run
                    deleted_runs += ", "
                deleted_runs = deleted_runs[:-2]
                flash("Runs deleted: {}".format(deleted_runs))
                return redirect(url_for('summary_table'))
                
            # Download excel table
            elif form.export.data:
                prepare_summary_table(runs)
                return redirect(url_for('download', filename="output/summary_table.xlsx"))

            # deletes mostly unnecessary data generated by the app over time
            elif form.cleanup.data:
                cleanup()
                flash("Cleanup completed")

    return render_template('summary_table.html', form=form, runs=runs)


# Prepares a downloadable zip, whose content the user decides
# (e.g. the runs to include, the files to download for those runs)
def prepare_summary_zip(form, to_process):
    print t.blue("Preparing summary zip")
    os.chdir(DATA_FOLDER)
    # If the zip file already exists, delete it
    if 'download.zip' in os.listdir('.'):
        os.remove('download.zip')

    # Put everything into download.zip
    for folder in to_process:
        if form.singles.data:
            process = subprocess.Popen("zip download.zip -r %s/Singletons" % folder, shell=True)
            process.wait()
        if form.core_genome.data:
            process = subprocess.Popen("zip download.zip %s/Core_Genome/CoreGenome.usearch" % folder, shell=True)
            process.wait()
        if form.pan_genome.data:
            process = subprocess.Popen("zip download.zip %s/Pan_Genome/PanGenome.usearch" % folder, shell=True)
            process.wait()

    # take the zip (in the data folder) and move to output folder
    shutil.move("download.zip", "{}/download.zip".format(OUTPUT_FOLDER))


# Prepares a summary table, which is an almost replica of the table
# rendered on the summary_table page
def prepare_summary_table(runs):
    print t.blue("Preparing summary table")
    df = pd.DataFrame.from_dict(runs)
    df = df.transpose()
    col_list = list(df)

    # reorder columns to be identical to the HTML table
    (col_list[7], col_list[0], col_list[6], col_list[5],
        col_list[3], col_list[10], col_list[8], col_list[4],
        col_list[1], col_list[9], col_list[11], col_list[2]) = col_list
    df = df[col_list]

    # prettify column names
    df.columns = ['Date', 'Run Name', 'Time Elapsed', '% Identity', '% ID/Length', 'Genomes', 
        'Genes/genome', 'Core Genes', 'Pan Genes', 'Singleton Genes', 'FAA SNPs', 'FFN SNPs']

    # write to output.xlsx
    os.chdir(OUTPUT_FOLDER)
    writer = pd.ExcelWriter('summary_table.xlsx')
    df.to_excel(writer,'Sheet1', index=False)
    writer.save()


# Deletes files related to runs to be deleted (to_process)
def delete_runs(to_process):
    os.chdir(DATA_FOLDER)
    # delete the actual run folders in /data
    for folder in to_process:
        shutil.rmtree(folder)

    update_runs_txt()

    # remove the run entries from runs.json and current_references.json
    with open('runs.json') as file:
        data = json.load(file)
    for folder in to_process:
        if folder in data.keys():
            data.pop(folder)
            print t.red('{} removed from runs.json'.format(folder))
    with open('runs.json', 'w') as out:
        json.dump(data, out, indent=4)

    with open('current_references.json') as file:
        data = json.load(file)
    for folder in to_process:
        if folder in data.keys():
            data.pop(folder)
            print t.red('{} removed from current_references.json'.format(folder))
    with open('current_references.json', 'w') as out:
        json.dump(data, out)

    # delete .json files used for pan-visualize and pan-components
    os.chdir(VISUALS_FOLDER)
    for folder in to_process:
        json_file = "{}.json".format(folder)
        json_flare_file = "{}_flare.json".format(folder)
        if json_file in os.listdir('.'):
            os.remove(json_file)
        if json_flare_file in os.listdir('.'):
            os.remove(json_flare_file)


# delete rarefaction graph
def delete_rarefaction(folder_name):
    os.chdir(RF_FOLDER)
    if '{}.jpg'.format(folder_name) in os.listdir('.'):
        os.remove('{}.jpg'.format(folder_name))
        print t.red('{}.jpg deleted'.format(folder_name))


# delete accessory files in the data folder that aren't always needed
def cleanup():
    # delete downloadable files; the app does not need them
    os.chdir(OUTPUT_FOLDER)
    all_files = os.listdir('.')
    for file in all_files:
        os.remove(file)
        print t.red(file + " deleted")
    
    # delete graph-related files (snps/aa changes)
    os.chdir(GRAPHS_FOLDER)
    folders = os.listdir('.')
    for folder in folders:
        shutil.rmtree(folder)
        print t.red(folder + " deleted")
    
    # check for lingering invalid run folders
    check_incomplete_runs()


# sees if there are any folders that are not represented in runs.txt
# runs.txt would only contain completed runs; therefore should delete
# invalid runs 
def check_incomplete_runs():
    os.chdir(DATA_FOLDER)
    # get the list of completed runs
    try:
        with open('runs.txt') as f:
            runs = [line.strip() for line in f]

        # get the list of folders
        folders = next(os.walk('.'))[1]
        if 'Tree' in folders:
            folders.remove('Tree')

        runs.sort()
        folders.sort()
        # because runs and folders are sorted, we can linearly
        # determine if a folder in folders is found in runs    
        i, j = 0, 0
        while i < len(folders) and j < len(runs):
            if folders[i] == runs[j]:
                i += 1
                j += 1
            else:
                print t.red("Deleting {}".format(folders[i]))
                delete_incomplete_run(folders[i])
                i += 1
        while i < len(folders):
            print t.red("Deleting {}".format(folders[i]))
            delete_incomplete_run(folders[i])
            i += 1
    except IOError as e:
        print t.red("I/O error({0}): {1}".format(e.errno, e.strerror))
        print t.red("Unable to find runs.txt. Quitting")
#================================== CORE_SUBSET =======================================
@app.route('/core_subset', methods=['GET', 'POST'])
@app.route('/core_subset/<folder_name>', methods=['GET', 'POST'])
def core_subset(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)

    strains = []
    run_name = ''
    # keys: genomes, values: list of core subset genes that 
    # correspond at the ith position, e.g.
    # {A: [1, 4, 12],               => genes A1, B12, and C86 are part of the same family
    #  B: [12, 148, 62],
    #  C: [86, 59, 43]}
    subset_genes = {}
    # a simple list of predicted functions for each core subset gene
    predicted_functions = []
    count = 0       # number of core subset genes
    form = None
    download_form = None

    # Render the form that allows the user to select a run
    select_list = get_select_list()

    # render table
    if folder_name:
        run_name = get_run_name(folder_name)
        form = CoreSubsetForm()
        # read STRAINS
        goto_folder(folder_name)
        with open(os.path.join(DATA_FOLDER, folder_name, 'STRAINS')) as f:
            strains = [line.strip() for line in f]

        # only render a download form if core subset analysis was already done
        if 'CORE_SUBSET.core.subset' in os.listdir('.'):
            # if core subset is empty, don't render the download form
            stats = os.stat('CORE_SUBSET.core.subset')
            if stats.st_size != 0:
                download_form = DownloadSubsetForm()
                subset_genes, count = complete_subset_genes(folder_name)
                predicted_functions = complete_predicted_functions(strains, subset_genes)

    if request.method == 'POST':
        # Submitting a list of strains for finding the core subset
        if form.submit.data:
            subset = request.form.getlist("check")
            if len(subset) == 0:
                flash("You must select at least one strain.")
            else:
                try:
                    # run the core subset analysis
                    create_core_subset_file(folder_name, subset)
                    run_core_subset(folder_name)
                    flash('The core subset analysis has successfully finished.')
                except ValueError:
                    flash('There appears to be no core genes for the subset you have provided.')
            return redirect(url_for('core_subset', folder_name=folder_name))
        
        # Downloading the subset_genes table in excel format
        elif download_form.download.data:
            prepare_core_subset_xlsx(subset_genes, predicted_functions)
            return redirect(url_for('download', filename="output/core_subset.xlsx"))

        # Download alignments for each core subset gene
        elif download_form.alignments.data:
            prepare_core_subset_alns(folder_name, subset_genes)
            return redirect(url_for('download',
                filename="output/core_subset.zip".format(folder_name)))

    return render_template('core_subset.html', form=form, select_list=select_list,
        download_form=download_form, strains=strains, subset_genes=subset_genes,
        predicted_functions=predicted_functions, count=count, run_name=run_name)


# returns a dictionary with genomes as keys and core subset genes as values (list)
# also returns the number of core subset genes
def complete_subset_genes(folder_name):
    count = 0
    subset_genes = {}

    goto_folder(folder_name)
    with open('CORE_SUBSET.core.subset') as file:
        for line in file:
            count += 1
            # find all patterns that look like '>Genome.Name_00001[...]'
            p = re.compile(r'>([^_]+)_(\d+)\[.*\]')
            matches = re.findall(r'>[^\t\n]*', line)
            for match in matches:
                strain, gene = p.match(match).groups()
                if strain in subset_genes.keys():
                    subset_genes[strain].append(gene)
                else:
                    subset_genes[strain] = [gene]
    
    return subset_genes, count


# returns a list of prokka predicted functions for each subset gene
def complete_predicted_functions(strains, subset_genes):
    predicted_functions = []
    # pick some genome, does not matter
    some_strain = subset_genes.keys()[0]
    # used as placeholder for the original order of subset genes
    pf_holder = {}
    
    with open('PROKKA_backups/{}.orig.faa'.format(some_strain)) as file:
        # matching everything to the right of a space gives us the predicted function
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                p = re.compile(r'[^ ]+_(\d+) (.+)')
                # gene: gene number; function: its pred. function
                gene, function = p.match(line).groups()
                pf_holder[gene] = function 
        predicted_functions = [pf_holder[gene] for gene in subset_genes[some_strain]]
    
    return predicted_functions


# prepares a zip of core subset gene alignments
# subset_genes is a dictionary, with keys as strains in your subset
# and values as lists of subset genes
def prepare_core_subset_alns(folder_name, subset_genes):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'ffn'))
    FASTA_HEADER = r'[^ ]+_(\d+) (.+)'

    # the first genome will server as our 'reference' point
    ref_genome = sorted(subset_genes.keys())[0]
    ref_genes = subset_genes[ref_genome]

    to_write = {ref_genome+'_'+gene+'.fasta': [] for gene in ref_genes}

    for genome, genes in subset_genes.items():
        # track each subset gene's sequence, for each genome
        seqs = {gene: '' for gene in genes} 
        sorted_genes = sorted(genes)  # a sorted list of prokka id's 
        write = False  # if True, then record sequence/header data (?)

        # go through our entire ffn file, look for genes that are SUBSET
        # genes (as they are in sorted_genes), and update our 'seqs' dictionary

        # since the genes in the ffn file are also sorted from top to bottom,
        # we can simply look at the 0th element in sorted_genes, compare it to the
        # current gene in the ffn, and if found, we update seqs and pop the 0th element
        # in sorted_genes.

        with open("{}.ffn".format(genome)) as f:
            for line in f:
                # end of previous sequence, possibly start of a new one
                if line.startswith('>'):
                    if write:
                        write = False
                        # get rid of all genes that have already been accounted for
                        sorted_genes.pop(0)
                        # don't need to read the rest of file if we covered all subset genes
                        if len(sorted_genes) == 0:
                            break
                    # let's consider if the Prokka number is in our sorted_genes
                    p = re.compile(FASTA_HEADER)
                    m = p.match(line)
                    gene_num = m.groups()[0]
                    if gene_num == sorted_genes[0]:
                        # we found a core subset gene
                        write = True
                        # get rid of gene number/annotation;
                        # causes issues with tree builder
                        line = line.split('_')[0] + "\n"
                        seqs[sorted_genes[0]] += line
                else:
                    if write:
                        seqs[sorted_genes[0]] += line 

        for prokka_id, seq in seqs.items():
            # convert prokka id to an index that gets you the prokka id of
            # the "reference" genome 
            i = genes.index(prokka_id)
            aln_filename = ref_genome + "_" + ref_genes[i] + ".fasta"
            to_write[aln_filename].append(seq)

    for l in to_write.values():
        l.sort()

    os.chdir(OUTPUT_FOLDER)
    if "core_subset.zip" in os.listdir('.'):
        os.remove("core_subset.zip")

    # perform the alignment
    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    p.map(worker_subset, to_write.items())

    for filename in to_write.keys():
        process = subprocess.Popen("zip core_subset.zip {}".format(filename), shell=True)
        process.wait()
        os.remove(filename)



def worker_subset(item):
    # item[0] is a two-tuple
    filename, contents = item
    with open(filename, 'w') as out:
        for seq in contents:
            out.write(seq)
    filename_new = "{}_out".format(filename)
    process = subprocess.Popen("muscle -in {0} -out {1}".format(filename, filename_new), shell=True)
    process.wait()

    os.remove(filename)
    os.rename(filename_new, filename)


# columns are strains (and predicted function), rows are core subset genes
def prepare_core_subset_xlsx(subset_genes, predicted_functions):
    os.chdir(OUTPUT_FOLDER)
    df = pd.DataFrame(columns=[strain for strain in subset_genes.keys()] + ['Predicted function'])

    for i in range(len(predicted_functions)):
        loc = [genes[i] for strain,genes in subset_genes.items()] + [predicted_functions[i]]
        df.loc[i] = loc

    writer = pd.ExcelWriter('core_subset.xlsx')
    df.to_excel(writer, 'Sheet1')
    writer.save()


# Creates a CORE_SUBSET file based on the subset provided by the user
def create_core_subset_file(folder_name, subset):
    goto_folder(folder_name)
    with open('CORE_SUBSET', 'w') as out:
        for strain in subset:
            out.write(strain)
            out.write('\n')


# Runs the core_subset analysis. Throws ValueError if there is no subset
# (although technically this is not an error) 
def run_core_subset(folder_name):
    goto_folder(folder_name)
    folder_dir = os.path.join(DATA_FOLDER, folder_name)
    process = subprocess.Popen(("HardCORE_usearch.pl -core_subset CORE_SUBSET "
        "-dir {} -strains STRAINS").format(folder_dir), shell=True)
    process.wait()
    # Check if the CORE_SUBSET file is empty or not
    with open('CORE_SUBSET.core.subset') as file:
        contents = file.read()
        if contents == '':
            raise ValueError("The core subset is empty")
# ================================ SNP_VIEWER ==========================================
@app.route('/core_genome', methods=['GET', 'POST'])
@app.route('/core_genome/<folder_name>', methods=['GET', 'POST'])
def core_genome(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)

    table = {}
    snp_form = None
    curr_ref = None
    run_name = ''

    select_list = get_select_list()

    if folder_name:
        run_name = get_run_name(folder_name)
        snp_form = SNPForm()
        os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
        unzip_fas(folder_name)
        
        # load .json with table contents to quickly render table in future
        if 'table.json' in os.listdir('.'):
            with open('table.json') as file:
                table = json.load(file)
        else:
            # create the initial table.json (might take a while)
            table = create_snp_json(folder_name)

        # update reference genome selection dropdown
        strains_location = os.path.join(DATA_FOLDER, folder_name, 'STRAINS')
        with open(strains_location) as file:
            strains = [line.strip() for line in file]
        snp_form.select_reference.choices = [(strain, strain) for strain in strains]

        # determine the current ref genome based on the run (folder_name)
        # also, if there is no reference, set ref to the first genome alphabetically
        with open(os.path.join(DATA_FOLDER, 'current_references.json')) as j:
            all_refs = json.load(j)
        if folder_name in all_refs:
            curr_ref = all_refs[folder_name]
        else:
            curr_ref = strains[0]
            update_reference(folder_name, curr_ref, all_refs)

    if request.method == 'POST':
        # Download the core gene alignment as NT sequences
        if 'aln_ffn' in request.form.keys():
            index = request.form['aln_ffn'][:-4]
            # convert to an actual gene filename
            gene = table[index]['name'] + ".fas"
            aln_file = 'data/{}/Core_Genome/usearch_SNPanalysis_ffn_files/{}'.format(folder_name, gene)
            return redirect(url_for('download', filename=aln_file))
        
        # Download alignment as AA sequences
        elif 'aln_faa' in request.form.keys():
            index = request.form['aln_faa'][:-4]
            gene = table[index]['name'] + ".fas"
            aln_file = 'data/{}/Core_Genome/usearch_SNPanalysis_faa_files/{}'.format(folder_name, gene)
            return redirect(url_for('download', filename=aln_file))

        # Download SNP .vcfs and .fas for selected genes
        elif snp_form.download_selected.data:
            try:
                # False -> do not download all
                prepare_snp_data(folder_name, table, False)  
                return redirect(url_for('download', filename="output/snp_data.zip".format(folder_name)))
            except ValueError:
                flash("You must select at least one gene for downloading!")
        
        # Download vcf and fas for every core gene
        elif snp_form.download_all.data:
            # True -> download all
            prepare_snp_data(folder_name, table, True)
            return redirect(url_for('download', filename="output/snp_data.zip".format(folder_name)))

        # Update the data by changing the reference genome
        elif snp_form.update_reference.data:
            new_ref = snp_form.select_reference.data
            if new_ref != curr_ref:
                update_reference(folder_name, new_ref, all_refs)
                flash('The reference has been updated to {}'.format(new_ref))
            else:
                flash('You have selected the same reference as before')
            return redirect(url_for('core_genome', folder_name=folder_name))
        
        # Download SNP view table in xlsx format
        elif snp_form.download_xlsx.data:
            prepare_snp_table(folder_name, table)
            return redirect(url_for('download', filename='output/snp_table.xlsx'.format(folder_name)))

    return render_template('core_genome.html', snp_form=snp_form, select_list=select_list, table=table,
        curr_ref=curr_ref, folder_name=folder_name, run_name=run_name)


# if the archives containing the .fas files haven't been unzipped yet, do so
def unzip_fas(folder_name):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    if not 'usearch_SNPanalysis_ffn_files' in os.listdir('.'):
        process = subprocess.Popen('unzip usearch_SNPanalysis_ffn_files.zip', shell=True)
        process.wait()
    if not 'usearch_SNPanalysis_faa_files' in os.listdir('.'):
        process = subprocess.Popen('unzip usearch_SNPanalysis_faa_files.zip', shell=True)
        process.wait()


# Returns a dictionary with keys as indexes, and values as core genes
# where for each core gene, statistics are recorded based on sequence length,
# snps/mutations, predicted function, and so on
def create_snp_json(folder_name):
    table = {}
    goto_folder(folder_name)
    os.chdir('Core_Genome')

    ffn_vcf = 'usearch_CORE_ffn.All_Genes.vcf'
    faa_vcf = 'usearch_CORE_faa.All_Genes.vcf'
    # --------------------------------------- Inner functions -----------------------------------------------
    def fill_in_basics(datatype):
        if datatype == 'nt':
            folder = 'usearch_SNPanalysis_ffn_files'
            init_entry = True
        elif datatype == 'aa':
            folder = 'usearch_SNPanalysis_faa_files'
            init_entry = False
        i = 0 # index, used as key

        # go into one of the fas folders, and read the sequence length
        # for each core gene
        os.chdir(folder)
        genes = sorted(glob.glob('*.fas'))
        for gene in genes:
            with open(gene) as file:
                file.readline()
                seq = ''
                line = file.readline().strip()
                while not line.startswith('>'):
                    seq += line
                    line = file.readline().strip()
            length = len(seq)
            if init_entry:
                table[str(i)] = {'length_nt': length,
                                 'snp_count': 0,
                                 'name': gene[:-4],
                                 'mutations': 0,
                                 'nt_sim': 100.0,
                                 'aa_sim': 100.0}
            else:
                table[str(i)]['length_aa'] = length
            i += 1
        os.chdir('../')
        # returns the list of globbed .fas files, for convenience
        return genes

    def fill_in_snps(datatype):
        if datatype == 'nt':
            vcf = ffn_vcf
            counter = "snp_count"
            length = "length_nt"
            similar = "nt_sim"
        elif datatype == 'aa':
            vcf = faa_vcf
            counter = "mutations"
            length = "length_aa"
            similar = "aa_sim"
        i = -1
        
        with open(vcf) as file:
            for line in file:
                line = line.strip()
                # found a snp/mut at our current gene
                if line.startswith("1"):
                    table[str(i)][counter] += 1
                elif line.startswith('>'):  # new core gene
                    if i > -1:
                        # calculate % diff within a core-gene
                        percent_diff = round(table[str(i)][counter] / float(table[str(i)][length]), 4) * 100
                        table[str(i)][similar] = 100 - percent_diff
                    i += 1
                    p = re.compile('>([^ ]*) (.*)')
                    m = p.match(line)
                    gene_name = m.groups()[0]

                    # this while loop is needed to find the next core gene to start counting snps
                    # this is because not all core genes will have snps
                    while table[str(i)]['name'] != gene_name:
                        i += 1
                        if i > -1:
                            # since we know genes that are skipped in the vcf
                            # don't have any snps
                            table[str(i)][similar] = 100.0
                        p = re.compile('>([^ ]*) (.*)')
                        m = p.match(line)
                        gene_name = m.groups()[0]
        
        percent_diff = round(table[str(i)][counter] / float(table[str(i)][length]), 4) * 100
        table[str(i)][similar] = 100 - percent_diff
    # --------------------------------------------------------------------------------------------------------
    # fill in basic information concerning the core genes
    genes = fill_in_basics('nt')
    fill_in_basics('aa')

    # open AllFAA.all to extract predicted functions for each core gene
    # if you run into a core gene, read its predicted function (but can
    # ignore the sequences themselves)
    with open('AllFAA.all') as file:
        i = 0
        length = len(table.keys())
        for line in file:
            # we're done; we included all core genes
            if i >= length:
                break
            # check if the gene is a core gene, if so, add it to
            # the 'table' dictionary
            if line.startswith('>'):
                if genes[i][:-4] in line:
                    line = line.strip()
                    p = re.compile(">([^ ]*) (.*)")
                    m = p.match(line)
                    gene, predicted_function = m.groups()
                    table[str(i)]['predicted_function'] = predicted_function
                    i += 1

    # calculate the number of snps/aa mutations for each core gene
    fill_in_snps('nt')
    fill_in_snps('aa')

    # write to a json file that the page will load to save time
    with open('table.json', 'w') as out:
        json.dump(table, out, indent=4)
    return table


def prepare_snp_data(folder_name, table, download_all):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    # Download all the core gene alignments, as well as their snps
    if download_all:
        process = subprocess.Popen("zip snp_data.zip -r usearch_SNPanalysis_ffn_files", shell=True)
        process.wait()
        process = subprocess.Popen("zip snp_data.zip -r usearch_SNPanalysis_faa_files", shell=True)
        process.wait()
        process = subprocess.Popen("zip snp_data.zip usearch_CORE_ffn.All_Genes.vcf", shell=True)
        process.wait()
        process = subprocess.Popen("zip snp_data.zip usearch_CORE_faa.All_Genes.vcf", shell=True)
        process.wait()
    # only download selected alns and snps
    else:
        to_process = request.form.getlist("check")
        if len(to_process) == 0:
            raise ValueError("You must specify at least one gene")
        genes = [table[index]['name'] for index in to_process]

        # write a subset of all_genes.vcf
        inputs = ['usearch_CORE_faa.All_Genes.vcf', 'usearch_CORE_ffn.All_Genes.vcf']
        outputs = ['usearch_CORE_faa.Selected_Genes.vcf', 'usearch_CORE_ffn.Selected_Genes.vcf']

        for i in range(2):
            with open(inputs[i]) as f:
                with open(outputs[i], 'w') as out:
                    for line in f:
                        line = line.lstrip()
                        # only lines that contain our genes of interest will be copied
                        if line.startswith('>'):
                            p = re.compile('>([^ ]*) (.*)')
                            m = p.match(line)
                            gene_name = m.groups()[0]
                            if gene_name in genes:
                                out.write(line)
                                copy_lines = True
                            else:
                                copy_lines = False
                        # actual snp data
                        elif line.startswith('1') and copy_lines:
                            out.write(line)
                        # column headers
                        elif line.startswith('#') and copy_lines:
                            out.write(line)

        # include raw .fas files for ease of view via SeaView
        for gene in genes:
            process = subprocess.Popen("zip snp_data.zip -r usearch_SNPanalysis_ffn_files/{}.fas".format(gene), shell=True)
            process.wait()
            process = subprocess.Popen("zip snp_data.zip -r usearch_SNPanalysis_faa_files/{}.fas".format(gene), shell=True)
            process.wait()
        process = subprocess.Popen("zip snp_data.zip usearch_CORE_ffn.Selected_Genes.vcf", shell=True)
        process.wait()
        process = subprocess.Popen("zip snp_data.zip usearch_CORE_faa.Selected_Genes.vcf", shell=True)
        process.wait()

    # move snp_data.zip to the output folder
    shutil.move("snp_data.zip", os.path.join(OUTPUT_FOLDER, "snp_data.zip"))


# Prepares the snp_table as an .xlsx for downloading
def prepare_snp_table(folder_name, table):
    # temporarily convert table's keys from str->int
    # (this will help with creating the data frame)
    table = {int(k):v for k,v in table.items()}
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    df = pd.DataFrame.from_dict(table)

    # cols are entries (ie. core genes), and rows are attributes (such as
    # name, predicted function); need to transpose and sort
    df = df.transpose()
    df = df.sort_index()

    # swap columns, new order should be:
    # gene name, predicted function, snp_count, length_nt, nt_sim, mutations, length_aa, aa_sim
    col_list = list(df)
    (col_list[7], col_list[6], col_list[3], col_list[5], 
    col_list[0], col_list[4], col_list[1], col_list[2]) = col_list
    df = df[col_list]

    # Prettify the column names
    df.columns = ['Gene name', 'Predicted function', '# SNPs', 'Gene length (bp)',
                  '% identity (nt)', '# A.A. mutations', 'Protein length (aa)',
                  '% identity (aa)']

    wb = xlsxwriter.Workbook('snp_table.xlsx')
    ws = wb.add_worksheet()

    # write column headers
    for i, colname in enumerate(df.columns):
        ws.write(0, i, colname)

    # keep track of longest predicted function, and
    # the length of the (first) gene name
    length_gene = 0
    length_pf = 0

    # write contents
    for index, row in df.iterrows():
        if len(row['Predicted function']) > length_pf:
            length_pf = len(row['Predicted function'])
        if len(row['Gene name']) > length_gene:
            length_gene = len(row['Gene name'])
        ws.write(index+1, 0, row['Gene name'])
        ws.write(index+1, 1, row['Predicted function'])
        ws.write(index+1, 2, row['# SNPs'])
        ws.write(index+1, 3, row['Gene length (bp)'])
        ws.write(index+1, 4, row['% identity (nt)'])
        ws.write(index+1, 5, row['# A.A. mutations'])
        ws.write(index+1, 6, row['Protein length (aa)'])
        ws.write(index+1, 7, row['% identity (aa)'])

    # my best attempt at an autofit
    ws.set_column(0, 0, length_gene)
    ws.set_column(1, 1, length_pf)
    wb.close()

    # transfer the file over to the 'outputs' folder
    shutil.move('snp_table.xlsx', os.path.join(OUTPUT_FOLDER, 'snp_table.xlsx'))


# updates the .fas (both ffn and faa) so that the new ref appears at the top
# also updates the .vcf files so that snps are relative to the new ref
# current_references.json keeps track of the ref for each run
def update_reference(folder_name, new_ref, all_refs):
    update_fas_from_reference(folder_name, new_ref)
    update_vcf_from_reference(folder_name, new_ref)
    all_refs[folder_name] = new_ref
    with open(os.path.join(DATA_FOLDER, 'current_references.json'), 'w') as out:
        json.dump(all_refs, out)


def update_fas_from_reference(folder_name, new_ref):
    # start with the faa files
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    # list of 3-tuples, consisting of filenames, datatype, and new reference genome name
    faa_files = [(f , 'aa', new_ref) for f in os.listdir('usearch_SNPanalysis_faa_files')]
    ffn_files = [(f , 'nt', new_ref) for f in os.listdir('usearch_SNPanalysis_ffn_files')]
    files = faa_files + ffn_files

    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    p.map(worker_update_fas, files)


def worker_update_fas(file):
    filename = file[0]
    datatype = file[1]
    new_ref = file[2]

    if datatype == 'aa':
        file_dir = os.path.join("usearch_SNPanalysis_faa_files", filename)
    elif datatype == 'nt':
        file_dir = os.path.join("usearch_SNPanalysis_ffn_files", filename)

    seqs = {}
    name = None

    with open(file_dir) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = seq
                name = line.rstrip()
                seq = ''
            else:
                seq += line
        # accounting for the last entry
        seqs[name] = seq

    # overwrite the existing faa_file    
    with open(file_dir, 'w') as out:
        # start with the reference genome and its sequence
        out.write('>' + new_ref + '\n')
        out.write(seqs['>' + new_ref])
        # write the other genomes' sequences, in alphabetical order
        for genome in sorted(seqs.keys()):
            # skip the reference; already done
            if genome == '>' + new_ref:
                continue
            out.write(genome + '\n')
            out.write(seqs[genome])


def update_vcf_from_reference(folder_name, new_ref):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))

    def rewrite_vcf(datatype):
        if datatype == 'aa':
            filename = 'usearch_CORE_faa.All_Genes.vcf'
        elif datatype == 'nt':
            filename = 'usearch_CORE_ffn.All_Genes.vcf'

        with open(filename) as file:
            with open(filename + '.TEMP', 'w') as out:
                for line in file:
                    if line.startswith('#'):
                        # line_arr is a list of elements separated by tabs
                        line_arr = line.rstrip().split('\t')
                        old_order = line_arr[9:]  # old_order contains the former order of genomes
                        new_order = line_arr[9:]  # new_order contains the new order (alphabetical order, excluding reference)
                        new_order.remove(new_ref)
                        new_order.sort()
                        # line_arr_out is also a list, but will be joined at the end so it can be written onto the new .vcf
                        line_arr_out = line_arr[:9]  # slice off the line array to exclude all the genome names
                        line_arr_out.append(new_ref)
                        # append the rest of the genomes in alphabetical order
                        for genome in new_order:
                            line_arr_out.append(genome)
                        # join back into a string
                        line = '\t'.join(line_arr_out) + '\n'
                    elif line.startswith('1'):
                        line_arr = line.rstrip().split('\t')
                        bases = {}
                        # populate locations of old bases
                        for i in range(len(old_order)):
                            if line_arr[9+i] == '.':  # same as the old reference
                                bases[old_order[i]] = line_arr[3]
                            else:
                                bases[old_order[i]] = line_arr[9+i]
                        line_arr_out = line_arr[:9]
                        line_arr_out[3] = bases[new_ref]  # change the reference column to match the new ref genome
                        line_arr_out.append('.')  # the 9th column is always the reference
                        for genome in new_order:
                            if bases[genome] == bases[new_ref]:
                                line_arr_out.append('.')
                            else:
                                line_arr_out.append(bases[genome])
                        # updating the ALT column; use a set to avoid duplicates
                        alt = set()
                        for genome, base in bases.items():
                            if base != bases[new_ref]:
                                alt.add(base)
                        alt = ','.join(list(alt))
                        line_arr_out[4] = alt
                        line = '\t'.join(line_arr_out) + '\n'
                    out.write(line)
        # Delete original vcf, and rename the new ones
        os.remove(filename)
        os.rename(filename + '.TEMP', filename)

    rewrite_vcf('nt')
    rewrite_vcf('aa')
# =====================================================================================
# SNP Mutations - displays table of snps/aa mutations, as well as buttons for viewing
# sortable bar graphs
@app.route('/snp_mutations/<folder_name>/<gene_name>', methods=['GET', 'POST'])
def snp_mutations(folder_name, gene_name):
    form = DownloadForm()

    # get the reference genome
    os.chdir(DATA_FOLDER)
    with open('current_references.json') as j:
        all_refs = json.load(j)
    ref_genome = all_refs[folder_name]

    # make folders for run date & ref genome
    if not os.path.isdir(GRAPHS_FOLDER):
        os.mkdir(GRAPHS_FOLDER)
    os.chdir(GRAPHS_FOLDER)
    if folder_name not in os.listdir('.'):
        os.mkdir(folder_name)
    os.chdir(folder_name)
    if ref_genome not in os.listdir('.'):
        os.mkdir(ref_genome)
    os.chdir(ref_genome)

    # read mutation data in [gene name].json
    if '{}.json'.format(gene_name) in os.listdir('.'):
        nt_mutations, aa_mutations, alts = read_snps_json(gene_name)
    # calculate mutations, and save to .tsv (for d3) and .json (for the app)
    else:
        # Find transition/transversion counts
        nt_mutations = count_snps(folder_name, gene_name, ref_genome)
        os.chdir(os.path.join(GRAPHS_FOLDER, folder_name, ref_genome))

        # Find conservative/non-conservative counts
        # need to use blast for aa, because blast has '+' and ' ', which
        # roughly indicates # of cons/non-cons mutations
        aa_mutations = count_aa_mutations(folder_name, gene_name, ref_genome)
        os.chdir(os.path.join(GRAPHS_FOLDER, folder_name, ref_genome))

        # create tsvs required by d3
        write_d3_tsvs(gene_name, nt_mutations, aa_mutations)
        # create json of snps/aa mutations
        write_snps_json(gene_name, nt_mutations, aa_mutations)

        # list of non-ref genomes in sorted order
        alts = sorted([alt for alt in nt_mutations.keys()]) 

    if request.method == 'POST':
        # download the excel table
        if form.download.data:
            prepare_mutation_table(folder_name, ref_genome, gene_name, nt_mutations, aa_mutations)
            return redirect(url_for('download', filename="output/{0}.xlsx".format(gene_name)))

    return render_template('snp_mutations.html', form=form, nt_mutations=nt_mutations, aa_mutations=aa_mutations,
        gene_name=gene_name, alts=alts, ref_genome=ref_genome)


# creates an excel table using pandas, and moves it to the output folder for serving to client
def prepare_mutation_table(folder_name, ref_genome, gene_name, nt_muts, aa_muts):
    df_nt = pd.DataFrame.from_dict(nt_muts)
    df_aa = pd.DataFrame.from_dict(aa_muts)
    df_nt = df_nt.transpose()
    df_aa = df_aa.transpose()

    # Add new column "Plot label"
    # len(df_nt) refers to the length of the index (i.e. how many alt genomes)
    df_nt['Plot Label'] = pd.Series(range(1,len(df_nt)+1), index=df_nt.index)
    df_nt.reset_index(level=0, inplace=True)
    # Prettify column names
    df_nt = df_nt.rename(columns = {'transitions': 'Transitions', 'transversions': 'Transversions', 'index': 'Genome'})
    # Reorder columns
    df_nt = df_nt[['Plot Label', 'Genome', 'Transitions', 'Transversions']]

    df_aa.reset_index(level=0, inplace=True)
    df_aa = df_aa.rename(columns = {'conservative': 'Conservative', 'non_conservative': 'Non-conservative'})
    df_aa = df_aa[['Conservative', 'Non-conservative']]

    # join df_nt and df_aa together
    df = df_nt.join(df_aa)

    # Save as excel
    wb = xlsxwriter.Workbook('{0}.xlsx'.format(gene_name))
    ws = wb.add_worksheet()

    # write column headers
    ws.write(0, 0, "Core Gene")
    ws.write(0, 1, "Reference Genome")
    ws.write(1, 0, gene_name)
    ws.write(1, 1, ref_genome)

    for i, colname in enumerate(df.columns):
        ws.write(3, i, colname)

    # write contents
    for index, row in df.iterrows():
        ws.write(index+4, 0, row['Plot Label'])
        ws.write(index+4, 1, row['Genome'])
        ws.write(index+4, 2, row['Transitions'])
        ws.write(index+4, 3, row['Transversions'])
        ws.write(index+4, 4, row['Conservative'])
        ws.write(index+4, 5, row['Non-conservative'])

    wb.close()
    shutil.move('{0}.xlsx'.format(gene_name), os.path.join(OUTPUT_FOLDER, '{0}.xlsx'.format(gene_name)))


# counts snps at a core gene, for each non-reference genome
# returns a dict with keys: genome name, values: snp counts (transition/transversion)
def count_snps(folder_name, gene_name, ref_genome):
    goto_folder(folder_name)
    os.chdir('Core_Genome/usearch_SNPanalysis_ffn_files')
    mutations = {}
    alts = []

    # do we want to do a whole core-genome analysis?
    if gene_name == 'all':
        all_genes = glob.glob('*.fas')
    else:
        all_genes = [gene_name + '.fas']

    for i, gene_file in enumerate(all_genes):
        # create a dictionary
        seqs = get_core_gene_sequences(gene_file)
        if i == 0:
            # initialize alts and mutations
            alts = sorted([genome for genome in seqs.keys() if genome != ref_genome])
            mutations = {
                alt: {'transitions': 0, 
                      'transversions': 0} 
                for alt in alts}
        ref_seq = seqs[ref_genome]

        # Compare the reference seq to every other seq
        # note: does not account for gaps
        for j, ref_base in enumerate(ref_seq):
            ref_base_type = base_type(ref_base)
            for alt in alts:
                alt_base = seqs[alt][j] 
                if alt_base != ref_base:  # SNP found
                    alt_base_type = base_type(alt_base)
                    # compare alt_base type and ref_base_type
                    # they both have to be one of the 4 recognized bases; i.e.
                    # '-' or 'N' will not be considered when counting the type of snp
                    if ref_base_type and alt_base_type:
                        if ref_base_type == alt_base_type:
                            mutation_type = 'transitions'
                        else:
                            mutation_type = 'transversions'
                        mutations[alt][mutation_type] += 1
    return mutations


# for one .fas file, fills out a dictionary with keys as genomes names
# and values as the core gene sequence for that genome
def get_core_gene_sequences(gene_file):
    curr_genome = ''
    seqs = {}
    with open(gene_file) as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('>'):  # new sequence
                if curr_genome:
                    seqs[curr_genome] = seq  # add previous sequence to seqs
                curr_genome = line[1:]  # ignore '>'
                seq = ''
            else:
            #elif line[0] in '-ACGT':  # is a base (or no base, if '-')
                seq += line
    seqs[curr_genome] = seq  # account for last genome
    return seqs


# similar to count_snps()
# gene_name can be 'all', or a particular core gene
def count_aa_mutations(folder_name, gene_name, ref_genome):
    # use the blosum substitution matrix
    os.chdir(APP_FOLDER)
    with open('matrix.json') as j:
        matrix = json.load(j)

    goto_folder(folder_name)
    os.chdir('Core_Genome/usearch_SNPanalysis_faa_files')

    aa_mutations = {}

    if gene_name == 'all':
        all_genes = glob.glob('*.fas')
    else:
        all_genes = [gene_name + '.fas']

    for i, gene_file in enumerate(all_genes):
        seqs = get_core_gene_sequences(gene_file)
        if i == 0:
            # initialize alts and mutations
            alts = sorted([genome for genome in seqs.keys() if genome != ref_genome])
            aa_mutations = {
                alt: {'conservative': 0, 
                      'non_conservative': 0} 
                for alt in alts}
        ref_seq = seqs[ref_genome]

        for j, ref_base in enumerate(ref_seq):
            if ref_base == '-':
                continue
            for alt in alts:
                alt_base = seqs[alt][j] 
                if alt_base == '-':
                    continue
                if alt_base != ref_base:
                    # look up using the matrix dictionary
                    sub = ref_base + alt_base
                    # if the key doesn't exist, try swapping aa
                    if sub in matrix.keys():
                        mutation_type = matrix[sub]
                    else:
                        mutation_type = matrix[sub[1] + sub[0]]

                    if mutation_type == 1:
                        aa_mutations[alt]['conservative'] += 1
                    elif mutation_type == -1:
                        aa_mutations[alt]['non_conservative'] += 1

    return aa_mutations


def write_d3_tsvs(gene_name, nt_mutations, aa_mutations):
    out1 = open("{}_transitions.tsv".format(gene_name), 'w')
    out2 = open("{}_transversions.tsv".format(gene_name), 'w')
    out3 = open("{}_cons.tsv".format(gene_name), 'w')
    out4 = open("{}_noncons.tsv".format(gene_name), 'w')

    out1.write("Genome\tMutations\n")
    out2.write("Genome\tMutations\n")
    out3.write("Genome\tMutations\n")
    out4.write("Genome\tMutations\n")

    sorted_keys = sorted(nt_mutations.keys())

    for i, genome in enumerate(sorted_keys):
        i = str(i+1)
        counts_nt = nt_mutations[genome]
        counts_aa = aa_mutations[genome]
        transitions, transversions = counts_nt['transitions'], counts_nt['transversions']
        out1.write(i + '\t' + str(transitions) + '\n')
        out2.write(i + '\t' + str(transversions) + '\n')
        cons, noncons = counts_aa['conservative'], counts_aa['non_conservative']
        out3.write(i + '\t' + str(cons) + '\n')
        out4.write(i + '\t' + str(noncons) + '\n')


def write_snps_json(gene_name, nt_mutations, aa_mutations):
    combined = {}
    for genome in aa_mutations.keys():
        combined[genome] = {'transitions': nt_mutations[genome]['transitions'],
                            'transversions': nt_mutations[genome]['transversions'],
                            'conservative': aa_mutations[genome]['conservative'],
                            'non_conservative': aa_mutations[genome]['non_conservative']}
    with open('{}.json'.format(gene_name), 'w') as out:
        json.dump(combined, out, indent=4)


# Determines if the nucleotide is a pyrimidine, purine, or None (-)
def base_type(base):
    if base == 'A' or base == 'G':
        return 'Purine'
    elif base == 'C' or base == 'T':
        return 'Pyrimidine'
    else:
        return None
# =====================================================================================
# SNP Graph: sortable bar graphs of SNPs or amino acid mutations, powered by the d3 library
@app.route('/snp_graph/<folder_name>/<gene_name>/<datatype>')
def snp_graph(folder_name, gene_name, datatype):
    # get the reference genome
    os.chdir(DATA_FOLDER)
    with open('current_references.json') as j:
        all_refs = json.load(j)
    ref_genome = all_refs[folder_name]

    # read either data1.tsv and data2.tsv, or data3.tsv and data4.tsv
    os.chdir(os.path.join(GRAPHS_FOLDER, folder_name, ref_genome))

    path = "graphs/%s/%s" % (folder_name, gene_name)
    tsv_files = ["{}_transitions.tsv".format(gene_name),
                 "{}_transversions.tsv".format(gene_name),
                 "{}_cons.tsv".format(gene_name),
                 "{}_noncons.tsv".format(gene_name)]
    filenames = ["graphs/%s/%s/%s" % (folder_name, ref_genome, tsv) for tsv in tsv_files]

    # get the alt_genome count by opening an arbitrary .tsv
    count = -1  # accounting for file header
    with open(tsv_files[0]) as f:
        for line in f:
            count += 1

    # take a subset of classes
    if datatype == 'nucleotide':
        classes = [".one", ".two"]
    elif datatype == 'amino_acid':
        classes = [".three", ".four"]

    # open up the json, and 
    nt_mutations, aa_mutations, alts = read_snps_json(gene_name)

    return render_template('snp_graph.html', filenames=filenames, datatype=datatype,
        count=count, classes=classes, alts=alts, ref_genome=ref_genome,
        gene_name=gene_name)


def read_snps_json(gene_name):
    with open('{}.json'.format(gene_name)) as f:
        data = json.load(f)
    nt_mutations = {}
    aa_mutations = {}
    alts = []
    for genome in data.keys():
        nt_mutations[genome] = {            
            'transitions': data[genome]['transitions'],
            'transversions': data[genome]['transversions']
        }
        aa_mutations[genome] = {
            'conservative': data[genome]['conservative'],
            'non_conservative': data[genome]['non_conservative']
        }
        alts.append(genome)
    alts.sort()
    return nt_mutations, aa_mutations, alts 
# =====================================================================================
# Rarefaction - running the rarefaction analysis to produce a curve that asymptotes at
# the 'true' core genome size, as the subset size is increased incrementally
@app.route('/rarefaction/<folder_name>', methods=['GET', 'POST'])
@app.route('/rarefaction', methods=['GET', 'POST'])
def rarefaction(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)

    rf_form = None
    run_name = ''

    # Initialize the form to select a run
    select_list = get_select_list()

    if folder_name:
        run_name = get_run_name(folder_name)
        rf_form = RarefactionForm()
        count = get_num_strains(folder_name)

    if request.method == 'POST':

        if rf_form.download_strains.data:
            return redirect(url_for('download', filename='data/{}/STRAINS'.format(folder_name)))
        elif rf_form.start_run.data:
            reps = rf_form.reps.data
            step = rf_form.step.data

            if rf_form.validate_on_submit():
                # make sure that step is smaller than count (# strains)
                if step > count:
                    flash ("step must be at most {}".format(count))
                    return redirect(url_for('rarefaction', folder_name=folder_name))

                strains_file = rf_form.strains.data
                goto_folder(folder_name)
                # check user-uploaded strains file
                if strains_file:
                    strains_file.save("STRAINS_LABELLED")
                    try:
                        check_labelled_strains(folder_name)
                    except ValueError as e:
                        flash(str(e))
                        os.remove("STRAINS_LABELLED")
                        return redirect(url_for('rarefaction', folder_name=folder_name))
                # create a default strains file
                else:
                    create_labelled_strains(folder_name)
                start_rarefaction_analysis(folder_name, reps, step)
                flash("The rarefaction analysis has been successfully completed!")
            else:
                if not reps or reps < 1 or reps > 10:
                    flash("reps must be between 1 and 10")
                if not step or step < 2:
                    flash("step must be at least 2")

    return render_template('rarefaction.html', rf_form=rf_form, select_list=select_list,
        folder_name=folder_name, run_name=run_name)


def create_labelled_strains(folder_name):
    with open('STRAINS') as file:
        with open ('STRAINS_LABELLED', 'w') as out:
            for line in file:
                line = line.replace('\n', '\t[A]\n')
                out.write(line)


def check_labelled_strains(folder_name):
    print t.green("Checking STRAINS_LABELLED")

    goto_folder(folder_name)
    with open('STRAINS') as file:
        STRAINS = [line.strip() for line in file]

    # check each line to verify that the genome names (and order)
    with open('STRAINS_LABELLED') as file:
        i = 0
        for line in file:
            if i >= len(STRAINS):
                raise ValueError("Labelled strains is too long")
            line = line.strip()
            # check that the line follows this pattern
            p = re.compile('(\S*)(\s+)(\[.+\])')
            m = p.match(line)
            if m == None:
                raise ValueError("Something is wrong with the labelled strains file")
            strain_name = m.groups()[0]
            # check strain name against the ith line of STRAINS
            if strain_name != STRAINS[i]:
                raise ValueError("Strain names don't match up")
            i += 1

    # make sure the two strains files have the same # lines
    if i != len(STRAINS):
        raise ValueError("Labelled strains is too short")

    print t.green("Check complete. The modified strains file is valid")


# Delete the 'improper' strains_labelled file
def delete_labelled_strains(folder_name):
    goto_folder(folder_name)
    if 'STRAINS_LABELLED' in os.listdir('.'):
        os.remove('STRAINS_LABELLED')


# fires off a bunch of consecutive hardcore analyses in a stepwise fashion,
# while recording the # of core genes in a tsv;
# this depends on HardCORE_Subset.py being in your path
def start_rarefaction_analysis(folder_name, reps, step):
    goto_folder(folder_name)
    # delete the old rarefaction results
    if 'CoreGenomeCount.tsv' in os.listdir('.'):
        os.remove('CoreGenomeCount.tsv')
    hardcore_subset_py = os.path.join(SCRIPTS_FOLDER, "HardCORE_Subset.py")
    process = subprocess.Popen("python {} --reps {} --step {}".format(hardcore_subset_py,
        reps, step), shell=True)
    process.wait()
# =====================================================================================
# Rarefaction Results - view the results of the rarefaction analysis, including a basic
# graph produced by matplotlib
@app.route('/rarefaction_results/<folder_name>', methods=['GET', 'POST'])
@app.route('/rarefaction_results', methods=['GET', 'POST'])
def rarefaction_results(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)

    rf_data = None
    filename = None
    download_form = None
    select_list = get_select_list()
    run_name = ''

    if folder_name:
        goto_folder(folder_name)
        # see if rarefaction was completed by looking for CGC.tsv
        if 'CoreGeneCount.tsv' in os.listdir('.'):
            # if so, then check if the rf_graph, table, etc. were already made.
            # if so, render; otherwise, make from new
            run_name = get_run_name(folder_name)
            # parse CoreGeneCount.tsv
            rf_data = read_core_gene_count(folder_name)
            # if the graph hasn't been generated yet, creat it using matploblib
            if '{0}.jpg'.format(folder_name) not in os.listdir('.'):
                draw_rarefaction_graph(folder_name, rf_data)
            # location of the rarefaction img
            filename = 'data/{0}/{0}.jpg'.format(folder_name)
            download_form = DownloadForm()
        else:
            flash(("It appears you haven't yet run the rarefaction. "
                "Visit the Rarefaction page "
                "and perform the analysis for this run, then come back!"))
            return redirect(url_for('rarefaction_results'))

    if request.method == 'POST':
        if download_form.download.data:
            # package the graph and excel data into a .zip
            prepare_rf_zip(folder_name, rf_data)
            return redirect(url_for('download', filename='output/rarefaction.zip'))

    return render_template('rarefaction_results.html', folder_name=folder_name, rf_data=rf_data,
        select_list=select_list, filename=filename, download_form=download_form, run_name=run_name)


# if CoreGeneCount.tsv exists, then parse it for visualization
# otherwise, throw an exception
def read_core_gene_count(folder_name):
    rf_data = {}
    goto_folder(folder_name)
    if 'CoreGeneCount.tsv' not in os.listdir('.'):
        raise ValueError("CoreGeneCount.tsv not found")
    with open('CoreGeneCount.tsv') as file:
        for line in file:
            line = line.strip()
            line_split = line.split('\t')
            size = int(line_split[0])  # the "step" size
            data = [int(val) for val in line_split[1:]]  #  measurements
            rf_data[size] = {'data': data}
            # round mean/stdev to 2 decimal places
            rf_data[size]['mean'] = round(np.mean(data), 2)
            rf_data[size]['stdev'] = round(np.std(data), 2)
    return rf_data


# uses the matplotlib library to draw a rarefaction graph
# uses standard deviation for error bars
def draw_rarefaction_graph(folder_name, rf_data):
    goto_folder(folder_name)
    x = []
    y = []
    y_err = []

    for size in sorted(rf_data.keys()):
        x.append(size)
        y.append(rf_data[size]['mean'])
        y_err.append(rf_data[size]['stdev'])

    plt.scatter(x, y)
    fig, ax = plt.subplots()
    ax.errorbar(x, y, yerr=y_err)

    # Labels and titles
    plt.xlabel('Number of Genomes')
    plt.ylabel('Number of Core Genes')
    plt.title('Rarefaction curve')

    plt.savefig('{0}.jpg'.format(folder_name))


# prepares a zip of the graph and results (in table format)
def prepare_rf_zip(folder_name, rf_data):
    prepare_rf_table(folder_name, rf_data)
    process = subprocess.Popen("zip rarefaction.zip {0}.jpg rarefaction_table.xlsx".format(folder_name), shell=True)
    process.wait()
    os.remove("rarefaction_table.xlsx")
    shutil.move("rarefaction.zip", os.path.join(OUTPUT_FOLDER, 'rarefaction.zip'))


def prepare_rf_table(folder_name, rf_data):
    goto_folder(folder_name)
    # split up data nicely (data, sd, mean)
    rf_mean = {size: rf_data[size]['mean'] for size in rf_data.keys()}
    rf_stdev = {size: rf_data[size]['stdev'] for size in rf_data.keys()}
    rf_data = {size: rf_data[size]['data'] for size in rf_data.keys()}

    df1 = pd.DataFrame.from_dict(rf_data)
    df2 = pd.DataFrame({'Mean': rf_mean, 'SD': rf_stdev})

    df1 = df1.transpose()

    df = df1.join(df2)
    df.reset_index(level=0, inplace=True)

    writer = pd.ExcelWriter('rarefaction_table.xlsx')
    df.to_excel(writer, 'Sheet1', index=False)
    writer.save()
# =====================================================================================
# Core Duplicates - find all duplicates of core genes, based on thresholds specified by the user
@app.route('/core_duplicates/<folder_name>', methods=['GET', 'POST'])
@app.route('/core_duplicates', methods=['GET', 'POST'])
def core_duplicates(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)

    duplicates = {}
    run_date = ''
    max_cols = 0
    core_dup_form = None
    plen = 0
    ident = 0
    row_count = 0
    family_count = 0
    run_name = ''
    select_list = get_select_list()

    if folder_name:
        run_name = get_run_name(folder_name)
        core_dup_form = CoreDuplicateForm()
        core_genes = initialize_core_genes(folder_name)
        run_date = prettify_date(folder_name)

        # check if we have already done a core duplicate analysis
        # (if we have produced a .json file, then yes)
        os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
        if 'core_duplicates.json' in os.listdir('.'):
            from_json = core_duplicates_read_json(folder_name)
            plen = float(from_json['plen'])
            ident = float(from_json['ident'])
            duplicates = from_json['duplicates']
            max_cols = int(from_json['max_cols'])
            row_count = int(from_json['row_count'])
            family_count = int(from_json['family_count'])

    if request.method == 'POST':
        # updating ident/plen parameters
        if core_dup_form.submit.data:
            if core_dup_form.validate_on_submit():
                ident = core_dup_form.ident.data
                plen = core_dup_form.plen.data
                # finding core duplicates that meet the user-specified threshold
                duplicates = find_core_duplicates(folder_name, ident, plen, core_genes)
                max_cols = get_max_columns(duplicates)

                # save results to .json to load in future
                family_count = len(duplicates.keys())
                for k, v in duplicates.items():
                    row_count += len(v['core_gene'])
                to_json = {
                    'ident': str(ident),
                    'plen': str(plen),
                    'duplicates': duplicates,
                    'max_cols': max_cols,
                    'row_count': row_count,
                    'family_count': family_count
                }
                core_duplicates_write_json(folder_name, to_json)
            else:
                if core_dup_form.ident.data > 1.0 or core_dup_form.ident.data < 0.5:
                    print "X"
                    flash("ident must be between 0.5 and 1.0")
                if core_dup_form.plen.data > 1.0 or core_dup_form.plen.data < 0.5:
                    flash("plen must be between 0.5 and 1.0")

                return redirect(url_for('core_duplicates', folder_name=folder_name))

        elif core_dup_form.download.data:
            prepare_duplicates_table(folder_name, duplicates)
            return redirect(url_for('download', filename='output/duplicates.xlsx'))

    return render_template('core_duplicates.html', select_list=select_list,
        core_dup_form=core_dup_form, plen=plen*100, ident=ident*100,
        run_date=run_date, duplicates=duplicates, max_cols=max_cols,
        run_name=run_name, row_count=row_count, family_count=family_count)


# returns a dict where keys: core gene family number,
# values: dict with keys as genome names and values as prokka id for this core gene
def initialize_core_genes(folder_name):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    core_genes = {}
    with open('CoreGenome.usearch') as f:
        for i, line in enumerate(f):
            core_gene = {}
            core_gene_arr = line.strip().split('\t')
            for gene in core_gene_arr:
                genome, gene_number = find_genome_regexp(gene)
                core_gene[genome] = gene_number
            core_genes[i+1] = core_gene
    return core_genes

def find_core_duplicates(folder_name, ident, plen, core_genes):
    # First, parses the files in the Cluster folder to see if there are "duplicated" genes that have
    # been removed from the analysis due to being too similar to the core gene
    
    # initialize empty dictionary
    tmp_duplicates = {i: {} for i in range(1, len(core_genes)+1)}

    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Clusters'))
    cluster_files = glob.glob('*')
    for cluster_file in cluster_files:
        with open(cluster_file) as f:
            for i, line in enumerate(f):
                line_arr = line.split('\t')
                first, second = line_arr[0], line_arr[1]
                # use regex to parse the confusing strings to obtain gene id's
                p = re.compile(r'gi\|(.*)_(\d+)\|ref\|1234\|')
                m1 = p.match(first)
                m2 = p.match(second)
                genome = m1.groups()[0]
                query = m1.groups()[1]      # the duplicate core gene (possibly)
                target = m2.groups()[1]     # the core gene itself
                # convert % to decimal (e.g. 90.0 -> 0.9)
                align_percent = float(line_arr[2]) / 100

                # check if the two genes are similar enough
                if align_percent >= ident:
                    query_len = int(line_arr[7])
                    target_len = int(line_arr[9])
                    max_len = max(query_len, target_len)
                    align_plen = float(line_arr[3]) / max_len

                    # check if the genes overlap considerably
                    if align_plen >= plen:                    
                        for i in range(1, len(core_genes.keys())+1):
                            # now check if this duplicated gene is a core gene
                            if target in core_genes[i][genome]:
                                # check if we need to make a new list for the entry
                                if genome in tmp_duplicates[i].keys():
                                    # append a tuple e.g. (01054, 94.5) where 01054 is the query #, and 94.5 is the % similarity
                                    tmp_duplicates[i][genome].append((query, line_arr[2]))  
                                else:
                                    tmp_duplicates[i][genome] = [(query, line_arr[2])]

    # we have a bunch of keys that do not have values, we want to remove those
    for i, core_gene in tmp_duplicates.items():
        if not core_gene:
            tmp_duplicates.pop(i)

    for i, core_gene in tmp_duplicates.items():
        for genome in core_gene.keys():
            core_gene_number = core_genes[i][genome]
            dup_core_list = [core_gene_number]
            dup_core_list.extend(core_gene[genome])
            core_gene[genome] = dup_core_list

    duplicates = {}
    # adding predicted function to each core gene that has a duplicate
    for i, core_gene in tmp_duplicates.items():
        duplicates[i] = {}
        duplicates[i]['core_gene'] = core_gene

        # go into PROKKA_BACKUPS to find predicted functions
        os.chdir(os.path.join(DATA_FOLDER, folder_name, 'PROKKA_backups'))
        first_genome = sorted(tmp_duplicates[i].keys())[0]

        # look up the core gene number
        gene_number = tmp_duplicates[i][first_genome][0]
        with open(first_genome + '.orig.faa') as f:
            predicted_function = find_predicted_function(f, gene_number)
        duplicates[i]['predicted_function'] = predicted_function

    return duplicates


# by passing in a prokka faa file and the core gene number, get its predicted function
# example: find_predicted_function(f, '00006') -> outputs 'Uronate isomerase'
def find_predicted_function(f, gene_number):
    for line in f:
        if gene_number in line:
            # try to look for the pattern
            # >*********_##### predicted function
            # where the ##### is the gene number
            p = re.compile(r'>.*_(\d+) (.*)')
            m = p.match(line.strip())
            line_number = m.groups()[0]
            line_pf = m.groups()[1]
            # check if the gene number at this line matches
            if line_number == gene_number:
                return line_pf


# returns the max # of columns needed to display all duplicate genes
# (not counting the core gene column)
def get_max_columns(duplicates):
    max_cols = 1
    for i,core_gene in duplicates.items():
        for genome, genes in core_gene['core_gene'].items():
            if len(genes)-1 > max_cols:
                max_cols = len(genes)-1
    return max_cols


def core_duplicates_write_json(folder_name, to_json):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    with open('core_duplicates.json', 'w') as out:
        json.dump(to_json, out)


def core_duplicates_read_json(folder_name):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    with open('core_duplicates.json') as f:
        contents = json.load(f)
    return contents
   

# Writes the table of core duplicates to an xlsx
def prepare_duplicates_table(folder_name, duplicates):
    os.chdir(DATA_FOLDER)
   
    col1 = []           # for core gene family
    col2 = []           # for genomes w dup core genes
    col3 = []           # predicted function
    rest = [[], []]     # core gene + duplicated gene id's
    j = 0               # determining how many Nones to insert(?)

    gene_families = sorted([int(k) for k in duplicates.keys()])

    for k in gene_families:
        k = str(k)
        # look into how many genomes are part of this core gene family
        genomes = sorted(duplicates[k]['core_gene'].keys())
        num_genomes = len(genomes)

        col1.append(k)
        col3.append(duplicates[k]['predicted_function'])
        for i in range(num_genomes-1):
            col1.append(None)
            col3.append(None)
        col2.extend(genomes)

        # we don't know how many columns there should be for "rest"
        # so we need to look in each genome to find the longest entry
        for genome in genomes:
            entries = duplicates[k]['core_gene'][genome]

            # if there are more entries than there are number of columns, need to create
            # additional columns and prepend with an appropriate amount of Nones (empty cells)
            if len(entries) > len(rest):
                diff = len(entries) - len(rest)
                for i in range(diff):
                    # new lists with 'j' many Nones
                    rest.append([None for none in range(j)])  

            for i, entry in enumerate(entries):
                if i == 0:
                    # rest[0] is the core gene itself, so don't need a %
                    rest[i].append(entry)
                else:
                    # rest[1:] are duplicates; write % similarity
                    rest[i].append(entry[0] + " (" + entry[1] + "%)")

            if len(entries) < len(rest):
                diff = len(rest) - len(entries)
                for i in range(diff):
                    rest[len(entries)+i].append(None)

            j += 1

    df = pd.DataFrame(
      {'Core Gene Family': col1,
      'Genomes with duplicate core genes': col2,
      'Predicted Function': col3,
      'Core gene number': rest[0],
      'Duplicate gene number(s)': rest[1]
      })

    # rearrange columns
    df = df[['Core Gene Family', 'Genomes with duplicate core genes', 'Predicted Function',
        'Core gene number', 'Duplicate gene number(s)']]

    # should account for rest[2:] (if they exist)
    if len(rest) > 2:
        for i in range(2, len(rest)):
            series = pd.Series(rest[i])
            df = pd.concat([df, series], axis=1)

    writer = pd.ExcelWriter('duplicates.xlsx')
    df.to_excel(writer,'Sheet1', index=False)
    writer.save()

    shutil.move('duplicates.xlsx', os.path.join(OUTPUT_FOLDER, 'duplicates.xlsx'))
# =====================================================================================
# Tree Builder - builds phylogenetic trees from either fasta files or core genes (including
# the entire core genome)
@app.route('/core_tree_builder/<folder_name>', methods=['GET', 'POST'])
@app.route('/core_tree_builder', methods=['GET', 'POST'])
@app.route('/tree_builder', methods=['GET', 'POST'])
def tree_builder(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)

    tree_form = None
    table = {}  # same table as in snp_viewer
    tree_props = {}
    select_list = []
    run_name = ''
    # rule <- "tree_builder" or "core_tree_builder"
    rule = str(request.url_rule) 
    rule = rule.split('/')[1]  

    if rule == 'tree_builder':
        tree_form = TreeFormWithUpload()
        os.chdir(TREE_FOLDER)
        tree_props = check_for_tree()
    
    elif rule == 'core_tree_builder':
        # only want to select a run if we are doing core tree building
        select_list = get_select_list()
        if folder_name:
            run_name = get_run_name(folder_name)
            tree_form = TreeForm()
            unzip_fas(folder_name)  # only if we have to
            os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
            if 'table.json' in os.listdir('.'):
                with open('table.json') as file:
                    table = json.load(file)
            else:
                table = create_snp_json(folder_name)
            if 'Tree' in os.listdir('.'):
                os.chdir('Tree')
                tree_props = check_for_tree(folder_name)

    if request.method == 'POST':
        if tree_form.submit.data:
            if tree_form.validate_on_submit():
                datatype = tree_form.datatype.data
                N = tree_form.num_runs.data
                algorithm = tree_form.algorithm.data
                multithread = tree_form.multithread.data
                merge = tree_form.merge.data
                builder = tree_form.builder.data

                if rule == 'core_tree_builder':
                    to_process = request.form.getlist("check")
                    # list of core genes is empty
                    if not to_process:
                        flash("You must select at least one gene")
                        return redirect(url_for('tree_builder', folder_name=folder_name))
                    else:
                        gene_list = [table[str(i)]['name'] for i in to_process]

                # validate user parameters
                try:
                    check_tree_parameters(datatype, algorithm, N, builder)
                except ValueError as e:
                    flash(str(e))
                    return redirect(url_for('tree_builder', folder_name=folder_name))

                # set default values, if unspecified
                if not algorithm:
                    if datatype == 'nt':
                        algorithm = DEFAULT_NT_ALGORITHM
                    elif datatype == 'aa':
                        algorithm = DEFAULT_AA_ALGORITHM
                if not N:
                    N = 10

                # upload, check, and parse the sequence files
                if rule == 'tree_builder':
                    try:
                        #prepare_tree_input(tree_form, merge)
                        validate_tree_input(tree_form)
                        gene_list = get_genelist()
                        prepare_tree_input(merge)
                    except ValueError as e:
                        flash(str(e))
                        return redirect(url_for('tree_builder', folder_name=folder_name))
                # obtain core gene sequences from the Core_Genome folder
                elif rule == 'core_tree_builder':
                    try:
                        prepare_core_tree_input(folder_name, datatype, gene_list, merge)
                    except ValueError as e:
                        flash(str(e))
                        return redirect(url_for('tree_builder', folder_name=folder_name))

                try:
                    # start building the tree
                    build_tree(algorithm, N, builder, multithread, merge, folder_name, datatype)

                    # generate a dictionary of tree-related properties
                    generate_tree_props(N, datatype, gene_list, builder, merge, algorithm, folder_name)
                    flash("Tree construction complete!")
                except OSError:
                    flash("ERROR: Two or fewer species (possibly as a result of collapsing)")
                return redirect(url_for('tree_builder', folder_name=folder_name))

            else:
                # user never specified a datatype
                if tree_form.datatype.data == "None":
                    flash("Must specify one of Amino acid/Nucleotide")
                return redirect(url_for('tree_builder', folder_name=folder_name))

        elif tree_form.dendroscope.data:
            try:
                dendroscope_view(folder_name)
            except ValueError:
                flash("Unable to find tree. Please generate a tree using 'Build Tree', then try again.")

        elif tree_form.ete.data or tree_form.pdf.data:
            if tree_form.ete.data:
                out_type = 'view'
            elif tree_form.pdf.data:
                out_type = 'pdf'

            core_flag = True if (rule == 'core_tree_builder') else False
            prepare_ete_tree(tree_form.show_align.data, tree_form.show_bootstrap.data,
                core_flag, out_type, folder_name)

            if tree_form.pdf.data:
                return redirect(url_for('download', filename='output/tree.pdf'))

        elif tree_form.download.data:
            try:
                prepare_tree_zip(tree_props['merge'], folder_name)
                return redirect(url_for('download', filename='output/tree.zip'))
            except ValueError:
                flash("You must prepare a tree first")  

    return render_template('tree_builder.html', tree_form=tree_form, select_list=select_list,
        table=table, rule=rule, run_name=run_name, tree_props=tree_props)


# upload the fasta/zip, and check that they are valid for building the tree
def validate_tree_input(form):
    os.chdir(TREE_FOLDER)
    for f in os.listdir('.'):
        os.remove(f)

    # get/download the uploaded file
    file = form.upload.data
    filename = secure_filename(file.filename)
    os.chdir(TREE_FOLDER)
    file.save(os.path.join(TREE_FOLDER, filename))

    # if user uploaded a zip, unzip contents
    if filename.endswith(".zip"):
        process = subprocess.Popen('unzip {}'.format(filename), shell=True)
        process.wait()
        os.remove(filename)

    # check that they are all fasta files
    files = glob.glob("*")
    for file in files:
        if not file.endswith(".fasta"):
            raise ValueError("Non-fasta file found")

    datatype = form.datatype.data
    try:
        check_tree_fastas(files, datatype)
    except ValueError as e:
        # tree fasta not validated, should delete the contents
        for f in os.listdir('.'):
            os.remove(f)
        raise e  # reraise the exception


# gets the list of genes from the Tree folder (for the generalized tree builder)
def get_genelist():
    os.chdir(TREE_FOLDER)
    gene_list = glob.glob("*.fasta")
    return gene_list


# modifies the tree files (e.g. concatenation, optionally merging nodes, getting rid of whitespace in headers)
def prepare_tree_input(merge):
    os.chdir(TREE_FOLDER)
    files = glob.glob("*.fasta")
    # once validated, get rid of all whitespace, and rename
    # extension to .fas (for fasconcat)
    for file in files:
        to_out = ''
        with open(file) as f:
            for line in f:
                if line.startswith('>'):
                    # get the bit before the first space
                    line = line.split(' ')[0]
                    line += '\n'
                to_out += line
        #need to rename the file extensions to fas
        os.remove(file)
        file = file[:-6] + '.fas'
        with open(file, 'w') as out:
            out.write(to_out)

    # Each header must be IDENTICAL
    fas_concat = os.path.join(SCRIPTS_FOLDER, "FASconCAT_v1.0.pl")
    process = subprocess.Popen('perl {} -s'.format(fas_concat), shell=True)
    process.wait()
    os.rename('FcC_smatrix.fas', 'concat_genes.fas')

    # should we clade collapse?
    if merge:
        # generates concat_genes.fas.trimmed.phy
        tree_node_merger = os.path.join(SCRIPTS_FOLDER, "TreeNodeMerger.py")
        process = subprocess.Popen('python {} concat_genes.fas -p'.format(tree_node_merger),
            shell=True)
    else:
        # generates concat_genes.phy
        fasta_2_phylip = os.path.join(SCRIPTS_FOLDER, "Fasta2Phylip.pl")
        process = subprocess.Popen('perl {} concat_genes.fas concat_genes.phy'.format(fasta_2_phylip),
            shell=True)
    process.wait()


def prepare_core_tree_input(folder_name, datatype, gene_list, merge):
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    if 'Tree' in os.listdir('.'):
        shutil.rmtree('Tree')
    os.mkdir('Tree')

    # take the faa alignments and put them together
    if datatype == 'aa':
        s = "usearch_SNPanalysis_faa_files"
    elif datatype == 'nt':
        s = "usearch_SNPanalysis_ffn_files"

    for gene in gene_list:
        shutil.copy('{0}/{1}.fas'.format(s, gene), 'Tree/{}.fas'.format(gene))
        shutil.copy('{0}/{1}.fas'.format(s, gene), 'Tree/{}.fas'.format(gene))
    os.chdir('Tree')

    # concat the fastas to make a single alignment
    fas_concat = os.path.join(SCRIPTS_FOLDER, "FASconCAT_v1.0.pl")
    process = subprocess.Popen('perl {} -s'.format(fas_concat), shell=True)
    process.wait()
    os.rename('FcC_smatrix.fas', 'concat_genes.fas')

    if merge:
        # generates the new .phy
        tree_node_merger = os.path.join(SCRIPTS_FOLDER, "TreeNodeMerger.py")
        process = subprocess.Popen('python {} concat_genes.fas -p'.format(tree_node_merger),
            shell=True)
    else:
        # convert the plain faa.fas file to phylip format
        fasta_2_phylip = os.path.join(SCRIPTS_FOLDER, "Fasta2Phylip.pl")
        print t.red(fasta_2_phylip)
        process = subprocess.Popen('perl {} concat_genes.fas concat_genes.phy'.format(fasta_2_phylip),
            shell=True)
    process.wait()


def build_tree(algorithm, N, builder, multithread, merge, folder_name, datatype):
    if folder_name:  # core tree builder
        os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome/Tree'))
    else:  # regular tree builder
        os.chdir(TREE_FOLDER)

    if builder == 'raxml':
        if multithread:
            # check if the multithreaded raxml is available
            if which('raxmlHPC-PTHREADS-SSE3'):
                rax = 'raxmlHPC-PTHREADS-SSE3'
            else:
                rax = 'raxmlHPC'
        else:
            rax = 'raxmlHPC'
        if merge:
            process = subprocess.Popen('{} -m {} -s concat_genes.fas.trimmed.phy -n Tree -p 12345 -x 12345 -N {} -f a'.format(
                rax, algorithm, N), shell=True)
        else:
            process = subprocess.Popen('{} -m {} -s concat_genes.phy -n Tree -p 12345 -x 12345 -N {} -f a'.format(
                rax, algorithm, N), shell=True)

        process.wait()
        os.rename('RAxML_bipartitions.Tree', 'out.Tree')

    elif builder == 'fasttree':
        # supply additional arguments if building a NT tree
        add_args = ''
        if datatype == 'nt':
            add_args = '-nt -gtr'
        if merge:
            process = subprocess.Popen('fasttree {} < concat_genes.fas.trimmed.fasta > out.Tree'.format(add_args), shell=True)
        else:
            process = subprocess.Popen('fasttree {} < concat_genes.fas > out.Tree'.format(add_args), shell=True)
        process.wait()

    # delete leftover fas files
    for file in os.listdir('.'):
        if file.endswith('.fas'):
            # don't delete usearch_CORE_ffn/faa.fas
            if not file.startswith('concat_genes'):
                os.remove(file)


def check_for_tree(folder_name=None):
    tree_props = {}
    if 'out.Tree' in os.listdir('.'):
        tree_props = retrieve_tree_props(folder_name)
    return tree_props
    

# keep track of genomes that were found in the first fasta file,
# and compare the contents of every other fasta to this list
def check_tree_fastas(files, datatype):
    first = []
    for i, file in enumerate(files):
        genomes_found = []
        with open(file) as f:
            for line in f:
                line = line.strip()
                # check that a genome is not repeated within the same file
                if line.startswith('>'):
                    genome = line
                    if i == 0:
                        if genome in first:
                            raise ValueError("The same genome was found twice")
                        first.append(genome)
                    genomes_found.append(genome)
                else:
                    # check that nt sequences have valid bases at each position
                    if datatype == 'nt':
                        for char in line:
                            if char not in 'ACGT-':
                                raise ValueError("Base is not one of ACGT")
        if i == 0:
            first.sort()
        genomes_found.sort()
        # make sure the genome contents of every other fasta is identical as the first
        if genomes_found != first:
            raise ValueError("Two files have different sets of genomes")


# Views the most recently generated tree in dendroscope
# specifically, opens up RAxML_bipartitions.FAA.Tree (or RAxML_bipartitions.FFN.Tree, whichever is available)
def dendroscope_view(folder_name=None):
    if folder_name:
        goto_folder(folder_name)
        os.chdir('Core_Genome/Tree')
    else:
        os.chdir(TREE_FOLDER)

    all_files = os.listdir('.')
    for file in all_files:
        if 'out.Tree' in file:
            process = subprocess.Popen(('dendroscope -x "open file=out.Tree; '
                'set drawer=RectangularPhylogram;expand direction=horizontal;'
                'expand direction=vertical;"'), shell=True)
            process.wait()


# prepares the tree.zip for download
def prepare_tree_zip(merge, folder_name=None):
    if folder_name:
        os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome/Tree'))
    else:
        os.chdir(TREE_FOLDER)

    if merge:
        process = subprocess.Popen(('zip tree.zip out.Tree '
            'concat_genes.fas.trimmed.phy concat_genes.fas.trimmed.fasta '
            'concat_genes.fas.legend.txt'), shell=True)
    else:
        process = subprocess.Popen(('zip tree.zip out.Tree '
            'concat_genes.fas concat_genes.phy'), shell=True)
    process.wait()

    # move zip to data folder
    shutil.move('tree.zip', os.path.join(OUTPUT_FOLDER, 'tree.zip'))


# This function retrieves/generates information related to a tree, e.g.
def generate_tree_props(N, datatype, gene_list, builder, merge, algorithm, folder_name=None):
    tree_props = {}
    if builder == 'raxml':
        tree_props['N'] = N
        tree_props['builder'] = 'raxml'
        tree_props['algorithm'] = algorithm
    elif builder == 'fasttree':
        tree_props['builder'] = 'fasttree'
        
    tree_props['datatype'] = datatype
    tree_props['gene_list'] = gene_list
    tree_props['merge'] = merge

    if folder_name:
        os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome/Tree'))
    else:
        os.chdir(TREE_FOLDER)

    with open('tree.json', 'w') as out:
        json.dump(tree_props, out, indent=4)


def retrieve_tree_props(folder_name=None):
    if folder_name:
        os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome/Tree'))
    else:
        os.chdir(TREE_FOLDER)
    with open('tree.json') as f:
        tree_props = json.load(f)
    return tree_props


# Checks if the tree parameters are good; this depends on two things:
# 1) Is this a core-gene tree builder, or the general tree builder?
# 2) Are you using raxml or fasttree?
def check_tree_parameters(datatype, algorithm, N, builder):
    # if using raxml, must check other arguments
    if builder == 'raxml':
        if N < 1 and N != None:  # N being None is okay
            print t.cyan(str(N))
            print t.red("N must be positive")
            raise ValueError("N must be positive")
        if datatype == 'nt':
            if algorithm not in ACCEPTED_NT_ALGORITHMS:
                print t.red("Invalid nucleotide algorithm")
                raise ValueError("Invalid nucleotide algorithm")
        elif datatype == 'aa':
            if algorithm not in ACCEPTED_AA_ALGORITHMS:
                print t.red("Invalid amino acid algorithm")
                raise ValueError("Invalid amino acid algorithm")

    print t.green("Tree parameters validated")


# Uses the python ete3 package to produce a visual tree, then
# serves it to the user for downloading
def prepare_ete_tree(show_align, show_bootstrap, core_flag, out_type, folder_name=None):
    if folder_name:
        os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome/Tree'))
    else:
        os.chdir(TREE_FOLDER)

    # checks if the fasta file was trimmed (ie. clade collapsing occurred)
    # by checking if a trimmed.fasta file exists
    is_trimmed = ''
    legend_file = ''
    for file in os.listdir('.'):
        if file.endswith('trimmed.fasta'):
            is_trimmed = '.trimmed.fasta'
        if file.endswith('.legend.txt'):
            legend_file = file

    b_param, i_param, L_param = '', '', ''
    
    # -b is used by TreeVisualizer -> enables bootstrap value display
    if show_bootstrap:
        b_param = '-b'
    
    # -i specifies the SNP alignment file
    if show_align:
        i_param = '-i concat_genes.fas{}'.format(is_trimmed)
    
    # -L enables displaying of legend; only valid if clade collapsing occurred
    # when the tree was generated
    if legend_file:
        L_param = '-L {}'.format(legend_file)

    # out_type: either 'view' or 'pdf':
    if out_type == 'pdf':
        out_param = "-out pdf"
    elif out_type == 'view':
        out_param = "-out view"

    tree_visualizer_py = os.path.join(SCRIPTS_FOLDER, "TreeVisualizer.py")
    process = subprocess.Popen(('python {} -t out.Tree '
        '{} {} {} {}').format(tree_visualizer_py, i_param, b_param, L_param, out_param), shell=True)
    process.wait()

    if out_type == 'pdf':
        shutil.move('tree.pdf', os.path.join(OUTPUT_FOLDER, 'tree.pdf'))
# =====================================================================================
# Gene Lookup - accepts a multi-fasta of either Amino Acid or Nucleotide sequences
# (not both at the same time), and searches against the database of all genomes in the run
@app.route('/gene_lookup/<folder_name>', methods=['GET', 'POST'])
@app.route('/gene_lookup', methods=['GET', 'POST'])
def gene_lookup(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)

    lookup_form = None
    blast_results = None
    core_families = None
    results_form = None
    run_name = ''
    select_list = get_select_list()

    if folder_name:
        run_name = get_run_name(folder_name)
        lookup_form = GeneLookupForm()

        # check if the database exists
        if 'blastdb' not in os.listdir(os.path.join(DATA_FOLDER, folder_name)):
            initialize_db(folder_name)
            flash("Usearch databases initialized")

        # if run was performed, render usearch results table
        if 'out.tsv' in os.listdir(os.path.join(DATA_FOLDER, folder_name, 'blastdb')):
            goto_folder(folder_name)
            os.chdir('blastdb')
            with open('search_results.json') as f:
                search_results = json.load(f)
            with open('search_params.json') as f:
                search_params = json.load(f)

            blast_results = search_results['blast_results']
            # core target genes are assigned numbers based on core family
            core_families = search_results['core_families']
            # list of target genes that aligned to each query gene
            gene_alignments = search_results['gene_alignments']

            # allows you to select & download alignment for a query gene
            results_form = LookupResultsForm()
            results_form.select.choices = [(k, k) for k in gene_alignments.keys()]

    if request.method == 'POST':
        # Run the analysis itself
        if lookup_form.submit.data:
            if lookup_form.validate_on_submit():
                datatype = request.form['datatype']  # aa/nt
                try:
                    upload_query_gene(folder_name, lookup_form)
                    ident = request.form['ident']
                    plen = request.form['plen']
                    perform_gene_lookup(folder_name, datatype, ident, plen)
                    create_search_results(folder_name)
                    create_search_params(folder_name, datatype, ident, plen)
                    flash("Gene search complete")
                except ValueError as e:
                    flash(str(e))
                return redirect(url_for('gene_lookup', folder_name=folder_name))  
            else:
                if not lookup_form.ident.data:
                    flash("Please type in a % Identity value between 0 and 1")
                elif lookup_form.ident.data > 1 or lookup_form.ident.data < 0:
                    flash("Please type in a % Identity value between 0 and 1")
                if not lookup_form.plen.data:
                    flash("Please type in a % Length value between 0 and 1")
                elif lookup_form.plen.data > 1 or lookup_form.plen.data < 0:
                    flash("Please type in a % Length value between 0 and 1")
                if not lookup_form.upload.data:
                    flash("Must specify a multi-fasta file")

        # download an alignment for a specified query gene
        elif results_form.download.data:
            selected_gene = results_form.select.data
            type_chosen = search_params['datatype']
            lookup_alignment(folder_name, gene_alignments, selected_gene, type_chosen)
            return redirect(url_for('download', filename="output/post_align.fasta"))

        elif results_form.genomes_included.data:
            prepare_genomes_included(folder_name, gene_alignments)
            return redirect(url_for('download', filename="output/genomes_included.tsv"))

        elif results_form.download_table.data:
            prepare_lookup_xlsx(folder_name, search_results)
            return redirect(url_for('download', filename="output/lookup_table.xlsx"))

    return render_template('gene_lookup.html', lookup_form=lookup_form, select_list=select_list,
        results_form=results_form, blast_results=blast_results, core_families=core_families,
        run_name=run_name)


def worker_db(file):
    if file.endswith('faa'):
        udb_filename = file[:-4] + '_{}.udb'.format('faa')
    elif file.endswith('ffn'):
        udb_filename = file[:-4] + '_{}.udb'.format('ffn')
    process = subprocess.Popen(('usearch -makeudb_usearch {0} '
        '-output {1} -quiet').format(file, udb_filename), shell=True)
    process.wait()


def initialize_db(folder_name):
    goto_folder(folder_name)
    os.mkdir('blastdb')

    # create a concatenated faa file
    all_faa = os.listdir('PROKKA_backups')
    for faa in all_faa:
        # we want to get rid of '.orig' in the filename
        shutil.copy('PROKKA_backups/{0}'.format(faa),
            'blastdb/{0}'.format(faa.replace('.orig', '')))
    # get rid of .orig in the all_faa list
    for i, faa in enumerate(all_faa):
        all_faa[i] = all_faa[i].replace('.orig', '')

    #create a concatenated ffn file
    all_ffn = os.listdir('ffn')
    for ffn in all_ffn:
        shutil.copy('ffn/{0}'.format(ffn), 'blastdb/{0}'.format(ffn))

    # create aa/nt databases
    os.chdir('blastdb')
    # usearch makedb function; uses multiple child processes
    all_files = all_faa + all_ffn
    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    p.map(worker_db, all_files)


def worker_usearch(to_process):
    file, extension, ident, plen = to_process

    udb_filename = file[:-4] + '_{}.udb'.format(extension)
    tsv_filename = file[:-4] + '_{}.tsv'.format(extension)

    process = subprocess.Popen(('usearch -usearch_global query.fasta -db {0} -id {1} '
        '-target_cov {2} -userout {3} -maxaccepts 0 -strand both '
        '-userfields query+target+id+tcov+ql+mism+gaps').format(
        udb_filename, ident, plen, tsv_filename), shell=True)
    process.wait()


def upload_query_gene(folder_name, form):
    # we need to delete the last entry!
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'blastdb'))
    to_delete = ["pre_align.fasta", "genomes_included.tsv", "search_results.json",
        "search_params.json", "out.tsv", "query.fasta"]
    for old_file in to_delete:
        if old_file in os.listdir('.'):
            os.remove(old_file)

    file = form.upload.data
    filename = secure_filename(file.filename)
    # check that the filename is a fasta
    if not filename.endswith(".fasta"):
        raise IOError("Only fasta files are allowed")
    file.save(filename)

    datatype = form.datatype.data
    # check the file itself (mostly we want to see that
    # if we chose nt, then we are getting nt sequences, and if we
    # chose aa, then we are getting aa sequences)
    with open(filename) as f:
        for line in f:
            if not line.startswith(">"):
                line = line.strip()
                for char in line:
                    if datatype == 'nt':
                        if char not in "ACGT":
                            os.remove(filename)
                            raise ValueError("Found invalid nucleotide: %s" % char)

    # if the file looks ok, rename it to query.fasta
    os.rename(filename, 'query.fasta')


# Run usearch on the query genes against the entire database of the specified run
# to find target genes that meet the ident/plen threshold
def perform_gene_lookup(folder_name, datatype, ident, plen):
    goto_folder(folder_name)
    os.chdir('blastdb')
    if datatype == 'nt':
        extension = 'ffn'
        all_files = glob.glob('*.ffn')
    elif datatype == 'aa':
        extension = 'faa'
        all_files = glob.glob('*.faa')

    # create an array of tuples for each file
    # due to how Pool works, it's easier to pass in parameters like ident
    # every time worker_usearch is called
    process_list = [(file, extension, ident, plen) for file in all_files]
    # for now, defaults to max cores
    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    p.map(worker_usearch, process_list)

    # get rid of previous gene-lookup results
    to_remove = [
        'out.tsv',
        'genomes_included.tsv',
        'search_params.json',
        'search_results.json'
    ]
    for elem in to_remove:
        if elem in os.listdir('.'):
            os.remove(elem)

    # put .tsv outputs into a single combined out.tsv
    all_tsv = glob.glob('*.tsv')

    combined_tsv = ''
    for tsv in all_tsv:
        with open(tsv) as f:
            contents = f.read()
        combined_tsv += contents
        os.remove(tsv)

    with open('out.tsv'.format(extension), 'w') as out:
        out.write(combined_tsv)


# writes to json the input parameters for the most recent lookup
def create_search_params(folder_name, datatype, ident, plen):
    summary = {
        'ident': ident,
        'plen': plen,
        'datatype': datatype
    }
    goto_folder(folder_name)
    os.chdir('blastdb')
    with open('search_params.json', 'w') as out:
        json.dump(summary, out, indent=4)


# parses the out.tsv from usearch; creates search_results.json, which
# stores information on whether a query gene hit a core gene,
# which genes hit with which core gene, etc.
def create_search_results(folder_name):
    with open(os.path.join(DATA_FOLDER, folder_name, 'blastdb/out.tsv')) as f:
        # each sublist is a usearch hit for a single target gene
        blast_results = [line.strip().split('\t') for line in f]

    # get rid of prokka annotations
    for i, result in enumerate(blast_results):
        result[1] = result[1].split(' ')[0]

    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'Core_Genome'))
    # each sublist is an entire core gene family
    core_genome = []
    with open('CoreGenome.usearch') as f:
        for line in f:
            line = line.strip()
            arr = line.split('\t')
            for i, elem in enumerate(arr):
                # get rid of square brackets and the '>'
                arr[i] = elem.split('[')[0][1:]
            core_genome.append(arr)

    # look for each target gene in core genome
    # if found, record the core gene family, otherwise keep as None
    core_families = {}

    for result in blast_results:
        blast_gene = result[1]
        core_families[blast_gene] = None 
        gene_family = []  # basically to copy over the core gene family
        index = -1  # placeholder; used to store the actual core gene family number

        # go through the core genome for every query gene;
        # perhaps can implement a better algorithm in the future
        for i, core_gene in enumerate(core_genome):
            if blast_gene in core_gene:
                gene_family = core_gene
                index = i+1  # since i starts at 0
                core_families[blast_gene] = index
                break

    # keys are query genes; values are all target genes associated with query gene
    # *regardless* of whether the target gene is core or not
    # useful for producing an alignment with query and all its targets
    gene_alignments = {}

    for result in blast_results:
        if result[0] not in gene_alignments.keys():
            gene_alignments[result[0]] = [result[1]]  # new list
        else:
            gene_alignments[result[0]].append(result[1])  # extend list

    # put all results into a single dictionary/object
    search_results = {
        'blast_results': blast_results,
        'core_families': core_families,
        'gene_alignments': gene_alignments
    }

    goto_folder(folder_name)
    os.chdir('blastdb')
    # write results as a json
    with open('search_results.json', 'w') as out:
        json.dump(search_results, out, indent=4)


# Looks up which target genes are associated with our query gene,
# so we can produce and view an alignment
def lookup_alignment(folder_name, gene_alignments, selected_gene, datatype):
    # get all the target genes for specified gene
    all_family_genes = gene_alignments[selected_gene]  
    # get all genome names by splitting to the left of underline
    all_genomes = sorted([ gene.split('_')[0] for gene in all_family_genes ])
    
    # convert datatype to file extension
    if datatype == 'nt':
        extension = 'ffn'
    elif datatype == 'aa':
        extension = 'faa'

    # look inside the blastdb database, and write only the target genes
    # included in subset_files
    subset_files = [ file+'.'+extension for file in all_genomes ]
    os.chdir(os.path.join(DATA_FOLDER, folder_name, 'blastdb'))
    
    # db_subset is the contents of the output alignment file
    db_subset = ''
    for file in subset_files:
        with open(file) as f:
            # if write is true, we "write" to the db_subset string
            write = False
            for line in f:
                if line.startswith('>'):
                    write = False
                    # ignore prokka annotation to the right
                    gene_only = line.split(' ')[0][1:]
                    # check if gene is found
                    if gene_only in all_family_genes:
                        write = True
                        db_subset += line
                else:
                    if write:
                        db_subset += line

    # do the same for the query genes; also, we are only searching for
    # one query gene, so break as soon as it is found
    query_subset = ''
    with open('query.fasta') as f:
        write = False
        for line in f:
            if line.startswith('>'):
                if write:
                    write = False
                    break
                # see if any of the family genes is contained
                gene_only = line.strip()[1:]
                if gene_only == selected_gene:
                    write = True
                    query_subset += line
            else:
                if write:
                    query_subset += line

    # write contents of query_subset and db_subset to pre_align.fasta
    # which must be aligned later by MUSCLE
    # the query gene is always the 'reference' (top of the alignment)
    with open('pre_align.fasta', 'w') as out:
        out.write(query_subset)
        out.write(db_subset)

    # run muscle
    process = subprocess.Popen("muscle -in pre_align.fasta -out tmp.fasta", shell=True)
    process.wait()
    # stable.py restores the original order of the genomes
    stable_py = os.path.join(SCRIPTS_FOLDER, "stable.py")
    process = subprocess.Popen(("python {} pre_align.fasta tmp.fasta > "
        "post_align.fasta".format(stable_py)), shell=True)
    process.wait()
    os.remove("tmp.fasta")

    # move post_align.fasta to output folder
    shutil.move("post_align.fasta", os.path.join(OUTPUT_FOLDER, "post_align.fasta"))


# create genomes_included.tsv, which indicates which genomes were
# present/absent for each query gene
def prepare_genomes_included(folder_name, gene_alignments):
    # get a complete list of genomes
    goto_folder(folder_name)
    strains = []
    with open('STRAINS') as f:
        for line in f:
            line = line.strip()
            strains.append(line)
    num_strains = len(strains)

    # keys: query genes, values: a dict with list of included genomes
    # and a list of excluded genomes
    incl_excl = {}  
    for q_gene, genomes in gene_alignments.items():
        # get rid of everything to the right of _
        included = sorted([genome.split('_')[0] for genome in genomes])
        incl_excl[q_gene] = {
            'included': included,
            'excluded': []
        }
        # if a genome from STRAINS is not found in included,
        # it is considered excluded
        i, j = 0,0
        num_included = len(included)
        while i < num_strains and j < num_included:
            if strains[i] == included[j]:
                i += 1
                j += 1
            else:
                incl_excl[q_gene]['excluded'].append(strains[i])
                i += 1
        while i < num_strains:
            incl_excl[q_gene]['excluded'].append(strains[i])
            i += 1

    # reorganize data in incl_excl so it is easy to write to .tsv
    table = []
    for q_gene, incl_dict in incl_excl.items():
        table.append(incl_dict['included'])
        table.append(incl_dict['excluded'])

    # table has sublists that are of different length, so
    # take the max length and iterate over that range
    max_len = 0
    for row in table:
        if len(row) > max_len:
            max_len = len(row)

    goto_folder(folder_name)
    os.chdir('blastdb')

    with open('genomes_included.tsv', 'w') as out:
        # write all query genes
        for q_gene in incl_excl.keys():
            out.write(q_gene + '\t\t\t')
        out.write('\n')
        # write headers
        for i in range(len(gene_alignments)):
            out.write('-Genomes Included-\t-Genomes Excluded-\t\t')
        out.write('\n')
        # write actual table contents
        for i in range(max_len):
            for j, row in enumerate(table):
                if i < len(row):
                    out.write(row[i])
                if j%2 == 1:
                    out.write('\t')
                out.write('\t')
            out.write('\n')

    shutil.move("genomes_included.tsv", os.path.join(OUTPUT_FOLDER, "genomes_included.tsv"))


def prepare_lookup_xlsx(folder_name, search_results):
    core_families = search_results['core_families']
    blast_results = search_results['blast_results']
    df = pd.DataFrame(blast_results)
    df2 = pd.DataFrame.from_dict(core_families.items())
    df = df.rename(index=str, columns={1:"index"})
    df2 = df2.rename(index=str, columns={0:"index"})

    df_combined = df.join(df2.set_index('index'), on='index')
    df_combined = df_combined.rename(index=str, columns={
        0: "Query Gene",
        'index': "Target Gene",
        2: "% Identity",
        3: "% Coverage",
        4: "Align Length",
        5: "Mismatches",
        6: "Open Gaps",
        1: "Core Gene"
    })

    cols = ['Core Gene', 'Query Gene', 'Target Gene', '% Identity', '% Coverage', 'Align Length', 'Mismatches', 'Open Gaps']
    df_combined.ix[:,cols]
    df_combined.sort_values(['Query Gene', 'Target Gene'])

    writer = pd.ExcelWriter('lookup_table.xlsx')
    df_combined.to_excel(writer, 'Sheet1', index=False)
    writer.save()
    shutil.move('lookup_table.xlsx', os.path.join(OUTPUT_FOLDER, 'lookup_table.xlsx'))
# =====================================================================================
# Migrate - create an importable zip of runs
# useful for transferring data from server to local copy of HardCORE
@app.route('/migrate', methods=['GET', 'POST'])
def migrate():
    form = MigrationForm()
    with open(os.path.join(DATA_FOLDER, 'runs.json')) as f:
        runs = json.load(f)

    if request.method == 'POST':
        if form.export_runs.data:
            # get the list of checkboxes
            to_process = request.form.getlist("check")
            if len(to_process) == 0:
                flash("You must select at least one run to export!")
            else:
                export_runs_to_zip(to_process)
                return redirect(url_for('download', filename="output/export.zip"))
        elif form.import_runs.data:
            import_runs_from_zip(form.upload.data)
            flash("Import successful")
            return redirect(url_for('migrate'))

    return render_template('migrate.html', form=form, runs=runs)


def export_runs_to_zip(runs):
    os.chdir(DATA_FOLDER)
    for folder_name in runs:
        process = subprocess.Popen("zip -r export.zip {0}".format(folder_name), shell=True)
        process.wait()
    to_upload = ['current_references.json', 'runs.json', 'runs.txt']
    for file in to_upload:
        process = subprocess.Popen("zip export.zip {0}".format(file), shell=True)
        process.wait()
    shutil.move("export.zip", os.path.join(OUTPUT_FOLDER, "export.zip"))


def import_runs_from_zip(file):
    # get/download the uploaded file
    filename = secure_filename(file.filename)
    os.chdir(DATA_FOLDER)
    os.mkdir("TMP")
    file.save(os.path.join(DATA_FOLDER, "TMP", filename))
    os.chdir("TMP")

    process = subprocess.Popen("unzip {0}".format(filename), shell=True)
    process.wait()

    folders = next(os.walk('.'))[1]

    with open('runs.json') as f:
        runs_json = json.load(f)
    runs_json = {k:v for k, v in runs_json.items() if k in folders}
    with open('current_references.json') as f:
        current_references = json.load(f)
    current_references = {k:v for k, v in current_references.items() if k in folders}

    os.chdir(DATA_FOLDER)
    for folder in folders:
        # delete old folder
        if folder in os.listdir('.'):
            shutil.rmtree(folder)
        # migrate the new folder over
        shutil.move("TMP/{0}".format(folder), './{0}'.format(folder))

    # update runs.json and current_references.json
    with open('runs.json') as f:
        runs_json_old = json.load(f)
    with open('current_references.json') as f:
        current_references_old = json.load(f)

    for k, v in runs_json.items():
        runs_json_old[k] = v
    for k,v in current_references.items():
        current_references_old[k] = v

    # write back to .json files
    with open('runs.json', 'w') as out:
        json.dump(runs_json_old, out, indent=4)
    with open('current_references.json', 'w') as out:
        json.dump(current_references_old, out, indent=4)

    shutil.rmtree("TMP")
    update_runs_txt()
# =====================================================================================
# Network - allows the user to view relationships between genomes based on
# how many gene families they share
@app.route('/network/<folder_name>', methods=['GET', 'POST'])
@app.route('/network', methods=['GET', 'POST'])
def network(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)

    select_list = get_select_list()

    if not os.path.isdir(VISUALS_FOLDER):
        os.mkdir(VISUALS_FOLDER)

    if folder_name:
        form = NetworkForm()
        # check for a post request
        if request.method == 'POST':
            if form.download.data:
                # serve STRAINS as a json file
                # return redirect(url_for('download', filename="data/{}/STRAINS".format(folder_name)))
                out_json = collections.OrderedDict()

                # use the strains file to get all genome names
                strains_dir = os.path.join(DATA_FOLDER, folder_name, "STRAINS")
                with open(strains_dir) as f:
                    k = 1
                    for line in f:
                        line = line.strip()
                        out_json[line] = k
                        k += 1

                labels_template_filename = "{}_labels_template.json".format(folder_name)
                with open(labels_template_filename, 'w') as out:
                    json.dump(out_json, out, indent=4)
                shutil.move(labels_template_filename, os.path.join(OUTPUT_FOLDER, labels_template_filename))
                return redirect(url_for('download', filename="output/{}".format(labels_template_filename)))

            elif form.submit.data:
                labels_file = form.upload.data
                if labels_file:
                    labels_filename = secure_filename(labels_file.filename)
                    os.chdir(VISUALS_FOLDER)
                    labels_file.save(labels_filename)

                    try:
                        categories = extract_categories(labels_filename)
                        os.rename(labels_filename, "{}_labels.json".format(folder_name))
                    except ValueError as e:
                        os.remove(labels_filename)

                return redirect(url_for('network', folder_name=folder_name))

            elif form.restore_default.data:
                labels_json = "{}_labels.json".format(folder_name)
                if labels_json.format(folder_name) in os.listdir(VISUALS_FOLDER):
                    os.remove(os.path.join(VISUALS_FOLDER, labels_json))

                return redirect(url_for('network', folder_name=folder_name))

        # create the pan_genome data json
        if "{}.json".format(folder_name) not in os.listdir(VISUALS_FOLDER):
            goto_folder(folder_name)
            os.chdir('Pan_Genome')
            script_call = os.path.join(SCRIPTS_FOLDER, 'create_pan_visualize_json.py {0}'.format(folder_name))
            process = subprocess.Popen("python {0}".format(script_call), shell=True)
            process.wait()
            shutil.move("data_pre.json", os.path.join(VISUALS_FOLDER, '{0}.json'.format(folder_name)))
        return render_template('network.html', folder_name=folder_name, form=form)
    else:
        return render_template('select_run.html', select_list=select_list, visual_type="network")


def extract_categories(filename):
    with open(filename) as f:
        categories = json.load(f)
    # do some error checking within categories, I guess?
    return categories

# =====================================================================================
# Sunburst - allows the user to view counts for gene families with a particular subset of genomes
@app.route('/sunburst/<folder_name>', methods=['GET', 'POST'])
@app.route('/sunburst', methods=['GET', 'POST'])
def sunburst(folder_name=None):
    # return 404 page if invalid folder_name
    if folder_name:
        try:
            validate_folder_name(folder_name)
        except ValueError:
            abort(404)
    
    if not os.path.isdir(VISUALS_FOLDER):
        os.mkdir(VISUALS_FOLDER)

    select_list = get_select_list()
    if folder_name:
        # check if the .json has not been made for this run yet
        if "{}_flare.json".format(folder_name) not in os.listdir(VISUALS_FOLDER):
            goto_folder(folder_name)
            script_call = os.path.join(SCRIPTS_FOLDER, 'create_pan_components_json.py {0}'.format(folder_name))
            process = subprocess.Popen("python {0}".format(script_call), shell=True)
            process.wait()

            shutil.move("Pan_Genome/flare.json", os.path.join(VISUALS_FOLDER, '{0}_flare.json'.format(folder_name)))
        return render_template('sunburst.html', folder_name=folder_name)
    else:
        return render_template('select_run.html', select_list=select_list, visual_type="sunburst")
# =====================================================================================
# About page
@app.route('/about')
def about():
    return render_template('about.html')
# =====================================================================================
# Download route
@app.route('/download/<path:filename>')
def download(filename):
    return send_from_directory(APP_FOLDER, filename, as_attachment=True)
# =====================================================================================
# 404 error page
@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OTHER METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# go into the folder specified by folder_name (a run)
def goto_folder(folder_name):
    if folder_name:
        os.chdir(os.path.join(DATA_FOLDER, folder_name))


#Input: a date string like 20170613_123529
#Output: a date string like Jun 13 2017, 12:35:29
def prettify_date(folder_name):
    dt_obj = datetime.datetime.strptime(folder_name, '%Y%m%d_%H%M%S')
    run_date = dt_obj.strftime("%b %d %Y, %H:%M:%S")
    return run_date


# # checks PanGenome.usearch and CoreGenome.usearch to ensure that the
# # run has been performed correctly
# # looks for two possible sources of errors:
# # 1) a genome should not appear twice or more in the same line
# # 2) the same gene should not appear twice
def check_pan_core_usearch(folder_name):
    to_process = [
        'Core_Genome/CoreGenome.usearch',
        'Pan_Genome/PanGenome.usearch'
    ]
    errors = {
        "Gene found twice in CoreGenome": False,
        "Genome appeared twice in the same line in PanGenome": False,
        "Genome appeared twice in the same line in CoreGenome": False
    }

    print t.cyan("Commencing Usearch Check for {0}".format(folder_name))
    goto_folder(folder_name)
    for usearch in to_process:
        found_genes = {}
        usearch_name = usearch.split('/')[-1]  # either CoreGenome.usearch or PanGenome.usearch
        print t.cyan("Checking %s..." % usearch_name)
        with open(usearch) as f:
            for line in f:
                line = line.strip().split('\t')
                genomes = []
                for gene in line:
                    genome, gene_number = find_genome_regexp(gene)
                    # keep track of what genomes were found
                    genomes.append(genome)
                    # check if genome/gene_number were previously found
                    if genome not in found_genes.keys():
                        found_genes[genome] = []
                    else:
                        if gene_number not in found_genes[genome]:
                            found_genes[genome].append(gene_number)
                        else:
                            # only for core genome we don't want to see the same gene twice anywhere
                            if usearch_name == 'CoreGenome.usearch':
                                errors["Gene found twice in CoreGenome"] = True

                # check if the same genome appears twice (or more)
                unique_genomes = set(genomes)
                if len(genomes) != len(unique_genomes):
                    if usearch_name == 'CoreGenome.usearch':
                        errors["Genome appeared twice in the same line in CoreGenome"] = True
                    else:
                        errors["Genome appeared twice in the same line in PanGenome"] = True

        print t.cyan("{0} check complete!".format(usearch_name))

    print t.cyan("Usearch Check Complete for {0}".format(folder_name))
    return errors


def find_genome_regexp(line):
    p = re.compile(r'>(.*)_(\d+)\[.*\]')
    m = p.match(line)
    genome = m.groups()[0]
    gene_number = m.groups()[1]
    return genome, gene_number


def get_num_strains(folder_name):
    goto_folder(folder_name)
    num = 0
    with open('STRAINS') as file:
        for line in file:
            num += 1
    return num


# determines if a program is inside PATH
# useful for checking if the user has the multithreaded raxml on their machine
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


# Gets the list of runs, where each entry is a 2-tuple of 
# the folder name, and its run name/prettified date
def get_select_list():
    choices = []
    with open(os.path.join(DATA_FOLDER, 'runs.json')) as f:
        summary = json.load(f)

    with open(os.path.join(DATA_FOLDER, 'runs.txt')) as runs:
        for line in runs:
            line = line.strip()
            # Converts YYYYMMDD_HHMMSS -> MMM DD YYYY, HH:MM:SS
            # Example: 20170710_181233 -> Jul 10 2017, 18:12:33
            dt_obj = datetime.datetime.strptime(line, '%Y%m%d_%H%M%S')
            readable_date = dt_obj.strftime("%b %d %Y, %H:%M:%S")

            # date_rep: what is shown to the user in the dropdown
            x = summary[line]['run_name']
            line_rep = x + " (" + readable_date + ")"
            choices.append((line, line_rep))

    return choices

# gets the run name that the user specified when they passed their files
# into the hardcore pipeline
def get_run_name(folder_name):
    with open(os.path.join(DATA_FOLDER, 'runs.json')) as f:
        summary = json.load(f)
    dt_obj = datetime.datetime.strptime(folder_name, '%Y%m%d_%H%M%S')
    readable_date = dt_obj.strftime("%b %d %Y, %H:%M:%S")
    return summary[folder_name]['run_name'] + " (" + readable_date + ")"


# throws an exception if the folder_name is invalid
def validate_folder_name(folder_name):
    os.chdir(DATA_FOLDER)
    folders = next(os.walk('.'))[1]
    if folder_name not in folders:
        raise ValueError("folder_name invalid")

