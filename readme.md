# readme.md

HardCORE is a combination of a set of Python and Perl scripts (the original HardCORE suite, written by Nicholas Petronella), and a full-stack web application that I (Brian Lee) have developed. It is designed mostly towards Canadian researchers in the federal government in response to emerging need to quickly and accurately compute the core and pan genome. This application will be useful for comparing two groups (such as genera or species), or determining the

Features include:
* Find the core genome of a set of strains that you are interested in, and view its constituent core genes
* Compare core genes across the different genomes
    * Sortable SNP bar charts that may possibly reveal patterns
* Perform subset analysis (which gene families are present in all members of a certain group, but in none of the other strains)
    * e.g. What genes do all of my E.coli have that Salmonella don't?
* Draw a rarefaction curve to obtain a representative core genome size
    * Important in publications to 'prove' that you have selected an appropriate number of strains for analysis
* Look for genes in your collection of genomes
* Create phylogenetic trees from either the core genome or aligned multi-fasta files
* Use interactive visualization tools

## Installation

There is no need to install anything if you are planning to use HardCORE on the BFSSI Bioinformatics Lab server. However, you will need Prokka-annotated genomes to perform the initial run. Please refer to the *Usage (for scientists)* section for instructions on how to get Prokka installed on your machine.
If you are planning on deployment onto the server, however, there are a number of bioinformatics tools and libraries that must be installed on your machine first.

* [Muscle](http://www.drive5.com/muscle/downloads.htm)
* [Usearch](http://www.drive5.com/usearch/download.html)
* [FastTree](http://www.microbesonline.org/fasttree)
* [ETE Toolkit](http://etetoolkit.org/download/)
   * Note: Please install ETE3 using the "Native" installation method
* [RAxML](https://github.com/stamatak/standard-RAxML)
* snp-sites (`$ sudo apt-get install snp-sites`)

### Special notes on installation
#### Usearch
After downloading the 32-bit executable (Linux version), rename it to just "usearch" and move it into /usr/bin using the following command:
```bash
sudo mv usearch /usr/bin
```

#### RAxML
Clone this GitHub repo and run the following lines
```bash
"make -f Makefile.gcc"
"rm *.o"
```
A multithreaded version can also be compiled, which is ~40% faster. To install this version of RAxML, run the following:
```bash
"make -f Makefile.SSE3.PTHREADS.gcc"
"rm *.o"
```

**Important**: make sure the standard-RAxML folder is in your PATH. To do so, append the following to the bottom of your .bashrc or .zshrc (depending on what shell you are using):
export PATH=/PATH/TO/standard-RAxML:$PATH

#### HardCORE Suite
The HardCORE Suite (credits to Nick, see below) is included in the application container. It can function on its own, however, as long as you run the HardCORE_usearch.pl script. Please refer to the source code in views.py for an example call using the perl script.

### Deployment
You may be interested in hosting a local copy of HardCORE, especially if you perform the main HardCORE run on the server, but want to move the results over to your local machine to perform further analysis. In this case, you will want to do the following steps.

1\. Open your Terminal and type in the following:
```
$ git clone https://github.com/hardcore/HardCORE.git
```
This will clone the repo into your current directory.

2\. There are a number of python requirements. To install them, go into the HardCORE folder using the terminal, and type:
```
pip install -r requirements.txt
```
This will use the requirements text file to install specific versions of Python dependencies, such as Flask, pandas, and matplotlib. Do note that using a virtual environment will likely not work (ETE doesn't like it in particular), so you will have to install the requirements globally. This may require superuser access (i.e. add the keyword 'sudo' before running the pip command).

3\. Go into the HardCORE directory using `cd HardCORE`. At this point, double-check that you have installed the dependencies listed in the "Installation" section. If you have everything needed to get HardCORE started, simply type
```
$ python run.py
```
and HardCORE will execute. You will want to see something similar to this:
```
 * Running on http://0.0.0.0:5000/ (Press CTRL+C to quit)
```
Note: Yes, HardCORE is powered by the default Flask development server. I have attempted to use Gunicorn as the WSGI server instead, but it was unfortunately causing problems related to multiple processes. For a project of this scale, which will not likely be deployed for use by the general public, I believe the development server should function well enough for your needs as a researcher.

4\. Type `ifconfig` on your terminal, and look the IP address of the host from where you will deploy HardCORE. Then, on a client machine that is on the same local network as the server machine, type into your browser (Chromium, Chrome, or Firefox should be fine):
```
http://(HOST IP ADDRESS):5000/
```
Congratulations! You are now running the HardCORE application using your own machine. Port 5000 is the default for serving the application, but feel free to change this in run.py in the HardCORE root directory.

Note: if you are using HardCORE on the same machine as the server, it will suffice to type:
```
http://0.0.0.0:5000/
```
or
```
http://localhost:5000/
```

## Usage

The HardCORE suite of Python/Perl scripts requires the user to provide fully annotated FFA and FFN gene sequences of their genomes, where the FAA files are in one zip, and the FFN files are in another. The annotation must be performed by Prokka, which you can download here: [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml) Note that Prokka has a number of dependencies, described under the "Requirements" header. Please ensure that they are installed first.

When you are on the index page, upload your FAA and FFN zips in their respective areas. You will notice a number of input boxes below the upload fields. These are where you will specify your percent identity (how similar in sequence should two genes be at the amino acid level to be considered part of the same gene family?) and percent length (over what percent of the sequence length should they express similarity?). The number of threads can usually be left blank to default to maximizing the number of CPU cores to be used, which will make the analysis run faster. You can optionally input a run name to organize your runs by the type of study (e.g. Brian-AMR-Summer17, Brian-AMR-Winter18, Lydia Heat Shock Vibrio); if none is provided, the run name is set to the current time.

<span style="color:red;">Please do not close your browser tab while running HardCORE!!</span>

Once the analysis is finished, you can begin to explore the other aspects of the HardCORE application. Whether you are interested in viewing the Core Genome (and its SNPs), building a phylogenetic tree from core genome sequences (including the whole core), or investigating genome relationships with the Network page, the world is your oyster! :fireworks:

## Additional Notes

There should be three text files in the /app/data/ folder. Don't delete these, even if they are empty. **runs.txt** is a newline-separated list of date-times, which represents the runs that reside in this app. **runs.json** is a json which stores properties about each run in more detail. If there are no runs, runs.jsons should still contain empty curly braces {} (as in, it cannot be completely empty, or that will cause problems). **current_references.json** is a list of currently selected reference genomes for each run in the Core Genome page. Again, even if there are no runs recorded in the json, the file should still contain {}.

The /app/output/ folder consists of files that have been served to clients. Feel free to use the *cleanup* button in the Summary Table page to safely remove these files.

The /app/static/ folder mostly contains css and js files. However, within the static folder, there is a graphs/ folder. It contains tsvs and jsons for rendering SNP tables at core genes (or the whole core), as well as the d3 bar charts. These can be safely deleted (and are indeed deleted using the cleanup button). There is also a visuals/ folder, which contains jsons for the Network and Sunburst pages. These can be deleted as well, but are currently not handled by the cleanup button. 

## Acknowledgements

Owen Neuber - for the source code of his BLISTR project. It served as a huge source of inspiration as to how to structure my code. [BLISTR](https://github.com/hc-BLISTR/BLISTR-main)
Nick Petronella - for his wisdom, intelligence, and mentorship.
Ana Pilar - for providing very-much needed feedback in regards to improving the UI and UX of the application
Sandeep, Robyn (BMH) - for taking the time to view our project during development, and offering suggestions such as the "Gene Lookup" feature
Mike Bostock - for his data visualization applets. View them here!
[Sortable Bar Chart](https://bl.ocks.org/mbostock/3885705)
[Force-Directed Graph](https://bl.ocks.org/mbostock/4062045)
[Bilevel Partition](https://bl.ocks.org/mbostock/5944371)

As well as the developers of the many programs and scripts that I have incorporated into our project, including ETE3, usearch/muscle, RAxML, FastTree, ... I am sorry if I have forgotten you!

## Contact

If you experience any bugs or difficulties in installing/using the HardCORE application, feel free to contact me (Brian) at brianhylee99@gmail.com . Nick at BFSSI may also be able to help; feel free to drop by his office or send an email at Nicholas.Petronella@hc-sc.gc.ca .

