from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import (BooleanField, SelectField, SubmitField, IntegerField, 
    DecimalField, RadioField, StringField, validators)
from wtforms.validators import DataRequired
from flask_uploads import UploadSet, ARCHIVES
from datetime import datetime
from app import app

arch = UploadSet('archives', ARCHIVES)
DATA_FOLDER = app.config['DATA_FOLDER']
# ====================================================================
# Core Duplicate Form, like IndexForm, needs ident/plen values to calculate what
# genes are considered "similar enough" to a core gene
class CoreDuplicateForm(FlaskForm):
    ident = DecimalField('ident', [validators.NumberRange(min=0.5, max=1)])
    plen = DecimalField('plen', [validators.NumberRange(min=0.5, max=1)])
    submit = SubmitField('Update Parameters')
    download = SubmitField('Download table as .xlsx')

# IndexForm is used on the main page to upload .faa and .ffn archives
# (currently only supports .zip and .tar)
# Optionally, the user can choose specify a run name, if not,
# the current date+time will be used instead
class HardCoreForm(FlaskForm):
    faa = FileField('faa', validators=[
        FileRequired(),
        FileAllowed(arch.extensions, message='Archives only!')
    ])
    ffn = FileField('ffn', validators=[
        FileRequired(),
        FileAllowed(arch.extensions, message='Archives only!')
    ])
    run_name = StringField('run_name')  # optional field
    ident = DecimalField('ident', [validators.NumberRange(min=0.5, max=1)])
    plen = DecimalField('plen', [validators.NumberRange(min=0.5, max=1)])
    threads = IntegerField('cores', [validators.NumberRange(min=1, max=None),
                                     validators.optional()])

# SummaryTableForm is used by the run_table page to view a
# summary of results for each run in a table format,
# also handles downloading of results
class SummaryTableForm(FlaskForm):
    singles = BooleanField('singles', default=False)
    pan_genome = BooleanField('pan_genome', default=False)
    core_genome = BooleanField('core_genome', default=False)
    download = SubmitField('Download')
    delete = SubmitField('Delete Entries')
    export = SubmitField('Download table as .xlsx')
    cleanup = SubmitField('Clean-up')


class CoreSubsetForm(FlaskForm):
    submit = SubmitField('Submit')

class DownloadExcelForm(FlaskForm):
    download = SubmitField('Download .xlsx')

class DownloadForm(FlaskForm):
    download = SubmitField('Download')

class DownloadSubsetForm(DownloadExcelForm):
    alignments = SubmitField('Download fasta alignments')

# SNPForm is used by the SNP Viewer page, which handles displaying of SNP data, as well as
# visualization of the SNPs through SeaView, updating the reference genome, and downloading data
class SNPForm(FlaskForm):
    select_reference = SelectField('select_reference', choices=[])
    update_reference = SubmitField('Update reference')
    download_selected = SubmitField('Download selected alignments')
    download_all = SubmitField('Download all alignments')
    download_xlsx = SubmitField('Download table as .xlsx')

class RarefactionForm(FlaskForm):
    step = IntegerField('step', [validators.required(), validators.NumberRange(min=2)])
    reps = IntegerField('reps', [validators.required(), validators.NumberRange(min=1, max=10)])
    strains = FileField('test', [validators.optional()])
    download_strains = SubmitField('Download Strains')
    start_run = SubmitField('Run')

class TreeForm(FlaskForm):
    # General parameters (note: fasttree does not use most of these)
    builder = SelectField('build', choices=[('raxml', 'RAxML'), ('fasttree', 'FastTree')])
    datatype = SelectField('datatype', choices=[('aa', 'Amino acid'), ('nt', 'Nucleotide')])
    algorithm = StringField('algorithm', [validators.optional()])
    num_runs = IntegerField('num_runs', [validators.optional()])
    f_param = StringField('f_param', [validators.optional()])
    multithread = BooleanField('multithread')
    merge = BooleanField('merge')  # ie. clade collapsing
    download = SubmitField('Download Newick ("out.Tree")')
    submit = SubmitField('Submit')

    # ETE3-specific parameters
    show_align = BooleanField('show_align')
    show_bootstrap = BooleanField('show_bootstrap')
    # Visualization
    ete = SubmitField('Visualize using ETE3')
    dendroscope = SubmitField('View Tree in Dendroscope')
    pdf = SubmitField('Download PDF Image')
    # Tree building algorithms
    raxml = SubmitField('Build (RAxML)')
    fasttree = SubmitField('Build (FastTree)')

# Extends upon the Tree Form to allow file uploads
class TreeFormWithUpload(TreeForm):
    # currently, upload only allows .zip files (even for a single .fasta)
    upload = FileField('upload', validators=[
        FileRequired(),
        FileAllowed(['zip', 'fasta'])
    ])

# Used in gene_lookup; user uploads a fasta file with query sequences, and searches
# against all genomes in the run using usearch
class GeneLookupForm(FlaskForm):
    upload = FileField('upload', validators=[FileRequired()])
    datatype = SelectField('datatype', choices=[('nt', 'Nucleotide'), ('aa', 'Amino Acid')])
    submit = SubmitField('Submit')
    ident = DecimalField('ident', [validators.NumberRange(min=0, max=1)])
    plen = DecimalField('plen', [validators.NumberRange(min=0, max=1)])
    
# This form has buttons that take you directly to a multi-alignment
# featuring the query gene, and all the 'hits'
class LookupResultsForm(FlaskForm):
    select = SelectField('select', choices=[])  # choices will be filled out
    download = SubmitField('Download alignment')
    genomes_included = SubmitField('Download list')
    download_table = SubmitField('Download')

class MigrationForm(FlaskForm):
    export_runs = SubmitField("Export")
    import_runs = SubmitField("Import")
    upload = FileField("upload", validators=[
        FileAllowed(['zip'])
    ])

class NetworkForm(FlaskForm):
    download = SubmitField("Download")
    upload = FileField("upload", validators=[ validators.optional(), FileAllowed(['txt']) ])
    submit = SubmitField("Submit")
    restore_default = SubmitField("Restore Default")