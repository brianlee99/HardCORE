import os

SECRET_KEY = '7d441f27d441f27567d441f2b6176a'

APP_FOLDER = os.path.join(os.getcwd(), 'app')
DATA_FOLDER = os.path.join(os.getcwd(), 'app/data/')
OUTPUT_FOLDER = os.path.join(os.getcwd(), 'app/output')
GRAPHS_FOLDER = os.path.join(os.getcwd(), 'app/static/graphs')
RF_FOLDER = os.path.join(os.getcwd(), 'app/static/rarefactions')
TREE_FOLDER = os.path.join(os.getcwd(), 'app/tree')
VISUALS_FOLDER = os.path.join(os.getcwd(), 'app/static/visuals')
HC_FOLDER = os.path.join(os.getcwd(), 'HardCORE_usearch_suite')
SCRIPTS_FOLDER = os.path.join(os.getcwd(), 'scripts')

#ALLOWED_EXTENSIONS = set(['zip', 'tar'])

# '' defaults to PROTGAMMALG
ACCEPTED_AA_ALGORITHMS = ['PROTGAMMABLOSUM62', 'PROTGAMMALG', 'PROTCATLG', 'PROTCATBLOSUM62', '']
# '' defaults to GTRGAMMA
ACCEPTED_NT_ALGORITHMS = ['GTRGAMMA', 'GTRCAT', '']

DEFAULT_NT_ALGORITHM = 'GTRGAMMA'
DEFAULT_AA_ALGORITHM = 'PROTGAMMALG'
