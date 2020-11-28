from smorfinder.model import run_model
from smorfinder.hmmsearch import run_hmmsearch_nooutput
from os.path import join, dirname, isfile
from os import makedirs
import wget
import sys

__version__ = '1.0.0'

# Check system OS, using Linux as catchall (Windows not explicitly supported)
platform = 'mac' if sys.platform == 'darwin' else 'linux'

PRODIGAL_PATH = join(dirname(__file__), 'bin/%s/prodigal' % platform)
HMMSEARCH_PATH = join(dirname(__file__), 'bin/%s/hmmsearch' % platform)
DSN1_MODEL_PATH = join(dirname(__file__), 'data/keras/dsn1_model.h5')
DSN2_MODEL_PATH = join(dirname(__file__), 'data/keras/dsn2_model.h5')
SMORFHMM_PATH = join(dirname(__file__), 'data/hmm/smorfams.hmm')

def download_data():
    if not isfile(SMORFHMM_PATH):
        print("Downloading HMM models...")
        makedirs(dirname(SMORFHMM_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/deepsmorfnet/downloads/smorfams.hmm', SMORFHMM_PATH)
        print()
    if not isfile(DSN1_MODEL_PATH):
        print("Downloading DSN1 model...")
        makedirs(dirname(DSN1_MODEL_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/deepsmorfnet/downloads/dsn1_model.h5', DSN1_MODEL_PATH)
        print()
    if not isfile(DSN2_MODEL_PATH):
        print("Downloading DSN2 model...")
        makedirs(dirname(DSN2_MODEL_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/deepsmorfnet/downloads/dsn2_model.h5', DSN2_MODEL_PATH)
        print()

def run_hmmsearch(faa):
    run_hmmsearch_nooutput(HMMSEARCH_PATH, SMORFHMM_PATH, faa)

def run_dsn1(infile):
    upstream, orf, downstream = [], [], []
    with open(infile) as infile:
        header = infile.readline().strip('\n').split('\t')
        for line in infile:
            line = line.strip('\n').split('\t')
            line = {header[i]:line[i] for i in range(len(header))}
            upstream.append(line['upstream'])
            orf.append(line['orf'])
            downstream.append(line['downstream'])
    model_preds = run_model(upstream, orf, downstream, DSN1_MODEL_PATH)
    return model_preds

def run_dsn2(infile):
    upstream, orf, downstream = [], [], []
    with open(infile) as infile:
        header = infile.readline().strip('\n').split('\t')
        for line in infile:
            line = line.strip('\n').split('\t')
            line = {header[i]:line[i] for i in range(len(header))}
            upstream.append(line['upstream'])
            orf.append(line['orf'])
            downstream.append(line['downstream'])
    model_preds = run_model(upstream, orf, downstream, DSN2_MODEL_PATH)
    return model_preds

download_data()
