import pandas as pd
from deepsmorfnet.model import run_model
from deepsmorfnet.hmmsearch import run_hmmsearch_nooutput
from os.path import join, dirname

PRODIGAL_PATH = join(dirname(__file__), 'bin/prodigal')
HMMSEARCH_PATH = join(dirname(__file__), 'bin/hmmsearch')
DSN_MODEL_PATH = join(dirname(__file__), 'data/keras/dsn_model.h5')
SMORFHMM_PATH = join(dirname(__file__), 'data/hmm/smorfams.hmm')

def run_hmmsearch(faa):
    run_hmmsearch_nooutput(HMMSEARCH_PATH, SMORFHMM_PATH, faa)


def run_dsn(infile):
    df = pd.read_csv(infile, sep='\t')
    model_preds = run_model(df.upstream, df.orf, df.downstream, DSN_MODEL_PATH)
    return model_preds