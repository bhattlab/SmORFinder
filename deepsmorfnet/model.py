import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import sys
stderr = sys.stderr
import numpy as np
from os.path import join, dirname


def run_model(fiveprime_seqs, orf_seqs, threeprime_seqs, sensitive):

    if sensitive == True:
        print("Using sensitive model...")
        MODEL_PATH = join(dirname(__file__), 'data/keras/sensitive_model.h5')
    else:
        print("Using precise model...")
        MODEL_PATH = join(dirname(__file__), 'data/keras/precise_model.h5')

    sys.stderr = open(os.devnull, 'w')
    from keras.models import load_model
    sys.stderr = stderr

    model = load_model(MODEL_PATH)

    X_fiveprime, X_orf, X_threeprime = prepare_dataset(fiveprime_seqs, orf_seqs, threeprime_seqs)

    predictions = model.predict([X_orf, X_fiveprime, X_threeprime])

    return predictions


def onehot_encode(seq, size, padding='threeprime'):
    vec = np.zeros((size, 4))
    if not isinstance(seq, str):
        seq = ''
    for i, c in enumerate(seq.upper()):
        if padding == 'fiveprime':
            i = i + size - len(seq)
        if c == 'A':
            vec[i, 0] = 1
        elif c == 'C':
            vec[i, 1] = 1
        elif c == 'G':
            vec[i, 2] = 1
        elif c == 'T':
            vec[i, 3] = 1

    return vec

def prepare_dataset(fiveprime, orf, threeprime):
    
    X_orf = np.array([onehot_encode(seq, size=153) for seq in orf])
    X_fiveprime = np.array([onehot_encode(seq, size=100, padding='left') for seq in fiveprime])
    X_threeprime = np.array([onehot_encode(seq, size=100, padding='right') for seq in threeprime])

    return X_fiveprime, X_orf, X_threeprime

def write_results_to_file(predictions, names, fiveprime, orf, threeprime, outfile):

    outfile = open(outfile, 'w')
    print('seqid', 'prob_smorf', '5p_seq', '3p_seq', 'orf_seq', sep='\t', file=outfile)
    for i, pred in enumerate(list(predictions[:,0])):
        print(names[i], pred, fiveprime[i], threeprime[i], orf[i], sep='\t', file=outfile)
    outfile.close()
