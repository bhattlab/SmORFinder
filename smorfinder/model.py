import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import sys
stderr = sys.stderr
import numpy as np
from tensorflow.keras import backend as K


def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall


def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision


def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2 * ((precision * recall) / (precision + recall + K.epsilon()))


def run_model(upstream_seqs, orf_seqs, downstream_seqs, dsn_model_path):

    sys.stderr = open(os.devnull, 'w')
    from tensorflow.keras.models import load_model
    sys.stderr = stderr

    model = load_model(dsn_model_path, custom_objects={'recall_m': recall_m, 'precision_m': precision_m, 'f1_m': f1_m})

    X_upstream, X_orf, X_downstream = prepare_dataset(upstream_seqs, orf_seqs, downstream_seqs)

    predictions = model.predict([X_orf, X_upstream, X_downstream], verbose=1)

    return predictions


def onehot_encode(seq, size, padding='downstream'):
    vec = np.zeros((size, 4))
    if not isinstance(seq, str):
        seq = ''
    for i, c in enumerate(seq.upper()):
        if padding == 'upstream':
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


def prepare_dataset(upstream, orf, downstream):
    
    X_orf = np.array([onehot_encode(seq, size=153) for seq in orf])
    X_upstream = np.array([onehot_encode(seq, size=100, padding='left') for seq in upstream])
    X_downstream = np.array([onehot_encode(seq, size=100, padding='right') for seq in downstream])

    return X_upstream, X_orf, X_downstream


def write_results_to_file(predictions_dsn1, predictions_dsn2, names, upstream, orf, downstream, outfile):

    outfile = open(outfile, 'w')
    print('seqid', 'dsn1_prob_smorf', 'dsn2_prob_smorf', '5p_seq', '3p_seq', 'orf_seq', sep='\t', file=outfile)
    for i, preds in enumerate(zip(list(predictions_dsn1[:,0]), list(predictions_dsn2[:,0]))):
        print(names[i], preds[0], preds[1], upstream[i], downstream[i], orf[i], sep='\t', file=outfile)
    outfile.close()
