from os.path import join, basename
from collections import defaultdict
from Bio import SeqIO
from os import makedirs
from shutil import move


exclude_smorfams = set(["smorfam03198", "smorfam03422", "smorfam03457", "smorfam03517", "smorfam03519", "smorfam03522", "smorfam03524", "smorfam03526", "smorfam03530", "smorfam03551", "smorfam03557", "smorfam03562", "smorfam03585", "smorfam03614", "smorfam03841"])

def _finalize(outdir, tmp_dir, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):

    final_prefix = basename(outdir)

    hmm_results = parse_hmm_results(join(tmp_dir, 'hmmsearch.tbl'))
    model_preds = parse_model_prediction(join(tmp_dir, 'model_predictions.tsv'))

    keep_ids = set()
    for seqid in model_preds:
        
        if seqid in hmm_results and hmm_results[seqid][0] in exclude_smorfams:
                continue

        phmm_evalue = None
        try:
            phmm_evalue = float(hmm_results[seqid][1])
        except:
            pass

        if meets_significance_cutoffs(float(model_preds[seqid]['dsn1_prob_smorf']), float(model_preds[seqid]['dsn2_prob_smorf']), phmm_evalue, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):
            keep_ids.add(seqid)

    keep_faa, keep_ffn = [], []
    for faa_rec, ffn_rec in zip(SeqIO.parse(join(tmp_dir, 'prodigal.small.faa'), 'fasta'),
                                SeqIO.parse(join(tmp_dir, 'prodigal.small.ffn'), 'fasta')):
        if faa_rec.id in keep_ids:
            keep_faa.append(faa_rec)
        if ffn_rec.id in keep_ids:
            keep_ffn.append(ffn_rec)

    SeqIO.write(keep_faa, join(outdir, final_prefix + '.faa'), 'fasta')
    SeqIO.write(keep_ffn, join(outdir, final_prefix + '.ffn'), 'fasta')

    final_table = []
    outgff = open(join(outdir, final_prefix + '.gff'), 'w')
    for rec in parse_gff(join(tmp_dir, 'prodigal.small.gff')):
        rec = rec.split('\t')
        rec_id = rec[-1].split(';')[0].split('=')[-1]

        if rec_id not in keep_ids:
            continue

        try:
            rec[-1] += 'smorfam=' + hmm_results[rec_id][0] + ';'
            rec[-1] += 'hmm_smorfam_evalue=' + str(hmm_results[rec_id][1]) + ';'
        except:
            rec[-1] += 'smorfam=None;hmm_smorfam_evalue=None;'

        rec[-1] += 'model_dsn1_prob_smorf=' + str(model_preds[rec_id]['dsn1_prob_smorf']) + ';'
        rec[-1] += 'model_dsn2_prob_smorf=' + str(model_preds[rec_id]['dsn2_prob_smorf']) + ';'
        rec[-1] += '5p_seq=' + model_preds[rec_id]['5p_seq'] + ';'
        rec[-1] += 'orf_seq=' + model_preds[rec_id]['orf_seq'] + ';'
        rec[-1] += '3p_seq=' + model_preds[rec_id]['3p_seq'] + ';'

        try:
            final_table.append(
                (rec_id, rec[0], rec[3], rec[4], rec[6], hmm_results[rec_id][0], str(hmm_results[rec_id][1]),
                 model_preds[rec_id]['dsn1_prob_smorf'], model_preds[rec_id]['dsn2_prob_smorf'], model_preds[rec_id]['5p_seq'], model_preds[rec_id]['orf_seq'],
                 model_preds[rec_id]['3p_seq'])
            )
        except:
            final_table.append(
                (rec_id, rec[0], rec[3], rec[4], rec[6], '', '',
                 model_preds[rec_id]['dsn1_prob_smorf'], model_preds[rec_id]['dsn2_prob_smorf'], model_preds[rec_id]['5p_seq'], model_preds[rec_id]['orf_seq'],
                 model_preds[rec_id]['3p_seq'])
            )

        print(*rec, sep='\t', file=outgff)
    outgff.close()

    with open(join(outdir, final_prefix + '.tsv'), 'w') as outfile:
        print('seqid', 'contig', 'start', 'end', 'orient', 'smorfam', 'hmm_smorfam_evalue', 'dsn1_prob_smorf', 'dsn2_prob_smorf',
              '5p_seq', 'orf', '3p_seq', sep='\t', file=outfile)
        for rec in final_table:
            rec = list(map(str, rec))
            print(*rec, sep='\t', file=outfile)


def meets_significance_cutoffs(dsn1_prob, dsn2_prob, phmm_evalue, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):

    
    if dsn1_prob > dsn1_indiv_cutoff or dsn2_prob > dsn2_indiv_cutoff:
        return True

    if phmm_evalue is not None and phmm_evalue < phmm_indiv_cutoff:
        return True

    if phmm_evalue is None:
        return False
    
    if dsn1_prob > dsn1_overlap_cutoff and dsn2_prob > dsn2_overlap_cutoff and phmm_evalue < phmm_overlap_cutoff:
        return True

    return False


def parse_gff(gff):
    with open(gff) as infile:
        for line in infile:
            if not line.startswith('#') and len(line) > 10:
                yield line.strip()

def parse_hmm_results(hmm_tbl):

    results = defaultdict(list)

    with open(hmm_tbl) as infile:
        for line in infile:
            if line.startswith('#'):
                continue

            line = line.strip().split()
            results[line[0]].append((line[2], float(line[4])))

    out = dict()
    for smorfam in results:
        res = sorted(results[smorfam], key=lambda x: x[1])
        smorfam_name = res[0][0]
        min_evalue = res[0][-1]

        out[smorfam] = (smorfam_name, min_evalue)

    return out

def parse_model_prediction(model_preds):

    model_preds_out = dict()
    with open(model_preds) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            line = {header[i]:line[i] for i in range(len(header))}
            model_preds_out[line['seqid']] = line

    return model_preds_out
