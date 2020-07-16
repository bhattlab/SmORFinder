import sys
import click
from os import makedirs
from os.path import join, isdir
import shutil
from Bio import SeqIO
from smorfinder.prodigal import *
from smorfinder.hmmsearch import run_hmmsearch
from smorfinder.model import *
from smorfinder.finalize import _finalize


def _run(fasta, outdir, threads, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff, mode):
    tmp_dir = join(outdir, 'tmp')
    if force and isdir(outdir):
        shutil.rmtree(outdir)
    try:
        makedirs(tmp_dir)
    except FileExistsError:
        click.echo("Output directory exists, please delete or overwrite with --force")
        sys.exit(1)

    click.echo("Running Prodigal...")
    if mode == 'meta':
        run_prodigal_multithread(prodigal_path, fasta, tmp_dir, threads)
    else:
        run_prodigal(prodigal_path, fasta, tmp_dir, meta=False)

    rename_proteins(join(tmp_dir, 'prodigal.faa'), join(tmp_dir, 'prodigal.ffn'), join(tmp_dir, 'prodigal.gff'))
    click.echo("Filtering to only predicted genes less than or equal to 50 aa in length...")
    filter_prodigal_small_genes(tmp_dir)
    click.echo("Running HMMSEARCH...")
    run_hmmsearch(hmmsearch_path, smorf_hmm_path, join(tmp_dir, 'prodigal.small.faa'), tmp_dir)
    click.echo("Extracting nucleotide sequences...")
    names, fiveprime, orf, threeprime = extract_sequences(fasta, join(tmp_dir, 'prodigal.small.gff'))
    click.echo("Running deep learning models on predicted smORFs...")
    dsn1_predictions = run_model(fiveprime, orf, threeprime, dsn1_model_path)
    dsn2_predictions = run_model(fiveprime, orf, threeprime, dsn2_model_path)

    write_results_to_file(dsn1_predictions, dsn2_predictions, names, fiveprime, orf, threeprime, join(tmp_dir, 'model_predictions.tsv'))
    click.echo("Finalizing results...")
    _finalize(outdir, tmp_dir, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff)

