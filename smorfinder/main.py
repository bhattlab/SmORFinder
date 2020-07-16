import click
from smorfinder import *
from smorfinder.help import CustomHelp
from smorfinder.run import _run

@click.group(cls=CustomHelp)
def cli():
    """Command-line tool to predict and annotate small protein sequences in genomic sequencing data"""
    pass


@cli.command(short_help='Run SmORFinder on a complete or draft genome assembly of a single species.', help_priority=1)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='smorf_output')
@click.option('--prodigal-path', '-pp', default=PRODIGAL_PATH, type=click.Path(exists=True))
@click.option('--dsn1-model-path', '-shp', default=DSN1_MODEL_PATH, type=click.Path(exists=True))
@click.option('--dsn2-model-path', '-shp', default=DSN2_MODEL_PATH, type=click.Path(exists=True))
@click.option('--smorf-hmm-path', '-shp', default=SMORFHMM_PATH, type=click.Path(exists=True))
@click.option('--hmmsearch-path', '-hp', default=HMMSEARCH_PATH, type=click.Path(exists=True))
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--dsn1-indiv-cutoff', '-idsn1', default=0.9999, help='Minimum cutoff necessary to keep prediction based on DSN1 significance cutoff alone. Between 0 and 1, default=0.9999')
@click.option('--dsn2-indiv-cutoff', '-idsn2', default=0.9999, help='Minimum cutoff necessary to keep prediction based on DSN2 significance cutoff alone. Between 0 and 1, default=0.9999')
@click.option('--phmm-indiv-cutoff', '-iphmm', default=1e-6, help='Minimum cutoff necessary to keep prediction based on pHMM significance cutoff alone. Between 0 and 1, default=1e-6')
@click.option('--dsn1-overlap-cutoff', '-odsn1', default=0.5, help='Minimum cutoff necessary to keep prediction based on DSN1 significance if both other models meet their respective cutoffs. Between 0 and 1, default=0.5')
@click.option('--dsn2-overlap-cutoff', '-odsn2', default=0.5, help='Minimum cutoff necessary to keep prediction based on DSN2 significance if both other models meet their respective cutoffs. Between 0 and 1, default=0.5')
@click.option('--phmm-overlap-cutoff', '-ophmm', default=1, help='Minimum cutoff necessary to keep prediction based on pHMM significance if both other models meet their respective cutoffs. Between 0 and 1, default=1')
def single(fasta, outdir, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(command='run', fasta=fasta, outdir=outdir, prodigal_path=prodigal_path, dsn1_model_path=dsn1_model_path,
               dsn2_model_path=dsn2_model_path,
               smorf_hmm_path=smorf_hmm_path, hmmsearch_path=hmmsearch_path, force=force, dsn1_indiv_cutoff=dsn1_indiv_cutoff, dsn2_indiv_cutoff=dsn2_indiv_cutoff, phmm_indiv_cutoff=phmm_indiv_cutoff, dsn1_overlap_cutoff=dsn1_overlap_cutoff, dsn2_overlap_cutoff=dsn2_overlap_cutoff, phmm_overlap_cutoff=phmm_overlap_cutoff)

    _run(fasta, outdir, 1, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff, mode='single')


@cli.command(short_help='Run SmORFinder on a metagenomic assembly.', help_priority=2)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='smorf_output')
@click.option('--threads', '-t', default=1)
@click.option('--prodigal-path', '-pp', default=PRODIGAL_PATH, type=click.Path(exists=True))
@click.option('--dsn1-model-path', '-shp', default=DSN1_MODEL_PATH, type=click.Path(exists=True))
@click.option('--dsn2-model-path', '-shp', default=DSN2_MODEL_PATH, type=click.Path(exists=True))
@click.option('--smorf-hmm-path', '-shp', default=SMORFHMM_PATH, type=click.Path(exists=True))
@click.option('--hmmsearch-path', '-hp', default=HMMSEARCH_PATH, type=click.Path(exists=True))
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--dsn1-indiv-cutoff', '-idsn1', default=0.9999, help='Minimum cutoff necessary to keep prediction based on DSN1 significance cutoff alone. Between 0 and 1, default=0.9999')
@click.option('--dsn2-indiv-cutoff', '-idsn2', default=0.9999, help='Minimum cutoff necessary to keep prediction based on DSN2 significance cutoff alone. Between 0 and 1, default=0.9999')
@click.option('--phmm-indiv-cutoff', '-iphmm', default=1e-6, help='Minimum cutoff necessary to keep prediction based on pHMM significance cutoff alone. Between 0 and 1, default=1e-6')
@click.option('--dsn1-overlap-cutoff', '-odsn1', default=0.5, help='Minimum cutoff necessary to keep prediction based on DSN1 significance if both other models meet their respective cutoffs. Between 0 and 1, default=0.5')
@click.option('--dsn2-overlap-cutoff', '-odsn2', default=0.5, help='Minimum cutoff necessary to keep prediction based on DSN2 significance if both other models meet their respective cutoffs. Between 0 and 1, default=0.5')
@click.option('--phmm-overlap-cutoff', '-ophmm', default=1, help='Minimum cutoff necessary to keep prediction based on pHMM significance if both other models meet their respective cutoffs. Between 0 and 1, default=1')
def meta(fasta, outdir, threads, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(command='run', fasta=fasta, outdir=outdir, threads=threads, prodigal_path=prodigal_path,
               dsn1_model_path=dsn1_model_path, dsn2_model_path=dsn2_model_path,
               smorf_hmm_path=smorf_hmm_path, hmmsearch_path=hmmsearch_path, force=force, dsn1_indiv_cutoff=dsn1_indiv_cutoff, dsn2_indiv_cutoff=dsn2_indiv_cutoff, phmm_indiv_cutoff=phmm_indiv_cutoff, dsn1_overlap_cutoff=dsn1_overlap_cutoff, dsn2_overlap_cutoff=dsn2_overlap_cutoff, phmm_overlap_cutoff=phmm_overlap_cutoff)

    _run(fasta, outdir, threads, prodigal_path, dsn1_model_path, dsn2_model_path, smorf_hmm_path, hmmsearch_path, force, dsn1_indiv_cutoff, dsn2_indiv_cutoff, phmm_indiv_cutoff, dsn1_overlap_cutoff, dsn2_overlap_cutoff, phmm_overlap_cutoff, mode='meta')

def log_params(**kwargs):
    click.echo("#### PARAMETERS ####")
    click.echo('\n'.join(list(map(lambda x: ': '.join(list(map(str, x))), kwargs.items()))))
    click.echo("####################")

if __name__ == '__main__':

    cli()
