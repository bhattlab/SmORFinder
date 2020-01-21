import click

from deepsmorfnet.help import CustomHelp
from deepsmorfnet.run import _run

from os.path import join, dirname


PRODIGAL_PATH = join(dirname(__file__), '../bin/prodigal')
HMMSEARCH_PATH = join(dirname(__file__), '../bin/hmmsearch')
SMORFHMM_PATH = join(dirname(__file__), '../data/hmm/smorfams.hmm')

@click.group(cls=CustomHelp)
def cli():
    """Command-line tool to predict and annotate small protein sequences in genomic sequencing data"""
    pass


@cli.command(short_help='Run deepsmorfnet on a complete or draft sequence of a single species.', help_priority=1)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='dsn_output')
@click.option('--prodigal-path', '-pp', default=PRODIGAL_PATH, type=click.Path(exists=True))
@click.option('--smorf-hmm-path', '-shp', default=SMORFHMM_PATH, type=click.Path(exists=True))
@click.option('--hmmsearch-path', '-hp', default=HMMSEARCH_PATH, type=click.Path(exists=True))
@click.option('--sensitive/--precise', default=False, help="Use a sensitive or precise model. default=precise")
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
def single(fasta, outdir, prodigal_path, smorf_hmm_path, hmmsearch_path, sensitive, force):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(command='run', fasta=fasta, outdir=outdir, prodigal_path=prodigal_path,
               smorf_hmm_path=smorf_hmm_path, hmmsearch_path=hmmsearch_path, sensitive=sensitive, force=force)

    _run(fasta, outdir, prodigal_path, smorf_hmm_path, hmmsearch_path, sensitive, force, mode='single')


@cli.command(short_help='Run deepsmorfnet on a metagenomic assembly.', help_priority=2)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='dsn_output')
@click.option('--threads', '-t', default=1)
@click.option('--prodigal-path', '-pp', default=PRODIGAL_PATH, type=click.Path(exists=True))
@click.option('--smorf-hmm-path', '-shp', default=SMORFHMM_PATH, type=click.Path(exists=True))
@click.option('--hmmsearch-path', '-hp', default=HMMSEARCH_PATH, type=click.Path(exists=True))
@click.option('--sensitive/--precise', default=False, help="Use a sensitive or precise model. default=precise")
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
def meta(fasta, outdir, threads, prodigal_path, smorf_hmm_path, hmmsearch_path, sensitive, force):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(command='run', fasta=fasta, outdir=outdir, threads=threads, prodigal_path=prodigal_path,
               smorf_hmm_path=smorf_hmm_path, hmmsearch_path=hmmsearch_path, sensitive=sensitive, force=force)

    _run(fasta, outdir, threads, prodigal_path, smorf_hmm_path, hmmsearch_path, sensitive, force, mode='meta')

def log_params(**kwargs):
    click.echo("#### PARAMETERS ####")
    click.echo('\n'.join(list(map(lambda x: ': '.join(list(map(str, x))), kwargs.items()))))
    click.echo("####################")

if __name__ == '__main__':

    cli()