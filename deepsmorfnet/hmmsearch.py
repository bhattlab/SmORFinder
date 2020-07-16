import subprocess
from os.path import join


def run_hmmsearch(hmmsearch_path, hmmfile_path, infile, outdir):
    bashCommand = '{hmmsearch} -E 1.0 -o {out} --tblout {outtbl} {hmmfile} {seqdb}'.format(
        hmmsearch=hmmsearch_path, out=join(outdir, 'hmmsearch.out'), outtbl=join(outdir, 'hmmsearch.tbl'),
        hmmfile=hmmfile_path, seqdb=infile
    )
    print('hmmsearch command:', bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

def run_hmmsearch_nooutput(hmmsearch_path, hmmfile_path, infile):
    bashCommand = '{hmmsearch} -E 1.0 {hmmfile} {seqdb}'.format(
        hmmsearch=hmmsearch_path,
        hmmfile=hmmfile_path, seqdb=infile
    )
    print('hmmsearch command:', bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

def filter_hmmsearch_table(tbl, outtbl):

    with open(tbl) as infile:
        for line in infile:
            pass
