import subprocess
from os.path import join
from os import remove, devnull
from Bio import SeqIO
from random import choices
import string
import gzip
import math
from multiprocessing import Pool
from itertools import cycle


def run_prodigal(prodigal_path, infile, outdir, meta):
    bashCommand = '{prodigal} -c -f gff -o {gff} -a {faa} -d {ffn} -i {input}'.format(
        prodigal=prodigal_path, gff=join(outdir, 'prodigal.gff'), faa=join(outdir, 'prodigal.faa'),
        ffn=join(outdir, 'prodigal.ffn'), input=infile
    )
    if meta == True:
        bashCommand += ' -p meta'
    print('prodigal command:', bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

def run_prodigal_simple(prodigal_path, infile, outprefix):
    FNULL = open(devnull, 'w')
    bashCommand = '{prodigal} -p meta -c -f gff -o {gff} -a {faa} -d {ffn} -i {input}'.format(
        prodigal=prodigal_path, gff=outprefix + '.gff', faa=outprefix + '.faa',
        ffn=outprefix + '.ffn', input=infile
    )

    print('prodigal command:', bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=FNULL, stderr=FNULL)
    output, error = process.communicate()


def run_prodigal_multithread(prodigal_path, infile, outdir, threads):

    print("Counting records in FASTA file...")
    record_count = count_records_in_fasta(infile)
    print("The FASTA file contains %d records..." % record_count)

    print("Writing FASTA file to batches for multithreading..." )
    record_count_per_file = math.ceil(record_count / threads)
    filecount = 0
    outrecs = []
    outfiles = []

    if infile.endswith('.gz'):
        infile = gzip.open(infile, "rt")

    for record in SeqIO.parse(infile, 'fasta'):
        outrecs.append(record)
        if len(outrecs) == record_count_per_file:
            filecount += 1
            outfile = join(outdir, 'input%d.fna' % filecount)
            SeqIO.write(outrecs, outfile, 'fasta')
            outfiles.append(outfile)
            outrecs = []

    if len(outrecs) > 0:
        filecount += 1
        outfile = join(outdir, 'input%d.fna' % filecount)
        SeqIO.write(outrecs, outfile, 'fasta')
        outfiles.append(outfile)
    del outrecs

    prodigal_files = [join(outdir, 'prodigal%d' % (i+1)) for i in range(len(outfiles))]
    with Pool(processes=threads) as pool:
        pool.starmap(run_prodigal_simple, zip(cycle([prodigal_path]), outfiles, prodigal_files))

    combine_files(outfiles, join(outdir, 'input.fna'))
    combine_files([f+'.faa' for f in prodigal_files], join(outdir, 'prodigal.faa'))
    combine_files([f+'.ffn' for f in prodigal_files], join(outdir, 'prodigal.ffn'))
    combine_files([f+'.gff' for f in prodigal_files], join(outdir, 'prodigal.gff'), ['#'])


def combine_files(files, outfile, exclude_startswith=[]):
    with open(outfile, 'w') as out:

        for f in files:

            with open(f) as infile:
                for line in infile:

                    skip = False
                    for exclude in exclude_startswith:
                        if line.startswith(exclude):
                            skip = True
                            break
                    if skip is True:
                        continue
                    out.write(line)

    for f in files:
        remove(f)



def revcomp(seq):
    out = ''
    for c in seq[::-1]:
        if c == 'A':
            out += 'T'
        elif c == 'T':
            out += 'A'
        elif c == 'C':
            out += 'G'
        elif c == 'G':
            out += 'C'
        else:
            out += c
    return out

def count_records_in_fasta(fasta):
    records = 0
    if fasta.endswith('.gz'):
        with gzip.open(fasta, "rt") as infile:
            for line in infile:
                if line.startswith('>'):
                    records += 1
    else:
        with open(fasta) as infile:
            for line in infile:
                if line.startswith('>'):
                    records += 1
    return records

def parse_gff(gff):
    with open(gff) as infile:
        for line in infile:
            if not line.startswith('#') and len(line) > 10:
                yield line.strip()


def rename_proteins(faa, ffn, gff):
    out_faa, out_ffn, out_gff = [], [], []
    prefix = ''.join(list(choices(string.ascii_uppercase, k=6)))
    num = 0
    for faa_rec, ffn_rec, gff_rec in zip(SeqIO.parse(faa, 'fasta'), SeqIO.parse(ffn, 'fasta'), parse_gff(gff)):
        num += 1
        new_id = prefix+'_'+str(num)

        faa_rec.id = new_id
        faa_rec.name = new_id
        desc = faa_rec.description.split()
        desc[-1] = 'ID='+new_id+';'+';'.join(desc[-1].split(';')[1:])
        faa_rec.description = ' '.join(desc)
        out_faa.append(faa_rec)

        ffn_rec.id = new_id
        ffn_rec.name = new_id
        desc = ffn_rec.description.split()
        desc[-1] = 'ID=' + new_id + ';' + ';'.join(desc[-1].split(';')[1:])
        ffn_rec.description = ' '.join(desc)
        out_ffn.append(ffn_rec)

        gff_rec = gff_rec.split('\t')
        gff_rec[-1] = 'ID=' + new_id + ';' + ';'.join(gff_rec[-1].split(';')[1:])
        out_gff.append('\t'.join(gff_rec))

    SeqIO.write(out_faa, faa, 'fasta')
    SeqIO.write(out_ffn, ffn, 'fasta')
    with open(gff, 'w') as outfile:
        print(*out_gff, sep='\n', file=outfile)


def filter_prodigal_small_genes(outdir):

    ffn_recs = []
    for rec in SeqIO.parse(join(outdir, 'prodigal.ffn'), 'fasta'):
        if len(rec.seq) <= 153:
            ffn_recs.append(rec)
    SeqIO.write(ffn_recs, join(outdir, 'prodigal.small.ffn'), 'fasta')

    faa_recs = []
    for rec in SeqIO.parse(join(outdir, 'prodigal.faa'), 'fasta'):
        if len(rec.seq) <= 51:
            faa_recs.append(rec)
    SeqIO.write(faa_recs, join(outdir, 'prodigal.small.faa'), 'fasta')

    outgff = open(join(outdir, 'prodigal.small.gff'), 'w')
    with open(join(outdir, 'prodigal.gff')) as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')

            if line[2] != 'CDS':
                continue

            if int(line[4]) - int(line[3]) + 1 <= 153:
                print('\t'.join(line), file=outgff)
    outgff.close()


def extract_sequences(fasta, gff):

    names, all_fiveprime, all_orfs, all_threeprime = [], [], [], []

    if fasta.endswith('.gz'):
        fasta_file = gzip.open(fasta, "rt")
        fasta = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, 'fasta')}
        fasta_file.close()
    else:
        fasta = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta, 'fasta')}

    with open(gff) as infile:
        for line in infile:


            line = line.strip().split('\t')
            contig, start, end, orient, seqname = line[0], int(line[3])-1, int(line[4]), line[6], line[-1].split(';')[0].split('=')[-1]
            orf = fasta[contig][start:end].upper()

            left_start = start - 100
            if left_start < 0:
                left_start = 0
            right_end = end + 100

            left_seq = fasta[contig][left_start:start].upper()
            right_seq = fasta[contig][end:right_end].upper()

            if orient == '-':
                orf = revcomp(orf)
                left_seq = revcomp(left_seq)
                right_seq = revcomp(right_seq)
                tmp = left_seq
                left_seq = right_seq
                right_seq = tmp

            names.append(seqname)
            all_fiveprime.append(left_seq)
            all_orfs.append(orf)
            all_threeprime.append(right_seq)

    return names, all_fiveprime, all_orfs, all_threeprime
