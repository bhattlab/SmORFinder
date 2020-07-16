# DeepSmORFNet
A command line tool to identify and annotate small proteins in microbial sequencing datasets.

# Installation

*DeepSmORFNet* can be installed with pip:

  pip install deepsmorfnet

You can then download the necessary data files just by running the command

  dsn

# Quick start

You can run *DeepSmORFnet* on a single isolate genome with the command:

  dsn single myGenome.fna
  
And you can run it on a metagenome, with multiple threads, with the command:

  dsn meta myMetagenome.fna

# DBsmORF

This is a web server that allows you to download pre-computed smORF annotations for hundreds of thousands of RefSeq genomes and HMP metagenomes. It can also annotate your own genomes and annotations through

## [Go to DBsMORF Now!](http://104.154.134.205:3838/DBsmORF/)
