# SmORFinder
A command line tool to identify and annotate small proteins in genomes and metagenomes.

# Installation

*SmORFinder* can be installed with pip:

    pip install smorfinder

You can then download the necessary data files just by running the command

    smorf

# Quick start

You can run *SmORFinder* on a single isolate genome with the command:

    smorf single myGenome.fna
  
And you can run it on a metagenome, with multiple threads, with the command:

    smorf meta myMetagenome.fna

# DBsmORF

This is a web server that allows you to download pre-computed smORF annotations for hundreds of thousands of RefSeq genomes and HMP metagenomes. It can also annotate your own genomes and annotations through a web form with *SmORFinder*.

## [Go to DBsMORF Now!](http://104.154.134.205:3838/DBsmORF/)
