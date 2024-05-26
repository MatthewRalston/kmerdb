---
title: Quickstart
layout: post
description: 'Everything you need to get started'
image: assets/images/dna3.gif
nav-menu: true
toc: true
---


* [Introduction](#introduction)
* [Install](#install)
* [Usage](#usage)
    * [kmerdb profile](#kmerdb-profile)
    * [kmerdb view](#kmerdb-view)
    * [kmerdb matrix](#kmerdb-matrix)
	* [kmerdb distance](#kmerdb-distance)
    * [kmerdb kmeans](#kmerdb-kmeans)
    * [kmerdb hierarchical](#kmerdb-hierarchical)
* [Documentation](#documentation)
* [Development](#development)
    * [Unit testing](#unit-testing)
* [API](#api)
    * [Importing](#importing)
    * [Reading](#reading)
    * [Writing](#writing)
* [License](#license)
* [Acknowledgements](#acknowledgements)


# Introduction

## What is `kmerdb`?

The K-mer database profile `kmerdb` is a file format (.kdb) and command-line utility (CLI) for creating k-mer count vector profiles from sequencing data. kmerdb views k-mers and their De Bruijn graphs as fundamental tools for biological sequence abundances and sequence similarities. K-mer methods have been noted as fast and are used in areas from NGS quality control to phylogenomic investigations.

kmerdb is based on the block GNU-zip file (bgzf) standard. Each .kdb file has a metadata header section, much like .bam files. It is essentially a tab-delimited format with a YAML header. Input file checksums, unique and total k-mer counts, and nullomer counts are stored in the metadata block at the top of the file.

Please visit the [Install](#/install) page for details on installation.

See the commands in the [Usage](#/usage) section for an idea of what functionality is built in to kdb.

See the original blog post on the concept [here](https://matthewralston.github.io/blog/kmer-database-format-part-1).


The kdb project was designed to facilitate conversation between heavily optimized legacy codebases without much public attention, like Jellyfish, regarding the utility of standardizing k-mer frameworks. These frameworks are used throughout assembly and alignment hashing/seed-matching strategies. The primary goal of this project is documenting data shapes, compression strategies (which of course related to efficiency of storage, transmission, rapid access, etc.), and anticipating UI possibilities with the increases in read/write speeds afforded by improving SSD technologies and utilization of more channels of more rapid interfaces for data transmission (i.e. m2, NVMe, PCIx). 

## `kmerdb` makes k-mer count vectors

The k-mer database count vector profile file format is rather simple. It contains a YAML metadata header section, followed by a tab delimited format. Both sections are combined in a compressed layer facilitated by Biopython's `bio.bgzf` module. 

Each file can be inspected with the 'view' and 'header' commands detailed in the [section below](#kdb-view).


# Install

[![PyPI version](https://img.shields.io/pypi/v/kdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kdb)
[![Coveralls code coverage](https://img.shields.io/coveralls/MatthewRalston/kdb/master.svg)][Coveralls]
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]

[pip]: https://pypi.org/project/kmerdb
[Pythons]: https://pypi.org/project/kmerdb/
[Coveralls]: https://coveralls.io/r/MatthewRalston/kmerdb?branch=master
[RTD]: https://kmerdb.readthedocs.io/en/latest/

The current version on PyPI is shown above.

```bash
pip install kmerdb
```

See the [install](/installation) page for development install instructions.


# Usage


## How to view --help

Use '-h' to view detailed usage information about the subcommands

```bash
# Usage    --help option    --debug mode
kmerdb --help # [+ --debug mode]
kmerdb usage graph

**** 
 o-O      |||
o---O     |||             [|[          kmerdb           ]|]
O---o     |||
 O-o      |||        version :     v0.8.2
  O       |||
 o-O      |||        GitHub  : https://github.com/MatthewRalston/kmerdb/issues
o---O     |||         PyPI   : https://pypi.org/project/kmerdb/
O---o     |||      Website   : https://matthewralston.github.io/kmerdb
 O-o      |||
                                                                       lang :         python
                                                                          v :      >= v3.7.4

                      package manger : pip
                        version      : >= 24.0
        package root : /home/user/.local/share/virtualenvs/kdb-venv/lib/python3.12/site-packages/kmerdb
        exe file     : /home/user/.local/share/virtualenvs/kdb-venv/lib/python3.12/site-packages/kmerdb/__init__.py

                      required packages : 9
                   development packages : 9

           ARGV : ['/home/user/.local/share/virtualenvs/kdb-venv/bin/kmerdb', 'usage', 'graph']
        
O---o
 O-o
  O
 o-O
o---O
O---o
 O-o
  O
 o-O
o---O
O---o
 O-o
  O
 o-O
o---O




Beginning program...




                          [ name ] :         graph

                   description : create a edge list in (block) .gz format from .fasta|.fa or .fastq format.




   :     4 column output : [ row idx | k-mer id node #1 | k-mer id node #2 | edge weight (adjacency count) ]

   :  make a deBruijn graph, count the number of k-mer adjacencies,  printing the edge list to STDOUT




                  +=============+====================+====================+=================================+
                  <    row idx  |  k-mer id node #1  |  k-mer id node #2  |  edge weight (adjacency count)  >
                  |             |                    |                    |                                 |
                  |             +
                  |
                  |
                  |
                  |
                  |





--------------------------


                    kmerdb graph -k 12 input_1.fa [example_2.fastq] output.12.kdbg

                    [-]    inputs : 

                           Input file can .fastq (or .fa).   - gzip.  Output is a weighted edge list in .kdb format (gzipped .csv with YAML header)

                    [-]    parameters : 

                           uses < -k > for k-mer size, --quiet to reduce runtime, -v, -vv to control logging. --



                    [-]    [ usage ]  :  kmerdb graph -k $K --quiet <input_1.fa.gz> [input_2.fq.gz] <output_edge_list_file.12.kdbg>















name: arguments
type: array
items:
- name: k
  type: int
  value: choice of k-mer size
- name: quiet
  type: flag
  value: Write additional debug level information to stderr?




name: inputs
type: array
items:
- name: <.fasta|.fastq>
  type: array
  value: gzipped or uncompressed input .fasta or .fastq file(s)
- name: .kdbg
  type: file
  value: Output edge-list filepath.




name: features
type: array
items:
- name: k-mer count arrays, linear, produced as file is read through sliding window.
    (Un)compressed support for .fa/.fq.
  shortname: parallel faux-OP sliding window k-mer shredding
  description: Sequential k-mers from the input .fq|.fa files are added to the De
    Bruijn graph. In the case of secondary+ sequences in the .fa or considering NGS
    (.fq) data, non-adjacent k-mers are pruned with a warning. Summary statistics
    for the entire file are given for each file read, + a transparent data structure.
- name: k-mer neighbors assessed and tallied, creates a unsorted edge list, with weights
  shortname: weighted undirected graph
  description: an edge list of a De Bruijn graph is generated from all k-mers in the
    forward direction of .fa/.fq sequences/reads. i.e. only truly neighboring k-mers
    in the sequence data are added to the tally of the k-mer nodes of the de Bruijn
    graph and the edges provided by the data.

...



#   +

# [ 3 main features: ]     k-mer counts (kmerdb profile -k 12 <input.fa|.fq> [<input.fa|.fq>])    'De Bruijn' graph (kmerdb graph)         [matrices, distances, and clustering!]

# Create a [composite] profile of k-mer counts from sequence files. (.fasta|.fastq|.fa.gz|.fq.gz)
kmerdb profile -k 8 example_1.fq.gz [example_2.fq.gz] profile_1.8.kdb

# Build a weighted edge list (+ node ids/counts = De Bruijn graph)
kmerdb graph -k 12 example_1.fq.gz example_2.fq.gz edges_1.kdbg

# View k-mer count vector
kmerdb view profile_1.8.kdb # -H for full header

# Note: zlib compatibility
#zcat profile_1.8.kdb

# View header (config.py[kdb_metadata_schema#L84])
kmerdb header profile_1.8.kdb

## Optional normalization, dim reduction, and distance matrix features:

# K-mer count matrix - Cython Pearson coefficient of correlation [ ssxy/sqrt(ssxx*ssyy) ]
kmerdb matrix pass *.8.kdb | kmerdb distance pearson STDIN
# 
# kmerdb matrix DESeq2 *.8.kdb
# kmerdb matrix PCA *.8.kdb
# kmerdb matrix tSNE *.8.kdb
#   # <pass> just makes a k-mer count matrix from k-mer count vectors.
# 

# Distances on count matrices [ SciPy ]  pdists + [ Cython ] Pearson correlation, scipy Spearman and scipy correlation pdist calculations are available ]
kmerdb distance -h
# 
#
# usage: kmerdb distance [-h] [-v] [--debug] [-l LOG_FILE] [--output-delimiter OUTPUT_DELIMITER] [-p PARALLEL] [--column-names COLUMN_NAMES] [--delimiter DELIMITER] [-k K]
#                       {braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,matching,minkowski,pearson,rogerstanimotorusselrao,seuclidean,sokalmichener,sokalsneath,spearman,sqeuclidean,yule} [<kdbfile1 kdbfile2 ...|input.tsv|STDIN> ...]

# +

#    Kmeans (sklearn, BioPython)
kmerdb kmeans -k 4 -i dist.tsv
#    BioPython Phylip tree + upgma
kmerdb hierarchical -i dist.tsv

```

## kmredb profile

A typical workflow first requires the generation of k-mer count vector profiles. The following command will generate multiple profiles at `$K`-mer resolution simultaneously.

```bash
parallel 'kmerdb profile -k $K {} {.}.$K.kdb' ::: $(/bin/ls test/data/*.fasta.gz)
```

NOTE: the test datasets include with the kmerdb git repository may be used to practice creating profiles, count matrices, distance matrices, and clustering output from small microbial genomic datasets.

## kmerdb

As mentioned before under [KDB format](#kdb-is-a-file-format), the kdb file consists of a header or metadata section, followed by data blocks until the end of the file. The header is YAML formatted and the data blocks are formatted as tab-separated value files (.tsv), with the last/right-most column being a JSON formatted metadata column. For developers, the YAML schema can be found in the config.py file.
```bash
# This should display the entire header of most files
>zcat example.8.kdb | head -n 30 
# This will also display just the header
>kmerdb header example.8.kdb
# The -H flag includes the header in the uncompressed output
>kmerdb view -H example.8.kdb
version: 0.8.3
metadata_blocks: 1
k: 8
total_kmers: 4132866
unique_kmers: 64103
unique_nullomers: 1433
metadata: false
sorted: false
kmer_ids_dtype: uint64
profile_dtype: uint64
count_dtype: uint64
frequencies_dtype: float64
tags: []
files:
- filename: test/data/Cacetobutylicum_ATCC824.fasta.gz
  md5: 0a0f73e1c8b8285703e29279bafaabef
  sha256: f9081291b62ff3387f1ca6ee2484669c849ed1840fdf2dd9dc3a0c93e9e87951
  total_reads: 2
  total_kmers: 4132866
  unique_kmers: 64103
  nullomers: 1433
header_offset: 324


========================

0	0	572	0.00013840274521361204
1	1	454	0.00010985112994227251
2	2	716	0.0001732453943582976
3	3	1256	0.00030390532865086844

     ...

```


## kmerdb matrix

Input to the matrix command is either a space separated list (positional arguments) of `.kdb` files, OR a input .tsv file OR otherwise, a .tsv can be read from STDIN by using 'STDIN' as the final positional argument.

The matrix subcommand generates the count matrix either un-normalized, normalized (via DESeq2), or with PCA or t-SNE dimensionality reduction applied. Note that default behavior of PCA if -n is not specified is to generate an elbow graph for the user to pick the appropriate choice of principal components for downstream analyses. The -n parameter is passed to the n_components parameter of sklearn.decomposition.PCA, which is commonly used for PCA in Python.



```bash
>kmerdb matrix -h
usage: kmerdb matrix [-h] [-v] [--debug] [-l LOG_FILE] [-nl {20,50,100,200}]
                     [-p PARALLEL] [--with-index]
                     [--column-names COLUMN_NAMES] [--delimiter DELIMITER]
                     [--output-delimiter OUTPUT_DELIMITER]
                     [--no-normalized-ints] [-k K] [-n N]
                     [--perplexity PERPLEXITY]
                     {PCA,tSNE,DESeq2,pass,Frequency}
                     [<kdbfile1 kdbfile2 ...|input.tsv|STDIN> ...]

positional arguments:
  {PCA,tSNE,DESeq2,pass,Frequency}
                        Choice of dimensionality reduction, normalization
                        method (DESeq2), or pass (no action)
  <kdbfile1 kdbfile2 ...|input.tsv|STDIN>
                        Two or more .kdb files, or another count matrix in
                        tsv/csv

options:
  -h, --help            show this help message and exit
  -v, --verbose         Prints warnings to the console by default
  --debug               Debug mode. Do not format errors and condense log
  -l LOG_FILE, --log-file LOG_FILE
                        Destination path to log file
  -nl {20,50,100,200}, --num-log-lines {20,50,100,200}
                        Number of logged lines to print to stderr. Default: 50
  -p PARALLEL, --parallel PARALLEL
                        Read files in parallel
  --with-index          Print the row indices as well
  --column-names COLUMN_NAMES
                        A filepath to a plaintext flat file of column names.
  --delimiter DELIMITER
                        The choice of delimiter to parse the input .tsv with.
                        DEFAULT: ' '
  --output-delimiter OUTPUT_DELIMITER
                        The output delimiter of the final csv/tsv to write.
                        DEFAULT: ' '
  --no-normalized-ints  Don't round normalized counts to the nearest integer
  -k K                  The k-dimension that the files have in common
  -n N                  The number of dimensions to reduce with PCA or t-SNE.
                        DEFAULT: an elbow graph will be generated if -n is not
                        provided to help the user choose -n
  --perplexity PERPLEXITY
                        The choice of the perplexity for t-SNE based
                        dimensionality reduction
>kmerdb matrix -n 3 PCA test/data/*.$K.kdb # Will use principal components analysis/SVD to reduce the dimensions of the count matrix
```

## kmerdb distance

Similar to the matrix subcommand, arguments to 'kmerdb distance' are multiple .kdb files, a input.tsv file, or the keyword 'STDIN' to read from standard input.

Suppose you want a distance matrix between profiles; this is made easy with the distance command, and prints the distance matrix as .tsv to STDOUT

```bash
>kmerdb distance correlation test/data/*.$K.kdb
```

The result is a symmetric matrix in tsv format with column headers formed from the filenames with the extensions/suffixes removed.



## kmerdb kmeans

Now you would like to cluster your count matrix, perhaps after reducing dimensionality of the dataset. tSNE is a recommended choice for using the count matrix to differentiate between strains or species. PCA is equally useful in understanding the differences between species and reducing the dimensionality to something reasonable is an important first step before clustering.


```bash
# You may supply any arbitrary count tsv to cluster, and it will run that directly.
>kmerdb kmeans -k 3 -i example.tsv

# Alternatively, you may pipe the result from kdb matrix directly to kdb cluster.
# The matrix command will not produce PCA reduced dimension matrix unless the -n parameter is passed as well.
# The choice of dimensions, n, can be identified from the elbow graph produced if -n is not passed.
>kmerdb matrix PCA -n 3 test/data/*.$K.kdb | kmerdb kmeans -k 3 -i STDIN

# Alternatively, t-SNE may be used to project the data into 2 dimensions for visualization.
>kmerdb matrix tSNE -n 2 test/data/*.$K.kdb | kmerdb cluster -k 3 -i STDIN
```

## kmerdb hierarchical

Usage for the 'hierarchical' subcommand is very similar to 'kmeans'. Be sure to use -i (or STDIN) to specify input, as its not positional in this command.



# Documentation

Documentation for the functions, parameters, positional arguments, and steps used by the subcommands may be obtained from the argparse -h,--help arguments, the 'usage' command (e.g. `kmerdb usage profile`), the README.md file, or the module's docstrings in the ReadTheDocs documentation.




# Development


[![PyPI version](https://img.shields.io/pypi/v/kmerdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kmerdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kmerdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kmerdb)
[![Coveralls code coverage](https://img.shields.io/coveralls/MatthewRalston/kmerdb/master.svg)][Coveralls]
[![ReadTheDocs status](https://readthedocs.org/projects/kmerdb/badge/?version=stable&style=flat)][RTD]


[pip]: https://pypi.org/project/kmerdb/
[Pythons]: https://pypi.org/project/kmerdb/
[Coveralls]: https://coveralls.io/r/MatthewRalston/kmerdb?branch=master
[RTD]: https://kmerdb.readthedocs.io/en/latest/

## Unit testing

The repository features Travis-CI pytest unit tests as the primary unit testing methodology for the functions, but not for all of the more complex IO-related operations, like those found in `bin/kdb`.

The method for installation and unit tests in a new development environment can be seen in '.travis.yml'. Currently, only unit tests are available for the suite. Acceptance testing has not been implemented. Until then, all acceptance testing is done manually prior to commits, squashes, rebases, and merges. Unit tests may be run with the following:

```bash
python setup.py test
```









# API

## Importing

The following section is designed to give developers a flavor of how to use the kdb modules for development purposes. The first topic of the API is importing.

```py
```python
# Simple full module import

# specific modules
from kmerdb import distance
from kmerdb import fileutil, graph
```
````

## Reading

The next obvious topic is reading kdb files. It's very easy to read them line by line, so you should have no problems inventing your own incremental and custom python distance metrics for a reader that is capable of compression but uses a generator style. Note that the lines are parsed raw for speed, and they don't form any memory intensive objects. You are free to make the processing of each line as complex or simple as you want.

```python

import yaml
from kdb import fileutil

files = ["foo.kdb", "bar.kdb", "baz.kdb"]

profiles = [fileutil.open(f, 'r', slurp=True) for f in files]

print(yaml.dump(profiles[0].metadata)) # Prints the YAML header

foo_profile = profiles[0]
# foo.metadata # YAML metadata
# foo.kmer_ids # kmer_ids
# foo.counts # a NumPy array of the counts

>for i, c in enumerate(o.counts):
    print(i, c) # Calculate something with each k-mer count
0    56
1    24


Check out the readthedocs or `__init__` source file for more details on module usage.

# License


Created by Matthew Ralston - [Scientist, Programmer, Musician](http://matthewralston.github.io) - [Email](mailto:mrals89@gmail.com)

Distributed under the Apache v2.0 license. See 'LICENSE.txt' for the copy distributed with this project. Open source software is not for everyone.


# Acknowledgements


Thank you to the authors of kPAL and Jellyfish for the early inspiration. And thank you to others for the encouragement along the way, who shall remain nameless. I wanted this library to be a good strategy for assessing these k-mer profiles, in a way that is both cost aware of the analytical tasks at play, capable of storing the exact profiles in sync with the current assemblies, and then updating the kmer databases only when needed to generate enough spectral signature information.

The intention is that more developers would want to add functionality to the codebase or even just utilize things downstream, but to build out directly with numpy and scipy/scikit as needed to suggest the basic infrastructure for the ML problems and modeling approaches that could be applied to such datasets. This project has begun under GPL v3.0 and hopefully could gain some interest.

Thanks to free-icon for the download icon used on the home page. It is also cited at the bottom of each page.


Thanks to my former mentors BC, MR, IN, CR, and my newer bosses PJ and KL.
Thanks of course to Liftco Gymdustries for the facelift.
Thanks to the Pap lab and the DOE for the first dataset that I continue to use.
Thank you to Ryan for the food and stuff.
Thanks to Blahah for tolerating someone snooping and imitating his Ruby style.
Thanks to Erin for getting my feet wet in this new field.
Thanks to Rachel for the good memories and friendship.
Thanks to Yasmeen for the usual banter.
Thanks to Max, Robin, and Robert for the halfway decent memories in St. Louis.
And thanks to my family and friends.
Go Blue Hens 2021.
