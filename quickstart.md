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
    * [kmerdb view](#kmerdb-view-and-kmerdb-header)
    * [kmerdb matrix](#kmerdb-matrix)
	* [kmerdb distance](#kmerdb-distance)
	* [kmerdb graph](#kmerdb-graph)
    * [kmerdb kmeans](#kmerdb-kmeans)
	* [kmerdb hierarchical](#kmerdb-hierarchical)
    * [kmerdb probability](#kmerdb-probability)
* [Documentation](#documentation)
* [Development](#development)
* [API](#api)
    * [Importing](#importing)
    * [Reading](#reading)
    * [Writing](#writing)
* [License](#license)
* [Acknowledgements](#acknowledgements)


# Introduction

## What is `.kdb`?

The K-mer database (.kdb) is a file format, a minimal python module (kmerdb) and command-line utility (CLI) packaged together in a pip module, for accessing precomputed k-mers from sequencing data. .kdb views k-mer graph databases as a fundamental technology to look at biological sequence identities and abundances. K-mer methods have been noted as fast and are used in areas from NGS quality control to phylogenomic investigations.

.kdb is based on the block GNU-zip file (bgzf) standard. Each .kdb file has a header or metadata section, much like .bam files. It is essentially a tab-delimited format with the last column unstructured for k-mer specific metadata. Input files and total k-mer counts are stored in the metadata block at the top of the file

Please visit the [Install](#/install) page for details on installation.

Visit the [Format](#/format) page for details on the `.kdb` format specification.

See the commands in the [Usage](#/usage) section for an idea of what functionality is built in to kdb.

See the original blog post on the concept [here](https://matthewralston.github.io/blog/kmer-database-format-part-1).


## What are k-mers?

Please see the [Wikipedia article](https://en.wikipedia.org/wiki/K-mer) on the subject for more details.

In brief, genomic/genetic, transcriptomics/dynamics, and perhaps even metagenomic/metagenetic spectra can be inferred directly from the k-mer spectra. And by keeping associations with reads, we can get to the spectral nature before any deeper or more performant methods need to be done on the datasets to generate other statistical calculations with a degree of precision perhaps not feasible with a pure k-mer frametwork at this point. But we can provide coverage information and pseudo-assembly information embedded in the graph. We next need a graph database layer and queries to perform certain de Brujin graph collapse methods in real time for the biologists to see.

## kdb is a file format

The k-mer database format is rather simple. It contains a metadata section, followed by a tab delimited format with an unstructured JSON as the final column. Both sections are combined in a compressed layer facilitated by Biopython's bio.bgzf module. 

Each file can be inspected with the view and header commands detailed in the [section below](#kdb-view).

## kdb features a bgzf graph database

Currently, the feature does not exist to transform the entire database to RDF for exploration with graph database visualization tools. But in the future, graphs could be loaded into Neptune databases for active app deployment. However, the fundamental relationships between the k-mers and each other, and the k-mers versus the dataset are all catloged with the `--all-metadata`


# Install

[![PyPI version](https://img.shields.io/pypi/v/kmerdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kmerdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kmerdb)
[![Coveralls code coverage](https://img.shields.io/coveralls/MatthewRalston/kmerdb/master.svg)][Coveralls]
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]

[pip]: https://pypi.org/project/kmerdb/
[Pythons]: https://pypi.org/project/kmerdb/
[Coveralls]: https://coveralls.io/r/MatthewRalston/kmerdb?branch=master
[RTD]: https://kdb.readthedocs.io/en/latest/

The current version on PyPI is shown above.

```bash
pip install kmerdb
```

See the [install](/installation) page for development install instructions.


# Usage



## How to view --help

Use '-h' to view detailed usage information about the subcommands

```bash
>kmerdb -h
usage: kmerdb [-h] {profile,header,graph,view,matrix,kmeans,hierarchical,distance,index,shuf,probability,citation} ...

positional arguments:
  {profile,header,graph,view,matrix,kmeans,hierarchical,distance,index,shuf,probability,citation}
                        Use --help with sub-commands
    profile             Parse data into the database from one or more sequence files
    header              Print the YAML header of the .kdb file and exit
    graph               Generate an adjacency list from .fa/.fq files
    view                View the contents of the .kdb file
    matrix              Generate a reduced-dimensionality matrix of the 4^k * n (k-mers x samples) data matrix.
    kmeans              Cluster the files according to their k-mer profile
    hierarchical        Use scipy.cluster.hierarchy to generate a dendrogram from a distance matrix
    distance            Calculate various distance metrics between profiles
    index               Create a index file that can be held in memory
    shuf                Create a shuffled .kdb file
    probability         Calculate the log-odds ratio of the Markov probability of a given sequence from the product (pi) of the transition probabilities(aij) times the frequency of the first k-mer (P(X1)), given the entire k-mer profile of a species. See https://matthewralston.github.io/quickstart#kmerdb-
                        probability for more details. 1. Durbin, R., Eddy, S.R., Krogh, A. and Mitchison, G., 1998. Biological sequence analysis: probabilistic models of proteins and nucleic acids. Cambridge university press.
    citation            Silence the citation notice on further runs

options:
  -h, --help            show this help message and exit
```

## kmerdb profile

A typical workflow first requires the generation of k-mer profiles. Complete metadata for each k-mer can be saved to the same database with `--all-metadata`. Note that this could cause significant increases in file size depending on the total k-mer coverage and the sequencing complexity. It is not recommended to experiment with `--all-metadata` at this time. Instead, we focus our attention on the numbers rather than the graph structure. Note that while individual profiles may be composite (i.e. you could mimic your own metagenomic compositions with downsampled fastq files to adjust proportions), the counts are stored in aggregate. All k-mer counts are stored in the header metadata, per-file.

```bash
usage: kmerdb profile [-h] [-v] [-k K] [-p PARALLEL] [--batch-size BATCH_SIZE] [-b FASTQ_BLOCK_SIZE] [-n N] [--both-strands] [--all-metadata] [--sorted] <.fasta|.fastq> [<.fasta|.fastq> ...] kdb

positional arguments:
  <.fasta|.fastq>       Fasta or fastq files
  kdb                   Kdb file

options:
  -h, --help            show this help message and exit
  -v, --verbose         Prints warnings to the console by default
  -k K                  Choose k-mer size (Default: 12)
  -p PARALLEL, --parallel PARALLEL
                        Shred k-mers from reads in parallel
  --batch-size BATCH_SIZE
                        Number of updates to issue per batch to PostgreSQL while counting
  -b FASTQ_BLOCK_SIZE, --fastq-block-size FASTQ_BLOCK_SIZE
                        Number of reads to load in memory at once for processing
  -n N                  Number of k-mer metadata records to keep in memory at once before transactions are submitted, this is a space limitation parameter after the initial block of reads is parsed. And during on-disk database generation
  --both-strands        Retain k-mers from the forward strand of the fast(a|q) file only
  --all-metadata        Include read-level k-mer metadata in the .kdb
  --sorted              Sort the output kdb file by count

```


The following command will generate multiple profiles at `$K`-mer resolution simultaneously. We have included 11 bacterial genomes under test/data for the user to explore. The following command, though in parallel, is applied to a single file, and thus represents the k-mer profile of the single file. 

```bash
parallel 'kmerdb profile -k $K {} {.}.$K.kdb' ::: $(/bin/ls test/data/*.fasta.gz)
```

## kmerdb view and kmerdb header

As mentioned before under [KDB format](#kdb-is-a-file-format), the kdb file consists of a header or metadata section, followed by data blocks until the end of the file. The header is YAML formatted and the data blocks are formatted as tab-separated value files (.tsv), with the last/right-most column being a JSON formatted metadata column. For developers, the YAML schema can be found in the config.py file.

It is worth noting that the 'files' attribute contains an array of file metadata for each file used to construct the profile, including total and unique k-mers per file *and* in aggregate. Note that you cannot sum unique k-mers across all files to retrieve the total unique count, rather an extra pass is required after the counts are summed to acquire the true nullomer/unique k-mer counts.

The view command also works on `.kdbg` files

```bash
# This should display the entire header of most files
>zcat test/data/Cdiff.8.kdb | head -n 30 
# This will also display just the header
>kmerdb header test/data/Cdiff.8.kdb
# The -H flag includes the header in the uncompressed output
# otherwise the command below produces just the tabular info.
>kmerdb view -H test/data/Cdiff.8.kdb
count_dtype: uint64
files:
- filename: test/data/Cdifficile_R3.fasta
  md5: ac539883c5104301f3d21de0c3813ad5
  mononucleotides:
    A: 1444649
    C: 601738
    G: 590294
    T: 1456466
  nullomers: 2337
  sha256: a96df44a1661bfa6db54b976da2ab0880e842c790188f7f18bb2b8b85f9caafc
  total_kmers: 4093140
  total_reads: 1
  unique_kmers: 63199
frequencies_dtype: float64
header_offset: 362
k: 8
kmer_ids_dtype: uint64
metadata: false
metadata_blocks: 1
profile_dtype: uint64
sorted: true
tags: []
total_kmers: 4093140
unique_kmers: 63199
unique_nullomers: 2337
version: 0.6.3


========================
0	1126	64835	0	0.0
1	1379	64140	0	0.0
2	1382	65205	0	0.0
3	1386	65501	0	0.0
4	1390	64989	0	0.0

...
```

## kmerdb matrix

The kmerdb matrix command generates the count matrix either un-normalized, normalized (via DESeq2), or with PCA or t-SNE dimensionality reduction applied. Note that default behavior of PCA if -n is not specified is to generate an elbow graph for the user to pick the appropriate choice of principal components for downstream analyses. The -n parameter is passed to the n_components parameter of sklearn.decomposition.PCA, which is commonly used for PCA in Python.

```bash
>kmerdb matrix -h
usage: kmerdb matrix [-h] [-v] [--with-index] [-d DELIMITER] [--output-delimiter OUTPUT_DELIMITER] [--no-normalized-ints] [-k K] [-n N] [--perplexity PERPLEXITY] {PCA,tSNE,Normalized,Unnormalized} [<kdbfile1 kdbfile2 ...|input.tsv|STDIN> ...]

positional arguments:
  {PCA,tSNE,Normalized,Unnormalized}
                        Choice of distance metric between two profiles
  <kdbfile1 kdbfile2 ...|input.tsv|STDIN>
                        Two or more .kdb files, or another count matrix in tsv/csv

options:
  -h, --help            show this help message and exit
  -v, --verbose         Prints warnings to the console by default
  --with-index          Print the row indices as well
  -d DELIMITER, --delimiter DELIMITER
                        The choice of delimiter to parse the input .tsv with. DEFAULT: ' '
  --output-delimiter OUTPUT_DELIMITER
                        The output delimiter of the final csv/tsv to write. DEFAULT: ' '
  --no-normalized-ints  Don't round normalized counts to the nearest integer
  -k K                  The k-dimension that the files have in common
  -n N                  The number of dimensions to reduce with PCA or t-SNE. DEFAULT: an elbow graph will be generated if -n is not provided to help the user choose -n
  --perplexity PERPLEXITY
                        The choice of the perplexity for t-SNE based dimensionality reduction

################################

# E x a m p l e s

################################


>kmerdb matrix -n 3 PCA test/data/*.$K.kdb
>kmerdb matrix Normalized test/data/*.$k.kdb
```



## kmerdb distance

Suppose you want a distance matrix between profiles; this is made easy with the distance command. The distance command supports all distance metrics used by `scipy.spatial.distance.pdist` to create the distance matrix/DataFrame and print it to STDOUT.

```bash
>kmerdb distance -h
usage: kmerdb distance [-h] [-v] [--output-delimiter OUTPUT_DELIMITER]
                    [-d DELIMITER] [-k K]
                    {braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,matching,minkowski,rogerstanimotorusselrao,seuclidean,sokalmichener,sokalsneath,spearman,sqeuclidean,yule}
                    <kdbfile1 kdbfile2 ...> [<kdbfile1 kdbfile2 ...> ...]
															
positional arguments:
  {braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,matching,minkowski,rogerstanimotorusselrao,seuclidean,sokalmichener,sokalsneath,spearman,sqeuclidean,yule}
                            Choice of distance metric between two profiles
  <kdbfile1 kdbfile2 ...>
                            Two or more .kdb files
																												
optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Prints warnings to the console by default
  --output-delimiter OUTPUT_DELIMITER
                        The output delimiter of the final csv/tsv to write.
  -d DELIMITER, --delimiter DELIMITER
                        The delimiter to use when printing the csv.
  -k K                  The k-dimension that the files have in common

>kmerdb distance spearman test/data/*.$K.kdb
>kmerdb distance correlation test/data/*.$K.kdb # Actually the Pearson correlation coefficient
```

The result is a symmetric matrix in tsv format with column headers formed from the filenames minus their extensions. It is presumed that to properly analyze the distance matrix, you would name the files after their sample name or their species, or some other identifying features. This naming convention holds for the `kdb matrix` command as well.

## kmerdb graph

Alternatively, you may want to see the edge list generated from one or more .fa/.fq files. This is the edge list generated by adjacencies supported by the provided data. It does not include all theoretically possible connections between all k-mers (i.e. the adjacency matrix is sparse). In fact, the `.kdbg` format is a (skinny) list of edges and their weights. Edge weights are calculated as the number of times the edge/k-mer pair was observed in the provided file(s). In the future, support for assembly may be available.



```bash
usage: kmerdb graph [-h] [-v] [-k K] [-p PARALLEL] [-b FASTQ_BLOCK_SIZE] [--both-strands] [--quiet] <.fasta|.fastq> [<.fasta|.fastq> ...] kdbg

positional arguments:
  <.fasta|.fastq>       Fasta or fastq files
  kdbg                  .kdbg file

options:
  -h, --help            show this help message and exit
  -v, --verbose         Prints warnings to the console by default
  -k K                  Choose k-mer size (Default: 12)
  -p PARALLEL, --parallel PARALLEL
                        Shred k-mers from reads in parallel
  -b FASTQ_BLOCK_SIZE, --fastq-block-size FASTQ_BLOCK_SIZE
                        Number of reads to load in memory at once for processing
  --both-strands        Retain k-mers from the forward strand of the fast(a|q) file only
  --quiet               Do not list all edges and neighboring relationships to stderr

```



## kmerdb kmeans

Now you would like to cluster your count matrix, perhaps after reducing dimensionality of the dataset. tSNE is a recommended choice for using the count matrix to differentiate between strains or species. PCA is equally useful in understanding the differences between species and reducing the dimensionality to something reasonable is an important first step before clustering.


```bash
# You may supply any arbitrary count tsv to cluster, and it will run that directly.
>kmerdb kmeans -k 3 sklearn example.tsv

# Alternatively, you may pipe the result from kdb matrix directly to kdb cluster.
# The matrix command will not produce PCA reduced dimension matrix unless the -n parameter
# is identified from the elbow graph it will produce as a side effect.
>kmerdb matrix PCA -n 3 test/data/*.$K.kdb | kmerdb kmeans -k 3 sklearn

# Alternatively, t-SNE may be used to project the data into 2 dimensions for visualization.
>kmerdb matrix tSNE -n 2 test/data/*.$K.kdb | kmerdb kmeans -k 3 Biopython
```

<!--
## kdb rarefy (deprecated)

Suppose you want supply a normalized matrix directly to [ecopy](https://github.com/Auerilas/ecopy) for rarefaction analysis. The rarefaction command requires a tsv and may be piped directly from 'kdb matrix'.


```bash
# You may supply any arbitrary count tsv to rarefy, and it will run that directly.
>./bin/kdb rarefy example.tsv

# Alternatively, you may pipe the result from kdb matrix directly to kdb rarefy.
>./bin/kdb matrix Unnormalized test/data/*.$K.kdb | ./bin/kdb rarefy
```
-->


## kmerdb hierarchical

Now you might want to use hierarchical clustering on a distance matrix for example.

```bash
# You may run the distance command with either direct kdb input or a count matrix
kmerdb distance spearman test/data/*.$K.kdb | kmerdb hierarchical
kmerdb matrix Normalized test/data/*.%K.kdb | kmerdb distance spearman [ /dev/stdin | STDIN ] | kmerdb hierarchical
```

## kmerdb probability

Calculate Markov chain probabilities from 

The following is equation 3.2 from p. 48 of reference 1.

<img src="https://render.githubusercontent.com/render/math?math=p(x) = P(X_{1})\prod_{i=2}^{N-k} a_{X_{i-1}X_{i}}">

where x is the full sequence, X<sub>1</sub> is the first k-mer subsequence of x, and a is the transition probability from X<sub>i-1</sub> to X<sub>i</sub> specified by the extension of equation 3.1 and 3.3 in p. 48-50.

<img src="https://render.githubusercontent.com/render/math?math=a_{st} = \frac{q_{t}}{\sum_{c=1}^{4} q_{s_{c}}}">

where c is one of the four possible suffixes prior to transition, q<sub>t</sub> is the frequency of the suffix transitioned to, q<sub>s</sub> are frequencies of each possible prefix.


I also must note the following equation I have had on my noteboard for about 12 months. The note at the bottom of the chalkboard says "take the log of both sides, since both represent the odds ratio."

<img src="https://render.githubusercontent.com/render/math?math=\frac{P(x|N)}{P(x|R)} = \frac{ p(x_{1}|N)\prod_{i=2}^{N-k} a_{X_{i-1}X_{i}} }{ p(x_{1}|R)\prod_{i=2}^{N-k} a_{ij} }">

Note that the a<sub>ij</sub> may be estimated from their maximum likelihood estimators for a first order markov model with a uniform random prior.

More specifically, p(x1\|R) is 1/4<sup>k</sup>. The transition probabilities, however, are estimated from the mononucleotide frequencies.

Then we log-transform the whole thing for a log-odds ratio test.

<img src="https://render.githubusercontent.com/render/math?math=\log10 \frac{P(x|N)}{P(x|R)} = \frac{ \log10 p(x_{1}|N) %2B \sum_{i=2}^{N-k} \log10 (a_{X_{i-1}X{i}}) }{ \log10 p(x_{1}|R) %2B \sum_{i=2}^{N-k} \log10 a_{ij} }">

The a<sub>ij</sub> will be further specified as follows. Each a<sub>ij</sub> is the frequency of the letter representing the transition in a first order Markov model

<img src="https://render.githubusercontent.com/render/math?math=P(x|R) = p(x_{1})\prod_{i=2}^{N-k} a_{ij}">

<img src="https://render.githubusercontent.com/render/math?math=\log10 P(x|R) = \log10 p(x_{1}) %2B \sum_{i=2}^{N-k} \log10 a_{ij}">

The meaning of three transition probability notations is as follows: a<sub>ij</sub> and a<sub>st</sub> and a<sub>Xi-1,Xi</sub> 

1. a<sub>st</sub>

The more generic form of a transition probability in general, the idea of the ratio of a specific frequency of an output sequence compared to all of the other 4 possible output states/nodes or whatever you want to refer to the states as in the Markov model, with respect to the prefix k-1mer that q<sub>t</sub> shares with each q<sub>sc</sub>, where c represents one other possible output states, as a single base on the end of the growing Markov chain.

2. a<sub>ij</sub>

A specific transition probability from k-mer i to k-mer j out of the L-k+1 k-mers of the input sequence, as modeled by the NULL model:

In this case, q<sub>t</sub> is the uniform random probability of a k-mer, namely the number of k-mers (proportional library size) divided by the number of possible k-mers (4<sup>k</sup>).

Then, for the actual model of the choice of state change, we add  the frequencies of the possible output sequences, which in the simplest first-order null model are the mononucleotide frequencies.

3. a<sub>Xi-1,Xi</sub>

The same transition probability as above, but under the ALTERNATIVE model. In our model, we use the frequency of the 1st/0th observed k-mer, P(X1\|R) = q<sub>X1</sub> to initialize the Markov chain.

Then, we use the observed frequencies of each outgoing state change as the denominator as before, to represent the first order state changes.



```bash
# This will create a tsv of probabilities and log-odds ratios.
kmerdb probability input.fasta example.kdb example.kdbi
```



1. Durbin, R., Eddy, S.R., Krogh, A. and Mitchison, G., 1998. Biological sequence analysis: probabilistic models of proteins and nucleic acids. Cambridge university press.



# Documentation


Documentation for the (sub)module can be found here: [https://kdb.readthedocs.io/en/latest](https://kdb.readthedocs.io/en/latest)

Additionally, running the distance, matrix, kmeans, or rarefy commands (which are arguably more complex in their usage) should be run with the DEBUG (-vv) verbosity setting on. This will yield additional information about the expected pipeline usage. That statement is echoed here.


```bash
The workflow is roughly as follows:

# # # # # # # # # #  #
# profile generation #
# # # # # # # # # #  #
# I have included -p and -b parameters that influence the rate of profile generation.
# The block size (-b) is primarily for the number of .fastq(.gz) records to read at a time.
# The -p parallel feature works within fastq parsing to process k-mer shredding/counting.
#
# -k $K is the parameter that specifies the k-mer resolution
#
# This command uses SQLite3 behind the scenes for on-disk k-mer counting
# since memory is rate limiting for profile generation when dealing 
# with biologically conventional choices of k (20 < k < 35).
parallel 'kdb profile -k $K {{}} {{.}}.$K.kdb' ::: $(/bin/ls test/data/*.fasta.gz)



# # # # # # #
# analysis  #
# # # # # # #
################################
# W A R N I N G :  M E M O R Y #
################################
# The first step of either rarefaction or clustering is to generate the k-mer profile matrix
# The matrix is not a new concept, it is just samples as columns of profiles. 
# Without metadata or graphical information embedded in the format, 
# we can create a matrix and do data science with the information.
# However, because the profiles are exponential in k, a linear increase in k 
# hypothetically results in an improvement of resolution that is at least superlinear in k.
# Therefore, the amount of memory can be calculated by the integer size times the profile resolution
# times the number of samples.
#
#
##################
# PCA + k-means
##################
# The first step ('kdb matrix') generates one from different profiles with the same choice of k.
# This command uses ecopy to normalize between sample k-mer total counts before PCA/tSNE.
# -n $N is the dimensionality of either PCA or tSNE. A good choice for tSNE is 2.
# If the command is run with PCA without a selected dimensionality, an elbow graph
# will be produced named 'PCA_variance_accumulation.png'. Please use this graph to select
# the number of principal components to use.
# The pipeline will not continue until -n $N is selected by the user.
# It is not recommended to feed Unnormalized or Normalized matrices directly to 'kdb kmeans'
# 
# The PCA/tSNE matrix will be dimReduced ($N) * N, where N is the number of samples/files/profiles.
#
# And finally, a k-means clustering will be done on the reduced dimensionality dataset
# Please note the randomness parameter 'random_state=42' for sklearn's kmeans is fixed at 42.
# Note here that the -k $K is not related to the choice of substring length 'k' for profile generation.
# The 'kdb kmeans' command produces two figures, first is an elbow graph looking at up to N clusters.
# This elbow graph will be written to 'kmeans_elbow_graph.png'.
# The second is the more typical scatterplot of the first two reduced dimensions
# and the k-means clustering labels shown over the scatter.
# This file will be written to 'kmeans_clustering.png'.
kdb matrix [-n $N] [ PCA | tSNE ] test/data/*.$K.kdb | kdb kmeans -k $K sklearn
kdb matrix [-n $N] [ PCA | tSNE ] test/data/*.$K.kdb | kdb kmeans -k $K --distance e Biopython
#
# If you wanted to save a matrix from kdb matrix for use on your own
# it is recommended that you consider gzip compressing it if it is the Normalized or Unnormalized matrix
# which we will see is used downstream in the rarefaction and hierarchical analytics pathways.
#


##################
# Hierarchical
##################
#
# The Normalized matrix goes to the distance subcommand, which can use any of scipy's pdist distances
# to form the m x m distance matrix.
# The third step (kdb hierarchical)  is to build a dendrogram with scipy.cluster.hierarchy.
# This final step produces a plot in addition to the tsvs produced in the prior steps,
# which can be captured as independent steps or with tee in a pipeline.
# The hierarchical clustering figure is written out to 'dendrogram.png'
kdb matrix [ Unnormalized | Normalized ] test/data/*.$K.kdb | kdb distance spearman | kdb hiearchical
```


# Development


[![PyPI version](https://img.shields.io/pypi/v/kmerdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kmerdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kmerdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kmerdb)
[![Coveralls code coverage](https://coveralls.io/repos/github/MatthewRalston/kmerdb/badge.svg?branch=master)](https://coveralls.io/github/MatthewRalston/kmerdb?branch=master)
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]


[pip]: https://pypi.org/project/kmerdb/
[Pythons]: https://pypi.org/project/kmerdb/
[RTD]: https://kdb.readthedocs.io/en/latest/

## Unit testing

The repository features Travis-CI pytest unit tests as the primary unit testing methodology for the functions, but not for all of the more complex IO-related operations, like those found in `bin/kdb`.

The method for installation and unit tests in a new development environment can be seen in '.travis.yml'. Currently, only unit tests are available for the suite. Acceptance testing has not been implemented. Until then, all acceptance testing is done manually prior to commits, squashes, rebases, and merges. Unit tests may be run with the following:

```bash
python setup.py test
```





# API

## Importing

The following section is designed to give developers a flavor of how to use the kdb modules for development purposes. The first topic of the API is importing.

```python
# Simple full module import
import kdb
# specific modules
from kdb import distance
from kdb import fileutil, index
```


## Reading

The next obvious topic is reading kdb files. It's very easy to read them line by line, so you should have no problems inventing your own incremental and custom python distance metrics for a reader that is capable of compression but uses a generator style. Note that the lines are parsed raw for speed, and they don't form any memory intensive objects. You are free to make the processing of each line as complex or simple as you want.

```python

from kdb import fileutil

files = ["foo.kdb", "bar.kdb", "baz.kdb"]

objects = [fileutil.open(f, 'r') for f in files]

#full_matrix = [o.slurp() for o in objects]

>for o in objects:
.  for line in o:
.    print(line) # Calculate something with each k-mer count
0    56
1    24

```
## Writing

Then we have the topic of writing kdb files. It's very easy again with the built in open method creating and returning instances of class KDBWriter (in this case). It's just a matter of iterating over your data and supplying it to the write method, which forwards to the Bio.bgzf write method through my wrapper class KDBWriter.


```python

from kdb import fileutil

with fileutil.open(f, 'w', header) as kdbfile:
  for x in list:
    kdbfile.write("\t".join(x) + "\n")

```

# License


Created by Matthew Ralston - [Scientist, Programmer, Musician](http://matthewralston.github.io) - [Email](mailto:mrals89@gmail.com)

Distributed under the Apache license. See 'LICENSE.txt' for the copy distributed with this project. Open source software is not for everyone, but for those of us starting out and trying to put the ecosystem ahead of ego, we march into the information age with this ethos.


# Acknowledgements


Thank you to the authors of kPAL and Jellyfish for the early inspiration. And thank you to others for the encouragement along the way, who shall remain nameless. I wanted this library to be a good strategy for assessing these k-mer profiles, in a way that is both cost aware of the analytical tasks at play, capable of storing the exact profiles in sync with the current assemblies, and then updating the kmer databases only when needed to generate enough spectral signature information.

The intention is that more developers would want to add functionality to the codebase or even just utilize things downstream, but to build out directly with numpy and scipy/scikit as needed to suggest the basic infrastructure for the ML problems and modeling approaches that could be applied to such datasets. This project has begun under GPL v3.0 and hopefully could gain some interest.

Also thank you to patelvivek (github/viviensio) for the Github ribbon I previously used
Thanks to free-icon for the download icon used on the home page. It is also cited at the bottom of each page.


Thanks to my former mentors BC, MR, IN, CR, and my newer bosses PJ and KL.
Thanks of course to Liftco Gymdustries for the facelift.
Thanks to the Pap lab and the DOE for the first dataset that I continue to use.
Thank you to Ryan for the food and stuff.
Thanks to Blahah for tolerating someone snooping and imitating his Ruby style.
Thanks to Erin for getting my feet wet in this new field.
Thanks to Rachel for the good memories and friendship.
Thanks to Yasmeen for the usual banter.
<!-- Thanks to Max, Robin, and Robert for the halfway decent memories in St. Louis. -->
And thanks to my family and friends.
Go Blue Hens 2021.
