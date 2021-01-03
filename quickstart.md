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
    * [kdb profile](#kdb-profile)
    * [kdb view](#kdb-view)
    * [kdb distance](#kdb-distance)
    * [kdb matrix](#kdb-matrix)
    * [kdb kmeans](#kdb-kmeans)
    * [kdb rarefy](#kdb-rarefy)
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

## What is `.kdb`?

The K-mer database (.kdb) is a file format and command-line utility (CLI) for accessing precomputed k-mers from sequencing data. .kdb views k-mer graph databases as a fundamental technology to look at biological sequence identities and abundances. K-mer methods have been noted as fast and are used in areas from NGS quality control to phylogenomic investigations.

.kdb is based on the block GNU-zip file (bgzf) standard. Each .kdb file has a header or metadata section, much like .bam files. It is essentially a tab-delimited format with the last column unstructured for k-mer specific metadata. Input files and total k-mer counts are stored in the metadata block at the top of the file

Please visit the [Install](#/install) page for details on installation.

See the commands in the [Usage](#/usage) section for an idea of what functionality is built in to kdb.

See the original blog post on the concept [here](https://matthewralston.github.io/blog/kmer-database-format-part-1).


The kdb project was designed to facilitate conversation between heavily optimized legacy codebases without much public attention, like Jellyfish, regarding the utility of standardizing k-mer frameworks. These frameworks are used throughout assembly and alignment hashing/seed-matching strategies. The primary goal of this project is documenting data shapes, compression strategies (which of course related to efficiency of storage, transmission, rapid access, etc.), and anticipating UI possibilities with the increases in read/write speeds afforded by improving SSD technologies and utilization of more channels of more rapid interfaces for data transmission (i.e. m2, NVMe, PCIx). 

## kdb is a file format

The k-mer database format is rather simple. It contains a metadata section, followed by a tab delimited format with an unstructured JSON as the final column. Both sections are combined in a compressed layer facilitated by Biopython's bio.bgzf module. 

Each file can be inspected with the view and header commands detailed in the [section below](#kdb-view).


# Install

[![PyPI version](https://img.shields.io/pypi/v/kdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kdb)
[![Coveralls code coverage](https://img.shields.io/coveralls/MatthewRalston/kdb/master.svg)][Coveralls]
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]

[pip]: https://pypi.org/project/kdb/
[Pythons]: https://pypi.org/project/kdb/
[Coveralls]: https://coveralls.io/r/MatthewRalston/kdb?branch=master
[RTD]: https://kdb.readthedocs.io/en/latest/

The current version on PyPI is shown above.

```bash
pip install kdb
```

See the [install](/installation) page for development install instructions.


# Usage


## How to view --help

Use '-h' to view detailed usage information about the subcommands

```bash
>./bin/kdb -h
usage: kdb [-h] {profile,header,view,matrix,rarefy,cluster,distance,index} ...

positional arguments:
  {profile,header,view,matrix,rarefy,cluster,distance,index}
                        Use --help with sub-commands
    profile             Parse data into the database from one or more sequence
                        files
    header              Print the YAML header of the .kdb file and exit
    view                View the contents of the .kdb file
    matrix              Generate a reduced-dimensionality matrix of the n *
                        4^k (sample x k-mer) data matrix.
    rarefy              Generate rarefaction information using
                        ecopy.diversity.rarefy for the supplied .kdb files
    cluster             Cluster the files according to their k-mer profile
    distance            Calculate various distance metrics between profiles
    index               Create a index file that can be held in memory

optional arguments:
  -h, --help            show this help message and exit
```

## kdb profile

A typical workflow first requires the generation of k-mer profiles. The following command will generate multiple profiles at `$K`-mer resolution simultaneously.

```bash
parallel './bin/kdb profile -k $K {} {.}.$K.kdb' ::: $(/bin/ls test/data/*.fasta.gz)
```

## kdb view

As mentioned before under [KDB format](#kdb-is-a-file-format), the kdb file consists of a header or metadata section, followed by data blocks until the end of the file. The header is YAML formatted and the data blocks are formatted as tab-separated value files (.tsv), with the last/right-most column being a JSON formatted metadata column. For developers, the YAML schema can be found in the config.py file.
```bash
# This should display the entire header of most files
>zcat test/data/foo.12.kdb | head -n 30 
# This will also display just the header
>./bin/kdb header test/data/foo.12.kdb
# The -H flag includes the header in the uncompressed output
>./bin/kdb view -H test/data/foo.12.kdb
version: 0.0.2
metadata_blocks: 1
k: 12
metadata: false
tags: []
files:
- filename: test/data/Cacetobutylicum_ATCC824.fasta.gz
  md5: 919357d5173cfa372e1e9a0f2b89c996
  mononucleotides:
	  A: 1427820
	  C: 637998
	  G: 640100
	  T: 1426962
  nullomers: 12880136
  sha256: b98be262a9904a3c2f84caa455679b7cebab7b2e9e15ca3105c69e001595abd6
  total_kmers: 8265716
  total_reads: 2
  unique_kmers: 3897080

0       4
1       1
2       0
```

## kdb distance

Suppose you want a distance matrix between profiles; this is made easy with the distance command

```bash
>./bin/kdb distance correlation test/data/*.$K.kdb
```

The result is a symmetric matrix in tsv format with column headers formed from the filenames minus their extensions. It is presumed that to properly analyze the distance matrix, you would name the files after their sample name or their species, or some other identifying features.



## kdb matrix

The kdb matrix command generates the count matrix either un-normalized, normalized (via ecopy), or with PCA or t-SNE dimensionality reduction applied. Note that default behavior of PCA if -n is not specified is to generate an elbow graph for the user to pick the appropriate choice of principal components for downstream analyses. The -n parameter is passed to the n_components parameter of sklearn.decomposition.PCA, which is commonly used for PCA in Python.

```bash
>./bin/kdb matrix -h
usage: kdb matrix [-h] [-v] [-k K] [-n N] [-d DELIMITER]
                  [--perplexity PERPLEXITY]
                  {PCA,tSNE,Normalized,Unnormalized} <.kdb> [<.kdb> ...]

positional arguments:
  {PCA,tSNE,Normalized,Unnormalized}
                        Choice of distance metric between two profiles
  <.kdb>                Two or more .kdb files

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Prints warnings to the console by default
  -k K                  The k-dimension that the files have in common
  -n N                  The number of dimensions to reduce with PCA or t-SNE.
                        DEFAULT: an elbow graph will be generated if -n is not
                        provided to help the user choose -n
  -d DELIMITER, --delimiter DELIMITER
                        The choice of delimiter to parse the DataFrame with
  --perplexity PERPLEXITY
                        The choice of the perplexity for t-SNE based
                        dimensionality reduction

>./bin/kdb matrix -n 3 PCA test/data/*.$K.kdb
```


## kdb kmeans

Now you would like to cluster your count matrix, perhaps after reducing dimensionality of the dataset. tSNE is a recommended choice for using the count matrix to differentiate between strains or species. PCA is equally useful in understanding the differences between species and reducing the dimensionality to something reasonable is an important first step before clustering.


```bash
# You may supply any arbitrary count tsv to cluster, and it will run that directly.
>./bin/kdb cluster -k 3 example.tsv

# Alternatively, you may pipe the result from kdb matrix directly to kdb cluster.
# The matrix command will not produce PCA reduced dimension matrix unless the -n parameter
# is identified from the elbow graph it will produce as a side effect.
>./bin/kdb matrix PCA -n 3 test/data/*.$K.kdb | ./bin/kdb cluster -k 3

# Alternatively, t-SNE may be used to project the data into 2 dimensions for visualization.
>./bin/kdb matrix tSNE -n 2 test/data/*.$K.kdb | ./bin/kdb cluster -k 3
```

## kdb rarefy

Suppose you want supply a normalized matrix directly to [ecopy](https://github.com/Auerilas/ecopy) for rarefaction analysis. The rarefaction command requires a tsv and may be piped directly from 'kdb matrix'.


```bash
# You may supply any arbitrary count tsv to rarefy, and it will run that directly.
>./bin/kdb rarefy example.tsv

# Alternatively, you may pipe the result from kdb matrix directly to kdb rarefy.
>./bin/kdb matrix Unnormalized test/data/*.$K.kdb | ./bin/kdb rarefy
```


# Documentation


Documentation for the (sub)module can be found here: [https://kdb.readthedocs.io/en/stable](https://kdb.readthedocs.io/en/stable)

Additionally, running the matrix, kmeans, or rarefy commands (which are arguably more complex in their usage) should be run with the DEBUG (-vv) verbosity setting on. This will yield additional information about the expected pipeline usage. That statement is echoed here.


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
parallel 'kdb profile -k $K {} {.}.$K.kdb' ::: $(/bin/ls test/data/*.fasta.gz)



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
# Cluster analysis
##################
# The first step ('kdb matrix') generates one from different profiles with the same choice of k.
# This command uses ecopy to normalize between sample k-mer total counts before PCA/tSNE.
# -n $N is the dimensionality of either PCA or tSNE. A good choice for tSNE is 2.
# If the command is run with PCA without a selected dimensionality, an elbow graph
# will be produced named '{0}'. Please use this graph to select
# the number of principal components to use.
# The pipeline will not continue until -n $N is selected by the user.
# It is not recommended to feed Unnormalized or Normalized matrices directly to 'kdb kmeans'
# 
# The PCA/tSNE matrix will be dimReduced ($N) * N, where N is the number of samples/files/profiles.
#
# And finally, a k-means clustering will be done on the reduced dimensionality dataset
# Note here that the -k $K is not related to the choice of substring length 'k' for profile generation.
# The 'kdb kmeans' command produces two figures, first is an elbow graph looking at up to N clusters.
# This elbow graph will be written to '{1}'.
# The second is the more typical scatterplot of the first two reduced dimensions
# and the k-means clustering labels shown over the scatter.
# This file will be written to '{2}'.
kdb matrix [-n $N] [ PCA | tSNE ] test/data/*.$K.kdb | kdb kmeans -k $K

#
# If you wanted to save a matrix from kdb matrix for use on your own
# it is recommended that you consider gzip compressing it if it is the Normalized or Unnormalized matrix
# which we will see is used downstream in the rarefaction analytical pathway.
#
##################
# Rarefaction
##################
# The Unnormalized and Normalized matrices go to ecopy's rarefy function to produce a rarefaction plot
# This plot will be written to '{3}'.
kdb matrix [ Unnormalized | Normalized ] test/data/*.$K.kdb | kdb rarefy
```


# Development


[![PyPI version](https://img.shields.io/pypi/v/kdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kdb)
[![Coveralls code coverage](https://img.shields.io/coveralls/MatthewRalston/kdb/master.svg)][Coveralls]
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]


[pip]: https://pypi.org/project/kdb/
[Pythons]: https://pypi.org/project/kdb/
[Coveralls]: https://coveralls.io/r/MatthewRalston/kdb?branch=master
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
Thanks to Max, Robin, and Robert for the halfway decent memories in St. Louis.
And thanks to my family and friends.
Go Blue Hens 2021.
