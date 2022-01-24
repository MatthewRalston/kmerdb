---
title: Developing
layout: post
description: 'Everything you need to get started developing with kmerdb'
image: assets/images/Python-transparent-greyscale.png
nav-menu: true
toc: true
---

*[Introduction](#introduction)
*[Install](#install)
*[Usage](#usage)
  * [fileutil](#fileutil)
  * [read file](#read)
	* [One file](#onefile)
	* [Multiple files in parallel](#multiple)
  * [distance](#distance)
  * [PCA](#PCA)
  * [t-SNE](#tSNE)
  * [UMAP](#umap)





# Introduction

## What is `.kdb` or what are k-mers?

Please refer to the [quickstart](quickstart) for basic guide on the command-line usage.

## Developing

This guide is meant to encourage other developers to leverage the capabilities of this `.kdb` bgzf spec.


The key file to watch is `config.py`, which keeps the version number as metadata for the metadata. This project makes no warantees about the stability of the format specification. It is unclear if the project goal is to store tabular count matrices, or to develop ingestion for graph databases stored as CYPHER text. At this point I can't decide on which features to develop since there is so little attention to this project. For this reason, metadata formats are destined to change. Perhaps someone would like to contribute to the jsonschema definitions?

This guide shows you how to do basic operations with the `kmerdb` submodules.

# Install

```bash
pip install kmerdb
```

```bash
git clone 
git clone git@github.com:MatthewRalston/kmerdb.git
cd kdb
# I know the following may be prone to fail, but I need to know if this is an issue, I always install from setup.py.
pip install --no-cache-dir . # python setup.py install 
```


# Usage

After installation is complete, you can import the module directly or focus on specific submodules you'd like to import.

## fileutil

The main module you will be using to interact with `.kdb` format files is `kmerdb.fileutil`, which provides the open method, and the `KDBReader`/`KDBWriter` classes.

We will use this in the section labelled [read](#read) to discuss the lazy loading behaviors of `KDBReader` objects.

## read

First we will import the `fileutil` module to access individual columns from bgzf compressed `.kdb` format files.

```python
from kmerdb import fileutil
```

### onefile

Here we can view the contents of one kmerdb file in Python as NumPy arrays. Default settings as of 0.6.5 are `uint64` arrays for counts and indices, and `float64` for frequencies.

```python
from kmerdb import fileutil
with fileutil.open("example.kdb", 'r', slurp=False) as kdb:
	# Check the metadata
	assert kdb.metadata["k"] == k, "Assertion failed to verify chosen k"
    kdb.slurp() # Actually read the data from the disk as needed.
	print(kdb.profile)
	print(kdb.kmer_ids)
	print(kdb.counts)
	print(kdb.frequencies)
```




### multiple

```python
from multiprocessing import pool
from kmerdb import fileutil

# Do some sanity checks
files = list(map(lambda f: fileutil.open(f, 'r', slurp=False), myFiles))
assert all(kdbrdr.k == k for kdbrdr in files), "Couldn't validate a uniform choice of k"

# Read in parallel
file_reader = fileutil.FileReader()
if cores > 1: # read in parallel
    with Pool(processses=cores) as pool:
    	files = pool.map(file_reader.load_file, myFiles)
else:
    files = list(map(lambda f: fileutil.open(f, 'r', slurp=True), myFiles))

data = np.array([kdbrdr.slurp() for kdbrdr in files], dtype=dtype)
```


## distance

## PCA

## tSNE

## umap

