---
category: 'What is .kdb?'
title: 'Goal of the project'
path: '/goal'

layout: nil
---

## How did this project start?

This project arose independently from kmer.js as I needed a Python or Ruby based version of a basic k-mer counting program. I decided to go with Python as I needed access to matrix functions on integers or shorts. 

KDB is a Python library designed for bioinformatics applications. It addresses the ['k-mer' problem](https://en.wikipedia.org/wiki/K-mer) (substrings of length k) in a simple and performant manner. It stores the k-mer counts/abundances and total counts. A per-kmer metadata feature is planned for the future. The k-mer spectrum of the fasta or fastq sequencing data is stored in the `.kdb` format spec, a bgzf file similar to `.bam`.

The principle goal of the library is k-mer statistics and rapid access to specific k-mers and associated abundances with a Python CLI and API. Other goals include access to the k-mer count distribution, k-mer transition probabilities, and more by leveraging the bgzf specification. Another low-hanging fruit could be approximating demultiplexing coefficients for artificial metagenomes.

## What is the goal of this project

The goal of this project is to explore the concept of k-mer signatures of sequences as having both chroma and abundance/intensity. While alignment remains the primary method of assigning sequence identity, k-mer methods and elaborate pre-hashing may solve another problem we currently face: memory. Although bam files do retain all information important for ensuring sequencing quality, each dense record of 4 lines in the file is used only to increment count, somewhere down the pipeline. Storing subprofiles/abundances directly may be a more efficient method of storing both alignment identities as well as abundance in deep sequencing datasets.


## What do I need to know about your programming style?

You need to understand that everything revolves around my logging statements. Default log level is warning and error, so you might not see what the program is doing until you activate logging statements with '-v' or the double v '-vv' for info and debug logging, respectively.

You also need to know that some of the newer commands 'distance' and 'matrix' in particular have a mixture of positional multiplexed, positional single channel, and implicitly required keyword arguments. To explain a little further, an argument that is multiplexed has more than one argument supplied in a single space at the end of the command, normally where the positional arguments sit. It's best to have multiple positional singles and to have few if any multiplexed channels, but it is possible in Python's argparse to have a singleton positional argument (not a multiplexed channel) followed by a multiplexed channel. Specifically, the following shows a positional argument 'PCA' followed by a list of files.

```bash
>./bin/kdb matrix [[ -n $N ]] PCA test1.kdb test2.kdb test3.kdb
# or the last 3 arguments can be contracted in the shell with a *.kdb.
```

The last term I defined earlier was implicitly required keyword arguments. Some keyword arguments may be required, like '-n' in the case of reduced dimensionality matrices from the 'kdb matrix' command. Also, -k is a required argument for 'kdb cluster', which (currently) uses kmeans to cluster the data.


