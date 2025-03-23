# README - kmerdb
> Python CLI and module for k-mer profiles, similarities, and graph files

NOTE: Beta-stage `.bgzf` and `zlib` compatible k-mer count vectors and DeBruijn graph edge-list formats.

## Development Status
![GitHub License](https://img.shields.io/github/license/matthewralston/kmerdb)
[![PyPI version](https://img.shields.io/pypi/v/kmerdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kmerdb.svg)][Pythons]
[![CircleCI](https://dl.circleci.com/status-badge/img/circleci/TxK1S2m7siJCSY89Dc6s4A/Dm3xDervECRDhDYKUkgUJN/tree/main.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/circleci/TxK1S2m7siJCSY89Dc6s4A/Dm3xDervECRDhDYKUkgUJN/tree/main)
[![codecov](https://codecov.io/gh/MatthewRalston/kmerdb/graph/badge.svg?token=8VB0RVRTSV)](https://codecov.io/gh/MatthewRalston/kmerdb)
[![ReadTheDocs status](https://readthedocs.org/projects/kmerdb/badge/?version=latest&style=flat)][RTD]
[![Downloads](https://static.pepy.tech/personalized-badge/kmerdb?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads)](https://pypi.org/project/kmerdb)
![PyPI - Downloads](https://img.shields.io/pypi/dm/kmerdb)
![GitHub Downloads (all assets, latest release)](https://img.shields.io/github/downloads/matthewralston/kmerdb/latest/total)
<!--
[![GitHub Downloads](https://img.shields.io/github/downloads/MatthewRalston/kmerdb/total.svg?style=social&logo=github&label=Download)](https://github.com/MatthewRalston/kmerdb/releases)
-->


[pip]: https://pypi.org/project/kmerdb/
[Pythons]: https://pypi.org/project/kmerdb/
[RTD]: https://kmerdb.readthedocs.io/en/latest/
[/]: https://kdb.info
[-]: ...



## Summary 

`kmerdb` is a Python CLI designed for k-mer counting and k-mer graph edge-lists. It addresses the ['k-mer' problem](https://en.wikipedia.org/wiki/K-mer) (substrings of length k) in a simple and performant manner. It stores the k-mer counts in a columnar format (input checksums, total and unique k-mer counts, nullomers, mononucleotide counts) with a YAML formatted metadata header in the first block of a `bgzf` formatted file. 

- [ x ] [Homepage:](https://matthewralston.github.io/kmerdb)
- [ x ] [Quick Start guide](https://matthewralston.github.io/kmerdb/quickstart)
- [ x ] `kmerdb usage subcommand_name`
  - `profile` - Make k-mer count vectors/profiles, calculate unique k-mer counts, total k-mer counts, nullomer counts. Import to read/write NumPy arrays from profile object attributes.
  - `graph` - Make a weighted edge list of kmer-to-kmer relationships, akin to a De Bruijn graph.
  - `minimizers` - Generate minimizers in a plain-text format from input sequences.
  - `alignment` - Generate Smith-Waterman alignment (not yet functional)
  - `usage` - Display verbose input file/parameter and algorithm details of subcommands.
  - `help` - Display verbose input file/parameter and algorithm details of subcommands.
  - `view` - View .tsv count/frequency vectors with/without preamble.
  - `header` - View YAML formatted header and aggregate counts
  - `matrix` - Collate multiple profiles into a count matrix for dimensionality reduction, etc.
  - `kmeans` - k-means clustering on a distance matrix via Scikit-learn or BioPython with kcluster distances
  - `hierarchical` - hierarchical clustering on a distance matrix via BioPython with linkage choices
  - `distance` - Distance matrices (from kmer count matrices) including SciPy distances, a Pearson correlation coefficient implemented in Cython, and Spearman rank correlation included as additional distances.
  - `index` - Create an index file for the kmer profile (Delayed:)
  - `shuf` - Shuffle a k-mer count vector/profile (Delayed:)
  - `version` - Display kmerdb version number
  - `citation` - Silence citation suggestion
- [ x ] `kmerdb subcommand -h|--help`


k-mer counts from .fa(.gz)/.fq(.gz) sequence data can be computed and stored for access to metadata and count aggregation faculties. For those familiar with `.bam`, a `view` and `header` functions are provided. This file is compatible with `zlib`.

Install with `pip install kmerdb`



Please see the [Quickstart guide](https://matthewralston.github.io/kmerdb/quickstart) for more information about the format, the library, and the project.


## Usage

```bash
# Usage    --help option    --debug mode
kmerdb --help # [+ --debug mode]
kmerdb usage profile


#   +

# [ 3 main features: ]     [ 1.   -    k-mer counts  ]

# Create a [composite] profile of k-mer counts from sequence files. (.fasta|.fastq|.fa.gz|.fq.gz)
kmerdb profile -vv -k 8 --output-name sample_1 sample_1_rep1.fq.gz [sample_1_rep2.fq.gz]
# Creates k-mer count vector/profile in sample_1.8.kdb. This is the input to other steps, including count matrix aggregation. --minK and --maxK options can be specified to create multiple k-mer profiles at once.
<!-- # Alternatively, can also take a plain-text samplesheet.txt with one filepath on each line. -->

#          De Bruijn graphs (not a main feature yet, delayed)
# Build a weighted edge list (+ node ids/counts = De Bruijn graph)
kmerdb graph -vv -k 12 example_1.fq.gz example_2.fq.gz edges_1.kdbg

# View k-mer count vector
kmerdb view profile_1.8.kdb # -H for full header

# Note: zlib compatibility
#zcat profile_1.8.kdb

# View header (config.py[kdb_metadata_schema#L84])
kmerdb header profile_1.8.kdb

## [ 3 main features: ]   [ 2. Optional normalization, PCA/tSNE, and distance metrics ]

# K-mer count matrix - Cython Pearson coefficient of correlation [ ssxy/sqrt(ssxx*ssyy) ]
kmerdb matrix from *.8.kdb | kmerdb distance pearson STDIN
# 
# kmerdb matrix -vv DESeq2 *.8.kdb
# kmerdb matrix -vv PCA *.8.kdb
# kmerdb matrix -vv tSNE *.8.kdb
#   # <from> just makes a k-mer count matrix from k-mer count vectors.
# 

# Distances on count matrices [ SciPy ]  pdists + [ Cython ] Pearson correlation, scipy Spearman and scipy correlation pdist calculations are available ]
kmerdb distance -h
# 
#
# usage: kmerdb distance [-h] [-v] [--debug] [-l LOG_FILE] [--output-delimiter OUTPUT_DELIMITER] [-p PARALLEL] [--column-names COLUMN_NAMES] [--delimiter DELIMITER] [-k K]
#                       {braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,matching,minkowski,pearson,rogerstanimotorusselrao,seuclidean,sokalmichener,sokalsneath,spearman,sqeuclidean,yule} [<kdbfile1 kdbfile2 ...|input.tsv|STDIN> ...]

# [ 3 main features: ]      [ 3. Clustering: k-means and hierarchical with matplotlib ]

#    Kmeans (sklearn, BioPython)
kmerdb kmeans -vv -k 4 -i dist.tsv
#    BioPython Phylip tree + upgma
kmerdb hierarchical -vv -i dist.tsv


```




## Usage example

![kmerdb.gif](kmerdb.gif)




## Installation

### OSX and Linux release:

```sh
pip install kmerdb
```

#### Optional DESeq2 normalization

DESeq2 is an optional R dependency for rpy2-mediated normalization. Make sure development libraries are installed from the repository.


```
pip install -r requirements-dev.txt
```
Next, install DESeq2 via bioconductor.
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```


## IUPAC support:

```bash
kmerdb profile -k $k -o output input.fa # This simply discards non-IUPAC characters.
```
IUPAC residues (ATCG+RYSWKM+BDHV) are kept throughout the k-mer counting. But non-IUPAC residues (N) and characters are trimmed from the sequences prior to k-mer counting. Non-standard IUPAC residues are counted as doublets or triplets.


## Documentation

Check out the [main webpage](https://matthewralston.github.io/kmerdb) and the [Readthedocs documentation](https://kdb.readthedocs.io/en/stable/), with examples and descriptions of the module usage.

Important features to usage that may be important may not be fully documented as the project is in beta.

For example, the IUPAC treatment is largely custom, and does the sensible thing when ambiguous bases are found in fasta files, but it could use some polishing. For example, the '`N`' residue rejection creates gaps in the k-mer profile from the real dataset by admittedly ommitting certain k-mer counts.
This is one method for counting k-mers and handling ambiguity. Fork it and play with it a bit.

Also, the parallel handling may not always be smooth, if you're trying to load dozens of 12+ mer profiles into memory. This would especially matter in the matrix command, before the matrix is generated. You can use single-core if your machine can't collate that much into main memory at once, depending on how deep the fastq dataset is. Even when handling small-ish k-mer profiles, you may bump into memory overheads rather quickly. 

Besides that, I'd suggest reading the source, the differente elements of the [main page](https://matthewralston.github.io/kmerdb) or the [RTD documentation](https://kdb.readthedocs.io/en/stable/).



## Minimum-viable product Documentation

### Problem statement: 

Calculate relevant metadata from k-mer profiles, manipulate count matrices and distance matrices on inter-species inter-profile distances, perform exploartory analysis using tools such as UPGMA hierarchical clustering, k-means clustering, PCA, and others.

Currently unsupported but targeted features include regression modeling, strassen multiplication, NumPy CPU multipication option, least squares and associated regression code.

Currently in-progress is the application note. I wouldn't mind ironing out a few more components of the program for the application note.

### Target audience:

Biologists and researchers using metagenomics and genomics tools for the assessment of inter-genomic distances, k-mer counts, nullomer counts, unique counts, and other metadata.

### Feature prioritization:

Features currently implemented in `kmerdb` include:


- [x] k-mer counting and profile generation
- [x] multiplexed profiles
- [x] De Bruijn graph structures
- [x] Count aggergation and normalization (DESeq2)
- [x] PCA
- [x] t-Stochastic Neighbor Embedding
- [x] Cython Pearson correlation coefficient
- [x] Other correlation coefficients (Spearman, Pearson via scipy)
- [x] hierarchical clustering
- [x] k-means clustering
- [x] data pipelining (Unix pipelines)

Features currently under consideration include: 

- [ ] regression modeling for metagenomnic populations
- [ ] alternative k-mer count vector distances
- [ ] alignment using k-mers as seed regions
- [ ] De Bruijn graph traversals and contig generation using DBA (De Bruijn graph) assembly

### Success Criteria:




### Technical specifications:

### User stories:

### Development timeline:

- 2/24 - 6/24 v0.7.7+ - appmap, De Bruijn parser, other updates

- 10/23 - 

### Testing and quality assurance plan:

### Feedback collection strategy:




## Development

https://matthewralston.github.io/kmerdb/developing

```bash
python setup.py test
```

## License

Created by Matthew Ralston - [Scientist, Programmer, Musician](http://matthewralston.github.io) - [Email](mailto:mralston.development@gmail.com)

Distributed under the Apache license. See `LICENSE.txt` for the copy distributed with this project. Open source software is not for everyone, and im the author and maintainer. cheers, on me. You may use and distribute this software, gratis, so long as the original LICENSE.txt is distributed along with the software. This software is distributed AS IS and provides no warranties of any kind.

```
   Copyright 2020 Matthew Ralston

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
```

## Contributing

1. Fork it (<https://github.com/MatthewRalston/kmerdb/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

## Acknowledgements

Thanks mom and dad and my hometown, awesome hs, and the University of Delaware faculty for support and encouragement. Thanks to my former mentors, bosses, and coworkers. It's been really enjoyable anticipating what the metagenomics community might want from a tool that can handle microbial k-mers well.

Thank you to the authors of kPAL and Jellyfish for the inspiration and bit shifting trick. And thank you to others for the encouragement along the way, who shall remain nameless. 

The intention is that more developers would want to add functionality to the codebase or even just utilize things downstream, but to build out directly with numpy and scipy/scikit as needed to suggest the basic infrastructure for the ML problems and modeling approaches that could be applied to such datasets. This project began under GPL v3.0 and was relicensed with Apache v2. Hopefully this project could gain some interest. I have so much fun working on this project. There's more to it than meets the eye. I'm working on a preprint, and the draft is included in some of the latest versions of the codebase, specifically .Rmd files.

More on the flip-side. It's so complex with technology these days...

<!--
Thanks of course to my fans (and haters). Yeah i see you.... but i dont.
Thanks to my former mentors BC, MR, IN, CR, and my newer bosses PJ and KL.
Thanks to the Pap lab for the first dataset that I continue to use.
Thank you to Ryan for the food and stuff. I actually made this project specifically so you and I could converse...
Thanks to Blahah for tolerating someone snooping and imitating his Ruby style.
Thanks to Erin for getting my feet wet in this new field. You are my mvp.
Thanks to Rachel for the good memories and friendship. And Sophie too. veggies n' R love.
Thanks to Yasmeen for the usual banter.
Thanks to A for the newer banter.
Thanks to Max, Robin, and Robert for the good memories in St. Louis. What's new?
Thanks to Fred for the good memories. Hope you're on soon.
Thanks to Nichole for the cookies and good memories. And your cute furballs too! Hope you're well
Thanks to S for the lessons, convos, and even embarassing moments. You're kind of awesome to me.
Thanks to a few friends I met in 2023 that reminded me I have a lot to learn about friendship, dating, and street smarts.
Thanks to them even more now that I got it xd up.

Thanks to the people of NCC for the Doordash money. It might not be much but I don't have it twisted (I do.)



Thanks to D from BCCS.
Thanks to C4H&H. I'm 'healthier' now, but I really still think I need more support than just BCCS. it's urgent.
Thanks to CT and family. Your love and support means the world to me.
Thanks to AS and family. Your support made a difference. Praying for better employment and opportunities.

And thanks to my family and friends.
Go Blue Hens
-->
