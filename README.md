# README - kmerdb
> Python CLI and module for k-mer profiles, similarities, and graph databases

NOTE: Beta-stage `.bgzf` and `zlib` compatible k-mer count and DeBruijn graph formats.

## Development Status
[![Downloads](https://static.pepy.tech/personalized-badge/kmerdb?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads)](https://pypi.org/project/kmerdb)
![PyPI - Downloads](https://img.shields.io/pypi/dm/kmerdb)
[![GitHub Downloads](https://img.shields.io/github/downloads/MatthewRalston/kdb/total.svg?style=social&logo=github&label=Download)](https://github.com/MatthewRalston/kmerdb/releases)
[![PyPI version](https://img.shields.io/pypi/v/kmerdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kmerdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kmerdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kmerdb)
[![Coveralls code coverage](https://coveralls.io/repos/github/MatthewRalston/kmerdb/badge.svg?branch=master)](https://coveralls.io/github/MatthewRalston/kmerdb?branch=master)
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]


[pip]: https://pypi.org/project/kmerdb/
[Pythons]: https://pypi.org/project/kmerdb/
[RTD]: https://kdb.readthedocs.io/en/latest/

## Summary 

k-mer counts or De Bruijn graph of .fa(.gz)/.fq(.gz) sequence data can be converted to `.kdb` (+new! `.kdbg` De Bruijn graph) format, a bgzf file similar to `.bam`. For those familiar with `.bam`, a `view` and `header` functions are provided. The output file is compatible with `zlib`.

`pip install kmerdb` is a Python CLI designed for k-mer counting and De Bruijn graphs. It addresses the ['k-mer' problem](https://en.wikipedia.org/wiki/K-mer) (substrings of length k) in a simple and performant manner. It stores the k-mer counts in a columnar format (input checksums, total and unique k-mer counts, nullomers, mononucleotide counts) with a YAML formatted metadata header in the first block of the `bgzf` formatted `.kdb` (k-mer database) file. 

Please see the [Quickstart guide](https://matthewralston.github.io/kmerdb/quickstart) for more information about the format, the library, and the project.


## Usage

```bash
# Usage    --help option    --debug mode
kmerdb --help # [+ --debug mode]
kmerdb usage -m graph

# Output:
 o-O      |||
o---O     |||             [|[          kmerdb           ]|]
O---o     |||
 O-o      |||        version :     v0.8.0
  O       |||
 o-O      |||        GitHub  : https://github.com/MatthewRalston/kmerdb/issues
o---O     |||         PyPI   : https://pypi.org/project/kmerdb/
O---o     |||      Website   : https://matthewralston.github.io/kmerdb
 O-o      |||
                                                                       lang :         python
                                                                          v :      >= v3.7.4

                      package manger : pip
                        version      : >= 24.0
        package root : /path2py/lib/python3.12/site-packages/kmerdb
        exe file     : /path2py/lib/python3.12/site-packages/kmerdb/__init__.py

                      required packages : 8
                   development packages : 14

           ARGV : ['/path2py/bin/kmerdb', 'usage', '-m', 'graph']
        
...

Beginning program...

                          name : graph
                   description : create edge nodes in (block) .gz compressed format from .fasta or .fq format.
        

* features

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




* steps

- name: read input file(s) from filesystem into k-mer arrays
  shortname: shred inputs into k-mer count arrays
  description: shred input sequences into k-mer count vector
  
- name: merge k-mer arrays and aggregate metadata
  shortname: merge k-mer count arrays for aggregate metadata (header)
  description: merge counts of nullomers, unique kmers, and total kmers.
  
- name: collate the weighted edge lists after reading multiple files. Output data
    consists of a edge_list, analogous metadata-header as YAML, kmer_counts, and nullomer_ids.
  shortname: extract undirected weighted graph
  description: consists of info from a .kdb file and a .kdbg file. The node IDs, the
    edges, and the number of times the pair was observed from forward sequences in
    the provided dataset
	
- name: print 'table' Final stats and close output file
  shortname: metrics and shutdown
  description: print final statistics, typically metadata values, and ensure file
    is closed.




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
kmerdb matrix pass *.8.kdb | kmerdb distance pearson -i STDIN
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


## Installation

### OSX and Linux release:

```sh
pip install --python-version 3.7.4 --pre kmerdb
```

#### Optional DESeq2 normalization

DESeq2 is an optional R dependency for rpy2-mediated normalization.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```


## IUPAC support:

```bash
kmerdb profile -k $k input.fa output.kdb # This may discard non-IUPAC characters, this feature lacks documentation!
```
IUPAC residues (ATCG+RYSWKM+BDHV) are kept throughout the k-mer counting. But non-IUPAC residues (N) and characters are trimmed from the sequences prior to k-mer counting. Non-standard IUPAC residues are counted as doublets or triplets.


## Documentation

Check out the [main webpage](https://matthewralston.github.io/kmerdb) and the [Readthedocs documentation](https://kdb.readthedocs.io/en/stable/), with examples and descriptions of the module usage.

Important features to usage that may be important may not be fully documented as the project is in beta.

For example, the IUPAC treatment is largely custom, and does the sensible thing when ambiguous bases are found in fasta files, but it could use some polishing. For example, the '`N`' residue rejection creates gaps in the k-mer profile from the real dataset by admittedly ommitting certain k-mer counts.
This is one method for counting k-mers and handling ambiguity. Fork it and play with it a bit.

Also, the parallel handling may not always be smooth, if you're trying to load dozens of 12+ mer profiles into memory. This would especially matter in the matrix command, before the matrix is generated. You can use single-core if your machine can't collate that much into main memory at once, depending on how deep the fastq dataset is, and the `--block-size` parameter in `kmerdb profile` is likely going to facilitate your memory overhead by reading chunks of `--block-size` reads into memory at once while accumulating the k-mer counts in a `uint64` array. Even when handling small-ish k-mer profiles, you may bump into memory overheads rather quickly. 

Besides that, I'd suggest reading the source, the differente elements of the [main page](https://matthewralston.github.io/kmerdb) or the [RTD documentation](https://kdb.readthedocs.io/en/stable/).




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
Thanks to Fred for the good memories.
Thanks to Nichole for the cookies and good memories. And your cute furballs too! Hope you're well
Thanks to S for the lessons, convos, and even embarassing moments. You're kind of awesome to me.
Thanks to a few friends I met in 2023 that reminded me I have a lot to learn about friendship, dating, and street smarts.
Thanks to them even more now that I got it xd up.

Thanks to the people of NCC for the Doordash money. It might not be much but I don't have it twisted (I do.)



Thanks to D from BCCS.
Thanks to C4H&H. I'm 'healthier' now, but I really still think I need more support than just BCCS. it's urgent.

And thanks to my family and friends.
Go Blue Hens
-->
