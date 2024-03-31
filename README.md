# README - kmerdb
> A Python CLI and module for k-mer profiles, similarities, and graph databases

NOTE: This project is in beta stage. Development is ongoing. But feel free to clone the repository and play with the code for yourself.

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

KDB is a Python library designed for bioinformatics applications. It addresses the ['k-mer' problem](https://en.wikipedia.org/wiki/K-mer) (substrings of length k) in a simple and performant manner. It stores the k-mer counts/abundances in a columnar format, with input file checksums, total counts, nullomers, and mononucleotide counts in a YAML formatted header in the first block of the `bgzf` formatted `.kdb` file. One restriction is that k-mers with unspecified sequence residues 'N' create gaps in the k-mer to sequence relationship space, and are excluded. That said, non-standard IUPAC residues are supported.


Please see the [Quickstart guide](https://matthewralston.github.io/kmerdb/quickstart) for more information about the format, the library, and the project.

The k-mer spectrum of the fasta or fastq sequencing data is stored in the `.kdb` format spec, a bgzf file similar to `.bam`. For those familiar with `.bam`, a `view` and `header` functions are provided to decompress a `.kdb` file into a standard output stream. The output file is compatible with `zlib`.



## Installation


### Dependencies

DESeq2 is required as a R dependency for rpy2-mediated normalization.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

All other dependencies are managed directly by pip. 


### OSX and Linux release:

```sh
pip install --python-version 3.7.4 --pre kmerdb
```



### Development installation:

```sh
git clone https://github.com/MatthewRalston/kmerdb.git
pip install .
```

## Usage Example

NOTE: Usage in detail can be found on the [quickstart page](https://matthewralston.github.io/kmerdb/quickstart#usage)

<!-- ## NOTE: Temporary usage pattern:
Migrating the repo from setup.py to the PEP606 standard PyProject.toml is borking my current invocation pattern. Sorry for the inconveniece... it's happening right now. -->

<!-- ```bash
python -m kmerdb [cmd] [options]
```

See `python -m kmerdb -h` for details.
-->

## CLI Usage

```bash
kmerdb --help
# Build a [composite] profile to a new .kdb file
kmerdb profile -k 8 example1.fq.gz example2.fq.gz profile.8.kdb

# Note: zlib compatibility
zcat profile.8.kdb

# Build a weighted edge list
kmerdb graph -k 12 example1.fq.gz example2.fq.gz edges.kdbg

# View the raw data
kmerdb view profile.8.kdb # -H for full header

# View the header
kmerdb header profile.8.kdb

# Collate the files. See 'kmerdb matrix -h' for more information.
# Note: the 'pass' subcommand passes the int counts through collation, without normalization.
# In this case the shell interprets '*.8.kdb' as all 8-mer profiles in the current working directory.
# The k-mer profiles are read in parallel (-p $cores), and collated into one Pandas dataframe, which is printed to STDOUT.
# Other options include DESeq2 normalization, frequency matrix, or PCA|tSNE based dimensionality reduction techniques.
kmerdb matrix -p $cores pass *.8.kdb > kmer_count_dataframe.tsv

# Calculate similarity between two (or more) profiles
# The correlation distance from Numpy is used on one or more profiles, or piped output from 'kmerdb matrix'.
kmerdb distance correlation profile1.kdb profile2.kdb (...) > distance.tsv

# A condensed, one-line invocation of the matrix and distance command using the bash shell's pipe mechanism is as follows.
kmerdb matrix pass *.8.kdb | kmerdb distance correlation STDIN > distance.tsv
```

## IUPAC support:

```bash
kmerdb profile -k $k input.fa output.kdb # This may discard non-IUPAC characters, this feature lacks documentation!
```
IUPAC residues (ATCG+RYSWKM+BDHV) are kept throughout the k-mer counting. But non-IUPAC residues (N) and characters are trimmed from the sequences prior to k-mer counting.



## Documentation

Check out the [main webpage](https://matthewralston.github.io/kmerdb) and the [Readthedocs documentation](https://kdb.readthedocs.io/en/stable/), with examples and descriptions of the module usage.

Important features to usage that may be important may not be fully documented.

For example, the IUPAC treatment is largely custom, and does the sensible thing when ambiguous bases are found in fasta files, but it could use some polishing.

In addition, the '`N`' residue rejection creates gaps in the k-mer profile from the real dataset by admittedly ommitting certain k-mer counts.
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

Distributed under the Apache license. See `LICENSE.txt` for the copy distributed with this project. Open source software is not for everyone, but we march into the information age with this ethos. I have the patent rights to this software. You may use and distribute this software, gratis, so long as the original LICENSE.txt is distributed along with the software. This software is distributed AS IS and provides no warranties of any kind.

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

Thank you to the authors of kPAL and Jellyfish for the early inspiration. And thank you to others for the encouragement along the way, who shall remain nameless. I wanted this library to be a good strategy for assessing these k-mer profiles, in a way that is both cost aware of the analytical tasks at play, capable of storing the exact profiles in sync with the current assemblies, and then updating the kmer databases only when needed to generate enough spectral signature information.

The intention is that more developers would want to add functionality to the codebase or even just utilize things downstream, but to build out directly with numpy and scipy/scikit as needed to suggest the basic infrastructure for the ML problems and modeling approaches that could be applied to such datasets. This project began under GPL v3.0 and was relicensed with Apache v2. Hopefully this project could gain some interest. I have so much fun working on just this one project. There's more to it than meets the eye. I'm working on a preprint, and the draft is included in some of the latest versions of the codebase, specifically .Rmd files.

More on the flip-side of this file. Literally. And figuratively. It's so complex with technology these days.

<!--
Thanks of course to that French girl from the children's series. 
Thanks to my former mentors BC, MR, IN, CR, and my newer bosses PJ and KL.
Thanks to the Pap lab for the first dataset that I continue to use.
Thank you to Ryan for the food and stuff. I actually made this project specifically so you and I could converse...
Thanks to Blahah for tolerating someone snooping and imitating his Ruby style.
Thanks to Erin for getting my feet wet in this new field. You are my mvp.
Thanks to Rachel for the good memories and friendship. And Sophie too. veggies n' R love.
Thanks to Yasmeen for the usual banter.
Thanks to Max, Robin, and Robert for the good memories in St. Louis.
Thanks to Freddy Miller for the good memories.
Thanks to Nichole for the cookies and good memories. And your cute furballs too!
Thanks to Stace for the lessons, convos, and even embarassing moments. You're kind of awesome to me.
Thanks to a few friends I met in 2023 that reminded me I have a lot to learn about friendship, dating, and street smarts.
And thanks to my family and friends.
Go Blue Hens
-->
