# README - kmerdb
> A Python CLI and module for k-mer profiles, similarities, and graph databases

NOTE: This project is pre-alpha, all of the badge links are broken and are just placeholders at the moment. Development is ongoing. But feel free to clone the repository and play with the code for yourself!

## Development Status

[![PyPI version](https://img.shields.io/pypi/v/kmerdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kmerdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kmerdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kmerdb)
[![Coveralls code coverage](https://coveralls.io/repos/github/MatthewRalston/kmerdb/badge.svg?branch=master)](https://coveralls.io/github/MatthewRalston/kmerdb?branch=master)
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]


[pip]: https://pypi.org/project/kmerdb/
[Pythons]: https://pypi.org/project/kmerdb/
[RTD]: https://kdb.readthedocs.io/en/latest/

## Summary 

KDB is a Python library designed for bioinformatics applications. It addresses the ['k-mer' problem](https://en.wikipedia.org/wiki/K-mer) (substrings of length k) in a simple and performant manner. It stores the k-mer counts/abundances and total counts. A per-kmer metadata feature is planned for the future. The k-mer spectrum of the fasta or fastq sequencing data is stored in the `.kdb` format spec, a bgzf file similar to `.bam`.

The principle goal of the library is k-mer statistics and rapid access to specific k-mers and associated abundances with a Python CLI and API. Other goals include access to the k-mer count distribution, k-mer transition probabilities, and more by leveraging the bgzf specification. Another low-hanging fruit could be approximating demultiplexing coefficients for artificial metagenomes.


## Installation

OS X and Linux release:

```sh
pip install kmerdb
```

Development installation:

```sh
git clone https://github.com/MatthewRalston/kmerdb.git
pip install -e .
```

## Usage Example

Usage in detail can be found on the [quickstart page](https://matthewralston.github.io/kdb/quickstart.html#usage)

CLI Usage

```bash
kmerdb --help
kmerdb summary --help
# Build a [composite] profile to a new or existing .kdb file
kmerdb profile example1.fq.gz example2.fq.gz profile.kdb
# Calculate similarity between two (or more) profiles
kmerdb distance correlation profile1.kdb profile2.kdb (...)
```


## Documentation

Check out the [main webpage](https://matthewralston.github.io/kmerdb) and the [Readthedocs documentation](https://kdb.readthedocs.io/en/stable/), with examples and descriptions of the module usage.

## Development

```bash
python setup.py test
```

## License

Created by Matthew Ralston - [Scientist, Programmer, Musician](http://matthewralston.github.io) - [Email](mailto:mrals89@gmail.com)

Distributed under the Apache license. See `LICENSE.txt` for the copy distributed with this project. Open source software is not for everyone, but for those of us starting out and trying to put the ecosystem ahead of ego, we march into the information age with this ethos.

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

1. Fork it (<https://github.com/MatthewRalston/kdb/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

## Acknowledgements

Thank you to the authors of kPAL and Jellyfish for the early inspiration. And thank you to others for the encouragement along the way, who shall remain nameless. I wanted this library to be a good strategy for assessing these k-mer profiles, in a way that is both cost aware of the analytical tasks at play, capable of storing the exact profiles in sync with the current assemblies, and then updating the kmer databases only when needed to generate enough spectral signature information.

The intention is that more developers would want to add functionality to the codebase or even just utilize things downstream, but to build out directly with numpy and scipy/scikit as needed to suggest the basic infrastructure for the ML problems and modeling approaches that could be applied to such datasets. This project has begun under GPL v3.0 and hopefully could gain some interest.

More on the flip-side of this file. Literally. And figuratively. It's so complex with technology these days.

Also thank you to patelvivek (github/viensio) for the Github forking ribbon on the kdb website.

<!--
Thanks to my former mentors BC, MR, IN, CR, and my newer bosses PJ and KL.
Thanks to the Pap lab for the first dataset that I continue to use.
Thank you to Ryan for the food and stuff.
Thanks to Blahah for tolerating someone snooping and imitating his Ruby style.
Thanks to Erin for getting my feet wet in this new field.
Thanks to Rachel for the good memories and friendship.
Thanks to Yasmeen for the usual banter.
Thanks to Max, Robin, and Robert for the halfway decent memories in St. Louis.
And thanks to my family and friends.
Go Blue Hens 2021.
-->
