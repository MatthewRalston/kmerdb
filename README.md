# README - KDB
> A Python CLI and module for k-mer profiles, similarities, and graph databases

NOTE: This project is pre-alpha, all of the badge links are broken and are just placeholders at the moment. Development is ongoing. But feel free to clone the repository and play with the code for yourself!

## Development Status

[![PyPI version](https://img.shields.io/pypi/v/kdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kdb.svg)][Pythons]
[![Travis build status](https://travis-ci.com/MatthewRalston/kdb.svg?branch=master)][TravisCI]
[![Coveralls code coverage](https://img.shields.io/coveralls/MatthewRalston/kdb/master.svg)][Coveralls]
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]


[pip]: https://pypi.org/project/kdb/
[Pythons]: https://pypi.org/project/kdb/
[TravisCI]: https://travis-ci.com/MatthewRalston/kdb
[Coveralls]: https://coveralls.io/r/MatthewRalston/kdb?branch=master
[RTD]: https://kdb.readthedocs.io/en/stable/

## Summary 

KDB is a Python library designed for bioinformatics applications. It addresses the ['k-mer' problem](https://en.wikipedia.org/wiki/K-mer) (substrings of length k) in a simple and performant manner. It generates a [De Brujin graph](https://en.wikipedia.org/wiki/De_Bruijn_graph) from the k-mer spectrum of fasta or fastq sequencing data and stores the graph and spectrum to the `.kdb` format spec, a bgzf file similar to BAM. 

It could include utilities for exporting the De Brujin graph to graph databases. It could include basic operations about graph properties, even exporting stats about node degree and high-level properties of graph structure could be useful metrics in describing sequence space complexity.

The reason for even including those metrics this early in the project before user interests are investigated is that it would facilitate the conversations about the sequence spaces under consideration. k-mer partitioning and phylogenetic considerations would be high-level abstractions of lower level efforts with sufficient optimization needs, to anticipate physical storage limitations and in memory efficiencies possible to retrieve those statistics, and would shed light on the types of queries common to highly pure and highly multiplexed sample states under consideration during different phases of data cleaning (read cleaning, trimming, data aggregation and merging processes. 

Artificial metagenomes would be a first dataset to be simulated, and the current distance metrics aren't useful in/with. But if you could partition reads, you could calculate profile distances of the partitions to the mass profile, which should be related to compositions after normalization.

The real question at my level is why would you need to give properties at this stage of the project to transient data projects that haven't even been entered on the issues. It's scope creep certainly, but I wanted to give readers an idea of where the project could have gone if there was more interest and development investment.

And profiling resource usage precludes an understanding or statistical investigation into the basic properties of what the software does at this point. If that's my approach to software.

## Installation

OS X and Linux release:

```sh
pip install kdb
```

Development installation:

```sh
git clone https://github.com/MatthewRalston/kdb.git
pip install requirements.txt#requirements-dev.txt
PYTHONPATH=$PYTHONPATH:$(pwd)
```

## Usage Example

CLI Usage

```bash
kdb --help
kdb summary --help
# Build a [composite] profile to a new or existing .kdb file
kdb profile example1.fq.gz example2.fq.gz profile.kdb
# Calculate similarity between two (or more) profiles
kdb similarity profile1.kdb profile2.kdb (...)
```

API usage

```python
from kdb import fileutil, kmer_util, profile

# Read a kdb file
kdb_rdr = fileutil.KDBReader(open("example.kdb", 'rb'))
kdb_rdr.read_profile()

# Print a profile
for c in kdb_rdr.profile:
  print(c)

# ... do something with the counts in the profile

# Save a kdb file
kdb_wrtr = fileutil.KDBWriter(open("example.kdb", 'wb'), kdb_rdr.get_header)
kdb_wrtr.write_profile(composite_profile, k)
```

## Documentation

Check out the [Readthedocs documentation](https://kdb.readthedocs.io/en/stable/), with examples and descriptions of the module usage.

## Development

```bash
pipenv run mamba test/*_spec.py
```

## License

Created by Matthew Ralston - [Scientist, Programmer, Musician](http://matthewralston.us) - [Email](mailto:mrals89@gmail.com)

Distributed under the GPL v3.0 license. See `LICENSE.txt` for the copy distributed with this project. Open source software is not for everyone, but for those of us starting out and trying to put the ecosystem ahead of ego, we march into the information age with this ethos.

## Contributing

1. Fork it (<https://github.com/MatthewRalston/kdb/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

