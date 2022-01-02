---
title: Install
nav-menu: true
layout: post
image: assets/images/installation-symbol.png
---


[![PyPI version](https://img.shields.io/pypi/v/kmerdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kmerdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kmerdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kmerdb)
[![Coveralls code coverage](https://coveralls.io/repos/github/MatthewRalston/kmerdb/badge.svg?branch=master)](https://coveralls.io/github/MatthewRalston/kmerdb?branch=master)
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]


[pip]: https://pypi.org/project/kmerdb/
[Pythons]: https://pypi.org/project/kmerdb/
[RTD]: https://kdb.readthedocs.io/en/latest/



## Optional R dependency

DESeq2 is required for rpy2-mediated normalization `kmerdb matrix Normalized *.kdb`

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

# Install

Install the latest version from PyPI

```bash
pip install kmerdb
```

Install locally for development

```bashn
git clone git@github.com:MatthewRalston/kmerdb.git
cd kdb
pip install --no-cache-dir . # python setup.py install
```

