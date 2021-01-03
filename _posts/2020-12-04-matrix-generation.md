---
category: Usage
title: 'kdb matrix'
path: '/kdb_matrix'

layout: nil
---

## How to generate normalized and reduced dimension count matrices.

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



