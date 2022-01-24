---
title: Data Analysis
nav-menu: true
layout: post
image: assets/images/southern_gel.jpg
---

<p style="font-size: 0.9rem;font-style: italic;"><img style="display: block;" src="https://live.staticflickr.com/68/205777138_995cea4253.jpg" alt="DNA - Blue"><a href="https://www.flickr.com/photos/99941535@N00/205777138">"DNA - Blue"</a><span> by <a href="https://www.flickr.com/photos/99941535@N00">Spanish Flea</a></span> is licensed under <a href="https://creativecommons.org/licenses/by-nc-nd/2.0/?ref=ccsearch&atype=html" style="margin-right: 5px;">CC BY-NC-ND 2.0</a><a href="https://creativecommons.org/licenses/by-nc-nd/2.0/?ref=ccsearch&atype=html" target="_blank" rel="noopener noreferrer" style="display: inline-block;white-space: none;margin-top: 2px;margin-left: 3px;height: 22px !important;"><img style="height: inherit;margin-right: 3px;display: inline-block;" src="https://search.creativecommons.org/static/img/cc_icon.svg?image_id=9b079eff-5fcd-41c8-80aa-6c74188f12e2" /><img style="height: inherit;margin-right: 3px;display: inline-block;" src="https://search.creativecommons.org/static/img/cc-by_icon.svg" /><img style="height: inherit;margin-right: 3px;display: inline-block;" src="https://search.creativecommons.org/static/img/cc-nc_icon.svg" /><img style="height: inherit;margin-right: 3px;display: inline-block;" src="https://search.creativecommons.org/static/img/cc-nd_icon.svg" /></a></p>

# How do I analyze k-mer profiles?

First, I would suggest that you begin with a solid review of the [wikipedia page](https://en.wikipedia.org/wiki/K-mer) for fundamentals. As I understand it, there are dozens of application spaces for k-mer in the bioinformatics application space. Without constraining the final application area or suggesting a particular use, we begin with what a k-mer profile can represent and then suggest some ways to generate matrices for further analysis.

# The fundamental matrix X

In this case we generate a matrix X that can be unnormalized, normalized, or with dimensionality reduction via PCA or t-SNE. Of course, additional transformations may be applied to the normalized or unnormalized X, by consuming the matrix X via the PyPI package `kmerdb` and then re-writing the transformed matrix X` to X.tsv.

```bash
#Check 'kmerdb matrix -h' for more details on usage

# Generating k-mer profiles
# Single profiles
kmerdb profile -k $k input1.fa input1.kdb # $k represents some common choice of k.
kmerdb profile -k $k input2.fa input2.kdb # Repeat via bash code as necessary to generate all inputs
# Compound profiles
#kmerdb profile -k $k -p 3 input1.fa input2.fa input3.fa first_three.kdb
#kmerdb profile -k $k -p 3 input4.fa input5.fa input6.fa second_three.kdb

# Unnormalized count matrix
kmerdb matrix -p $cores pass *.$k.kdb > X.tsv 
# In this case, X.tsv is a tsv for import into Pandas
```

Then implement a custom normalization, transformation, or processing step.
```python
import pandas as pd
df = pd.read_csv("X.tsv", sep="\t")

# Perform some normalization on X, and print as tsv
final_df.to_csv("X1.tsv", sep="\t", index=False)
```

Then, you can use the tsv as input to further exporatory steps such as clustering or distance matrix generation. All commands (matrix, distance, kmeans, hierarchical) are designed to mostly read and write tsv/csv through STDIN/STDOUT and are thus pipeable.

```bash
# Generate a Spearman correlation coefficient distance matrix
kmerdb distance spearman X.tsv
# Calculate the ssxy/sqrt(ssxx*ssyy) Pearson correlation coefficient
kmerdb distance pearson X.tsv
# Use SciPy to calculate correlation distance
kmerdb distance correlation X.tsv
#cat X.tsv | kmerdb distance spearman STDIN # or '/dev/stdin', has a bizarre syntax for reading from STDIN
# kmerdb distance -h

kmerdb kmeans -k 20 -i X.tsv sklearn
#cat X.tsv | kmerdb kmeans -k 20 Biopython
# kmerdb kmeans -h
kmerdb hierarchical -i X.tsv
#cat X.tsv | kmerdb hierarchical
# kmerdb hierarchical -h
```

# Transforming the data matrix with Unix pipes

```bash
# Create the initial data matrix X to use as tsv input to other commands
kmerdb matrix Unnormalized sample1.kdb sample2.kdb ... sampleN.kdb > X.tsv
# or, again with the awful syntax for reading from STDIN.
kmerdb matrix Unnormalized *.kdb | kmerdb matrix Normalized STDIN | kmerdb distance spearman STDIN | kmerdb kmeans -k 10 sklearn STDIN
```


I encourage you to check out the [CLI documentation](/quickstart#usage) for details on what functions are used for normalization, dimensionality reduction, and distance matrix generation.

Also, please take a look at [the example_report README.md](https://github.com/MatthewRalston/kmerdb/tree/master/examples/example_report) for more details about how to populate the report with metadata about an analysis of samples via kdb.

Then, please check out the template [index.Rmd](https://github.com/MatthewRalston/kmerdb/blob/master/examples/example_report/index.Rmd) for information about the statistical analyses performed and how these become the primary index.html page for the results.

And finally, please check out the [finished report](https://github.com/MatthewRalston/kmerdb/blob/master/examples/example_report/index.pdf).
