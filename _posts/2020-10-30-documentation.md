---
category: 'What is .kdb?'
title: 'Documentation'
path: '/documentation'

layout: nil
---

The kdb project was designed to facilitate conversation between heavily optimized legacy codebases without much public attention, like Jellyfish, regarding the utility of standardizing k-mer frameworks. These frameworks are used throughout assembly and alignment hashing/seed-matching strategies. The primary goal of this project is documenting data shapes, compression strategies (which of course related to efficiency of storage, transmission, rapid access, etc.), and anticipating UI possibilities with the increases in read/write speeds afforded by improving SSD technologies and utilization of more channels of more rapid interfaces for data transmission (i.e. m2, NVMe, PCIx). 


With this in mind, the documentation for the libraries themselves can be found here: https://kdb.readthedocs.io/en/stable

Additionally, running the matrix, cluster, or rarefy commands (which are arguably more complex in their usage) should be run with the DEBUG (-vv) verbosity setting on. This will yield additional information about the expected pipeline usage. That statement is echoed here.


```bash
The workflow is roughly as follows:

# # # # # # # # # #  #
# profile generation #
# # # # # # # # # #  #
# I have included -p and -b parameters that influence the rate of profile generation.
# The block size (-b) is primarily for the number of .fastq(.gz) records to read at a time.
# The -p parallel feature works within fastq parsing to process k-mer shredding/counting.
#
# -k $K is the parameter that specifies the k-mer resolution
#
# This command uses SQLite3 behind the scenes for on-disk k-mer counting
# since memory is rate limiting for profile generation when dealing 
# with biologically conventional choices of k (20 < k < 35).
parallel 'kdb profile -k $K {{}} {{.}}.$K.kdb' ::: $(/bin/ls test/data/*.fasta.gz)



# # # # # # #
# analysis  #
# # # # # # #
################################
# W A R N I N G :  M E M O R Y #
################################
# The first step of either rarefaction or clustering is to generate the k-mer profile matrix
# The matrix is not a new concept, it is just samples as columns of profiles. 
# Without metadata or graphical information embedded in the format, 
# we can create a matrix and do data science with the information.
# However, because the profiles are exponential in k, a linear increase in k 
# hypothetically results in an improvement of resolution that is at least superlinear in k.
# Therefore, the amount of memory can be calculated by the integer size times the profile resolution
# times the number of samples.
#
#
##################
# Cluster analysis
##################
# The first step ('kdb matrix') generates one from different profiles with the same choice of k.
# This command uses ecopy to normalize between sample k-mer total counts before PCA/tSNE.
# -n $N is the dimensionality of either PCA or tSNE. A good choice for tSNE is 2.
# If the command is run with PCA without a selected dimensionality, an elbow graph
# will be produced named '{0}'. Please use this graph to select
# the number of principal components to use.
# The pipeline will not continue until -n $N is selected by the user.
# It is not recommended to feed Unnormalized or Normalized matrices directly to 'kdb cluster'
# 
# The PCA/tSNE matrix will be dimReduced ($N) * N, where N is the number of samples/files/profiles.
#
# And finally, a k-means clustering will be done on the reduced dimensionality dataset
# Note here that the -k $K is not related to the choice of substring length 'k' for profile generation.
# The 'kdb cluster' command produces two figures, first is an elbow graph looking at up to N clusters.
# This elbow graph will be written to '{1}'.
# The second is the more typical scatterplot of the first two reduced dimensions
# and the k-means clustering labels shown over the scatter.
# This file will be written to '{2}'.
kdb matrix [-n $N] [ PCA | tSNE ] test/data/*.$K.kdb | kdb cluster -k $K

#
# If you wanted to save a matrix from kdb matrix for use on your own
# it is recommended that you consider gzip compressing it if it is the Normalized or Unnormalized matrix
# which we will see is used downstream in the rarefaction analytical pathway.
#
##################
# Rarefaction
##################
# The Unnormalized and Normalized matrices go to ecopy's rarefy function to produce a rarefaction plot
# This plot will be written to '{3}'.
kdb matrix [ Unnormalized | Normalized ] test/data/*.$K.kdb | kdb rarefy
```
