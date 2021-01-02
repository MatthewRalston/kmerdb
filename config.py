VERSION="0.0.2"

header_schema = {
    "type": "object",
    "properties": {
        "version": {"type": "string"},
        "metadata_blocks": {"type": "number"},
        "k": {"type": "number"},
        "metadata": {"type": "boolean"},
        "tags": {
            "type": "array",
            "items": {"type": "string"}
        },
        "files": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "filename": {"type": "string"},
                    "sha256": {
                        "type": "string",
                        "minLength": 64,
                        "maxLength": 64
                    },
                    "md5": {
                        "type": "string",
                        "minLength": 32,
                        "maxLength": 32
                    },
                    "total_reads": {"type": "number"},
                    "total_kmers": {"type": "number"},
                    "unique_kmers": {"type": "number"},
                    "nullomers": {"type": "number"},
                    "mononucleotides": {
                        "type": "object",
                        "properties": {
                            "A": {"type": "number"},
                            "C": {"type": "number"},
                            "G": {"type": "number"},
                            "T": {"type": "number"}
                        }
                    }
                },
                "required": ["filename", "sha256", "md5", "total_reads", "total_kmers", "unique_kmers", "nullomers"]
            }
        },
        "comments": {
            "type": "array",
            "items": {"type": "string"}
        }
    },
    "required": ["version", "metadata_blocks", "k", "tags", "files"]
}




























pca_variance_fig_filepath = "PCA_variance_accumulation.png"
kmeans_elbow_graph_fig_filepath = "kmeans_on_pca_elbow_graph.png"
kmeans_clustering_fig_filepath = "kmeans_clustering_of_kmer_profiles.png"
ecopy_rarefaction_fig_filepath = "ecopy_rarefaction_curve.png"
files = (pca_variance_fig_filepath, kmeans_elbow_graph_fig_filepath, kmeans_clustering_fig_filepath, ecopy_rarefaction_fig_filepath)

DEFAULT_MASTHEAD = """

==========================================
                 k d b
==========================================

Thank you for using kdb. Please feel free
to submit issues and requests to 
https://github.com/MatthewRalston/kdb

Copyright 2020 Matt Ralston (mrals89@gmail.com)

# First links
https://matthewralston.github.io/blog/kmer-database-format-part-1
https://github.com/MatthewRalston/kdb


Please cite my repository in your work!

Feel free to e-mail me or reach out!

"""
DONE = """

==========================================
----------------D O N E-------------------
==========================================

"""


MATRIX_MASTHEAD = """

==========================================
           ./bin/kdb matrix
==========================================

Matrix generation beginning!

Matrix will be written to STDOUT as this is the first step of the pipeline.

"""
for i in range(42, 1, -1):
    MATRIX_MASTHEAD += i*"=" + "\n"


CLUSTER_MASTHEAD = """

==========================================
           ./bin/kdb cluster
==========================================

K-means clustering beginning!

"""
for i in range(42, 1, -1):
    CLUSTER_MASTHEAD += i*"=" + "\n"

RAREFY_MASTHEAD = """

==========================================
           ./bin/kdb rarefy
==========================================


Rarefaction beginning!

"""
for i in range(42, 1, -1):
    RAREFY_MASTHEAD += i*"=" + "\n"





DEBUG_MASTHEAD = '''

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
'''.format(*files)

