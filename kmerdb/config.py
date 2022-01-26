'''
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

'''



VERSION="0.6.6"
header_delimiter = "\n" + ("="*24) + "\n"

metadata_schema = {
    "type": "object",
    "properties": {
        "version": {"type": "string"},
        "metadata_blocks": {"type": "number"},
        "k": {"type": "number"},
        "total_kmers": {"type": "number"},
        "unique_kmers": {"type": "number"},
        "unique_nullomers": {"type": "number"},
        "metadata": {"type": "boolean"},
        "sorted": {"type": "boolean"},
        "profile_dtype": {"type": "string"},
        "kmer_ids_dtype": {"type": "string"},
        "count_dtypes": {"type": "string"},
        "frequencies_dtype": {"type": "string"},
        "dtype": {"type": "string"},
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
                    "nullomers": {"type": "number"}
                },
                "required": ["filename", "sha256", "md5", "total_reads", "total_kmers", "unique_kmers", "nullomers"]
            }
        },
        "comments": {
            "type": "array",
            "items": {"type": "string"}
        }
    },
    "required": ["version", "metadata_blocks", "total_kmers", "unique_kmers", "unique_nullomers", "k", "tags", "files", "kmer_ids_dtype", "profile_dtype", "frequencies_dtype", "count_dtype"]
}




























pca_variance_fig_filepath = "PCA_variance_accumulation.png"
kmeans_elbow_graph_fig_filepath = "kmeans_elbow_graph.png"
kmeans_clustering_fig_filepath = "kmeans_clustering_of_kmer_profiles.png"
ecopy_rarefaction_fig_filepath = "ecopy_rarefaction_curve.png"
hierarchical_clustering_dendrogram_fig_filepath = "dendrogram.png"
spearman_upgma_tree_phy = "kdb_spearman_upgma_tree.phyloxml"
files = (pca_variance_fig_filepath, kmeans_elbow_graph_fig_filepath, kmeans_clustering_fig_filepath, ecopy_rarefaction_fig_filepath, hierarchical_clustering_dendrogram_fig_filepath)

DEFAULT_MASTHEAD = """

==========================================
               k m e r d b
==========================================

Thank you for using kmerdb. Please feel free
to submit issues and requests to 
https://github.com/MatthewRalston/kdb

Copyright 2020 Matt Ralston (mrals89@gmail.com)

# First links
==========================================-

|      Project     p a g e               |

===========================================

https://matthewralston.github.io/kmerdb

===========================================

|     P y P I                             |

===========================================
https://pypi.org/project/kmerdb/


===========================================

|          G i t h u b                    |

===========================================

https://github.com/MatthewRalston/kdb


============================================

|   N o t   Very   H u m e r u s           |

============================================
https://matthewralston.github.io/blog/kmer-database-format-part-1

Please cite my repository in your work!

Feel free to e-mail me or reach out!

|||Reddy47||| - Matt Ralston on Spotify
On top : what she sounds like when she's on top of the logging aesthetics.
Aston Martin Music : what she sounds like when she's on air in the logging aesthetics.
Pemex : A dichotomy between a wreckage and it's stubborn bgzf format. Just here for the noise, bass, and the smoke.
Sipping on Yak: What she sounds like when she's purring along in her new dtype specifications. Where did that come from... this stuff is kind of cool...
Death by Dishonor: What she sounds like when the fileutil slurps into memory. I have generally no idea how she sounds when things are running right, but this song give me the sense of urgency and mystery to keep listening.
Belial | Nit Grit : What she sounds like when she's churning through the lexical sort on the primary index. Distinguishing the count from the k-mer id is difficult to describe. 1:30 on, anyways. Checking in. And when the fileutil is barely dressing up the issues relating to absorption issues.
Who Dat? - Terror Reid, Getter

The soundtrack to my program. \s Also, hi mom.

"""
DONE = """

==========================================
----------------D O N E-------------------
==========================================

"""


DISTANCE_MASTHEAD = """

==========================================
           kmerdb distance
==========================================

Distance matrix generation beginning!

Distance matrix will be written to STDOUT as this is the first step of the pipeline.

"""
# for i in range(42, 1, -1):
#     DISTANCE_MASTHEAD += i*"=" + "\n"


MATRIX_MASTHEAD = """

==========================================
           kmerdb matrix
==========================================
kmerdb matrix <pass> test/data/*.8.kdb
kmerdb matrix PCA test/data/*.8.kdb
kmerdb matrix tSNE test/data/*.8.kdb
kmerdb matrix DESeq2 test/data/*.8.kdb
kmerdb matrix Frequency test/data/*.8.kdb

pass <==/==> don't normalize, adjust frequencies, perform any modification of k-mers.

This is what 95% of people are looking for: aggregating profiles into tabular .tsv format.

Matrix generation beginning!

Matrix will be written to STDOUT as this is the first step of the pipeline.

"""
# for i in range(42, 1, -1):
#     MATRIX_MASTHEAD += i*"=" + "\n"


KMEANS_MASTHEAD = """

==========================================
           kmerdb kmeans
==========================================

K-means clustering beginning!

"""
# for i in range(42, 1, -1):
#     KMEANS_MASTHEAD += i*"=" + "\n"


HIERARCHICAL_MASTHEAD = """

==========================================
           kmerdb hierarchical
==========================================

Hierarchical clustering beginning!

"""
    


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
# This command used to use PostgreSQL behind the scenes for on-disk k-mer counting
# since memory is limiting for profile generation when dealing 
# with biologically conventional choices of k (20 < k < 35).
#
# Now we're all in memory.
# I STRONGLY SUGGEST YOU START WITH MORE MODERATE CHOICES OF K (10 < k < 15)
parallel 'kmerdb profile -k $K {{}} {{.}}.$K.kdb' ::: $(/bin/ls test/data/*.fasta.gz)



# # # # # # #
# analysis  #
# # # # # # #

##################
# normalization
##################
# Use rpy2 and DESeq2 to normalize NB-distributed k-mer counts
# Graphical comparison can be made by comparing counts of unnormalized data to normalized
# kmerdb matrix [ Unnormalized | Normalized ] test/data/*.$K.kdb > (un)normalized_matrix.tsv
# 
#
# Install DESeq2 with the following, if not installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("DESeq2")

##################
# distance matrix
##################
# kmerdb distance spearman normalized_matrix.tsv > normalized_spearman_dist.tsv




################################
# W A R N I N G :  M E M O R Y #
################################
# The following is some 'basic' guesses about expected memory usage.
# I believe this is an oversimplification of the memory profile for technical language and representation issues.
# The profile loading function reads the array into memory from the file, assuming a certain data encoding.
# The coding is now set to 'uint64' by default, which in theory allows us to only be limited by numpy array size.
# Now that we have that out of the way.
# The integer depth we have is large.
#
#
#
##################
# dimensionality reduction + kmeans
##################
# The first step ('kdb matrix') generates one from different profiles with the same choice of k.
# This command uses ecopy to normalize between sample k-mer total counts before PCA/tSNE.
# -n $N is the dimensionality of either PCA or tSNE. A good choice for tSNE is 2.
# If the command is run with PCA without a selected dimensionality, an elbow graph
# will be produced named '{0}'. Please use this graph to select
# the number of principal components to use.
# The pipeline will not continue until -n $N is selected by the user.
# It is not recommended to feed Unnormalized or Normalized matrices directly to 'kmerdb kmeans'
# 
# The PCA/tSNE matrix will be dimReduced ($N) * N, where N is the number of samples/files/profiles.
#
# And finally, a k-means clustering will be done on the reduced dimensionality dataset
# Please note the randomness parameter 'random_state=42' for sklearn's kmeans is fixed at 42.
# Note here that the -k $K is not related to the choice of substring length 'k' for profile generation.
# The 'kmerdb kmeans' command produces two figures, first is an elbow graph looking at up to N clusters.
# This elbow graph will be written to '{1}'.
# The second is the more typical scatterplot of the first two reduced dimensions
# and the k-means clustering labels shown over the scatter.
# This file will be written to '{2}'.
kmerdb matrix [-n $N] [ PCA | tSNE ] normalized_matrix.tsv | kmerdb kmeans -k $K sklearn
kmerdb matrix [-n $N] [ PCA | tSNE ] normalized_matrix.tsv | kmerdb kmeans -k $K --distance e Biopython
#
# If you wanted to save a matrix from kdb matrix for use on your own
# it is recommended that you consider gzip compressing it if it is the Normalized or Unnormalized matrix
# which we will see is used downstream in the rarefaction and hierarchical analytics pathways.
#

##################
# Hierarchical
##################
#
# The Normalized matrix goes to the distance subcommand, which can use any of scipy's pdist distances
# to form the m x m distance matrix.
# The third step (kdb hierarchical)  is to build a dendrogram with scipy.cluster.hierarchy.
# This final step produces a plot in addition to the tsvs produced in the prior steps,
# which can be captured as independent steps or with tee in a pipeline.
kmerdb matrix [ Normalized ] test/data/*.$K.kdb | kmerdb distance spearman | kmerdb hiearchical

'''.format(*files)

