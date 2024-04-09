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



VERSION="0.7.8"
REQUIRES_PYTHON="3.7.4"
header_delimiter = "\n" + ("="*24) + "\n"


requirements_count = 8
requirements_dev_count = 14 # 4/9/24 there are some duplicates sure, but the requirements evolve faster than the dev do, are more essential for function, and I dont change the -dev file much. 

subcommands = ["usage", "help", "profile", "graph", "index", "shuf", "matrix", "distance"] # kmeans and hierarchical and probability commands deprecated



graph_schema = {
    "type": "object",
    "properties": {
        "version": {"type": "string"},
        "metadata_blocks": {"type": "number"},
        "k": {"type": "number"},
        "total_kmers": {"type": "number"},
        "unique_kmers": {"type": "number"},
        "total_nullomers": {"type": "number"},
        "sorted": {"type": "boolean"},
        "n1_dtype": {"type": "string"},
        "n2_dtype": {"type": "string"},
        "weights_dtype": {"type": "string"},
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
    "required": ["version", "metadata_blocks", "k", "tags", "files", "total_kmers", "unique_kmers", "unique_nullomers", "n1_dtype", "n2_dtype", "weights_dtype"]
}

kdb_metadata_schema = {
    "type": "object",
    "properties": {
        "version": {"type": "string"},
        "metadata_blocks": {"type": "number"},
        "k": {"type": "number"},
        "total_kmers": {"type": "number"},
        "unique_kmers": {"type": "number"},
        "total_nullomers": {"type": "number"},
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









statfile =               "NAME.stats.txt"
logfile =                "NAME.stderr.log" # Primary logfile
output_usage_note =                "NAME.stdout.kmerdb.SUBCOMMAND.usage.txt"
#                "NAME.stderr.ascii.log",
primary_output =                "NAME.K.kdbg"


# FORMAT IN 
# kdb_program_metadata = {
#     "type": "object",
#     "properties": {
#         "outputs": {
#             "type": "array",
#             "items": [
#                 output_usage_note,
#                 logfile,
#                 statfile,
#                 primary_output,
#                 "NAME.K.kdb"
#             ]
#         },
#         "output_dir":
#         "logfile"
#     },
#     required = []
# }

# kdbg_program_metadata = {
#     "type": "object",
#     "properties": {
#         "outputs": {
#             "type": "array",
#             "items": [
#                 output_usage_note,
#                 logfile,
#                 statfile,
#                 primary_output,
#                 "NAME.K.kdbg"
#             ]
#         }
#         "output_dir": {
#             "type": "string"
#         }
#     },
#     required = ["output_dir", "outputs"]
# }
















# pca_variance_fig_filepath = "PCA_variance_accumulation.png"
# kmeans_elbow_graph_fig_filepath = "kmeans_elbow_graph.png"
# kmeans_clustering_fig_filepath = "kmeans_clustering_of_kmer_profiles.png"
# #ecopy_rarefaction_fig_filepath = "ecopy_rarefaction_curve.png"
# hierarchical_clustering_dendrogram_fig_filepath = "dendrogram.png"
# spearman_upgma_tree_phy = "kdb_spearman_upgma_tree.phyloxml"
# files = (pca_variance_fig_filepath, kmeans_elbow_graph_fig_filepath, kmeans_clustering_fig_filepath, ecopy_rarefaction_fig_filepath, hierarchical_clustering_dendrogram_fig_filepath)

#######################################################

#          L o g o s

#######################################################

KMERDB_LOGO = """
 o-O      |||
o---O     |||             [|[          kmerdb           ]|]
O---o     |||
 O-o      |||        version :     v{0}
  O       |||
 o-O      |||        GitHub  : https://github.com/MatthewRalston/kmerdb/issues
o---O     |||         PyPI   : https://pypi.org/project/kmerdb/
O---o     |||      Website   : https://matthewralston.github.io/kmerdb
 O-o      |||
""".format(VERSION)

GITHUB_LOGO = """
.--------------------------------------------------.
|                 .mmMMMMMMMMMMMMMmm.              |
|             .mMMMMMMMMMMMMMMMMMMMMMMMm.          |
|          .mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm.       |
|        .MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM.     |
|      .MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM.   |
|     MMMMMMMM'  `"MMMMM"""""""MMMM""`  'MMMMMMMM  |
|    MMMMMMMMM                           MMMMMMMMM |
|   MMMMMMMMMM:                         :MMMMMMMMMM|
|  .MMMMMMMMMM                           MMMMMMMMMM.
|  MMMMMMMMM"                             "MMMMMMMMM
|  MMMMMMMMM                               MMMMMMMMM
|  MMMMMMMMM                               MMMMMMMMM
|  MMMMMMMMMM                             MMMMMMMMMM
|  `MMMMMMMMMM                           MMMMMMMMMMM
|   MMMMMMMMMMMM.                     .MMMMMMMMMMMMM
|    MMMMMM  MMMMMMMMMM         MMMMMMMMMMMMMMMMMMM|
|     MMMMMM  'MMMMMMM           MMMMMMMMMMMMMMMM` |
|      `MMMMMM  "MMMMM           MMMMMMMMMMMMMM`   |
|        `MMMMMm                 MMMMMMMMMMMM`     |
|          `"MMMMMMMMM           MMMMMMMMM"`       |
|             `"MMMMMM           MMMMMM"`          |
|                 `""M           M""`              |
'--------------------------------------------------'
"""




SPONGEBOB = """
⢀⣠⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣤⣄⡀
⣿⡋⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⢙⣿
⣿⡇⠀⠀⠀⣠⣴⠾⠿⠷⣦⡀⢀⣴⠾⠿⠷⣦⣄⠀⠀⠀⢸⣿
⢸⡇⠀⠀⢰⡟⠁⠀⣀⣀⠈⢿⡿⠁⣀⣀⠀⠈⢻⡆⠀⠀⢸⡇
⢸⣷⠀⠀⢸⡇⠀⠀⠿⠿⠁⣸⣇⠈⠿⠿⠀⠀⢸⡇⠀⠀⣾⡇
⠘⣿⠀⠀⠈⠻⣦⣄⣀⣤⣾⠛⠛⣷⣤⣀⣠⣴⠟⠁⠀⠀⣿⠃
⠀⣿⠀⠘⢷⣄⠀⠉⠉⠙⢿⠄⠠⡿⠋⠉⠉⠀⣠⡾⠃⠀⣿⠀
⠀⣿⠀⠀⠀⠙⠻⢶⣦⣤⣤⣤⣤⣤⣤⣴⡶⠟⠋⠀⠀⠀⣿⠀
⠀⣿⡄⠀⠀⠀⠀⠀⣿⣀⣹⡇⢸⣏⣀⣿⠀⠀⠀⠀⠀⢠⣿⠀
⠀⢿⡇⠀⠀⠀⠀⠀⠉⠉⠉⠁⠈⠉⠉⠉⠀⠀⠀⠀⠀⢸⡿⠀
⠀⢸⣿⣿⣿⣿⣟⠛⠛⠛⣻⡿⢿⣟⠛⠛⠛⣻⣿⣿⣿⣿⡇⠀
⠀⢸⣿⣿⣿⣿⣿⣷⣤⣾⣿⣶⣶⣿⣷⣤⣾⣿⣿⣿⣿⣿⡇⠀
⠀⠘⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠘⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀
⠀⠀⠉⠉⠉⠉⠉⠉⠉⠉⠻⣧⣼⠟⠉⠉⠉⠉⠉⠉⠉⠉⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀

"""




#######################################################

#          E n v i r o n m e n t

#######################################################






LANGUAGE = """


                                                                      lang :         python
                                                                         v :         v3.8.0
                                                                             
                                                                                      
                                                                                     
"""

# A language specific version number
#PYTHON_VERSION = """

                                                               
#"""

# A package manager version
PIP_VERSION = """
                      package manger : pip
                        version      : v24.0
"""



DEPENDENCY_COUNT = """
                       dependencies  : {0}
           development_dependencies  : {1}
""".format(requirements_count, requirements_dev_count)




#SUBCOMMAND_NAME
#VERSION




# OUTPUTS determined at runtime
# OUTPUT_DIR/PREFIX determined at runtime
# Logfile (in output_dir) determined at runtime
# Verbosity determined at runtime


#######################################################

#          b a n n e r s  / mini - h e a d e r s

#######################################################


GITHUB_PROJECT_INFO = """
=======================================================
                  ||      G i t H u b     ||
=======================================================
                         Repo: kmerdb
                 Pinned issue: #130
               Feature branch: graph_algo
-------------------------------------------------------
"""












#######################################################

#          s p a c e r s

#######################################################













FIVE_LINES = """





"""

DNA_SPACER_1 = """
=================================
=================================


O       o O       o O       o O
| O   o | | O   o | | O   o | | O
| | O | | | | O | | | | O | | | |
| o   O | | o   O | | o   O | | o
o       O o       O o       O O


=================================
"""

DNA_SPACER_lol = """
=================================
=================================
Carbon rules everything around me

O       o O       o O       o O
| O   o | | O   o | | O   o | | O
| | O | | | | O | | | | O | | | |
| o   O | | o   O | | o   O | | o
o       O o       O o       O O


=================================
"""


DNA_COLUMN_1 = """
O---o
 O-o
  O
 o-O
o---O
O---o
 O-o
  O
 o-O
o---O
O---o
 O-o
  O
 o-O
o---O
"""






"""
# "@->- -|>" ... (@->-)
# "..., -|> ?"
"""


#######################################################

#          F o o t e r s

#######################################################


FOOTER_COMMAND_SUMMARY = """
command : COMMAND
runtime : RUNTIME
logfile : LOGFILE
exit_status : EXITCODE


0---------------------------0

.kdb_config file : CONFIG_PATH
was_pipeline : IS_PIPELINE

-----

conda_environment_validated? : CONDA_READY

"""


















#DONE = + usage + issue_tracker
# Accept the citation notice with 'kmerdb citation'
DONE = """
DONE
"""












































#######################################################

#          O l d       M a s t h e a d s

#######################################################



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
kmerdb matrix <* pass *> test/data/*.8.kdb
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


GRAPH_MASTHEAD = """

==========================================
           kmerdb graph
==========================================
current interface v.0.7.8

kmerdb graph -k 12 example_seqs1.fa.gz example_seqs2.fa.gz         < ... microbe1.fa.gz microbe2.fa.gz > output.kdbg

[ note ]

kmerdb graph -h





-----

future interface, p. v0.8.0
kmerdb graph example.fa example.kdbg --prefix temp_output
-----


.
.
.
output.kdbg
output.kdb
output.kdbi
output.kdb.gi
error.log
debug.log
output.kdbg.stats.neighbors.txt
output.kdbg.stats.edges.txt
output.kdbg.stats.adjacency.txt
---jk
output.kdbg.stats_txt
---------------------
distribution.txt
quick_stats.txt
---------------------
 m e t a d a t a
---------------------
kmers-per-file (array ok)

kmers-total  :
kmers-avg   :
singletons :
doublets :
triplets :
---------------------
notes...
-----











...final metadata schema unspecified

k-mer uniqueness  -- as ( a decimal value between o an 1)    = unique_kmers / size of kspace / total_kmers     (( ...note : not included in header metadata as of v0.7.8 )
'meaningful_kmers = unique_kmers + total_kmers'  (per file)     :  {0}      =    ' {1} '  ' + ' ' {2} '   
grand_total_of_kmers     :    
--------------------------------------------------------------------------------------

(table structure

kmer_totals  :  (array okay)  (per file)
unique kmers :  (array okay)  (per file)
nullomers    :  (array okay)  (per file)


unique_nullomers

output.kdbg.stats.doublets
output.kdbg.stats.triplets
output.kdb.stats.average
output.k
output.kdbg.stats.paths.txt


"""



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
    

