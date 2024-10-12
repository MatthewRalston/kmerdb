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



VERSION="0.8.9"
REQUIRES_PYTHON="3.7.4"
header_delimiter = "\n" + ("="*24) + "\n"


requirements_count = 11
requirements_dev_count = 8

# 9/22/24 Kind of? I mean...
# 4/9/24 there are some duplicates sure, but the requirements evolve faster than the dev do, are more essential for function, and I dont change the -dev file much. 

# 12 subcommands, and associated function names in __init__.py
subcommands = ["profile", "graph", "matrix", "distance", "kmeans", "hierarchical", "view", "header", "index", "shuf", "usage", "help", "version", "composition"]
subcommand_functions = ["profile", "make_graph", "get_matrix", "distances", "kmeans", "hierarchical", "view", "header", "index_file", "shuf", "expanded_help", "expanded_help", "version", "composition"]

default_logline_choices = (20, 50, 100, 200)
KDB_COLUMN_NUMBER = 4



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
        "unique_nullomers": {"type": "number"},
        "kmer_coverage": {"type": "number"},
        "kmer_coverage_histogram_dtype": {"type": "number"},
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
    "required": ["version", "metadata_blocks", "total_kmers", "unique_kmers", "unique_nullomers", "k", "tags", "files", "kmer_ids_dtype", "profile_dtype", "frequencies_dtype", "count_dtype", "unique_nullomers", "kmer_coverage", "kmer_coverage_histogram_dtype"]
}


exit_summary_schema = {
    "type": "object",
    "properties": {
        "subcommand": {"type": "string"},
        "kmerdb-version": {"type": "string"},
        "python-version": {"type": "string"},
        "feature": {"type": "number"},
        "feature_name": {"type": "string"},
        "feature_shortname": {"type": "string"},
        "feature_description": {"type": "string"},
        "step": {"type": "number"},
        "step_name": {"type": "string"},
        "step_shortname": {"type": "string"},
        "step_description": {"type": "string"},
        "log_file": {"type": "string"},
        "traceback": {"type": "string"},
        "error_file_name": {"type": "string"},
        "error_line_number": {"type": "number"},
        "error": {"type": "string"},        
    },
    "required": ["subcommand", "kmerdb-version", "python-version", "feature", "feature_name", "feature_shortname", "feature_description", "step", "step_name", "step_shortname", "step_description", "traceback", "error_file_name", "error_line_number", "error", ]
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
















pca_variance_fig_filepath = "PCA_variance_accumulation.png"
kmeans_elbow_graph_fig_filepath = "kmeans_elbow_graph.png"
kmeans_clustering_fig_filepath = "kmeans_clustering_of_kmer_profiles.png"
#ecopy_rarefaction_fig_filepath = "ecopy_rarefaction_curve.png"
hierarchical_clustering_dendrogram_fig_filepath = "dendrogram.png"
upgma_tree_phylip = "kdb_upgma_tree.phyloxml"
# files = (pca_variance_fig_filepath, kmeans_elbow_graph_fig_filepath, kmeans_clustering_fig_filepath, ecopy_rarefaction_fig_filepath, hierarchical_clustering_dendrogram_fig_filepath)

#######################################################

#          L o g o s

#######################################################



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









thanks = "\n\n\n" + "="*40 + """\n

  kmerdb v{0}    -     view -h|--help (or usage subcommand) for detailed information.


  please report any issues to https://github.com/MatthewRalston/kmerdb/issues and thanks.
""".format(VERSION)








#DONE = + usage + issue_tracker
# Accept the citation notice with 'kmerdb citation'
DONE = """
DONE\n
"""

































