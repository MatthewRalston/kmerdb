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
import sys
import os

import yaml
from collections import OrderedDict






from kmerdb import config, util


yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)

PROGRAM_BANNER = """
**** 
 o-O      |||
o---O     |||             [|[          kmerdb           ]|]
O---o     |||
 O-o      |||        version :     v{0}
  O       |||
 o-O      |||        GitHub  : https://github.com/MatthewRalston/kmerdb/issues
o---O     |||         PyPI   : https://pypi.org/project/kmerdb/
O---o     |||      Website   : https://matthewralston.github.io/kmerdb
 O-o      |||
""".format(config.VERSION)
INTERPRETER = "                                                                       lang :         python\n"
# hardcoded

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

THREE_LINES = """



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

GITHUB_PROJECT_BANNER = """
=======================================================
                  ||      G i t H u b     ||
=======================================================
                         Repo: kmerdb
               Feature branch: graph_algo
-------------------------------------------------------
"""



PINNED_ISSUE = """
                 Pinned issue: #130
"""



class kmerdb_appmap:



        


    
    def __init__(self, argv):

        self.MODULE_ROOT = os.path.join("..", os.path.dirname(__file__))
        self.COMMAND_FILE = os.path.join(self.MODULE_ROOT, "__init__.py")
        self.PACKAGE_MANAGER = """
                      package manger : pip
                        version      : v24.0
        package root : {0}
        exe file     : {1}

                      required packages : {2}
                   development packages : {3}

           ARGV : {4}
        """.format(self.MODULE_ROOT, self.COMMAND_FILE, config.requirements_count, config.requirements_dev_count, argv)
        self.REQUIRES_PYTHON = config.REQUIRES_PYTHON
        self.VERSION_HARDCODED = "                                                                          v :         v{0}\n".format(config.REQUIRES_PYTHON)

        #
        # loaded_modules
        #

        #
        # dependencies
        #
        #      required:
        #      optional:

        # usage_notes.txt


        self.GRAPH_BANNER = """
                          name : graph
                   description : create edge nodes in (block) .gz compressed format from .fasta or .fq format.
        """

        self.GRAPH_PARAMS = OrderedDict({
            "name": "arguments",
            "values": [
                {
                    "name"  : "k",
                    "type"  : "int",
                    "value" : "choice of k-mer size"
                    },
                {
                    "name"  : "quiet",
                    "type"  : "flag",
                    "value" : "Write additional debug level information to stderr?"
                }
            ]
        })
        
        

        self.GRAPH_INPUTS = OrderedDict({
            "name": "inputs",
            "values": [
                {
                    "name"  : "<.fasta|.fastq>",
                    "type"  : "array",
                    "value" : "gzipped or uncompressed .fasta or .fastq file(s)"
                },
                {
                    "name"  : ".kdbg",
                    "type"  : "file",
                    "value" : "edge list."
                }
            ]
        })


        self.GRAPH_FEATURES = OrderedDict({
            "name": "features",
            "type": "array",
            "items": [
                OrderedDict({
                    "title": "k-mer id arrays organized linearly as file is read through sliding window. (Un)compressed support for .fa/.fq.",
                    "shortname": "parallel OOP sliding window k-mer shredding",
                    "description": "a Reader class is instantiated, and invoked on every file, supporting message passing and output collation. The data are read sequentially, so the possible edges for consideration are available by the identity of the k-mer and necessarily its 8 nearest neighbors. The appropriate pair observed in the dataset, and this pair is considered an 'edge' in the edge space constrained under k. It is added to the edge list by virtue of proximity in the file's base vector. In the case of secondary, tertiary, etc. sequences in the .fa or under massively parallel sequencing conditions (such as that by virtue of sequencing-by-synthesis) the offending edge is removed for the beginning and end of each sequence, and a warning is given to the user. Summary statistics for the entire file are given for each file read, as well as a cohesive structure, provided to the user before edge generation begins"
                }),
                OrderedDict({
                    "title": "k-mer neighbors assessed and tallied, creates a unsorted edge list, weight weights",
                    "shortname": "weighted undirected graph",
                    "description": "an edge list is generated from all k-mers in the forward direction of .fa/.fq sequences/reads. i.e. only truly neighboring k-mers in the sequence data are added to the tally of the k-mer nodes of the de Bruijn graph and the edges provided by the data."
                })
            ]
        })


        self.GRAPH_STEPS = OrderedDict({
            "name": "steps",
            "type": "array",
            "items": [
                OrderedDict({
                    "name": "read input file(s) from filesystem into k-mer arrays",
                    "shortname": "shred inputs into k-mer count arrays",
                    "description": "shred input sequences into k-mer count vector",
                }),
                OrderedDict({
                    "name": "merge k-mer arrays and aggregate metadata",
                    "shortname": "merge k-mer count arrays for aggregate metadata (header)",
                    "description": "merge counts of nullomers, unique kmers, and total kmers."
                }),
                OrderedDict({
                    "name": "collate the weighted edge list from graph.parsefile's Parser object during parallel file reading. This returns a edge_list, header, counts, and a nullomer_array.",
                    "shortname": "extract undirected weighted graph",
                    "description": "consists of info from a .kdb file and a .kdbg file. The node IDs, the edges, and the number of times the pair was observed from forward sequences in the provided dataset"
                }),
                OrderedDict({
                    "name": "print 'table' Final stats and close output file",
                    "shortname": "metrics and shutdown",
                    "description": "print final statistics, typically metadata values, and ensure file is closed."
                })
                
            ]
        })
        

        
        self.PROFILE_BANNER = """
                          name : profile
                   description : create a k-mer count vector. g-zipped and tallied.
        """

        self.PROFILE_PARAMS = OrderedDict({
            "name": "arguments",
            "values": [
                {
                    "name"  : "k",
                    "type"  : "parameter",
                    "value" : "choice of k-mer size parameter 'k'",
                },
                {
                    "name" : "sorted",
                    "type" : "flag",
                    "value": "descending in k-mer count"
                },
                {
                    "name"  : "quiet",
                    "type"  : "flag",
                    "value" : "Write additional debug level information to stderr?"
                }
            ]
        })

        self.PROFILE_INPUTS = OrderedDict({
            "name": "inputs",
            "values": [
                {
                    "name"  : "<.fasta|.fastq>",
                    "type"  : "array",
                    "value" : "gzipped or uncompressed .fasta or .fastq file(s)"
                },
                {
                    "name"  : ".kdb",
                    "type"  : "file",
                    "value" : "4 column table (row number, index number, k-mer id, count)."
                }
            ]
        })


        self.PROFILE_FEATURES = OrderedDict({
            "name": "features",
            "type": "array",
            "items": [
                OrderedDict({
                    "title": "k-mer id arrays organized linearly as file is read through sliding window. (Un)compressed support for .fa/.fq.",
                    "shortname": "parallel OOP sliding window k-mer shredding",
                    "description": "a Reader class is instantiated, and invoked on every file, supporting message passing and output collation. The data are read sequentially, so a length N = 4^k array is populated from the k-mer counts in the file. A k-mer id and count are recorded in the .kdb file."
                }),
                OrderedDict({
                    "title": "k-mers are tallied, and metadata merged from the input files",
                    "shortname": "merge counts, metadata from across inputs",
                    "description": "an final metadata structure and output metrics are collected for display to the user."
                })

            ]
        })


        self.PROFILE_STEPS = OrderedDict({
            "name": "steps",
            "type": "array",
            "items": [
                OrderedDict({
                    "name": "read input file(s) from filesystem into k-mer arrays",
                    "shortname": "shred inputs into k-mer count arrays",
                    "description": "shred input sequences into k-mer count vector",
                }),
                OrderedDict({
                    "name": "merge k-mer arrays and aggregate metadata",
                    "shortname": "merge k-mer count arrays for aggregate metadata (header)",
                    "description": "merge counts of nullomers, unique kmers, and total kmers."
                }),
                OrderedDict({
                    "name": "collate the count array",
                    "shortname": "collate the count array",
                    "description": "aggregate counts from different input files"
                }),
                OrderedDict({
                    "name": "print 'table' Final stats and close output file",
                    "shortname": "metrics and shutdown",
                    "description": "print final statistics, typically metadata values, and ensure file is closed."
                })

            ]

        })


        
        self.MATRIX_BANNER = """
                          name : matrix
                   description : create a k-mer array as a matrix/data-frame. tsv to stdout and/or output file.
        """
        
        self.MATRIX_PARAMS = OrderedDict({
            "name": "arguments",
            "type": "array",
            "items": [
                {
                    "name" : "function",
                    "type": "parameter",
                    "value" : "dimensionality reduction, DESeq2 normalization (via rpy2: assumes you have DESeq2 installed via Bioconductor [OPTIONAL DEPENDENCY]), (pass) no-transformation, or a k-mer Frequency matrix. [PCA, tSNE, DESeq2, pass, Frequency]",
                },
                {
                    "name" : "n",
                    "type":  "int",
                    "value": "dimension choise for scipy PCA or tSNE"
                },
                {
                    "name": "perplexity",
                    "type": "int",
                    "value" : "perplexity parameter for scipy tSNE"
                },
                {
                    "name": "no-normalized-ints",
                    "type": "flag",
                    "value": "do not round counts to nearest integer when importing from RPy2 from DESeq2 normalization"
                }
            ]
        })


        self.MATRIX_INPUTS = OrderedDict({
            "name": "inputs",
            "type": "array",
            "items": [
                {
                    "name": "kdbfiles",
                    "type": "file(s)",
                    "value": ".kdb count vector files produced from kmerdb profile"
                },
                {
                    "name": "input.tsv",
                    "type": "file",
                    "value": "A table of csv/tsv file counts, possibly for dimensionality reduction."
                },
                {
                    "name": "STDIN",
                    "type": "flag",
                    "value": "Read from STDIN, possibly in a unix pipe."
                }
                
            ]
        })

        self.MATRIX_FEATURES = OrderedDict({
            "name": "features",
            "type": "array",
            "items": [
                OrderedDict({
                    "name": "read/import datasets",
                    "shortname": "read/import multiple count vectors",
                    "description" : "count vectors in columns of csv/tsv input over 'STDIN', .tsv/csv file (with --delimiter option), or multiple .kdb files."
                }),
                OrderedDict({
                    "name": "Pandas library .csv sanitation",
                    "shortname": "validate input",
                    "description": "sanitize and structure validation on inputs (TODO: 4/5/24 not done yet. Uses a pandas csv loader on tsv or STDIN input. Otherwise collates multiple .kdb count vectors and coerces the input into a NumPy array otherwise failure.)"
                }),
                OrderedDict({
                    "name": "transformation on data matrix",
                    "shortname": "transformation on data matrix",
                    "description": "Apply 'transformation' to input matrix. Not a distance calculation, but a restructuring of the data. Options are PCA (+scipy), t-stochastic neighbor embedding(scipy), 'DESeq2' median-of-ratios normalization, or retrieve frequencies for related applications."
                })
            ]
        })

        self.MATRIX_STEPS = OrderedDict({
                "name": "steps",
                "type": "array",
                "items": [
                    OrderedDict({
                        "name": "read input file(s) from filesystem into k-mer arrays",
                        "shortname": "read .kdb files into k-mer count arrays",
                        "description": "read multiple .kdb files into k-mer count matrix",
                    }),
                    OrderedDict({
                        "name": "parse .csv into Pandas DataFrame",
                        "shortname": "Pandas .csv import from input file or STDIN",
                        "description": "merge counts of nullomers, unique kmers, and total kmers."
                    }),
                    OrderedDict({
                        "name": "Apply PCA dimensionality reduction with provided dimensionality n",
                        "shortname": "PCA dim reduction",
                        "description": "Uses scipy PCA to apply dimensionality reduction on input matrix"
                    }),
                    OrderedDict({
                        "name": "Apply t-SNE dimensionality reduction with provided dimensionality n",
                        "shortname": "t-SNE dim reduction",
                        "description": "Uses scipy t-SNE to apply dimensionality reduction on input matrix"
                    }),
                    OrderedDict({
                        "name": "Apply DESeq-2 'median-of-ratios' count normalization",
                        "shortname": "DESeq-2 normalize",
                        "description": "Uses rpy2 to call the R server and pass the data matrix for DESeq-2 normalization"

                    })

                ]

        })

    
    
        self.DISTANCE_BANNER = """
                          name : distance
                   description : create a distance matrix, using scipy distances, or even my Cython Pearson corr coeff (pearson).
                                  Spearman rank corr coeff also availble (vis scipy). Non-cython correlation coefficient availalble
                                  as 'correlation' distance in scipy.
        """

        self.DISTANCE_PARAMS = OrderedDict({
            "name": "arguments",
            "type": "array",
            "items": [
                {
                    "name" : "distance method",
                    "type": "parameter",
                    "value" : "braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,matching,minkowski,pearson,rogerstanimotorusselrao,seuclidean,sokalmichener,sokalsneath,spearman,sqeuclidean,yule distances. spearman and correlation distances provided by scipy. Pearson is a custom Cython correlation coefficient. all other provided by scipy",
                },
            ]
        })

        
        self.DISTANCE_INPUTS = OrderedDict({
            "name": "inputs",
            "type": "array",
            "items": [
                {
                    "name": "kdbfiles",
                    "type": "file(s)",
                    "value": ".kdb count vector files produced from kmerdb profile"
                },
                {
                    "name": "input.tsv",
                    "type": "file",
                    "value": "A table of csv/tsv file counts, possibly for dimensionality reduction."
                },
                {
                    "name": "STDIN",
                    "type": "flag",
                    "value": "Read from STDIN, possibly in a unix pipe."
                }
                
            ]
        })

        self.DISTANCE_FEATURES = OrderedDict({
            "name": "features",
            "type": "array",
            "items": [
                OrderedDict({
                    "name": "read/import datasets",
                    "shortname": "read/import multiple count vectors",
                    "description" : "count vectors in columns of csv/tsv input over 'STDIN', .tsv/csv file (with --delimiter option), or multiple .kdb files."
                }),
                OrderedDict({
                    "name": "Pandas library .csv sanitation",
                    "shortname": "validate input",
                    "description": "sanitize and structure validation on inputs (TODO: 4/5/24 not done yet. Uses a pandas csv loader on tsv or STDIN input. Otherwise collates multiple .kdb count vectors and coerces the input into a NumPy array otherwise failure.)"
                }),
                OrderedDict({
                    "name": "generates distance matrix",
                    "shortname": "generates distance matrix",
                    "description": "Choose 'distance method' to apply to input matrix. SciPy distance methods plus a Cython Pearson coefficient of correlation. All distances availalable on -h|--help"
                })
            ]
        })
    
        self.DISTANCE_STEPS = OrderedDict({
            "name": "steps",
                "type": "array",
                "items": [
                    OrderedDict({
                        "name": "read input file(s) from filesystem into k-mer arrays",
                        "shortname": "read .kdb files into k-mer count arrays",
                        "description": "read multiple .kdb files into k-mer count matrix",
                    }),
                    OrderedDict({
                        "name": "parse .csv into Pandas DataFrame",
                        "shortname": "Pandas .csv import from input file or STDIN",
                        "description": "merge counts of nullomers, unique kmers, and total kmers."
                    }),
                    OrderedDict({
                        "name": "calculate distance matrix from input data matrix",
                        "shortname": "distance matrix generation",
                        "description": "Use input data matrix to create a distance matrix, produced by SciPy and a choice of distance method"
                    })

                ]

        })
        # KMEANS_BANNER = """
        #                   name : kmeans
        #            description : 
        # """

        # HIERARCHICAL_BANNER = """

        # """


        # config.py features
        #SUPPORTED_FEATURES = {}

        # ASCII + RICH.py checklist
        

    def print_program_header(self):
        sys.stderr.write(PROGRAM_BANNER)
        sys.stderr.write(INTERPRETER)
        sys.stderr.write(self.VERSION_HARDCODED)
        sys.stderr.write(self.PACKAGE_MANAGER)

        # Spacer
        sys.stderr.write(DNA_COLUMN_1)

    def print_graph_header(self):

        
        sys.stderr.write(self.GRAPH_BANNER)

        sys.stdedrr.write(THREE_LINES)
        
        sys.stderr.write(yaml.dump(self.GRAPH_PARAMS))
        
        sys.stdedrr.write(THREE_LINES)
        
        sys.stderr.write(yaml.dump(self.GRAPH_INPUTS))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.GRAPH_FEATURES))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.GRAPH_STEPS))

    def print_profile_header(self):
        
        sys.stderr.write(self.PROFILE_BANNER)

        sys.stdedrr.write(THREE_LINES)
            
        sys.stderr.write(yaml.dump(self.PROFILE_PARAMS))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.PROFILE_INPUTS))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.PROFILE_FEATURES))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.PROFILE_STEPS))

    def print_matrix_header(self):

        sys.stderr.write(self.MATRIX_BANNER)

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.MATRIX_PARAMS))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.MATRIX_INPUTS))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.MATRIX_FEATURES))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.MATRIX_STEPS))


    def print_distance_header(self):

        sys.stderr.write(self.DISTANCE_BANNER)

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.DISTANCE_PARAMS))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.DISTANCE_INPUTS))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.DISTANCE_FEATURES))

        sys.stdedrr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(self.DISTANCE_STEPS))






# #
# #          logger and selected log lines
# #


# # LOGGER BANNER + SPACER




# #               POST_LOG_HEADER
# #

# PROFILE_BANNER = """
# =======================================================
#                ||     p r o f i l e     ||
# =======================================================
#                     inputs  : 
#                     outputs :
#                  output_dir :
#                     logfile :
#               relevant_loc  :     # an hash of value: documentation strings and lines-of-code(loc) given a feature or error.
#            relevant_issues  : 



#                  parameters : {0}
# -------------------------------------------------------

# PARAMS_YAML:


# parameters : [
#  - param1 (DEFAULT: SOME_DEFAULT_VALUE, type=) : 
#  - param2 (some_functionality)  :


# -------------------------------------------------------

# """
# HEADER_BANNER =
# VIEW_BANNER =
# MATRIX_BANNER =
# DISTANCE_BANNER =
# GRAPH_BANNER = 
# KMEANS_BANNER =
# HIERARCHICAL_BANNER = 


# # From usage-notes (include description
# #
# # parameters:
# #      - (short), (long), type, description, examples, help 
# #      - ...
# #      -
# #


# # ISO-8601
# RUNTIME


# LOGFILE

# EXIT_CODE

# TRACEBACK

# # caught exception
# LOGGABLE_LINE


# # 20, 50, 100, 200 -n
# LAST_N_LINES_OF_LOG



# METADATA # [metadata] via YAML representer, key indices, key/values, last program stage
# METADATA_VALUES_DESCRIPTIONS # Just simple key-value descriptions


# METADATA_SCHEMA_TXT # Use YAML representer...
# METADATA_SCHEMA_DESCRIPTION # What the loops, maps, conditionals/branches, and data structure validation mean in this data structure, how it changes, where invalid scenarios may be involved, which test datasets to use and test command.
# SCHEMA_COUNT = "                      schema count: {0}".format(schemas_count)


# PROGRAM_STAGE_HEADER = """
# ========================================================

#            s t a g e : 1,2,3,...

# ========================================================

#              name     :
#          description  :

#                          LOREM IPSUM
# """

# TYPES_OF_ERRORS = """
#               name: error1
#       title   : something_descriptive_for_example
#           description :

#                  LOREM IPSUM

#           related_issues : (N/A) 133, 132, 101
#                  exit_codes : 1, 2, 3, ...etc.,
# """
















