


from collections import OrderedDict


#        yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)





from kmerdb import config


# pyproject.toml REQUIRED_PYTHON_VERSION
REQUIRES_PYTHON=config.


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
"""
INTERPRETER = "                                                                       lang :         python"
# hardcoded
VERSION_HARDCODED = "                                                                          v :         v3.8.0"
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


DNA_SPACER_1 = 
"""
=================================
=================================


O       o O       o O       o O
| O   o | | O   o | | O   o | | O
| | O | | | | O | | | | O | | | |
| o   O | | o   O | | o   O | | o
o       O o       O o       O O


=================================
"""

DNA_SPACER_lol =
"""
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



DNA_COLUMN_1 =
"""
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




def class kmerdb_appmap:

    def __init__(self, argv):









        MODULE_ROOT = path.join("..", path.dirname(__file__))
        COMMAND_FILE = path.join(MODULE_ROOT, "__init__.py")
        PACKAGE_MANAGER = """
                      package manger : pip
                        version      : v24.0
                        package root : {0}
                        exe file     : {1}


           ARGV : {2}
        """.format(MODULE_ROOT, COMMAND_FILE, argv)

        #
        # loaded_modules
        #

        #
        # dependencies
        #
        #      required:
        #      optional:

        # usage_notes.txt


        GRAPH_BANNER = """
                          name : graph
                   description : create edge nodes in (block) .gz compressed format from .fasta or .fq format.
        """

        GRAPH_PARAMS = OrderedDict({
            "name": "arguments",
            "values": [
                {
                    "name"  : "k",
                    "type"  :,
                    "value" : "choice of k-mer size"
                    },
                {
                    "name"  : "quiet",
                    "type"  : "flag",
                    "value" : "Write additional debug level information to stderr?"
                }
            ]
        })
        
        

        GRAPH_INPUTS = OrderedDict({
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


        GRAPH_FEATURES = OrderedDict({
            "name": "features",
            "type": "array",
            "items": [
                {
                    "title": "k-mer id arrays organized linearly as file is read through sliding window. (Un)compressed support for .fa/.fq."
                    "shortname": "parallel OOP sliding window k-mer shredding"
                    "description": "a Reader class is instantiated, and invoked on every file, supporting message passing and output collation. The data are read sequentially, so the possible edges for consideration are available by the identity of the k-mer and necessarily its 8 nearest neighbors. The appropriate pair observed in the dataset, and this pair is considered an 'edge' in the edge space constrained under k. It is added to the edge list by virtue of proximity in the file's base vector. In the case of secondary, tertiary, etc. sequences in the .fa or under massively parallel sequencing conditions (such as that by virtue of sequencing-by-synthesis) the offending edge is removed for the beginning and end of each sequence, and a warning is given to the user. Summary statistics for the entire file are given for each file read, as well as a cohesive structure, provided to the user before edge generation begins"
                },
                {
                    "title": "k-mer neighbors assessed and tallied, creates a unsorted edge list, weight weights",
                    "shortname": "weighted undirected graph",
                    "description": "an edge list is generated from all k-mers in the forward direction of .fa/.fq sequences/reads. i.e. only truly neighboring k-mers in the sequence data are added to the tally of the k-mer nodes of the de Bruijn graph and the edges provided by the data."
                }
            ]
        })


        GRAPH_STEPS = OrderedDict({
            "name": "steps",
            "type": "array",
            "items": [
                {
                    "name": "read input file(s) from filesystem into k-mer arrays",
                    "shortname": "shred inputs into k-mer count arrays",
                    "description": "shred input sequences into k-mer count vector",
                },
                {
                    "name": "merge k-mer arrays and aggregate metadata",
                    "shortname": "merge k-mer count arrays for aggregate metadata (header)",
                    "description": "merge counts of nullomers, unique kmers, and total kmers."
                },
                {
                    "name": "collate the weighted edge list from graph.parsefile's Parser object during parallel file reading. This returns a edge_list, header, counts, and a nullomer_array.",
                    "shortname": "extract undirected weighted graph",
                    "description": "consists of info from a .kdb file and a .kdbg file. The node IDs, the edges, and the number of times the pair was observed from forward sequences in the provided dataset"
                },
                {
                    "name": "print 'table' Final stats and close output file",
                    "shortname": "metrics and shutdown",
                    "description": "print final statistics, typically metadata values, and ensure file is closed."
                }
                
            ]
        })
        

        
        PROFILE_BANNER = """
                          name : profile
                   description : create a k-mer count vector. g-zipped and tallied.
        """

        PROFILE_PARAMS = OrderedDict({
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

        PROFILE_INPUTS = OrderedDict({
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


        PROFILE_FEATURES = OrderedDict({
            "name": "features",
            "type": "array",
            "items": [
                {
                    "title": "k-mer id arrays organized linearly as file is read through sliding window. (Un)compressed support for .fa/.fq."
                    "shortname": "parallel OOP sliding window k-mer shredding"
                    "description": "a Reader class is instantiated, and invoked on every file, supporting message passing and output collation. The data are read sequentially, so a length N = 4^k array is populated from the k-mer counts in the file. A k-mer id and count are recorded in the .kdb file."
                },
                {
                    "title": "k-mers are tallied, and metadata merged from the input files",
                    "shortname": "merge counts, metadata from across inputs",
                    "description": "an final metadata structure and output metrics are collected for display to the user."
                }

            ]
        })


        PROFILE_STEPS = OrderedDict({
            "name": "steps",
            "type": "array",
            "items": [
                {
                    "name": "read input file(s) from filesystem into k-mer arrays",
                    "shortname": "shred inputs into k-mer count arrays",
                    "description": "shred input sequences into k-mer count vector",
                },
                {
                    "name": "merge k-mer arrays and aggregate metadata",
                    "shortname": "merge k-mer count arrays for aggregate metadata (header)",
                    "description": "merge counts of nullomers, unique kmers, and total kmers."
                },
                {
                    "name": "collate the count array"
                    "shortname": "collate the count array",
                    "description": "aggregate counts from different input files"
                },
                {
                    "name": "print 'table' Final stats and close output file",
                    "shortname": "metrics and shutdown",
                    "description": "print final statistics, typically metadata values, and ensure file is closed."
                }

            ]

        })


        
        MATRIX_BANNER = """
                          name : matrix
                   description : create a k-mer array as a matrix/data-frame. tsv to stdout and/or output file.
        """
        
        MATRIX_PARAMS = OrderedDict({
            "name" "arguments",
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


        MATRIX_INPUTS = OrderedDict({
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

        MATRIX_FEATURES = OrderedDict({
            "name": "features",
            "type": "array",
            "items": [
                {
                    "name": "read/import datasets",
                    "shortname": "read/import multiple count vectors"m
                    "description" : "count vectors in columns of csv/tsv input over 'STDIN', .tsv/csv file (with --delimiter option), or multiple .kdb files."
                },
                {
                    "name": "Pandas library .csv sanitation",
                    "shortname": "validate input",
                    "description": "sanitize and structure validation on inputs (TODO: 4/5/24 not done yet. Uses a pandas csv loader on tsv or STDIN input. Otherwise collates multiple .kdb count vectors and coerces the input into a NumPy array otherwise failure.)"
                },
                {
                    "name": "transformation on data matrix",
                    "shortname": "transformation on data matrix",
                    "description": "Apply 'transformation' to input matrix. Not a distance calculation, but a restructuring of the data. Options are PCA (+scipy), t-stochastic neighbor embedding(scipy), 'DESeq2' median-of-ratios normalization, or retrieve frequencies for related applications."
                }
            ]
        })

        MATRIX_STEPS = OrderedDict({
            {
                "name": "steps",
                "type": "array",
                "items": [
                    {
                        "name": "",
                        "shortname": "",
                        "description": ""
                    },
                    {

                    }
                ]
            }
        })

    
    
        DISTANCE_BANNER = """
                          name : distance
                   description : create a distance matrix, using scipy distances, or even my Cython Pearson corr coeff (pearson).
                                  Spearman rank corr coeff also availble (vis scipy). Non-cython correlation coefficient availalble
                                  as 'correlation' distance in scipy.
        """

        DISTANCE_PARAMS = OrderedDict({
            "name" "arguments",
            "type": "array",
            "items": [
                {
                    "name" : "distance method",
                    "type": "parameter",
                    "value" : "braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,matching,minkowski,pearson,rogerstanimotorusselrao,seuclidean,sokalmichener,sokalsneath,spearman,sqeuclidean,yule distances. spearman and correlation distances provided by scipy. Pearson is a custom Cython correlation coefficient. all other provided by scipy",
                },
            ]
        })

        
        DISTANCE_INPUTS = OrderedDict({
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

        DISTANCE_FEATURES = OrderedDict({
            "name": "features",
            "type": "array",
            "items": [
                {
                    "name" 
                    "shortname"

                    "description"
            
        })
    
        # KMEANS_BANNER = """
        #                   name : kmeans
        #            description : 
        # """

        # HIERARCHICAL_BANNER = """

        # """


        # config.py features
        SUPPORTED_FEATURES = {}

        # ASCII + RICH.py checklist
        




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



#
#          logger and selected log lines
#


# LOGGER BANNER + SPACER




#               POST_LOG_HEADER
#

PROFILE_BANNER = """
=======================================================
               ||     p r o f i l e     ||
=======================================================
                    inputs  : 
                    outputs :
                 output_dir :
                    logfile :
              relevant_loc  :     # an hash of value: documentation strings and lines-of-code(loc) given a feature or error.
           relevant_issues  : 



                 parameters : {0}
-------------------------------------------------------

PARAMS_YAML:


parameters : [
 - param1 (DEFAULT: SOME_DEFAULT_VALUE, type=) : 
 - param2 (some_functionality)  :


-------------------------------------------------------

"""
HEADER_BANNER =
VIEW_BANNER =
MATRIX_BANNER =
DISTANCE_BANNER =
GRAPH_BANNER = 
KMEANS_BANNER =
HIERARCHICAL_BANNER = 


# From usage-notes (include description
#
# parameters:
#      - (short), (long), type, description, examples, help 
#      - ...
#      -
#


# ISO-8601
RUNTIME


LOGFILE

EXIT_CODE

TRACEBACK

# caught exception
LOGGABLE_LINE


# 20, 50, 100, 200 -n
LAST_N_LINES_OF_LOG



METADATA # [metadata] via YAML representer, key indices, key/values, last program stage
METADATA_VALUES_DESCRIPTIONS # Just simple key-value descriptions


METADATA_SCHEMA_TXT # Use YAML representer...
METADATA_SCHEMA_DESCRIPTION # What the loops, maps, conditionals/branches, and data structure validation mean in this data structure, how it changes, where invalid scenarios may be involved, which test datasets to use and test command.
SCHEMA_COUNT = "                      schema count: {0}".format(schemas_count)


PROGRAM_STAGE_HEADER = """
========================================================

           s t a g e : 1,2,3,...

========================================================

             name     :
         description  :

                         LOREM IPSUM
"""

TYPES_OF_ERRORS = """
              name: error1
      title   : something_descriptive_for_example
          description :

                 LOREM IPSUM

          related_issues : (N/A) 133, 132, 101
                 exit_codes : 1, 2, 3, ...etc.,
"""
















