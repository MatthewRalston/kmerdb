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

import jsonschema


from kmerdb import config, util



yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)

default_logline_choices = (20, 50, 100, 200)
PINNED_ISSUES = (132, 133, 137)

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
==============================================================
                  ||      G i t H u b     ||
==============================================================
                         Repo: kmerdb
               Feature branch: interface

Issue Tracker: https://github.com/MatthewRalston/kmerdb/issues
-------------------------------------------------------
"""



PINNED_ISSUE = """
                 Pinned issue: {0}
""".format(", ".join(list(map(str, PINNED_ISSUES))))





"""
===============================
COMMAND INFO
===============================
"""





command_1_name = "profile"
command_2_name = "make_graph"
command_3_name = "get_matrix"
command_4_name = "distance"
command_5_name = "view"
command_6_name = "header"
command_7_name = "kmeans"
command_8_name = "hierarchical"
command_9_name = "index"
command_10_name = "shuf"

COMMANDS = [
    command_1_name, #"profile",
    command_2_name, #"make_graph",
    command_3_name, #"get_matrix",
    command_4_name, #"distance",
    command_5_name, #"view",
    command_6_name, #"header",
    command_7_name, #"kmeans",
    command_8_name, #"hierarchical",
    command_9_name, #"index",
    command_10_name, #"shuf",
]

        
command_1_description = "create a k-mer count vector from fasta input. g-zipped/.bgzf   -  4 column table with YAML header."
command_1_description_long = """

      :           4 columns output table       :    [   row idx   |   k-mer id  |    count   |  frequency  ]

      :    print a k-mer count vector to STDOUT





           Finally, print a summary of k-mers processed.



>kmerdb profile -k 12 test/data/Cdifficile_R3.fa.gz output.12.kdb


...


-------------------------+
                         |
      Example:           |
                         |
-------------------------+

Final stats:


Total k-mers processed: 4093136
Unique nullomer count:   14509305
Unique 12-mer count:     2267911
Theoretical 12-mer number (4^12):     16777216
==============================



==============================
------------------------------

Output file:

==============================

version: 0.8.0
metadata_blocks: 1
k: 12
total_kmers: 4093136
unique_kmers: 2267911
unique_nullomers: 14509305
metadata: false
sorted: false
kmer_ids_dtype: uint64
profile_dtype: uint64
count_dtype: uint64
frequencies_dtype: float64
tags: []
files:
- filename: test/data/Cdifficile_R3.fasta.gz
  md5: ed09041ad0e991a8037edd0f05643571
  sha256: f2978ff92a19b82d569b10cfe3602eb37a613d6e48db22af82126799beefb1e6
  total_reads: 1
  total_kmers: 4093136
  unique_kmers: 2267911
  nullomers: 14509305

========================



                  +=============+====================+====================+=================================+
                  <    row idx  |  k-mer id          |  k-mer count       |  k-mer frequency        >                        
                  |             |                    |                    |                                 |
                  |             +
                  |
                  |
                  |
                  |
                  |




"""


command_1_parameters = "use < -k > for k-mer size, --quiet to stop STDOUT redundancy, -vv for debug level logging. -v for info."
command_1_inputs = "Input file can be < .fasta | .fastq | .fa.gz >. Last positional arg is output, a gzipped .tsv of k-mer counts and frequencies"
command_1_usage = "kmerdb profile -k $K --quiet <example_1.fa.gz> [example_2.fq.gz] <output_count_vector_file.12.kdb>"
COMMAND_1_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb profile -k 12 input_1.fa output.12.kdb

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]



""".format(command_1_name, command_1_description, command_1_description_long, command_1_inputs, command_1_parameters, command_1_usage)

COMMAND_1_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
        {
            "name"  : "k",
            "type"  : "parameter",
            "value" : "choice of k-mer size parameter 'k'.",
        },
        {
            "name" : "sorted",
            "type" : "flag",
            "value": "descending in k-mer count"
        },
        {
            "name"  : "quiet",
            "type"  : "flag",
            "value" : "write additional debug level information to stderr?"
        }
    ]
})

COMMAND_1_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name"  : "<.fasta|.fastq>",
            "type"  : "array",
            "value" : "gzipped or uncompressed input .fasta or .fastq file(s)"
        },
        {
            "name"  : "<output>.$K.kdb - a k-mer count vector in a 4 column output table. zlib compatible",
            "type"  : "file",
            "value" : "Output file. K-mer count vector in a 4 column table (index number, k-mer id, count, and frequency)."
        }
    ]
})


COMMAND_1_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "k-mer count array produced as sequences are read by sliding window approach. (Un)compressed support for .fa/.fq.",
            "shortname": "parallel OOP sliding window k-mer shredding",
            "description": "N = 4^k count-vector from one or more sequence files. (index, k-mer id, count, and frequency [as float64])"
        }),
        OrderedDict({
            "name": "k-mers are tallied, and metadata merged from the input files",
            "shortname": "merge counts, metadata from across inputs",
            "description": "a final metadata structure and output metrics are collected for display to the user."
        })

            ]
})


COMMAND_1_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "read input file(s) from filesystem into k-mer arrays",
            "shortname": "shred inputs into k-mer count arrays",
            "description": "shred input sequences into k-mer count vector",
        }),
        OrderedDict({
            "name": "collate the count array",
            "shortname": "collate the count array",
            "description": "aggregate counts from different input files"
        }),
        OrderedDict({
            "name": "merge k-mer arrays and aggregate metadata",
            "shortname": "merge k-mer count arrays for aggregate metadata (header)",
            "description": "merge counts of nullomers, unique kmers, and total kmers."
        }),
        OrderedDict({
            "name": "print 'table' Final stats and close output file",
            "shortname": "metrics and shutdown",
            "description": "print final statistics, typically metadata values, and ensure file is closed."
        })

    ]

})


command_2_description = "create a edge list in (block) .gz format from .fasta|.fa or .fastq format."
command_2_description_long = """


   :     4 column output : [ row idx | k-mer id node #1 | k-mer id node #2 | edge weight (adjacency count) ]

   :  make a deBruijn graph, count the number of k-mer adjacencies,  printing the edge list to STDOUT




                  +=============+====================+====================+=================================+
                  <    row idx  |  k-mer id node #1  |  k-mer id node #2  |  edge weight (adjacency count)  >
                  |             |                    |                    |                                 |
                  |             +
                  |
                  |
                  |
                  |
                  |
"""
        
command_2_parameters = "uses < -k > for k-mer size, --quiet to reduce runtime, -v, -vv to control logging. --"
command_2_inputs = "Input file can .fastq (or .fa).   - gzip.  Output is a weighted edge list in .kdb format (gzipped .csv with YAML header)"
command_2_usage = "kmerdb graph -k $K --quiet <input_1.fa.gz> [input_2.fq.gz] <output_edge_list_file.12.kdbg>" 


COMMAND_2_BANNER = """



                          [ name ] :         {0}

                   description : {1}

{2}




--------------------------


                    kmerdb {0} -k 12 input_1.fa [example_2.fastq] output.12.kdbg

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}











""".format(command_2_name, command_2_description, command_2_description_long, command_2_inputs, command_2_parameters, command_2_usage)

        

        
COMMAND_2_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
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
        
        

COMMAND_2_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name"  : "<.fasta|.fastq>",
            "type"  : "array",
            "value" : "gzipped or uncompressed input .fasta or .fastq file(s)"
        },
        {
            "name"  : ".kdbg",
            "type"  : "file",
            "value" : "Output edge-list filepath."
        }
    ]
})


COMMAND_2_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "k-mer count arrays, linear, produced as file is read through sliding window. (Un)compressed support for .fa/.fq.",
            "shortname": "parallel faux-OP sliding window k-mer shredding",
            "description": "Sequential k-mers from the input .fq|.fa files are added to the De Bruijn graph. In the case of secondary+ sequences in the .fa or considering NGS (.fq) data, non-adjacent k-mers are pruned with a warning. Summary statistics for the entire file are given for each file read, + a transparent data structure."
        }),
        OrderedDict({
            "name": "k-mer neighbors assessed and tallied, creates a unsorted edge list, with weights",
            "shortname": "weighted undirected graph",
            "description": "an edge list of a De Bruijn graph is generated from all k-mers in the forward direction of .fa/.fq sequences/reads. i.e. only truly neighboring k-mers in the sequence data are added to the tally of the k-mer nodes of the de Bruijn graph and the edges provided by the data."
        })
    ]
})


COMMAND_2_STEPS = OrderedDict({
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
            "name": "collate the weighted edge lists after reading multiple files. Output data consists of a edge_list, analogous metadata-header as YAML, kmer_counts, and nullomer_ids.",
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
        

command_3_description = "Aggregate k-mer count vectors into k-mer count matrix. Output is plain text .tsv"
command_3_description_long = """








...





"""
command_3_parameters = "Use PCA or tSNE for dim. reduction. Use DESeq2 for count vector sample normalization. Use 'pass' to just aggregate count vectors into .tsv"
command_3_inputs = "Several input files with the same [-k]. Or a k-mer count matrix as .tsv. Use 'STDIN' in place of an input file to read .tsv from STDIN in a Unix pipe. Prints .tsv to stdout."
command_3_usage = "kmerdb matrix <pass|PCA|tSNE|DESeq2> [input1.12.kdb input2.12.kdb ...] [OR] kmerdb matrix DESeq2 <input_count_matrix.tsv> [OR] kmerdb matrix PCA -k 5 STDIN"

        
COMMAND_3_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb matrix <pass|PCA|tSNE|DESeq2> input1.12.kdb input2.12.kdb ...

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_3_name, command_3_description, command_3_description_long, command_3_inputs, command_3_parameters, command_3_usage)
        
        
COMMAND_3_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
        {
            "name" : "function",
            "type": "parameter",
            "value" : "dimensionality reduction, DESeq2 normalization (via rpy2), or pass - no-transformation, or (deprecated) k-mer Frequency matrix. [PCA, tSNE, DESeq2, pass, Frequency]",
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


COMMAND_3_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "kdbfiles",
            "type": "file(s)",
            "value": ".kdb count vector files produced from command #1 in the profile <-> k-mer count matrix <-> count distance matrix pipeline"
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

COMMAND_3_FEATURES = OrderedDict({
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
            "description": "sanitize and structure validation on inputs (TODO: 4/5/24 not done yet. Uses a pandas csv loader on tsv/csv files or 'STDIN' to read a matrix/DataFrame input from a standard input stream. Otherwise collates multiple .kdb count vectors and coerces the input into a NumPy array otherwise failure.)"
        }),
        OrderedDict({
            "name": "transformation on data matrix",
            "shortname": "transformation on data matrix",
            "description": "Apply 'transformation' to input matrix. Not a distance calculation, but a restructuring of the data. Options are PCA (+scipy), t-stochastic neighbor embedding(scipy), 'DESeq2' median-of-ratios normalization, or retrieve frequencies for related applications."
        })
    ]
})

COMMAND_3_STEPS = OrderedDict({
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

                    }),
        OrderedDict({
            "name": "pass unnormalized counts",
            "shortname": "unnormalized data",
            "description": "Passes the unnormalized data in Data Frame tsv format to STDOUT/file"
        })

    ]

})



command_4_description = "Create distance matrix from input count matrix (.tsv). Output is plain text .tsv"
command_4_description_long = """








...





"""
command_4_parameters = "Choose a distance metric. 23 distances available. 'pearson' 'correlation' and 'spearman' are relevant for k-mer count vector similarity."
command_4_inputs = "Several input files with the same [-k]. Or a k-mer count matrix as .tsv. Use 'STDIN' in place of an input file to read .tsv from STDIN in a Unix pipe. Prints .tsv to stdout."
command_4_usage = "kmerdb distance  [input1.12.kdb input2.12.kdb ...] [OR] kmerdb distance correlation <input.tsv> [OR] kmerdb distance spearman STDIN"

        
COMMAND_4_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb distance pearson [input1.12.kdb input2.12.kdb ...]  OR kmerdb distance pearson <input.tsv> or kmerdb distance pearson STDIN

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_4_name, command_4_description, command_4_description_long, command_4_inputs, command_4_parameters, command_4_usage)

        

COMMAND_4_PARAMS = OrderedDict({
    "name": "arguments",
            "type": "array",
    "items": [
        {
            "name" : "distance method",
            "type": "parameter",
            "value" : "braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,matching,minkowski,pearson,rogerstanimotorusselrao,seuclidean,sokalmichener,sokalsneath,spearman,sqeuclidean,yule distances. 'spearman' and 'correlation' distances provided by scipy. 'pearson' is a custom Cython correlation coefficient: ssxy/sqrt(ssxx^2 x ssyy^2) | all others provided by scipy",
        },
    ]
})

        
COMMAND_4_INPUTS = OrderedDict({
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

COMMAND_4_FEATURES = OrderedDict({
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
            "description": "Choose 'distance method' to apply to input matrix. Available SciPy distances include 'braycurtis,canberra,chebyshev,cityblock,correlation,cosine,dice,euclidean,hamming,jaccard,jensenshannon,kulsinski,mahalanobis,matching,minkowski,pearson,rogerstanimotorusselrao,seuclidean,sokalmichener,sokalsneath,spearman,sqeuclidean,yule' distances. 'spearman' and 'correlation' distances provided by scipy. 'pearson' is a custom Cython correlation coefficient: ssxy/sqrt(ssxx^2 x ssyy^2) | all others provided by scipy. All distances availalable on -h|--help"
        })
    ]
})
    
COMMAND_4_STEPS = OrderedDict({
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

















command_5_description = "Decompress and view contents of .kdb(g) file "
command_5_description_long = """








...





"""
command_5_parameters = "Use -H|--header to additionally print header."
command_5_inputs = "Input may be a v{0} count vector (.kdb) or edge list (.kdbg) file.".format(config.VERSION)
command_5_usage = "kmerdb view -H <input_1.12.kdb>"

        
COMMAND_5_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb view -H <input_1.12.kdb>

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_5_name, command_5_description, command_5_description_long, command_5_inputs, command_5_parameters, command_5_usage)

        

COMMAND_5_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
        {
            "name" : "header flag",
            "type": "parameter",
            "value" : "Use -H|--header to print YAML metadata header to STDOUT before table section"
        },
    ]
})

        
COMMAND_5_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "k-mer database file (.kdb)",
            "type": "file",
            "value": ".kdb count vector files produced from kmerdb profile"
        }
        
            ]
})

COMMAND_5_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "read inputs: parse and validate input file's header and table",
            "shortname": "read input data into count vectors",
            "description" : "parse input data, a custom format, and load data into memory"
        }),
        OrderedDict({
            "name": "print file contents",
            "shortname": "print to STDOUT",
            "description": ""
        })
    ]
})
    
COMMAND_5_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "read input file from filesystem into k-mer count vector",
            "shortname": "read .kdb file into k-mer count vector",
            "description": "read a single .kdb file into a k-mer count vector",
        }),
        OrderedDict({
            "name": "print file contents to STDOUT",
            "shortname": "decompress and print input file",
            "description": ""
        })

    ]

})




command_6_description = "Show YAML metadata header from .kdb(g) file"
command_6_description_long = ""
command_6_parameters = "N/A"
command_6_inputs = "Input may be a v{0} count vector (.kdb) or edge list (.kdbg) file.".format(config.VERSION)
command_6_usage = "kmerdb header <input_1.12.kdb>"

        
COMMAND_6_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb header <input_1.12.kdb>

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_6_name, command_6_description, command_6_description_long, command_6_inputs, command_6_parameters, command_6_usage)

        

COMMAND_6_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
        {
            "name" : "json",
            "type": "parameter",
            "value" : "Print results in JSON format"
        },
    ]
})

        
COMMAND_6_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "k-mer database file (.kdb)",
            "type": "file",
            "value": ".kdb count vector file produced from kmerdb profile"
        }
        
            ]
})

COMMAND_6_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "read inputs: parse and validate input file's header and table",
            "shortname": "read input data into count vectors",
            "description" : "parse input data, a custom format, and load data into memory"
        }),
        OrderedDict({
            "name": "print YAML metadata header",
            "shortname": "print to STDOUT",
            "description": ""
        })
    ]
})
    
COMMAND_6_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "read input file from filesystem into k-mer count vector",
            "shortname": "read .kdb file into k-mer count vector",
            "description": "read a single .kdb file into a k-mer count vector",
        }),
        OrderedDict({
            "name": "print file contents to STDOUT",
            "shortname": "decompress and print input file",
            "description": ""
        })

    ]

})


















command_7_description = "K-means clustering with biopython or scikit-learn"
command_7_description_long = "Produces eblow graph if k is not determined a-priori. Uses matplotlib for graphics."
command_7_parameters = "Use -k to control clustering with k-means. Choice of 'sklearn' or 'biopython' k-means clustering. Choice of distance metrics offered by kcluster"
command_7_inputs = "Input is a .tsv file. Use STDIN to read input from standard input".format(config.VERSION)
command_7_usage = "kmerdb kmeans -k 5 -i <count_or_distance_matrix.tsv> sklearn"

        
COMMAND_7_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb kmeans -k 5 <count_matrix.tsv> biopython

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_7_name, command_7_description, command_7_description_long, command_7_inputs, command_7_parameters, command_7_usage)

        

COMMAND_7_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
        {
            "name" : "k",
            "type": "parameter",
            "value" : "Controls k-means clustering"
        },
    ]
})

        
COMMAND_7_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "Input data matrix",
            "type": "file",
            "value": "Input may be a count matrix or distance/correlation matrix."
        }
        
            ]
})

COMMAND_7_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse input count matrix or distance/correlation matrix",
            "shortname": "read input data matrix",
            "description" : "load data into memory from input file"
        }),
        OrderedDict({
            "name": "Use biopython or scikit-learn for k-means clutsering",
            "shortname": "kmeans clustering",
            "description": "creates elbow graph if choice of k is ambiguous, and other output graphics from distance matrix"
        })
    ]
})
    
COMMAND_7_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse .csv file",
            "shortname": "Loads data into Pandas DataFrame",
            "description": "Infers datatype from .csv data",
        }),
        OrderedDict({
            "name": "Run k-means clustering",
            "shortname": "k-means clustering",
            "description": "Uses biopython or scikit-learn for k-means clustering, and matplotlib for graphics"
        })

    ]

})





command_8_description = "Hierarchical clustering with biopython"
command_8_description_long = "Uses matplotlib for graphics."
command_8_parameters = "-m|--method determines linkage fitting (see --help for details)"
command_8_inputs = "Input is a .tsv file. Use STDIN to read input from standard input".format(config.VERSION)
command_8_usage = "kmerdb hierarchical -m ward -i <count_or_distance_matrix.tsv>"

        
COMMAND_8_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb kmeans -k 5 <count_matrix.tsv> biopython

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_8_name, command_8_description, command_8_description_long, command_8_inputs, command_8_parameters, command_8_usage)

        

COMMAND_8_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
        {
            "name" : "method",
            "type": "parameter",
            "value" : "Linkage method for hierarchical clustering. Choices are single,complete,average,weighted,centroid,median,ward"
        },
    ]
})

        
COMMAND_8_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "Input data matrix",
            "type": "file",
            "value": "Input may be a count matrix or distance/correlation matrix."
        }
        
    ]
})

COMMAND_8_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse input count matrix or distance/correlation matrix",
            "shortname": "read input data matrix",
            "description" : "load data into memory from input file"
        }),
        OrderedDict({
            "name": "Use biopython for hierarchical clutsering",
            "shortname": "hierarchical clustering",
            "description": "Uses Bio.Cluster.treecluster for hierarchical clustering. Produces a phylip UPGMA tree."
        })
    ]
})
    
COMMAND_8_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse .csv file",
            "shortname": "Loads data into Pandas DataFrame",
            "description": "Infers datatype from .csv data",
        }),
        OrderedDict({
            "name": "Run hierarchical clustering",
            "shortname": "hierarchical clustering",
            "description": "Uses biopython hierarchical clustering. Produces a UPGMA tree in Phylip format"
        })

    ]

})























command_9_description = "Create an index file (deprecated)"
command_9_description_long = ""
command_9_parameters = "N/A"
command_9_inputs = "Input is a v{0} .kdb count vector file".format(config.VERSION)
command_9_usage = "kmerdb index <kmer_count_vector.kdb>"

        
COMMAND_9_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb index <input.kdb>

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_9_name, command_9_description, command_9_description_long, command_9_inputs, command_9_parameters, command_9_usage)

        

COMMAND_9_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
    ]
})

        
COMMAND_9_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "Input k-mer count vector (.kdb) file",
            "type": "file",
            "value": "File to index."
        }
        
            ]
})

COMMAND_9_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse input k-mer count vector .kdb file",
            "shortname": "read input file",
            "description" : "load data into memory from input file"
        }),
        OrderedDict({
            "name": "Create index file",
            "shortname": "index input file",
            "description": "Create a plaintext index file for on-disk operations when memory may be an issue. (Deprecated)"
        })
    ]
})
    
COMMAND_9_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse input k-mer db .kdb file",
            "shortname": "Loads data from input file",
            "description": "Validate and parse input data",
        }),
        OrderedDict({
            "name": "Create index file",
            "shortname": "Create index file",
            "description": "Create a plaintext index file (.kdbi) for on disk operations (Deprecated)"
        })

    ]

})


        

command_10_description = "Shuffle k-mer count vector (Deprecated)"
command_10_description_long = ""
command_10_parameters = "N/A"
command_10_inputs = "Input is a v{0} .kdb count vector file".format(config.VERSION)
command_10_usage = "kmerdb shuf <kmer_count_vector.kdb>"

        
COMMAND_10_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb shuf <input.kdb>

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_10_name, command_10_description, command_10_description_long, command_10_inputs, command_10_parameters, command_10_usage)

        

COMMAND_10_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
    ]
})

        
COMMAND_10_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "Input k-mer count vector (.kdb) file",
            "type": "file",
            "value": "File to shuffle."
        }
        
    ]
})

COMMAND_10_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse input k-mer count vector .kdb file",
            "shortname": "read input file",
            "description" : "load data into memory from input file"
        }),
        OrderedDict({
            "name": "Shuffle k-mer count vector",
            "shortname": "Shuffle k-mer counts",
            "description": "(Deprecated)"
        })
    ]
})
    
COMMAND_10_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse input k-mer db .kdb file",
            "shortname": "Loads data from input file",
            "description": "Validate and parse input data",
        }),
        OrderedDict({
            "name": "Shuffle k-mer count vector",
            "shortname": "Shuffle k-mer counts",
            "description": "(Deprecated)"
        })

    ]

})


class kmerdb_appmap:



        


    
    def __init__(self, argv, logger=None):


        if logger is not None:
            self.logger = logger
            self.logfile = logger.logfile
        else:
            self.logger = None
            self.logfile = None
        self._loggable = logger is not None
        
        
        self.MODULE_ROOT = os.path.join("..", os.path.dirname(__file__))
        self.COMMAND_FILE = os.path.join(self.MODULE_ROOT, "__init__.py")
        self.PACKAGE_MANAGER = """
                      package manger : pip
                        version      : >= 24.0
        package root : {0}
        exe file     : {1}

                      required packages : {2}
                   development packages : {3}

           ARGV : {4}
        """.format(self.MODULE_ROOT, self.COMMAND_FILE, config.requirements_count, config.requirements_dev_count, argv)
        self.REQUIRES_PYTHON = config.REQUIRES_PYTHON
        self.VERSION_HARDCODED = "                                                                          v :      >= v{0}\n".format(config.REQUIRES_PYTHON)


        
        #
        # loaded_modules
        #

        #
        # dependencies
        #
        #      required:
        #      optional:

        # usage_notes.txt





        


        self.ALL_PARAMS = {
            "profile": COMMAND_1_PARAMS["items"],
            "make_graph": COMMAND_2_PARAMS["items"],
            "get_matrix": COMMAND_3_PARAMS["items"],
            "distance": COMMAND_4_PARAMS["items"],
            "view": COMMAND_5_PARAMS["items"],
            "header": COMMAND_6_PARAMS["items"],
            "kmeans": COMMAND_7_PARAMS["items"],
            "hierarchical": COMMAND_8_PARAMS["items"],
            "index": COMMAND_9_PARAMS["items"],
            "shuf": COMMAND_10_PARAMS["items"]

        }

        ALL_INPUTS = {
            "profile": COMMAND_1_INPUTS["items"],
            "make_graph": COMMAND_2_INPUTS["items"],
            "get_matrix": COMMAND_3_INPUTS["items"],
            "distance": COMMAND_4_INPUTS["items"],
            "view": COMMAND_5_INPUTS["items"],
            "header": COMMAND_6_INPUTS["items"],
            "kmeans": COMMAND_7_INPUTS["items"],
            "hierarchical": COMMAND_8_INPUTS["items"],
            "index": COMMAND_9_INPUTS["items"],
            "shuf": COMMAND_10_INPUTS["items"]
        }

        ALL_FEATURES = {
            "profile": COMMAND_1_FEATURES["items"],
            "make_graph": COMMAND_2_FEATURES["items"],
            "get_matrix": COMMAND_3_FEATURES["items"],
            "distance": COMMAND_4_FEATURES["items"],
            "view": COMMAND_5_FEATURES["items"],
            "header": COMMAND_6_FEATURES["items"],
            "kmeans": COMMAND_7_FEATURES["items"],
            "hierarchical": COMMAND_8_FEATURES["items"],
            "index": COMMAND_9_FEATURES["items"],
            "shuf": COMMAND_10_FEATURES["items"]


        }


        ALL_STEPS = {

            "profile": COMMAND_1_STEPS["items"],
            "make_graph": COMMAND_2_STEPS["items"],
            "get_matrix": COMMAND_3_STEPS["items"],
            "distance": COMMAND_4_STEPS["items"],
            "view": COMMAND_5_STEPS["items"],
            "header": COMMAND_6_STEPS["items"],
            "kmeans": COMMAND_7_STEPS["items"],
            "hierarchical": COMMAND_8_STEPS["items"],
            "index": COMMAND_9_STEPS["items"],
            "shuf": COMMAND_10_STEPS["items"]
        }
        
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
        """
        Only used in 'usage' method. Spacer at bottom is permissible because this may be used externally to usage.
        """
        
        sys.stderr.write(PROGRAM_BANNER)
        sys.stderr.write(INTERPRETER)
        sys.stderr.write(self.VERSION_HARDCODED)
        sys.stderr.write(self.PACKAGE_MANAGER)

        # Spacer
        sys.stderr.write(DNA_COLUMN_1)

        sys.stderr.write(THREE_LINES)

    def print_graph_header(self):
        """
        'kmerdb usage graph'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.

        """
        
        sys.stderr.write(COMMAND_2_BANNER)

        sys.stderr.write(THREE_LINES)
        
        sys.stderr.write(yaml.dump(COMMAND_2_PARAMS))
        
        sys.stderr.write(THREE_LINES)
        
        sys.stderr.write(yaml.dump(COMMAND_2_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_2_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_2_STEPS))



    def print_profile_header(self):
        """
        'kmerdb usage profile'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.

        """
        sys.stderr.write(COMMAND_1_BANNER)

        sys.stderr.write(THREE_LINES)
            
        sys.stderr.write(yaml.dump(COMMAND_1_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_1_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_1_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_1_STEPS))



    def print_matrix_header(self):
        """
        'kmerdb usage matrix'


        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.
        """

        sys.stderr.write(COMMAND_3_BANNER)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_3_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_3_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_3_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_3_STEPS))


    def print_distance_header(self):
        """

        'kmerdb usage distance'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.
        """

        sys.stderr.write(COMMAND_4_BANNER)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_4_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_4_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_4_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_4_STEPS))





    def print_view_header(self):
        """

        

        'kmerdb usage header'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.


        """

        sys.stderr.write(COMMAND_5_BANNER)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_5_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_5_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_5_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_5_STEPS))


    def print_header_header(self):
        """

        

        'kmerdb usage header'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.


        """

        sys.stderr.write(COMMAND_6_BANNER)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_6_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_6_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_6_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_6_STEPS))



    def print_kmeans_header(self):
        """

        

        'kmerdb usage kmeans'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.


        """

        sys.stderr.write(COMMAND_7_BANNER)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_7_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_7_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_7_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_7_STEPS))


    def print_hierarchical_header(self):
        """

        

        'kmerdb usage hierarchical'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.


        """

        sys.stderr.write(COMMAND_8_BANNER)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_8_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_8_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_8_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_8_STEPS))



        
    def print_index_header(self):
        """

        

        'kmerdb usage index'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.


        """

        sys.stderr.write(COMMAND_9_BANNER)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_9_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_9_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_9_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_9_STEPS))




        
    def print_shuf_header(self):
        """

        

        'kmerdb usage shuf'

        Writes input files, parameters, steps, and features as verbose --help on the parameters effects on the runtime.


        """

        sys.stderr.write(COMMAND_10_BANNER)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_10_PARAMS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_10_INPUTS))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_10_FEATURES))

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(yaml.dump(COMMAND_10_STEPS))



        
        
    def print_github_block(self):
        """
        Github repo block. May be useful to some users when browsing logilfe.
        """


        sys.stderr.write(THREE_LINES)
        
        sys.stderr.write(DNA_SPACER_1)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(GITHUB_LOGO)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(GITHUB_PROJECT_BANNER)

        sys.stderr.write(PINNED_ISSUE)

        sys.stderr.write(THREE_LINES)

        sys.stderr.write(DNA_SPACER_1)

        sys.stderr.write(THREE_LINES)

        
        

    def exit_gracefully(self, e:Exception, subcommand:str=None, step:int=None, feature:int=None, logs:list=None, n_logs:int=None):
        """
        We need to handle exit gracefully. The 'step' and 'feature' categories/flags/ints are passed from __init__ or down its callstack to 
        """
        import traceback
        

        if e is None:
            raise ValueError("Need an error to exit")
        elif not isinstance(e, Exception):
            raise ValueError("Need an error to exit")

        
        if n_logs is None or type(n_logs) is not int:
            raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument n_logs to be a int")
        elif logs is None or type(logs) is not list:
            raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument logs to be a list")
        elif feature is not None and type(feature) is not int:
            raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument feature to be a int")
        elif step is not None and type(step) is not int:
            raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument step to be a int")
        elif subcommand is not None and type(subcommand) is not str:
            raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument subcommand to be a str")



        
        N = len(logs) 
        loggable_line = N
        assert subcommand in config.subcommands, "Unknown subcommand"




        

        # This is the "Error blocks" metadata
        exit_summary = OrderedDict({
            "subcommand": subcommand,
            "kmerdb-version": config.VERSION,
            "python-version": config.REQUIRES_PYTHON,
            "feature": feature,
            "feature_name": self.ALL_FEATURES[subcommand][feature]["name"],
            "feature_shortname": self.ALL_FEATURES[subcommand][feature]["shortname"],
            "feature_description": self.ALL_FEATURES[subcommand][feature]["description"],
            "step" : step,
            "step_name": self.ALL_STEPS[subcommand][step]["name"],
            "step_shortname": self.ALL_STEPS[subcommand][step]["shortname"],
            "step_description": self.ALL_STEPS[subcommand][step]["description"],
            # The *total* number of logged lines produced by the program and returned to the global 'logs' var in __init__.py
            "log_file": self.logfile,
            "traceback": str(traceback.extract_tb(e.__traceback__)),
            "last_logged_line": loggable_line, 
            "error": e.__str__(),
        })


        yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
        

        try:
            jsonschema.validate(instance=exit_summary, schema=config.exit_summary_schema)
        except jsonschema.ValidationError as e:
            sys.stderr.write("Failed to validate the exit summary. Internal Error.\n")
            raise e
                







        self.print_github_block()

        """
        Print last n lines of log
        """
        for i in range(n_logs):

            try:
                if self._loggable:
                    sys.stderr.write("{0} - last line of log\n".format(n_logs - i))
                    self.logger.log_it(logs[i], "ERROR")
                else:
                    sys.stderr.write("{0} - last line of log\n".format(n_logs - i))
                    sys.stderr.write(logs[i], "ERROR")
            except Exception as e:
                raise e
                
        sys.stderr.write("-" * 80 + "\n")
        if self._loggable:
            
            self.logger.log_it("...displaying last {0} lines of the log. Please see '{1}' for more details...".format(n_logs, self.logfile), "ERROR")
            self.logger.log_it(e.__str__())

            sys.stderr.write(THREE_LINES)

            sys.stderr.write(DNA_SPACER_1)

            sys.stderr.write(THREE_LINES)

            self.logger.log_it("="*40 + "\n", "ERROR")

            self.logger.log_it(" "*20 + "ERROR: Program exit summary:\n", "ERROR")
            
            self.logger.log_it("="*40 + "\n", "ERROR")
            
            self.logger.log_it("\n" + yaml.dump(exit_summary), "ERROR")

            self.logger.log_it("="*40 + "\n")

            sys.stderr.write(THREE_LINES)
            
            
        else:
            sys.stderr.write("...displaying last {0} lines of the log. Please see '{1}' for more details...\n".format(n_logs, self.logfile))
            sys.stderr.write(e.__str__())

            sys.stderr.write(THREE_LINES)

            sys.stderr.write(DNA_SPACER_1)

            sys.stderr.write(THREE_LINES)

            sys.stderr.write("="*40 + "\n\n")

            sys.stderr.write(" "*20 + "ERROR: Program exit summary:\n\n")
            
            sys.stderr.write("="*40 + "\n")
            
            sys.stderr.write("\n" + yaml.dump(exit_summary) + "\n")

            sys.stderr.write("="*40 + "\n")

            sys.stderr.write(THREE_LINES)



            
        return exit_summary


# #
# #          logger and selected log lines
# #

















