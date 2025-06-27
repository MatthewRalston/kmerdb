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
PINNED_ISSUES = (140, 143, 150, 152, 153,)

PROGRAM_BANNER = """


 o-O      |||
o---O     |||             [|[          kmerdb           ]|]
O---o     |||
 O-o      |||        version :     v{0}
  O       |||
 o-O      |||        GitHub  : https://github.com/MatthewRalston/kmerdb/issues
o---O     |||         PyPI   : https://pypi.org/project/kmerdb/
O---o     |||      Website   : https://matthewralston.github.io/kmerdb
 O-o      |||









# |||||||||||||||||||||||||||||||||||||||
#      [ Usage ] :        |||||||||||||||
# |||||||||||||||||||||||||||||||||||||||


#  Check test/data for example fasta files.


# -----------------------------
# Generate k-mer count profiles
# -----------------------------
kmerdb profile -k 12 -o profile_1 input_1.fa.gz [input_2.fq] ...
...


# -----------------------------
# Merge profiles
# -----------------------------
kmerdb matrix from profile_1.12.kdb profile_2.12.kdb profile_3.12.kdb ... > count_matrix.tsv

# -----------------------------
# Generate inter-profile distances
# -----------------------------
kmerdb distance pearson count_matrix.tsv


# -----------------------------
# Pipeline form
# -----------------------------

kmerdb matrix from ... | kmerdb distance pearson STDIN > pearson_correlation_matrix.tsv

# -----------------------------
# Okay, how about PCA, t-SNE?
# -----------------------------
kmerdb matrix PCA -n 25 ... > 25_pseudomer_count_profiles.tsv
# kmerdb matrix tSNE -n 25


# -----------------------------
# k-means clustering?
# -----------------------------
kmerdb kmeans <Biopython|sklearn> -k 4 --distance e -i input_distance_matrix.tsv # Using 'e' for Euclidean distance with Biopython. Check the source, Biopython RTD, and sklearn RTD.
# Produces

# -----------------------------
# okay, now straight-up hierarchical clustering:
# -----------------------------
kmerdb hierarchical -i input_distance_matrix.tsv --method complete # Uses complete linkage



""".format(config.VERSION)
INTERPRETER = "                                                                       lang :         python\n"
# hardcoded

# print_program_header
#        sys.stderr.write(PROGRAM_BANNER)
#        sys.stderr.write(INTERPRETER)
#        sys.stderr.write(self.VERSION_HARDCODED)
#        sys.stderr.write(self.PACKAGE_MANAGER)



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
               Feature branch: main

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





command_1_name, command_2_name, command_3_name, command_4_name, command_5_name, command_6_name, command_7_name, command_8_name, command_9_name, command_10_name, command_11_name, command_12_name, command_13_name, command_14_name, command_15_name, command_16_name, command_17_name = config.subcommands

COMMANDS = [
    command_1_name, #"profile",
    command_2_name, #"view",
    command_3_name, #"header",
    command_4_name, #"matrix",
    command_5_name, #"distance",
    command_6_name, #"kmeans",
    command_7_name, #"hierarchical",
    command_8_name, #"codons",
    command_9_name, #"CUB",
    command_10_name, #"graph",
    command_11_name, #"get_minimizers",
    command_12_name, #"get_alignments",
    command_13_name, #"index", #
    command_14_name, #"shuf",
    command_15_name, #"usage",
    command_16_name, #"help",
    command_17_name, #"version"]
]




        
command_1_description = "create a k-mer count vector from fasta/fastq input. g-zipped/.bgzf   -  4 column table with YAML header."
command_1_description_long = """

      :           4 columns output table       :    [   row idx   |   k-mer id  |    count   |  frequency  ]

      :    print a k-mer count vector to STDOUT





           Finally, print a summary of k-mers processed.



>kmerdb profile -k 12 -o profile_1 test/data/Cdifficile_R3.fa.gz


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
command_1_usage = "kmerdb profile -k $K --quiet --debug -vv -o output_profile1 <example_1.fa.gz> [example_2.fq.gz]"
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
            "description": "aggregate counts from different sequences within a single file"
        }),
        OrderedDict({
            "name": "merge k-mer arrays and aggregate metadata",
            "shortname": "merge k-mer count arrays and metadata (header) across all files",
            "description": "merge counts of nullomers, unique kmers, and total kmers across all files."
        }),
        OrderedDict({
            "name": "print 'table' Final stats and close output file",
            "shortname": "metrics and shutdown",
            "description": "print final statistics, typically metadata values, and ensure file is closed."
        })

    ]

})


command_2_description = "Decompress and view contents of .kdb(g) file "
command_2_description_long = """








...





"""
command_2_parameters = "Use -H|--header to additionally print header."
command_2_inputs = "Input may be a v{0} count vector (.kdb) or edge list (.kdbg) file.".format(config.VERSION)
command_2_usage = "kmerdb view -H <input_1.12.kdb>"

        
COMMAND_2_BANNER = """












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
""".format(command_2_name, command_2_description, command_2_description_long, command_2_inputs, command_2_parameters, command_2_usage)

        

COMMAND_2_PARAMS = OrderedDict({
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

        
COMMAND_2_INPUTS = OrderedDict({
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

COMMAND_2_FEATURES = OrderedDict({
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
    
COMMAND_2_STEPS = OrderedDict({
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


command_3_description = "Show YAML metadata header from .kdb(g) file"
command_3_description_long = ""
command_3_parameters = "N/A"
command_3_inputs = "Input may be a v{0} count vector (.kdb) or edge list (.kdbg) file.".format(config.VERSION)
command_3_usage = "kmerdb header <input_1.12.kdb>"

        
COMMAND_3_BANNER = """












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
""".format(command_3_name, command_3_description, command_3_description_long, command_3_inputs, command_3_parameters, command_3_usage)

        

COMMAND_3_PARAMS = OrderedDict({
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

        
COMMAND_3_INPUTS = OrderedDict({
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

COMMAND_3_FEATURES = OrderedDict({
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
    
COMMAND_3_STEPS = OrderedDict({
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


command_4_description = "Aggregate k-mer count vectors into k-mer count matrix. Output is plain text .tsv"
command_4_description_long = """








...





"""
command_4_parameters = "Use PCA or tSNE for dim. reduction. Use DESeq2 for count vector sample normalization. Use 'pass' to just aggregate count vectors into .tsv"
command_4_inputs = "Several input files with the same [-k]. Or a k-mer count matrix as .tsv. Use 'STDIN' in place of an input file to read .tsv from STDIN in a Unix pipe. Prints .tsv to stdout."
command_4_usage = "kmerdb matrix <from|PCA|tSNE|DESeq2> [input1.12.kdb input2.12.kdb ...] [OR] kmerdb matrix DESeq2 <input_count_matrix.tsv> [OR] kmerdb matrix PCA -k 5 STDIN"

        
COMMAND_4_BANNER = """












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
""".format(command_4_name, command_4_description, command_4_description_long, command_4_inputs, command_4_parameters, command_4_usage)
        
        
COMMAND_4_PARAMS = OrderedDict({
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


COMMAND_4_INPUTS = OrderedDict({
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
            "description": "sanitize and structure validation on inputs (TODO: 4/5/24 not done yet. Uses a pandas csv loader on tsv/csv files or 'STDIN' to read a matrix/DataFrame input from a standard input stream. Otherwise collates multiple .kdb count vectors and coerces the input into a NumPy array otherwise failure.)"
        }),
        OrderedDict({
            "name": "transformation on data matrix",
            "shortname": "transformation on data matrix",
            "description": "Apply 'transformation' to input matrix. Not a distance calculation, but a restructuring of the data. Options are PCA (+scipy), t-stochastic neighbor embedding(scipy), 'DESeq2' median-of-ratios normalization, or retrieve frequencies for related applications."
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



command_5_description = "Create distance matrix from input count matrix (.tsv). Output is plain text .tsv"
command_5_description_long = """








...





"""
command_5_parameters = "Choose a distance metric. 23 distances available. 'pearson' 'correlation' and 'spearman' are relevant for k-mer count vector similarity."
command_5_inputs = "Several input files with the same [-k]. Or a k-mer count matrix as .tsv. Use 'STDIN' in place of an input file to read .tsv from STDIN in a Unix pipe. Prints .tsv to stdout."
command_5_usage = "kmerdb distance  [input1.12.kdb input2.12.kdb ...] [OR] kmerdb distance correlation <input.tsv> [OR] kmerdb distance spearman STDIN"

        
COMMAND_5_BANNER = """












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
""".format(command_5_name, command_5_description, command_5_description_long, command_5_inputs, command_5_parameters, command_5_usage)

        

COMMAND_5_PARAMS = OrderedDict({
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

        
COMMAND_5_INPUTS = OrderedDict({
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

COMMAND_5_FEATURES = OrderedDict({
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
    
COMMAND_5_STEPS = OrderedDict({
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



command_6_description = "K-means clustering with biopython or scikit-learn"
command_6_description_long = "Produces eblow graph if k is not determined a-priori. Uses matplotlib for graphics. Prints {0} and {1} as primary outputs, in addition to the 'kmeans_elbow_graph.png' which is printed if no k is supplied.".format(config.pca_variance_fig_filepath, config.kmeans_clustering_fig_filepath)
command_6_parameters = "Use -k to control clustering with k-means. Choice of 'sklearn' or 'biopython' k-means clustering. Choice of distance metrics offered by kcluster"
command_6_inputs = "Input is a .tsv file. Use STDIN to read input from standard input".format(config.VERSION)
command_6_usage = "kmerdb kmeans -k 5 -i <count_or_distance_matrix.tsv> sklearn"

        
COMMAND_6_BANNER = """












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
""".format(command_6_name, command_6_description, command_6_description_long, command_6_inputs, command_6_parameters, command_6_usage)

        

COMMAND_6_PARAMS = OrderedDict({
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

        
COMMAND_6_INPUTS = OrderedDict({
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

COMMAND_6_FEATURES = OrderedDict({
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
    
COMMAND_6_STEPS = OrderedDict({
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





command_7_description = "Hierarchical clustering with biopython"
command_7_description_long = "Uses matplotlib for graphics. Creates {0} and {1} as primary outputs.".format(config.hierarchical_clustering_dendrogram_fig_filepath, config.upgma_tree_phylip)
command_7_parameters = "-m|--method determines linkage fitting (see --help for details)"
command_7_inputs = "Input is a .tsv file. Use STDIN to read input from standard input".format(config.VERSION)
command_7_usage = "kmerdb hierarchical -m ward -i <count_or_distance_matrix.tsv>"

        
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
            "name" : "method",
            "type": "parameter",
            "value" : "Linkage method for hierarchical clustering. Choices are single,complete,average,weighted,centroid,median,ward"
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
            "name": "Use biopython for hierarchical clutsering",
            "shortname": "hierarchical clustering",
            "description": "Uses Bio.Cluster.treecluster for hierarchical clustering. Produces a phylip UPGMA tree."
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
            "name": "Run hierarchical clustering",
            "shortname": "hierarchical clustering",
            "description": "Uses biopython hierarchical clustering. Produces a UPGMA tree in Phylip format"
        })

    ]

})

"""

Command 8 codons
"""
command_8_description = "Codon usage table"
command_8_description_long = "Calculate codon usage counts (or frequencies) for codons in a reference sequence population"
command_8_parameters = "--as-frequencies changes counts to frequencies within a synonymous amino-acid family of codons, --ignore-noncanonicals omits sequences that do not contain standard CDS definition (non-ATG start codon, non-TAA/TAG/TGA stop codons), --ignore-invalid-cds omits sequences that are considered invalid (length not divisible by 3, otherwise uses --ignore-noncanonicals definition). --dont-ignore-start-codons includes counts from start codons for downstream analysis, --dont-ignore-stop-codons includes counts from stop codons for downstream analysis"
command_8_inputs = "Input is a .fna fasta file of CDS sequences"
command_8_usage = "kmerdb codons -vv input.fna"

        
COMMAND_8_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb codons <input.fna>

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
            "name" : "--as-frequencies",
            "type": "parameter",
            "value" : "Output is frequencies of codon usage within an synonymous amino-acid family"
        },
        {
            "name": "--ignore-noncanonicals",
            "type": "parameter",
            "value": "Ignore (dont throw error) sequences with non-canonical start/stop codons (more permissive)",
        },
        {
            "name": "--ignore-invalid-cds",
            "type": "parameter",
            "value": "Ignore (dont throw error) sequences that do not fit the canonical definition of a CDS: length divisible by 3, standard start/stop codons (If you want the program to error out on sequences that are considered invalid then dont use this parameter. Otherwise this parameter will just ignore and omit those sequences in the final table)"
        },
        {
            "name": "--include-start-codons",
            "type": "parameter",
            "value": "Include start codon counts in the count/frequency table"
        },
        {
            "name": "--include-stop-codons",
            "type": "parameter",
            "value": "Include stop codon counts in the count/frequency table"
        },
        {
            "name": "--output-delimiter",
            "type": "parameter",
            "value": "The choice of output delimiter for the codon count/frequency table"
        }
    ]
})

        
COMMAND_8_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "Reference .fna fasta file",
            "type": "file",
            "value": "Input is a .fna file of CDS sequences."
        }
        
    ]
})

COMMAND_8_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Read fasta input",
            "shortname": "Parse fasta file",
            "description" : "Parse CDS sequences from a .fna fasta file"
        }),
        OrderedDict({
            "name": "Calculate codon counts",
            "shortname": "codon counts",
            "description": "Calculate counts of codons (3-mers) along the CDS sequence "
        }),
        OrderedDict({
            "name": "Convert to synonymous codon frequencies",
            "shortname": "synonymous frequencies",
            "description": "Will calculate the frequencies of a specific codons usage within a synonymous codon family for an amino-acid"
        }),
        OrderedDict({
            "name": "Print .tsv of codon counts/frequencies",
            "shortname": "Print table of codon counts",
            "description": "Prints a .tsv table of codon counts or synonymous usage frequencies"
        })
    ]
})
    
COMMAND_8_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse .fna fasta sequences of CDS",
            "shortname": "Parse fasta sequences",
            "description": "Parse each sequence of the .fna file to validate and calculate counts",
        }),
        OrderedDict({
            "name": "Validate each CDS",
            "shortname": "CDS validation",
            "description": "Checks if length of CDS is divisible by 3, and contains valid/canonical start/stop codons"
        }),
        OrderedDict({
            "name": "Calculate 3-mer/codon counts",
            "shortname": "count codons",
            "description": "Make codon counts along each sequence"
        }),
        OrderedDict({
            "name": "Table of codon counts",
            "shortname": "codon count table",
            "description": "Print nx64 count/frequency matrix of n sequences and their 64 codon counts"
        })

    ]

})

"""
Command 9
"""

command_9_description = "Codon usage bias"
command_9_description_long = "Use codon counts ('kmerdb codons') for Chi-Square tests of overrepresentation of specific codons within a synonymous codon family"
command_9_parameters = "--as-frequencies changes counts to frequencies within a synonymous amino-acid family of codons, --ignore-noncanonicals omits sequences that do not contain standard CDS definition (non-ATG start codon, non-TAA/TAG/TGA stop codons), --ignore-invalid-cds omits sequences that are considered invalid (length not divisible by 3, otherwise uses --ignore-noncanonicals definition). --dont-ignore-start-codons includes counts from start codons for downstream analysis, --dont-ignore-stop-codons includes counts from stop codons for downstream analysis"
command_9_inputs = "Input is a table of codon counts from a reference population and a .fna fasta file of CDS sequences to compare to the reference family's codon usage"
command_9_usage = "kmerdb CUB -vv <--sequences input.fna> <codon_usage.tsv>"

        
COMMAND_9_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb CUB -vv <--sequences input.fna> <codon_usage.tsv>

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
        {
            "name" : "--as-frequencies",
            "type": "parameter",
            "value" : "Output is frequencies of codon usage within an synonymous amino-acid family"
        },
        {
            "name": "--ignore-noncanonicals",
            "type": "parameter",
            "value": "Ignore (dont throw error) sequences with non-canonical start/stop codons (more permissive)",
        },
        {
            "name": "--ignore-invalid-cds",
            "type": "parameter",
            "value": "Ignore (dont throw error) sequences that do not fit the canonical definition of a CDS: length divisible by 3, standard start/stop codons (If you want the program to error out on sequences that are considered invalid then dont use this parameter. Otherwise this parameter will just ignore and omit those sequences in the final table)"
        },
        {
            "name": "--include-start-codons",
            "type": "parameter",
            "value": "Include start codon counts in the count/frequency table"
        },
        {
            "name": "--include-stop-codons",
            "type": "parameter",
            "value": "Include stop codon counts in the count/frequency table"
        },
        {
            "name": "--output-delimiter",
            "type": "parameter",
            "value": "The choice of output delimiter for the codon count/frequency table"
        }
    ]
})

        
COMMAND_9_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "--sequences",
            "type": "file",
            "value": "Input is a .fna file of CDS sequences."
        },
        {
            "name": "input",
            "type": "file",
            "value": "'input' is a codon_usage table of counts (not frequencies) from 'kmerdb codons'"
        }
        
    ]
})

COMMAND_9_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Read fasta input",
            "shortname": "Parse fasta file",
            "description" : "Parse CDS sequences from a .fna fasta file"
        }),
        OrderedDict({
            "name": "Read codon usage .tsv table",
            "shortname": "Read codon count table",
            "description": "Read codon count table as .tsv from 'kmerdb codons'"
        }),
        OrderedDict({
            "name": "Chi-square test",
            "shortname": "chi-square test",
            "description": "Calculate chi-square test statistics across each codon's expected and observed counts, based on the reference family of sequences given by 'kmerdb codons', within a synonymous amino-acid family. Performs 20 (20 different amino-acids) chi-square tests for each statistic. Working on the multiple-hypothesis correction framework."
        }),
        OrderedDict({
            "name": "Print .tsv of chi-square test results",
            "shortname": "Print table of chi-square results",
            "description": "Prints a .tsv table of chi-square values and p-values for each of the 20 synonymous usage families"
        })
    ]
})
    
COMMAND_9_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse the codon count (not frequency) table",
            "shortname": "parse codon counts",
            "description": "Parse the codon counts .tsv table from 'kmerdb codons'"
        }),
        OrderedDict({
            "name": "Parse .fna fasta sequences of CDS",
            "shortname": "Parse fasta sequences",
            "description": "Parse each sequence of the .fna file to validate and calculate counts",
        }),
        OrderedDict({
            "name": "Validate each CDS",
            "shortname": "CDS validation",
            "description": "Checks if length of CDS is divisible by 3, and contains valid/canonical start/stop codons"
        }),
        OrderedDict({
            "name": "Calculate Chi-square test statistics and p-values",
            "shortname": "chi-square tests",
            "description": "Calculate 20 different chi-square tests for each sequence (one for each amino acid)"
        }),
        OrderedDict({
            "name": "Print table of chi-square results",
            "shortname": "print chi-square results",
            "description": "Print nx20 chi-square result matrix of n sequences and their 20 amino acid results from chi-square"
        })

    ]

})







"""

Command 10
"""






command_10_description = "create a edge list in (block) .gz format from .fasta|.fa or .fastq format."
command_10_description_long = """


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
        
command_10_parameters = "uses < -k > for k-mer size, --quiet to reduce runtime, -v, -vv to control logging."
command_10_inputs = "Input file can .fastq (or .fa).   - gzip.  Output is a weighted edge list in .kdb format (gzipped .csv with YAML header)"
command_10_usage = "kmerdb graph -k $K --quiet <input_1.fa.gz> [input_2.fq.gz] <output_edge_list_file.12.kdbg>" 


COMMAND_10_BANNER = """



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











""".format(command_10_name, command_10_description, command_10_description_long, command_10_inputs, command_10_parameters, command_10_usage)

        

        
COMMAND_10_PARAMS = OrderedDict({
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
        
        

COMMAND_10_INPUTS = OrderedDict({
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


COMMAND_10_FEATURES = OrderedDict({
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


COMMAND_10_STEPS = OrderedDict({
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
        

"""

Command 11
"""




command_11_description = "Calculate minimizers"
command_11_description_long = """
    Calculate minimizers, outputs a binary array to associate with a index array.
    


"""
command_11_parameters = "Parameter of interest is window size"
command_11_inputs = "Input is a v{0} .kdb count vector file".format(config.VERSION)
command_11_usage = "kmerdb minimizers -vv --window-size W --debug input1.12.kdb > input1.12.kdb.kdbi"


        
COMMAND_11_BANNER = """












[--------------------------------------------------------------------------------------]








                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------

                    kmerdb minimizers --window-size W input1.12.kdb > input1.12.kdb.kdbi

                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_11_name, command_11_description, command_11_description_long, command_11_inputs, command_11_parameters, command_11_usage)

        

COMMAND_11_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
        {
            "name": "Reference .kdb file to create minimizers from",
            "type": "file",
            "value": "The k-mer count profile to calculate minimizers."
        }
    ]
})

        
COMMAND_11_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "",
            "type": "file",
            "value": "File to decompose the k-mer profile into its parts."
        }
        
    ]
})

COMMAND_11_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "linear regression for ",
            "shortname": "",
            "description" : ""
        }),
        OrderedDict({
            "name": "",
            "shortname": "",
            "description": "(Deprecated)"
        })
    ]
})
    
COMMAND_11_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "",
            "shortname": "",
            "description": "(uhhhh...)",
        }),
        OrderedDict({
            "name": "",
            "shortname": "Shuffle k-mer counts",
            "description": "(Deprecated)"
        })

    ]

})



command_12_description = "Sequence alignment"
command_12_description_long = """
        Perform Smith-Waterman alignment on a reference database using seed regions (minimizers) determined from the reference .fasta sequence and the associated .kdb file.
    


"""
command_12_parameters = "Parameters of note are match score, mismatch score, gap open/extend penalties"
command_12_inputs = "Input is a v{0} .kdbi minimizer index file and two fasta files, a query and a reference".format(config.VERSION)
command_12_usage = "kmerdb alignment -vv --debug query.fasta reference.fasta input1.8.kdbi"


        
COMMAND_12_BANNER = """












[--------------------------------------------------------------------------------------]


        kmerdb alignment query.fasta reference.fasta input1.8.kdbi





                        [  n a m e    ]         :  -   {0}

                   description : {1}

{2}




--------------------------



                    [-]    inputs : 

                           {3}

                    [-]    parameters : 

                           {4}



                    [-]    [ usage ]  :  {5}


















[--------------------------------------------------------------------------------------]
""".format(command_12_name, command_12_description, command_12_description_long, command_12_inputs, command_12_parameters, command_12_usage)

        

COMMAND_12_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
        {
            "name": "Window size",
            "type": "int",
            "value": "The window-size to use to count minimizers with."
        }
    ]
})

        
COMMAND_12_INPUTS = OrderedDict({
    "name": "inputs",
    "type": "array",
    "items": [
        {
            "name": "The reference file's .kdbi minimizer index file.",
            "type": "file",
            "value": "The reference k-mer minimizers."
        },
        {
            "name": "Reference fasta file.",
            "type": "file",
            "value": "The reference fasta file to perform alignment against."
        },
        {
            "name": "Query fasta file.",
            "type": "file",
            "value": "The query fasta file to calculate minimizers and perform alignment with."
        }

    ]
})

                                

COMMAND_12_FEATURES = OrderedDict({
    "name": "features",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Parse reference .fna fasta sequence file ",
            "shortname": "Parse references",
            "description" : "Load .fna reference fasta sequences"
        }),
        OrderedDict({
            "name": "Parse query .fna fasta sequence file ",
            "shortname": "Parse queries",
            "description" : "Load .fna query fasta sequences"
        }),
        OrderedDict({
            "name": "Parse reference minimizers .kdbi file",
            "shortname": "Parse minimizers",
            "description": "Parse .kdbi minimizer index file for reference sequences"
        }),
        OrderedDict({
            "name": "Perform Smith-Waterman alignment on queries against references",
            "shortname": "Smith-Waterman alignment",
            "description": "Use minimizers as possible seed regions to match queries against references and extend alignment using Smith-Waterman"
        })
    ]
})
    
COMMAND_12_STEPS = OrderedDict({
    "name": "steps",
    "type": "array",
    "items": [
        OrderedDict({
            "name": "Load queries and references",
            "shortname": "load fasta files",
            "description": "Read the .fna query/reference fasta sequences",
        }),
        OrderedDict({
            "name": "Load the minimizer index",
            "shortname": "load .kdbi minimizers",
            "description": "Read the .kdbi minimizer index file for seed regions"
        }),
        OrderedDict({
            "name": "Calculate seed matches of each query against references",
            "shortname": "Seed matches",
            "description": "Seed possible alignment matches using the matches of .kdbi minimizers of references sequences with seed regions on the query for Smith-Waterman extension"
        }),
        OrderedDict({
            "name": "Perform SW alignment",
            "shortname": "Smith-Waterman alignment",
            "description": "Perform full alignment/extension from seeded matches of query sequences against the reference sequences."
        })

    ]

})




command_13_description = "Create an index file (deprecated)"
command_13_description_long = ""
command_13_parameters = "N/A"
command_13_inputs = "Input is a v{0} .kdb count vector file".format(config.VERSION)
command_13_usage = "kmerdb index <kmer_count_vector.kdb>"

        
COMMAND_13_BANNER = """












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
""".format(command_13_name, command_13_description, command_13_description_long, command_13_inputs, command_13_parameters, command_13_usage)

        

COMMAND_13_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
    ]
})

        
COMMAND_13_INPUTS = OrderedDict({
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

COMMAND_13_FEATURES = OrderedDict({
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
    
COMMAND_13_STEPS = OrderedDict({
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


        

command_14_description = "Shuffle k-mer count vector (Deprecated)"
command_14_description_long = ""
command_14_parameters = "N/A"
command_14_inputs = "Input is a v{0} .kdb count vector file".format(config.VERSION)
command_14_usage = "kmerdb shuf <kmer_count_vector.kdb>"

        
COMMAND_14_BANNER = """












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
""".format(command_14_name, command_14_description, command_14_description_long, command_14_inputs, command_14_parameters, command_14_usage)

        

COMMAND_14_PARAMS = OrderedDict({
    "name": "arguments",
    "type": "array",
    "items": [
    ]
})

        
COMMAND_14_INPUTS = OrderedDict({
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

COMMAND_14_FEATURES = OrderedDict({
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
    
COMMAND_14_STEPS = OrderedDict({
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




# command_12_description = "K-mer Markov sequence probability feature  (Deprecated)"
# command_12_description_long = """
#  Uses conditional probabilities and multiplication rule along with Markov model of sequence to use likelihood/odds ratios to test likelihood of the sequence(s) given the inputs. Unadjusted for multiple-hypothesis testing.

# Conditional probability :

# P(X|Y) = P(XY)/P(Y)


# Multiplication rule :

# P(X) = P(an|an-1,an-2,...a1)     a { N  , X = an-1,an-2,...,a1




# """
# command_12_parameters = "inputs may be one or more fasta files, and the .kdb files needed for the model's output probability."
# command_12_inputs = "Input is a v{0} .kdb count vector file".format(config.VERSION)
# command_12_usage = "kmerdb prob --db <kmer_count_vector_1.kdb> [--db kmer_count_vector_2.kdb] <query_sequences_1.fasta.gz> [query_sequences_2.fasta.gz]"


        
# COMMAND_12_BANNER = """












# [--------------------------------------------------------------------------------------]








#                         [  n a m e    ]         :  -   {0}

#                    description : {1}

# {2}




# --------------------------

#                     kmerdb prob <--db input.kdb> sequences.fasta[.gz]

#                     [-]    inputs : 

#                            {3}

#                     [-]    parameters : 

#                            {4}



#                     [-]    [ usage ]  :  {5}


















# [--------------------------------------------------------------------------------------]
# """.format(command_12_name, command_12_description, command_12_description_long, command_12_inputs, command_12_parameters, command_12_usage)

        

# COMMAND_12_PARAMS = OrderedDict({
#     "name": "arguments",
#     "type": "array",
#     "items": [
#         {
#             "name": "K-mer database file.",
#             "type": "file",
#             "value": "[NOTE] : multiple may be specified. | The k-mer frequencies to use in the Markov model probability calculations, and likelihood/odds-ratio tests"
#         }
#     ]
# })

        
# COMMAND_12_INPUTS = OrderedDict({
#     "name": "inputs",
#     "type": "array",
#     "items": [
#         {
#             "name": "Input .fa|.fasta|.fa.gz files.",
#             "type": "array",
#             "value": "File(s) to query against the k-mer database."
#         }
        
#     ]
# })

# COMMAND_12_FEATURES = OrderedDict({
#     "name": "features",
#     "type": "array",
#     "items": [
#         OrderedDict({
#             "name": "[FIXME! 7/28/24 Hi fross, glhf!]",
#             "shortname": "",
#             "description" : ""
#         }),
#         OrderedDict({
#             "name": "",
#             "shortname": "",
#             "description": "(Deprecated)"
#         })
#     ]
# })
    
# COMMAND_12_STEPS = OrderedDict({
#     "name": "steps",
#     "type": "array",
#     "items": [
#         OrderedDict({
#             "name": "",
#             "shortname": "",
#             "description": "(uhhhh...)",
#         }),
#         OrderedDict({
#             "name": "",
#             "shortname": "Shuffle k-mer counts",
#             "description": "(Deprecated)"
#         })

#     ]

# })




# command_14_description = "Compositional analysis"
# command_14_description_long = """
#     Linear regression wth two inputs.
#     A 2D array 'A' with two dimensions, and non-necessarily square.

#     and a vector of data determining the linear constants in the regression.
    

#     least squares problems occur when b is not in the column space of A.
    
#     They are solvable if and only if the matrix A:
#       * is non-singular
#       * has independent rows/columns
#       * a non-zero determinant (i.e. it is *not* row identical to the identity matrix)
#       * and is full-rank. 

    
#     ############
#     # definition
#     ############
#     At(b-Ax) = 0 OR
#     AtAx = Atb

#     but, less developed...
#     Ax = b

#     These are the normal equations.


# """
# command_14_parameters = "inputs are the table of counts to use in the decomposition, and the composite/collated metagenomic k-mer profile (as .kdb) to decompose."
# command_14_inputs = "Input is a v{0} .kdb count vector file".format(config.VERSION)
# command_14_usage = "kmerdb composition -vv --debug input1.12.kdb 12_mer_profiles.tsv"


        
# COMMAND_14_BANNER = """












# [--------------------------------------------------------------------------------------]








#                         [  n a m e    ]         :  -   {0}

#                    description : {1}

# {2}




# --------------------------

#                     kmerdb composition composite1.12.kdb kmer_counts.12.tsv

#                     [-]    inputs : 

#                            {3}

#                     [-]    parameters : 

#                            {4}



#                     [-]    [ usage ]  :  {5}


















# [--------------------------------------------------------------------------------------]
# """.format(command_14_name, command_14_description, command_14_description_long, command_14_inputs, command_14_parameters, command_14_usage)

        

# COMMAND_14_PARAMS = OrderedDict({
#     "name": "arguments",
#     "type": "array",
#     "items": [
#         {
#             "name": "K-mer database file.",
#             "type": "file",
#             "value": "The k-mer count proifle of a composite/metagenomic population to decompose into percentages"
#         }
#     ]
# })

        
# COMMAND_14_INPUTS = OrderedDict({
#     "name": "inputs",
#     "type": "array",
#     "items": [
#         {
#             "name": "Input count matrix (.tsv)",
#             "type": "file",
#             "value": "File to decompose the k-mer profile into its parts."
#         }
        
#     ]
# })

# COMMAND_14_FEATURES = OrderedDict({
#     "name": "features",
#     "type": "array",
#     "items": [
#         OrderedDict({
#             "name": "linear regression for ",
#             "shortname": "",
#             "description" : ""
#         }),
#         OrderedDict({
#             "name": "",
#             "shortname": "",
#             "description": "(Deprecated)"
#         })
#     ]
# })
    
# COMMAND_14_STEPS = OrderedDict({
#     "name": "steps",
#     "type": "array",
#     "items": [
#         OrderedDict({
#             "name": "",
#             "shortname": "",
#             "description": "(uhhhh...)",
#         }),
#         OrderedDict({
#             "name": "",
#             "shortname": "Shuffle k-mer counts",
#             "description": "(Deprecated)"
#         })

#     ]

# })




###################################################

#            F i n a l     c o m m a n d    a g g r e g a t e

###################################################


ALL_PARAMS = {
    "profile": COMMAND_1_PARAMS["items"],
    "graph": COMMAND_2_PARAMS["items"],
    "matrix": COMMAND_3_PARAMS["items"],
    "distance": COMMAND_4_PARAMS["items"],
    "kmeans": COMMAND_7_PARAMS["items"],
    "hierarchical": COMMAND_8_PARAMS["items"],    
    "view": COMMAND_5_PARAMS["items"],
    "header": COMMAND_6_PARAMS["items"],
    "index": COMMAND_9_PARAMS["items"],
    "shuf": COMMAND_10_PARAMS["items"],
    "minimizers": COMMAND_11_PARAMS["items"],
    "alignment": COMMAND_12_PARAMS["items"],
    #"composition": COMMAND_14_PARAMS["items"],
}

ALL_INPUTS = {
    "profile": COMMAND_1_INPUTS["items"],
    "graph": COMMAND_2_INPUTS["items"],
    "matrix": COMMAND_3_INPUTS["items"],
    "distance": COMMAND_4_INPUTS["items"],
    "kmeans": COMMAND_7_INPUTS["items"],
    "hierarchical": COMMAND_8_INPUTS["items"],
    "view": COMMAND_5_INPUTS["items"],
    "header": COMMAND_6_INPUTS["items"],
    "index": COMMAND_9_INPUTS["items"],
    "shuf": COMMAND_10_INPUTS["items"],
    
    "minimizers": COMMAND_11_INPUTS["items"],
    "alignment": COMMAND_12_INPUTS["items"],
    #"composition": COMMAND_14_INPUTS["items"],
}

ALL_FEATURES = {
    "profile": COMMAND_1_FEATURES["items"],
    "graph": COMMAND_2_FEATURES["items"],
    "matrix": COMMAND_3_FEATURES["items"],
    "distance": COMMAND_4_FEATURES["items"],
    "kmeans": COMMAND_7_FEATURES["items"],
    "hierarchical": COMMAND_8_FEATURES["items"],
    "view": COMMAND_5_FEATURES["items"],
    "header": COMMAND_6_FEATURES["items"],
    "index": COMMAND_9_FEATURES["items"],
    "shuf": COMMAND_10_FEATURES["items"],
    "minimizers": COMMAND_11_FEATURES["items"],
    "alignment": COMMAND_12_FEATURES["items"],
    #"composition": COMMAND_14_FEATURES["items"],

}


ALL_STEPS = {

    "profile": COMMAND_1_STEPS["items"],
    "graph": COMMAND_2_STEPS["items"],
    "matrix": COMMAND_3_STEPS["items"],
    "distance": COMMAND_4_STEPS["items"],
    "kmeans": COMMAND_7_STEPS["items"],
    "hierarchical": COMMAND_8_STEPS["items"],
    "view": COMMAND_5_STEPS["items"],
    "header": COMMAND_6_STEPS["items"],
    "index": COMMAND_9_STEPS["items"],
    "shuf": COMMAND_10_STEPS["items"],
    "minimizers": COMMAND_11_FEATURES["items"],
    "alignment": COMMAND_12_FEATURES["items"],
    #"composition": COMMAND_14_STEPS["items"],
}




        
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





        


        
        # KMEANS_BANNER = """
        #                   name : kmeans
        #            description : 
        # """

        # HIERARCHICAL_BANNER = """

        # """


        # config.py features
        #SUPPORTED_FEATURES = {}

        # ASCII + RICH.py checklist
        
    def print_verbosity_header(self):
        """
        This prints the verbosity warnings and default logging behavior for the user.
        """

        sys.stderr.write("="*40 + "\n")
        sys.stderr.write("[ DEBUG ] : Default is warning-only logging. [-v,-vv, --debug, --quiet]\n")
        sys.stderr.write("[ INFO  ] : Default is warning-only logging. [-v,-vv, --debug, --quiet]\n")
        sys.stderr.write("="*40 + "\n")
        sys.stderr.write("WARNING. Some features are experimental. Note: you are tracking kmerdb/main. Visit the README.md header, usage section, or quickstart page for additl.\n")
        sys.stderr.write("-"*40 + "\n")

        #print("x..xx")
        #print("xd up ya brakes")

        #
        #
        #...


        #' ...#x'"' ''"'?"    # ... #xdupyabrakexsz '#xdupyabreakesz' '"'"''"""""""...'"'. #XduPyabrakx. '"' ....

        return




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
        
        sys.stderr.write(DNA_SPACER_lol)

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


        # if e is None:
        #     raise ValueError("Need an error to exit")
        # elif not isinstance(e, Exception):
        #     raise ValueError("Need an error to exit")

        
        # if n_logs is None or type(n_logs) is not int:
        #     raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument n_logs to be a int")
        # elif logs is None or type(logs) is not list:
        #     raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument logs to be a list")
        # elif feature is not None and type(feature) is not int:
        #     raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument feature to be a int")
        # elif step is not None and type(step) is not int:
        #     raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument step to be a int")
        # elif subcommand is not None and type(subcommand) is not str:
        #     raise TypeError("kmerdb.appmap.exit_gracefully expects the keyword argument subcommand to be a str")



        
        N = min(len(logs), n_logs)

        tb = traceback.extract_tb(e.__traceback__)
        last_traceback_FrameSummary = tb[-1]
        error_file_name = last_traceback_FrameSummary.filename
        error_line_number = last_traceback_FrameSummary.lineno

        
        #assert subcommand in config.subcommand_functions, "Unknown subcommand"

        if self._loggable:
            self.logger.log_it("Program error! Collecting error metadata and formatting...", "ERROR")
        # This is the "Error blocks" metadata


        sys.stderr.write("Aggregating program metadata, if this fails without error without the --debug flag, please report to the GitHub issue tracker with the title 'Error summary convenience function'.")
        sys.stderr(subcommand)
        sys.stderr(config.VERSION)
        sys.stderr(config.REQUIRES_PYTHON)
        sys.stderr(feature)
        sys.stderr(ALL_FEATURES[subcommand][feature]["name"])
        sys.stderr(ALL_FEATURES[subcommand][feature]["shortname"])
        sys.stderr(ALL_FEATURES[subcommand][feature]["description"])
        sys.stderr(step)
        sys.stderr(ALL_STEPS[subcommand][step]["name"])
        sys.stderr(ALL_STEPS[subcommand][step]["shortname"])
        sys.stderr(ALL_STEPS[subcommand][step]["description"])
            # The *total* number of logged lines produced by the program and returned to the global 'logs' var in __init__.py
        sys.stderr(self.logfile)
        sys.stderr(str(tb))
        sys.stderr(error_file_name)
        sys.stderr(error_line_number)
        sys.stderr(e.__str__())

        
        e_sum = {
            "subcommand": subcommand,
            "kmerdb-version": config.VERSION,
            "python-version": config.REQUIRES_PYTHON,
            "feature": feature,
            "feature_name": ALL_FEATURES[subcommand][feature]["name"],
            "feature_shortname": ALL_FEATURES[subcommand][feature]["shortname"],
            "feature_description": ALL_FEATURES[subcommand][feature]["description"],
            "step" : step,
            "step_name": ALL_STEPS[subcommand][step]["name"],
            "step_shortname": ALL_STEPS[subcommand][step]["shortname"],
            "step_description": ALL_STEPS[subcommand][step]["description"],
            # The *total* number of logged lines produced by the program and returned to the global 'logs' var in __init__.py
            "log_file": self.logfile,
            "traceback": str(tb),
            "error_file_name": error_file_name,
            "error_line_number": error_line_number,
            "error": e.__str__(),
        }

        
        exit_summary = OrderedDict(e_sum)

        
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



        
        for i in range(N):

            try:
                if self._loggable:
                    sys.stderr.write("{0} - last line of log\n".format(n_logs - i))
                    self.logger.log_it(logs[i], "ERROR")
                else:
                    sys.stderr.write("{0} - last line of log\n".format(n_logs - i))
                    sys.stderr.write(logs[i])
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


    








