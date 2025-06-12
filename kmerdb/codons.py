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
import logging



logger = logging.getLogger(__file__)



from kmerdb import kmer





codon_table = { # A 3-mer id hashed to the value, which is a amino acid letter or ASCII char code
    63: "F",
    61: "F",
    60: "L",
    62: "L", # Possible start codon
    31: "L",
    29: "L",
    28: "L",
    30: "L",
    15: "I",
    13: "I",
    12: "I",
    14: "M", # START codon
    46: "V", # possible start codon
    45: "V",
    44: "V",
    47: "V",
    52: "S",
    53: "S",
    54: "S",
    55: "S",
    20: "P",
    21: "P",
    22: "P",
    23: "P",
    24: "P",
    4: "T",
    5: "T",
    6: "T",
    7: "T",
    36: "A",
    37: "A",
    38: "A",
    39: "A",
    49: "Y",
    51: "Y",
    48: None, # STOP codon
    50: None, # STOP codon
    17: "H",
    19: "H",
    16: "Q",
    18: "Q",
    1: "N",
    3: "N",
    0: "K",
    2: "K",
    35: "D",
    33: "D",
    32: "E",
    34: "E",
    59: "C",
    57: "C",
    56: None, # STOP codon
    58: "W",
    24: "R",
    25: "R",
    26: "R",
    27: "R",
    11: "S",
    9: "S",
    40: "G",
    41: "G",
    42: "G",
    43: "G"
}


synonymous_usage = {
    "F": [61, 63],
    "L": [60, 62, 31, 29, 28, 30],
    "I": [13, 15, 12],
    "M": [14],
    "V": [44, 45, 46, 47],
    "S": [52, 53, 54, 55],
    "P": [20, 21, 22, 23],
    "T": [4, 5, 6, 7],
    "A": [36, 37, 38, 39]
    "Y": [49, 51],
    "H": [17, 19],
    "Q": [16, 18],
    "N": [1, 3],
    "K": [0, 2],
    "D": [33, 35],
    "E": [32, 34],
    "C": [57, 59],
    "W": [58],
    "R": [24, 25, 26, 27],
    "S": [9, 11]
    "G": [40, 41, 42, 43]
}
