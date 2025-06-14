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
CODON_LIST = list(codon_table.keys())

STOP_CODONS = [48, 50, 56]
POSSIBLE_START_CODONS = [14, 46, 62]

synonymous_usage = {
    "F": [61, 63],
    "L": [60, 62, 31, 29, 28, 30],
    "I": [13, 15, 12],
    "M": [14],
    "V": [44, 45, 46, 47],
    "S": [52, 53, 54, 55],
    "P": [20, 21, 22, 23],
    "T": [4, 5, 6, 7],
    "A": [36, 37, 38, 39],
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
    "S": [9, 11],
    "G": [40, 41, 42, 43]
}



def codon_frequency_table(seq):
    """
    Takes an input sequence 'seq' and produces 3-mer frequencies
    Returns a tuple of codon ids, and codon counts, and codon frequencies
    """

    if type(seq) is not str:
        raise TypeError("kmerdb.codons.codon_frequency_table expects a str as its first positional argument")

    
    #seq_len = len(seq) # Do not use seq_len as proxy for CDS length, as stop codons will be included
    if seq_len % 3 != 0:
        raise ValueError("kmerdb.codons.codon_frequency_table expects the sequence to be divisible by three")
    
    codons = get_codons_in_order(seq) # Returns a list of int codon ids, in order of translation frame
    codon_counts = []
    for cdn in CODON_LIST:
        codon_counts.append(codons.count(cdn))
    # Instead of using seq_len as CDS length, use length of codon list as it doesnt include STOP codons.
    return (CODON_LIST, codon_counts, list(map(lambda c: float(c)/len(codons), codon_counts)))

def get_codons_in_order(seq):
    if type(seq) is not str:
        raise TypeError("kmerdb.codons.codon_frequency_table expects a str as its first positional argument")
    elif len(seq) % 3 != 0:
        raise ValueError("kmerdb.codons.codon_frequency_table expects the sequence to be divisible by three")
    seq_len = len(seq)

    codons = []
    codon = ''
    
    i = 0
    for j, c in enumerate(seq):
        codon += c
        if i == 2:
            if kmer.kmer_to_id(codon) in STOP_CODONS:
                if j == seq_len:
                    break
                else:
                    raise ValueError("Internal stop codon identified: {0}\n\nOffending sequence: {1}".format(codon, seq))
            else:
                codons.append(kmer.kmer_to_id(codon)) # Convert to k-mer id before returning
                codon = ''
                i = 0
        else:
            i+=1
        
    return codons


