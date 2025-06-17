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


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
    8: "R",
    10: "R",
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


def is_sequence_cds(seq, ignore_noncanonicals:bool=False, ignore_invalid_cds:bool=False):
    if type(seq) is not str and not isinstance(seq, Seq) and not isinstance(seq, SeqRecord):
        raise TypeError("kmerdb.codons.is_sequence_cds() expects a sequence as its first positional argument")
    seqid, s  = (seq.id, seq.seq)#, seq.desc)
    seq = str(s)
    seq_len = len(seq)
    if seq_len % 3 != 0:
        logger.error("kmerdb.codons.is_sequence_cds() expects a valid coding sequence as its first positional argument. Length was {0} which is not evenly divisible by 3".format(seq_len))
        #raise ValueError("kmerdb.codons.is_sequence_cds() was given a CDS with length {0}, which is not evenly divisible by 3".format(seq_len))
        return False

    #logger.warning("Ignore invalid CDS? {0}".format(ignore_invalid_cds))
    codons = get_codons_in_order(seq, seq_id=seqid, ignore_invalid_cds=ignore_invalid_cds)
    #print(codons)
    #print(list(map(lambda c: kmer.id_to_kmer(c, 3), codons)))
    if codons is None:
        logger.warning("kmerdb.codons.is_sequence_cds() cannot mark this sequence as a CDS when containing an unknown residue")
        return False
    elif codons[0] not in POSSIBLE_START_CODONS:
        logger.warning("kmerdb.codons.is_sequence_cds() found a non-standard start codon '{0}' in sequence with id '{1}' ... ignoring non-canonical start codons".format(kmer.id_to_kmer(codons[0], 3), seqid))
        if ignore_noncanonicals is True:
            return False
        else:
            raise ValueError("kmerdb.codons.is_sequence_cds() expects the coding sequence {0} starting with {1} to be one of the start codons: {2} as the first positional argument".format(seqid, codons[0], ", ".join(list(map(lambda c: kmer.id_to_kmer(c, 3), POSSIBLE_START_CODONS)))))
    else:
        return True
    
    
def get_codons_in_order(seq, seq_id:str=None, ignore_invalid_cds:bool=False, include_stop_codons:bool=False):
    if type(seq) is not str:
        raise TypeError("kmerdb.codons.get_codons_in_order() expects a str as its first positional argument")
    elif seq_id is not None and type(seq_id) is not str:
        raise TypeError("kmerdb.codons.get_codons_in_order() expects the keyword argument 'seq_id' to be a bool")
    elif type(ignore_invalid_cds) is not bool:
        raise TypeError("kmerdb.codons.get_codons_in_order() expects the keyword argument 'ignore_invalid_cds' to be a bool")
    
    if not kmer.is_sequence_na(seq):
        raise ValueError("kmerdb.codons.get_codons_in_order() expects a nucleic acid sequence as its first positional argument")
    elif len(seq) % 3 != 0:
        raise ValueError("kmerdb.codons.get_codons_in_order() expects the sequence to be divisible by three")
    seq_len = len(seq)

    codons = []
    codon = ''
    
    i = 0

    for j, c in enumerate(seq):
        #print(seq_len, j, c)
        codon += c
        if i == 2:
            if j == seq_len - 1:
                if kmer.kmer_to_id(codon) not in STOP_CODONS:
                    if ignore_invalid_cds is True:
                        return None
                    else:
                        logger.warning("Ignore invalid CDS? {0}".format(ignore_invalid_cds))
                        msg = "Invalid stop codon identified: {0}. Expected to be one of: {1}\n\nOffending sequence: ID: {2}\n{3}".format(codon, ", ".join(list(map(lambda c: kmer.id_to_kmer(c, 3), STOP_CODONS))), seq_id, seq)
                        logger.error(msg)
                        raise ValueError(msg)
                elif kmer.kmer_to_id(codon) in STOP_CODONS:
                    #raise ValueError("Internal stop codon identified: {0}\n\nOffending sequence: {1}".format(codon, seq))
                    """
                    What should I do with a stop codon? Include it?
                    """
                    if include_stop_codons is True:
                        codons.append(kmer.kmer_to_id(codon)) # Comment if you don't like including stop codons
                    break
            else:
                codon_id = kmer.kmer_to_id(codon) # Convert to k-mer id before returning
                if codon_id is None:
                    logger.warn("kmerdb.codons.get_codons_in_order() found an unknown nucleic-acid residue {0} at position {1} of input sequence '{2}'...\ncannot translate this sequence\n\n".format(c, j, seq_id))
                    logger.warn("kmerdb.codons.get_codons_in_order() found an unknown amino-acid residue {0} in the input sequence".format(codon))
                    return None
                codons.append(codon_id)
                codon = ''
                i = 0
        else:
            i+=1
        
    return codons # Note, terminal codon is stop codon



    


def codon_frequency_table(seq, seqid, ignore_invalid_cds:bool=False):
    """
    Takes an input sequence 'seq' and produces 3-mer frequencies
    Returns a tuple of codon ids, and codon counts, and codon frequencies
    """

    if type(seq) is not str:
        raise TypeError("kmerdb.codons.codon_frequency_table expects a str as its first positional argument")
    elif type(seqid) is not str:
        raise TypeError("kmerdb.codons.codon_frequency_table expects a str as its second positional argument")
    elif type(ignore_invalid_cds) is not bool:
        raise TypeError("kmerdb.codons.codon_frequency_table expects the keyword argument 'ignore_invalid_cds' to be a bool")
    
    seq_len = len(seq) # Do not use seq_len as proxy for CDS length, as stop codons will be included
    if seq_len % 3 != 0:
        if ignore_invalid_cds is True:
            return (None, None, None)
        else:
            raise ValueError("kmerdb.codons.codon_frequency_table expects the sequence to be divisible by three")
    
    codons = get_codons_in_order(seq, seq_id=seqid, ignore_invalid_cds=ignore_invalid_cds) # Returns a list of int codon ids, in order of translation frame
    if codons is None:
        logger.warning("kmerdb.codons.codon_frequency_table() cannot return frequencies of invalid CDS")
        return (None, None, None)
    codon_counts = []
    for cdn in CODON_LIST:
        codon_counts.append(codons.count(cdn))
    # Instead of using seq_len as CDS length, use length of codon list as it doesnt include STOP codons.
    return (CODON_LIST, codon_counts, list(map(lambda c: float(c)/len(codons), codon_counts)))

