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
import json


logger = logging.getLogger(__file__)


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from kmerdb import kmer

import pandas as pd
import numpy as np



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
    48: None, # STOP codon, TAA/UAA, ochre
    50: None, # STOP codon, TAG/UAG, amber
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
    56: None, # STOP codon, TGA/UGA, opal
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
CODON_LIST = list(range(64))
#CODON_LIST = list(codon_table.keys())
CODON_IUPAC_CODES = list(codon_table.values())
#CODON_LIST = list(set(codon_table.values()))
#CODON_LIST.remove(None)

if len(codon_table.keys()) != 64:
    raise ValueError("kmerdb.codons : cannot have more than 64 codons")


STOP_CODONS = [48, 50, 56]
POSSIBLE_START_CODONS = [14, 46, 62]

synonymous_codons = {
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


def is_sequence_cds(seq, include_noncanonicals:bool=False):
    if type(seq) is not str and not isinstance(seq, Seq) and not isinstance(seq, SeqRecord):
        raise TypeError("kmerdb.codons.is_sequence_cds() expects a sequence as its first positional argument")
    seqid, s  = (seq.id, seq.seq)#, seq.desc)
    seq = str(s)
    seq_len = len(seq)
    if seq_len % 3 != 0:
        logger.error("kmerdb.codons.is_sequence_cds() expects a valid coding sequence as its first positional argument. Length was {0} which is not evenly divisible by 3".format(seq_len))
        return False


    codons, non_canonical = get_codons_in_order(seq, seq_id=seqid)
    #print(codons)
    #print(list(map(lambda c: kmer.id_to_kmer(c, 3), codons)))
    if codons is None:
        logger.warning("kmerdb.codons.is_sequence_cds() cannot mark this sequence as a CDS when containing an unknown residue")
        return False
    elif non_canonical is True:
        start_codons = ", ".join(list(map(lambda c: kmer.id_to_kmer(c, 3), POSSIBLE_START_CODONS)))
        stop_codons = ", ".join(list(map(lambda c: kmer.id_to_kmer(c, 3), STOP_CODONS)))
        noncanon_start = kmer.id_to_kmer(codons[0], 3)
        noncanon_stop = kmer.id_to_kmer(codons[-1], 3)

        msg = "kmerdb.codons.is_sequence_cds() expects the coding sequence '{0}' starting with '{1}' to be one of the start codons: {2} OR ending with '{3}' to be one of the stop codons {4}".format(seqid, noncanon_start, start_codons, noncanon_stop, stop_codons)
        if codons[0] not in POSSIBLE_START_CODONS:
            if include_noncanonicals is True:
                return True
            elif include_noncanonicals is False:
                logger.warning(msg)
                return False
    else:
        return True
    
    
def get_codons_in_order(seq:str, seq_id:str=None):
    """
    :param seq: a fasta nucleic-acid sequence as a string
    :type bool:
    :param seq_id: the sequence identifier for the fasta sequence
    :type str:
    :raises TypeError: When the seq 1st positional argument is not a str
    :raises TypeError: When the seq_id keyword argument is not a str
    :raises ValueError: When the sequence does not use IUPAC nucleic acid characters
    :raises ValueError: When the sequence has a length whos modulus of 3 is not  0 (L%3 != 0)
    :returns: (codon_list, is_non_canonical:bool)
    :rtype: tuple
    """
    if type(seq) is not str:
        raise TypeError("kmerdb.codons.get_codons_in_order() expects a str as its first positional argument")
    elif seq_id is None or type(seq_id) is not str:
        raise TypeError("kmerdb.codons.get_codons_in_order() expects the keyword argument 'seq_id' to be a bool")
    
    if not kmer.is_sequence_na(seq):
        raise ValueError("kmerdb.codons.get_codons_in_order() expects a nucleic acid sequence as its first positional argument")
    elif len(seq) % 3 != 0:
        raise ValueError("kmerdb.codons.get_codons_in_order() expects the sequence to be divisible by three")
    
    seq_len = len(seq)
    is_non_canonical = False
    codons = []
    codon = ''
    
    i = 0

    for j, c in enumerate(seq):
        codon += c
        if i == 2:
            codon_id = kmer.kmer_to_id(codon)
            if codon_id is None:
                """
                               h a n d l e    a    i n v a l i d    k m e r    i d
                """
                is_non_canonical = True
                logger.error("\t\tParameters causing error:\n\n")
                msg = "Sequence '{0}' contained an invalid codon at position {1}. The codon encountered at this point was '{2}'. Omitting this sequence from consideration".format(seq_id, j, codon)
                logger.error(msg)
                """
                If an invalid codon/3-mer id is encountered, reset and skip this codon count.
                """
                codon = ''
                i = 0
                continue
                #return None, None # This is icky behavior
            elif j == 2 and codon_id not in POSSIBLE_START_CODONS:
                is_non_canonical = True
            if j == seq_len - 1:
                if codon_id not in STOP_CODONS:
                    is_non_canonical = True
                # Do not append stop codon to count
            codons.append(codon_id)
            codon = ''
            i = 0
        else:
            i+=1
    return codons, is_non_canonical


def count_codons(codon_list:list, include_start_codons:bool=False, include_stop_codons:bool=False):
    """
    :param codon_list: List of codon/3-mer ids in the order found in the sequence
    :type list:
    :param include_start_codons: include the start codon in the counts?
    :type bool:
    :param include_stop_codons: include the stop codon in the counts?
    :type bool:
    :returns: Vector of 64 elements corresponding to 3-mer/codon ids.
    :rtype: np.ndarray
    """
    if type(codon_list) is not list or not all(type(c) is int for c in codon_list):
        raise ValueError("kmerdb.codons.count_codons() expects a list of codon ids as its first positional argument")
    elif type(include_start_codons) is not bool:
        raise TypeError("kmerdb.codons.count_codons() expects the keyword argument 'include_start_codons' to be a bool")
    elif type(include_stop_codons) is not bool:
        raise TypeError("kmerdb.codons.count_codons() expects the keyword argument 'include_stop_codons' to be a bool")


    codon_counts = np.zeros(64, dtype="uint32")

    lenseq = len(codon_list)
    for i, cdn in enumerate(codon_list):
        if i == 0:
            if include_start_codons is True:
                codon_counts[cdn] += 1
            else:
                continue # Ommit the start codon from the counts
        elif i == lenseq - 1:
            if include_stop_codons is True:
                codon_counts[cdn] += 1
            else:
                continue # Omit the stop codon from the counts
        else:
            codon_counts[cdn] += 1
    return codon_counts





def codon_frequency_table(seq, seqid, include_stop_codons:bool=False, include_start_codons:bool=False):
    """
    :param seq: A nucleic acid CDS sequence 
    :type str:
    :param seqid: A fasta identifier for the sequence
    :type str:
    :param include_stop_codons: do not include codon counts for stop codons
    :type bool:
    :param include_start_codons: do not include codon counts for start codons
    :type bool:
    :returns: A typle of (codon_ids, codon_counts, codon_frequencies_wrt_length, codon_frequencies_in_family)
    :rtype: tuple
    Takes an input sequence 'seq' and produces 3-mer frequencies
    Returns a tuple of codon ids, and codon counts, and codon frequencies (wrt length), and codon frequencies within an synonymous AA family
    """

    if type(seq) is not str:
        raise TypeError("kmerdb.codons.codon_frequency_table() expects a str as its first positional argument")
    elif type(seqid) is not str:
        raise TypeError("kmerdb.codons.codon_frequency_table() expects a str as its second positional argument")
    elif type(include_start_codons) is not bool:
        raise TypeError("kmerdb.codons.codon_frequency_table() expects the keyword argument 'include_start_codons' to be a bool")
    elif type(include_stop_codons) is not bool:
        raise TypeError("kmerdb.codons.codon_frequency_table() expects the keyword argument 'include_stop_codons' to be a bool")
    
    
    seq_len = len(seq) # Do not use seq_len as proxy for CDS length, as stop codons will be included
    if seq_len % 3 != 0:
        raise ValueError("kmerdb.codons.codon_frequency_table() expects the sequence length to be divisible by three")
        
    codons, non_canonical = get_codons_in_order(seq, seq_id=seqid)

    
    num_codons = len(codons)
    if num_codons != len(seq) / 3 and num_codons != (len(seq) / 3) - 1 and num_codons != (len(seq) / 3) - 2: # Remove start or stop codons or both
        msg = "For sequence '{0}', the number of codons ({1}) should be proportional to the length ({2}/3) of the sequence...".format(seqid, num_codons, len(seq))
        raise ValueError(msg)

    codon_counts = count_codons(codons, include_start_codons=include_start_codons, include_stop_codons=include_stop_codons)
    codon_frequencies_wrt_length = np.array(codon_counts) / np.sum(codon_counts)

    codon_frequencies_in_family = np.zeros(64)
    for kmer_id in range(64):
        #kmer_id = CODON_LIST[i]
        corresponding_aa = codon_table[kmer_id]

        if corresponding_aa is not None:
            cdn_cnts = codon_counts[synonymous_codons[corresponding_aa]]
            total_counts_for_family = cdn_cnts.sum()
            cdn_cnt = codon_counts[kmer_id]
            if total_counts_for_family != 0:
                codon_frequencies_in_family[kmer_id] = float(cdn_cnt) / total_counts_for_family
            #print("\n\n\n\nAmino acid: {0}\nCodon ids: {1}\nSpecific codon id: {2}\nCodon counts: {3}\nCorresponding count: {4}\nTotal counts: {5}\nFrequency: {6}\n\n\n\n".format(corresponding_aa, synonymous_codons[corresponding_aa], kmer_id, cdn_cnts, cdn_cnt, total_counts_for_family, codon_frequencies_in_family[kmer_id]))

    return (list(range(64)), codon_counts, codon_frequencies_wrt_length.tolist(), codon_frequencies_in_family.tolist())


def get_expected_codon_frequencies(L:int, codon_counts:pd.DataFrame):
    """
    Calculate the expected frequencies/counts of a sequence of length L based on the counts and their synonymous families in the pandas.DataFrame

    :param L: the length of the sequence in question
    :type L: int
    :param codon_counts: the codon frequency table (rows=genes, columns=codon counts)
    :type codon_counts: pandas.DataFrame
    :raises TypeError: if type of L (length of sequence) is not int
    :raises TypeError: if the codon_counds is not a pandas.DataFrame
    :returns: expected codon (counts, frequencies) from sequence of length L
    :rtype: tuple
    """
    if type(L) is not int:
        raise TypeError("kmerdb.codons.get_expected_codon_frequencies() expects L to be an int")
    elif not isinstance(codon_counts, pd.DataFrame):
        raise TypeError("kmerdb.codons.get_expected_codon_frequencies() expects codon_count_vector to be an Pandas DataFrame")
    elif codon_counts.shape[1] != 64:
        raise ValueError("kmerdb.codons.get_expected_codon_frequencies() expects a 64 column data frame as its second positional argument")

    final_frequencies = np.zeros(64)
    total_cdn_cnts = np.array(codon_counts.sum(axis=0), dtype="float64") # Column sum of codon_counts across all genes
    frac_cdn_cnts = total_cdn_cnts / int(total_cdn_cnts.sum())
    exp_cnts = np.array(frac_cdn_cnts * (L/3), dtype="uint32")
    return exp_cnts, exp_cnts / exp_cnts.sum()
