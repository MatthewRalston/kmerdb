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
CODON_LIST = list(codon_table.keys())
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
        if ignore_noncanonicals is True or ignore_invalid_cds is True:
            return False
        else:
            raise ValueError("kmerdb.codons.is_sequence_cds() expects the coding sequence {0} starting with {1} to be one of the start codons: {2} as the first positional argument".format(seqid, codons[0], ", ".join(list(map(lambda c: kmer.id_to_kmer(c, 3), POSSIBLE_START_CODONS)))))
    else:
        return True
    
    
def get_codons_in_order(seq, seq_id:str=None, ignore_invalid_cds:bool=False, ignore_stop_codons:bool=False, ignore_start_codons:bool=False):
    if type(seq) is not str:
        raise TypeError("kmerdb.codons.get_codons_in_order() expects a str as its first positional argument")
    elif seq_id is not None and type(seq_id) is not str:
        raise TypeError("kmerdb.codons.get_codons_in_order() expects the keyword argument 'seq_id' to be a bool")
    elif type(ignore_invalid_cds) is not bool:
        raise TypeError("kmerdb.codons.get_codons_in_order() expects the keyword argument 'ignore_invalid_cds' to be a bool")
    # elif type(ignore_stop_codons) is not bool:
    #     raise TypeError("kmerdb.codons.get_codons_in_order() expects the keyword argument 'include_stop_codons' to be a bool")
    # elif type(ignore_start_codons) is not bool:
    #     raise TypeError("kmerdb.codons.get_codons_in_order() expects the keyword argument 'include_start_codons' to be a bool")
    
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
            if j == 2 and i == 2:
                if ignore_start_codons is True and kmer.kmer_to_id(codon) in POSSIBLE_START_CODONS: # If we should ignore the start codon counts and it is a valid start codon, skip this codon count and continue
                    continue
                else:
                    codons.append(kmer.kmer_to_id(codon))
                    codon = ''
                    i = 0
            elif  j == seq_len - 1:
                if kmer.kmer_to_id(codon) in STOP_CODONS:
                    #raise ValueError("Internal stop codon identified: {0}\n\nOffending sequence: {1}".format(codon, seq))
                    """
                    What should I do with a stop codon? Include it?
                    """
                    if ignore_stop_codons is True: # If we should not include them, break in the last loop
                        break
                    else: # Otherwise, include the valid stop codon as the last codon in the sequence
                        codons.append(kmer.kmer_to_id(codon)) # Comment if you don't like including stop codons

                elif kmer.kmer_to_id(codon) not in STOP_CODONS:
                    if ignore_invalid_cds is True: # If we are just ignoring invalid CDS, return None and do nothing
                        return None
                    else: # Otherwise throw an error
                        logger.warning("Ignore invalid CDS? {0}".format(ignore_invalid_cds))
                        msg = "Invalid stop codon identified: {0}. Expected to be one of: {1}\n\nOffending sequence: ID: {2}\n{3}".format(codon, ", ".join(list(map(lambda c: kmer.id_to_kmer(c, 3), STOP_CODONS))), seq_id, seq)
                        logger.error(msg)
                        raise ValueError(msg)
            else: # If this is not the start or stop codon, then convert to codon kmer-id
                codon_id = kmer.kmer_to_id(codon) # Convert to k-mer id before returning
                if codon_id is None:
                    logger.warn("kmerdb.codons.get_codons_in_order() found an unknown nucleic-acid residue {0} at position {1} of input sequence '{2}'...\ncannot translate this sequence\n\n".format(c, j, seq_id))
                    logger.warn("kmerdb.codons.get_codons_in_order() found an unknown amino-acid residue {0} in the input sequence".format(codon))
                    return None
                codons.append(codon_id) # Then add the codon, reset codon to empty and i to 0
                codon = ''
                i = 0
        else:
            # Otherwise we are in the middle of a codon and do not have 3 consecutive bases yet.
            # The letter is appeneded to the codon string at the top of this loop
            i+=1
        
    return codons # Note, terminal codon is stop codon



    


def codon_frequency_table(seq, seqid, ignore_invalid_cds:bool=False, include_stop_codons:bool=False, include_start_codons:bool=False):
    """
    :param seq: A nucleic acid CDS sequence 
    :type str:
    :param seqid: A fasta identifier for the sequence
    :type str:
    :param ignore_invalid_cds: return (None, None, None) if the CDS is invalid
    :type bool:
    :param ignore_stop_codons: do not include codon counts for stop codons
    :type bool:
    :param ignore_start_codons: do not include codon counts for start codons
    :type bool:
    :returns: A typle of (codon_ids, codon_counts, codon_frequencies)
    :rtype: tuple
    Takes an input sequence 'seq' and produces 3-mer frequencies
    Returns a tuple of codon ids, and codon counts, and codon frequencies (wrt length), and codon frequencies within an synonymous AA family
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
            return (None, None, None, None)
        else:
            raise ValueError("kmerdb.codons.codon_frequency_table expects the sequence to be divisible by three")
    
    codons = get_codons_in_order(seq, seq_id=seqid, ignore_invalid_cds=ignore_invalid_cds) #ignore_stop_codons=ignore_stop_codons, ignore_start_codons=ignore_start_codons) # Returns a list of int codon ids, in order of translation frame
    if codons is None:
        logger.warning("kmerdb.codons.codon_frequency_table() cannot return frequencies of invalid CDS")
        return (None, None, None, None)
    else:
        num_codons = len(codons)
        if num_codons == len(seq) / 3  or num_codons != (len(seq) / 3) - 1 or num_codons != (len(seq) / 3) - 2: # Remove start or stop codons or both
            pass
        else:
            raise ValueError("The number of codons, {0}, should be proportional to the length ({1}/3) of sequence '{2}' ".format(num_codons, len(seq), seqid))
    #codon_counts = []
    codon_counts = np.zeros(len(CODON_LIST), dtype="uint32")
    for i, cdn in enumerate(codons):
        if i == 0:
            if include_start_codons is True:
                if cdn in POSSIBLE_START_CODONS:
                    # if cdn == 14 and include_start_codons is True and cdn in POSSIBLE_START_CODONS:
                    #     raise RuntimeError("This should work")
                    codon_counts[cdn] += 1
                elif cdn not in POSSIBLE_START_CODONS:
                    if ignore_invalid_cds is True:
                        logger.warn("Sequence '{0}' has an invalid start codon '{1}' and you chose to ignore invalid CDSes (those without a standard start codon (kmerdb.codons)). Omitting this sequence from the resulting table...".format(seqid, kmer.id_to_kmer(cdn, 3)))
                        return (None, None, None, None) # We are ignoring invalid cds so return Non
                    elif ignore_invalid_cds is False:
                        logger.warn("Sequence '{0}' has an invalid start codon '{1}' but you chose not to ignore invalid CDSes (those without a standard start codon (kmerdb.codons)). Throwing an exception to alert you to the presence of an invalid CDS".format(seqid, kmer.id_to_kmer(cdn, 3)))
                        raise ValueError("Internal Error: Found invalid situation (ignoring certain codons) with {0}th codon '{1}' in sequence '{2}'... ".format(i, kmer.id_to_kmer(cdn, 3), seqid))
            elif include_start_codons is False:
                if cdn in POSSIBLE_START_CODONS:
                    continue # Explicitily omitting the count of the valid start codon
                elif cdn not in POSSIBLE_START_CODONS:
                    if ignore_invalid_cds is True:
                        logger.warn("Sequence '{0}' has an invalid start codon '{1}' and you chose to ignore invalid CDSes (those without a standard start codon (kmerdb.codons)). Omitting this sequence from the resulting table".format(seqid, kmer.id_to_kmer(cdn, 3)))
                        return (None, None, None, None)
                    elif ignore_invalid_cds is False: #
                        logger.warn("Sequence '{0}' has an invalid start codon '{1}' but you chose not to ignore invalid CDSes (those without a standard start codon (kmerdb.codons)). Throwing an exception to alert you to the presence of an invalid CDS".format(seqid, kmer.id_to_kmer(cdn, 3)))
                        raise ValueError("Internal Error: Found invalid situation (ignoring certain codons) with {0}th codon '{1}' in sequence '{2}'... ".format(i, kmer.id_to_kmer(cdn, 3), seqid))
        elif i == num_codons - 1:
            if include_stop_codons is True:
                if cdn in STOP_CODONS:
                    codon_counts[cdn] += 1 # Explicitly including the valid stop codon counts in the table
                elif cdn not in STOP_CODONS:
                    if ignore_invalid_cds is True:
                        logger.warn("Sequence '{0}' has an invalid stop codon '{1}' and you chose to ignore invalid CDSes (those without a standard start codon (kmerdb.codons)). Omitting this sequence from the resulting table".format(seqid, kmer.id_to_kmer(cdn, 3)))
                        return (None, None, None, None) # We are ignoring invalid CDS anyways so return None
                    elif ignore_invalid_cds is False:
                        logger.warn("Sequence '{0}' has an invalid stop codon '{1}' but you chose not to ignore invalid CDSes (those without a standard start codon (kmerdb.codons)). Throwing an exception to alert you to the presence of an invalid CDS".format(seqid, kmer.id_to_kmer(cdn, 3)))
                        raise ValueError("Internal Error: Found invalid situation (ignoring certain codons) with {0}th codon '{1}' in sequence '{2}'... ".format(i, kmer.id_to_kmer(cdn, 3), seqid))
                        
            elif include_stop_codons is False:
                if cdn in STOP_CODONS:
                    continue # Explicitly omitting the count of the valid stop codon
                elif cdn not in STOP_CODONS:
                    if ignore_invalid_cds is True:
                        logger.warn("Sequence '{0}' has an invalid stop codon '{1}' and you chose to ignore invalid CDSes (those without a standard start codon (kmerdb.codons)). Omitting this sequence from the resulting table".format(seqid, kmer.id_to_kmer(cdn, 3)))
                        return (None, None, None, None)

                    elif ignore_invalid_cds is False:
                        logger.warn("Sequence '{0}' has an invalid stop codon '{1}' but you chose not to ignore invalid CDSes (those without a standard start codon (kmerdb.codons)). Throwing an exception to alert you to the presence of an invalid CDS".format(seqid, kmer.id_to_kmer(cdn, 3)))
                        raise ValueError("Internal Error: Found invalid situation (ignoring certain codons) with {0}th codon '{1}' in sequence '{2}'...".format(i, kmer.id_to_kmer(cdn, 3), seqid))
        elif i != 0 and i!= num_codons - 1:
            #codon_counts.append(codons.count(cdn))
            codon_counts[cdn] += 1
        else:
            raise ValueError("Internal Error: Found invalid situation (ignoring certain codons) with {0}th codon '{1}' in sequence '{2}'... ".format(i, kmer.id_to_kmer(cdn, 3), seqid))

    #codon_counts = np.array(codon_counts)
    # Instead of using seq_len as CDS length, use length of codon list as it doesnt include STOP codons.
    """
    Okay, so there are two types of codon frequencies:
    the codon frequencies wrt total sequence length (sum to 1)
    the codon frequencies in the codon family (all codon frequencies in family sum to 1, total of codon frequencies is larger than 1)

    # Codons is a list of codon ids corresponding to the codons in the sequence
    # codon_counts[i] is a single codons counts, len(codon_counts) = 64
    # CODON_LIST[i] is the codon id of the corresponding count in codon_counts
    # codon_table[CODON_LIST[i]] is the amino-acid IUPAC residue code corresponding to the ith codon
    # len(synonymous_codons[codon_table[CODON_LIST[i]]]) is the number of codons for the amino acid corresponding to CODON_LIST[i]
    # codon_counts[i] / len(synonymous_codons[codon_table[CODON_LIST[i]]]) is the frequency of the ith codon in its family.
    """
    codon_frequencies_wrt_length = np.array(codon_counts) / np.sum(codon_counts)
    #codon_frequencies_wrt_length = codon_frequencies_wrt_length.tolist()
    #codon_frequencies_wrt_length = list(map(lambda c: float(c)/np.sum(codon_counts), codon_counts))

    # print("Printing 1th codon kmerid...")
    # print(CODON_LIST[1])
    # print("Printing the 1th codon's amino acid")
    # print(codon_table[CODON_LIST[1]])
    # print("Printing the array corresponding to amino acid 'F'")
    # print(synonymous_codons[codon_table[CODON_LIST[1]]])
    # print("Printing the array of codon ids for the synonymous codons to codon 1")

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
        raise TypeError("kmerdb.codons.get_expected_codon_frequencies expects L to be an int")
    elif not isinstance(codon_counts, pd.DataFrame):
        raise TypeError("kmerdb.codons.get_expected_codon_frequencies expects codon_count_vector to be an Pandas DataFrame")
    elif codon_counts.shape[1] != 64:
        raise ValueError("kmerdb.codons.get_expected_codon_frequencies expects a 64 column data frame as its second positional argument")


    #cdn_column_names = codon_counts.columns.tolist()
    final_frequencies = np.zeros(64)

    """
    A row contains different codon counts for a single gene
    """
    total_cdn_cnts = np.array(codon_counts.sum(axis=0), dtype="float64") # Column sum of codon_counts across all genes
    frac_cdn_cnts = total_cdn_cnts / int(total_cdn_cnts.sum())
    exp_cnts = np.array(frac_cdn_cnts * (L/3), dtype="uint32")
    # print("Sequence length: {0} amino-acids: {1}".format(L, L/3))
    # print("Sum of fractions: {0}".format(frac_cdn_cnts.sum()))
    # print("Sum of counts: {0}".format(exp.sum()))
    # sys.exit(1)
    return exp_cnts, exp_cnts / exp_cnts.sum()
    # for aa, indices in synonymous_codons.items():
    #     cdn_cnts_for_aa = total_cdn_cnts[indices] # Change to .iloc
    #     family_aa_total_cnt = np.sum(cdn_cnts_for_aa)
    #     family_frequencies = cdn_cnts_for_aa / family_aa_total_cnt
    #     for i, idx in enumerate(indices):
    #         final_frequencies[idx] = family_frequencies[i]
    # return final_frequencies, np.array(np.floor(L * final_frequencies), dtype="uint32")
