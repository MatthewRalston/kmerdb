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
import os
import sys
import copy
import random
from itertools import chain, repeat, product

import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__file__)



#############################
#
# Dictionaries/maps
#
#############################

"""
a list of the unicode character encodings for the DNA alphabet
note 65 is A, 67 is C, 71 is G, 84 is T.
"""

letterToBinaryNA={ # Unicode UTF-8 byte codes for standard NA residues: ACGT NOTE: letterToBinaryNA
    65: 0,
    67: 1,
    71: 2,
    84: 3
}
binaryToLetterNA=['A', 'C', 'G', 'T'] # FIXME: binaryToLetterNA
standard_lettersNA=set("ACTG")


"""
IUPAC support mappings for the k-mer counter (NOT USED IN THE ASSEMBLER)
"""
IUPAC_NA_DOUBLETS = {
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"]
}
IUPAC_NA_DOUBLET_CHARS=set(IUPAC_NA_DOUBLETS.keys())
IUPAC_NA_TRIPLETS = {
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"]
}
IUPAC_NA_TRIPLET_CHARS=set(IUPAC_NA_TRIPLETS.keys())
# *extremely* necessary variable


permitted_NA_characters = IUPAC_NA_TRIPLET_CHARS.union(IUPAC_NA_DOUBLET_CHARS).union(standard_lettersNA) # set("ACTGRYSWKMBDHV")
permitted_NA_characters_with_N = permitted_NA_characters.union(set("N"))




"""
Amino acid dictionaries
"""

letterToBinaryAA = {
    65: 0, # A
    67: 1, # C
    68: 2, # D
    69: 3, # E
    70: 4, # F
    71: 5, # G
    72: 6, # H
    73: 7, # I
    75: 8, # K
    76: 9, # L
    77: 10,# M
    78: 11,# N
    80: 12,# P
    81: 13,# Q
    82: 14,# R
    83: 15,# S
    84: 16,# T
    86: 17,# V
    87: 18,# W
    88: 19, # X = Unknown amino acid
    89: 20,# Y
    33: 21, # ! = dummy
    34: 22, # " = dummy
    35: 23, # # = dummy
    36: 24, # $ = dummy
    37: 25, # % = dummy
    38: 26, # & = dummy
    39: 27, # ' = dummy
    40: 28, # (
    41: 29, # ) 
    42: 30, # *
    43: 31, # +
    44: 32, # ,
}


AMINO_ACID_IUPAC_CODES = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
binaryToLetterAA = AMINO_ACID_IUPAC_CODES + ["X", "!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ","] # original 20 elements + 12 dummy characters

standard_lettersAA = set(binaryToLetterAA)

IUPAC_AA_DOUBLETS = {
    "B": ["R", "N"], # Aspartate (R) or Aspartamine (N)
    "Z": ["E", "Q"], # Glutamate (E) or Glutamine (Q)
    "J": ["L", "I"]  # Leucine (L) or Isoleucine (I)
}
IUPAC_AA_DOUBLET_CHARS=set(IUPAC_AA_DOUBLETS)


permitted_AA_characters = standard_lettersAA.union(IUPAC_AA_DOUBLET_CHARS)
permitted_AA_characters_with_X = permitted_AA_characters.union(set("X"))


#############################
#
# Functions
#
#############################




def is_sequence_na(s:str):
    """
    :param s: The sequence to ask is the sequence an amino-acid fasta sequence
    :type s: str
    :raises TypeError: Raises a TypeError if the sequence is not str or Bio.SeqRecord
    :raises ValueError: Raises a ValueError if the sequence contains non-IUPAC letter codes for nucleic acid OR amino acid sequences
    :returns: Whether or not the sequence is an amino-acid using standard IUPAC letter codes
    :rtype: bool

    """
    if type(s) is str:
        letters = set(str(s))
    elif isinstance(s, SeqRecord):
        letters = set(str(s.seq))
    else:
        raise TypeError("kmer.is_sequence_na expects a str/Bio.SeqRecord as its only positional argument")

    if letters.difference(permitted_NA_characters_with_N) == set():
        return True
    elif letters.difference(permitted_AA_characters_extended) == set():
        return False
    else:
        sys.stderr.write("Letters in sequence: {0}\n".format(letters))
        raise ValueError("Unable to determine IUPAC nature of the following sequence\n{0}".format(s))


    
def is_sequence_aa(s:str):
    """
    :param s: The sequence to ask is the sequence an amino-acid fasta sequence
    :type s: str
    :raises TypeError: Raises a TypeError if the sequence is not str or Bio.SeqRecord
    :raises ValueError: Raises a ValueError if the sequence contains non-IUPAC letter codes for nucleic acid OR amino acid sequences
    :returns: Whether or not the sequence is an amino-acid using standard IUPAC letter codes
    :rtype: bool
    """

    if type(s) is str:
        letters = set(str(s))
    elif isinstance(s, SeqRecord):
        letters = set(str(s.seq))
    else:
        raise TypeError("kmer.is_sequence_na expects a str/Bio.SeqRecord as its only positional argument")
    if letters.difference(permitted_NA_characters_with_N) == set():
        return False
    elif letters.difference(permitted_AA_characters_extended) == set():
        return True
    else:
        raise ValueError("Unable to determine IUPAC nature of the following sequence\n{0}".format(s))

def left_neighbors(kmer:str):
    if type(kmer) is not str:
        raise TypeError("kmerdb.kmer.left_neighbors() expected a str as its argument")
    last_char_removed = kmer[:-1]
    return list(map(lambda c: kmer_to_id(c + last_char_removed), binaryToLetterNA))

def right_neighbors(kmer:str):
    if type(kmer) is not str:
        raise TypeError("kmerdb.kmer.right_neighbors() expected a str as its argument")
    first_char_removed = kmer[1:]
    return list(map(lambda c: kmer_to_id(first_char_removed + c), binaryToLetterNA))

def is_right_neighbor(kmer:str, right_neighbor):
    return kmer[1:] == right_neighbor[:-1]

def is_left_neighbor(kmer:str, left_neighbor):
    return left_neighbor[1:] == kmer[:-1]

def are_neighbors(kmer1:str, kmer2:str):
    if type(kmer1) is not str:
        raise TypeError("kmerdb.kmer.are_neighbors() expected a str as its first positional argument")
    elif type(kmer2) is not str:
        raise TypeError("kmerdb.kmer.are_neighbors() expected a str as its second positional argument")
    return is_right_neighbor(kmer1, kmer2) or is_left_neighbor(kmer1, kmer2)
    

def kmer_to_id(s, is_aa:bool=False):
    """Convert a fixed length k-mer string to the binary encoding parameterized upon that same k

    Note that the conversion of a k-mer string to an id integer
    is consistent regardless of k, 
    because the k is implicit in the k-mer string's size.

    Therefore, this method does not need to be wrapped in the k-mer class

    Acknowledgements for the 'idx = idx << 2' bit-shifting trick goes to the authors of kPAL.

    @article{anvar2014determining,
    title={Determining the quality and complexity of next-generation sequencing data without a reference genome},
    author={Anvar, Seyed Yahya and Khachatryan, Lusine and Vermaat, Martijn and van Galen, Michiel and Pulyakhina, Irina and Ariyurek, Yavuz and Kraaijeveld, Ken and den Dunnen, Johan T and de Knijff, Peter and ’t Hoen, Peter AC and others},
  journal={Genome biology},
    volume={15},
    pages={1--15},
    year={2014},
    publisher={Springer}
    }


    :param s: The input k-mer as string
    :type s: str
    :param is_aa: Is the input an amino acid?
    :type is_aa: bool
    :raises TypeError: if the input (first argument) is not a str
    :raises KeyError: if non-standard (ATCG) character detected
    :raises ValueError: if non-IUPAC single letter codes for nucleic-acids OR amino-acids are detected
    :returns: The kPal-inspired binary encoding (thanks!)
    :rtype: int

    """
        # and not isinstance(s, Bio.Seq.Seq)  # Shortened the checker as this is a common function.
    if type(s) is not str and not isinstance(s, SeqRecord) and not isinstance(s, Seq): # Typecheck the input k-mer
        logger.error("Seq: {0}, type: {1}".format(s, type(s)))
        raise TypeError("kmerdb.kmer.kmer_to_id expects a str or Bio.SeqRecord.SeqRecord as its positional argument")
    """
    Determine if sequence is NA or AA
    Offloaded to a bool
    
    Omitting the is_aa bool typecheck as this is a frequently used function
    """
    idx = 0
    if is_aa is True and s.find("X") != -1:
        return None # Explicitly return None if an unknown amino acid is found
    elif is_aa is False and s.find('N') != -1: # k-mer with 'N' do not have a binary encoding
        logger.debug(TypeError("kdb.kmer.kmer_to_id expects the letters to contain only nucleotide symbols ATCG"))
        return None

    # print("Is sequence NA? {0}".format(is_sequence_na(str(s))))
    # print("Is sequence AA? {0}".format(is_sequence_aa(str(s))))
        
    if is_aa is False and is_sequence_na(str(s)) is True: #and is_sequence_aa(str(s)) is False:
        #logger.debug("k-mer '{0}' is determined as an nucleic acid sequence".format(s))
        pass
    elif is_aa is True and is_sequence_aa(str(s)) is True: #and is_sequence_na(str(s)) is False:
        #logger.debug("k-mer '{0}' is determined as an amino acid sequence".format(s))
        pass
    else:
        raise ValueError("Could not determine whether sequence was amino acid or nucleic acid: \n{0}".format(str(s)))

    if isinstance(s, Seq) or isinstance(s, SeqRecord):
        s = str(s)        
    for c in bytes(s, "UTF-8"): # Use byteshifting for fast conversion to binary encoding
        #print("Full sequence: {0}".format(s))
        #print("Character: {0}".format(c))
        if is_aa is True:
            idx = idx << 5
            try:
                idx = idx | letterToBinaryAA[c]
            except KeyError as e:
                sys.stderr.write("Entire sequence: {0}".format(s) + "\n")
                sys.stderr.write("Problematic character: {0}".format(c) + "\n")
                raise e
        elif is_aa is False:
            idx = idx << 2
            try:
                idx = idx | letterToBinaryNA[c]
            except KeyError as e:
                sys.stderr.write("Entire sequence: {0}".format(s) + "\n")
                sys.stderr.write("Problematic character: {0}".format(c) + "\n")
                raise e
    if idx is None:
        raise ValueError("kmer.kmer_to_id produced an invalid k-mer id")
    return idx


def id_to_kmer(id, k, is_aa:bool=False):
    """
    Convert an id_to_kmer. I don't understand this docstring's purpose.

    :param id: The int id is the input to the id_to_kmer conversion
    :type id: int
    :param k: The int k is used to byte convert the id to the sequence of characters
    :type k: int
    :param is_aa: The sequence is expected (unchecked) to be a amino-acid sequence
    :raises TypeError: if input (first argument) kmer-id encoding is not an int
    :raises TypeError: if input (second argument) of choice of k is not an int
    :raises ValueError: Internal error: if the k-mer length from the id is erroneously 0.
    :raises ValueError: Internal error: if the k-mer length from the id is erroneously not equal to k.
    :returns: The k-mer as string
    :rtype: str
    """
    if type(id) is not int:
        sys.stderr.write("kmer.id_to_kmer was given the following type as an id\n")
        sys.stderr.write(str(type(id)) + "\n")
        raise TypeError("kmerdb.id_to_kmer expects an int as its first positional argument")
    elif type(k) is not int:
        sys.stderr.write(str(type(k)) + "\n")
        raise TypeError("kmerdb.id_to_kmer expects an int as its second positional argument")
    kmer = ""
    for i in range(k):
        if is_aa is False:
            kmer += binaryToLetterNA[id & 0x03]
            id = id >> 2
        elif is_aa is True:
            kmer += binaryToLetterAA[id & 0x1F] # 1F otherwise
            id = id >> 5
    kmer = list(kmer)
    kmer_len = len(kmer)
    if kmer_len == 0:
        raise ValueError("kmer.id_to_kmer returned an empty k-mer. The kmer-id was '{0}', the choice of k is {1}".format(id, k))
    elif kmer_len != k:
        raise ValueError("kmer.id_to_kmer encountered an inconsistency error. Expected lenght of k-mer sequence was {0}".format(k, kmer_len))
    #just_reversed = kmer.reverse()
    kmer.reverse()
    return ''.join(kmer)


def neighbors(kmer:str, kmer_id:int,  k:int, quiet:bool=True):
    """
    3/11/24 revived. given a k-mer of length k, give its neighbors. This is so ugly.

    But rock on :
    :param kmer: The sequence as a string that will be sliced for k-mers
    :type kmer: str
    :param kmer_id: The k-mer id for neighbor generation
    :type kmer_id: int
    :param k: The int k 
    :type k: int
    :raises TypeError: if the input (first argument) k-mer is not a str or Bio.SeqRecord
    :raises TypeError: if the input (second argument) kmer_id is not an int
    :raises TypeError: if the input (third argument) k is not an int
    :raises ValueError: if the input kmer has lenght not equal to k. Shouldn't happen
    :returns: A list of 'neighboring' or 'adjacent' k-mer ids from derived from the input k-mer
    :rtype: list
    """
    if not isinstance(kmer, str):
        raise TypeError("kmerdb.kmer.neighbors expects a Biopython Seq object as its first positional argument")
    elif type(kmer_id) is not int:
        raise TypeError("kmerdb.kmer.neighbors expects an int as its second positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.kmer.neighbors expects an int as its third positional argument")
    elif len(kmer) != k:
        raise ValueError("kmerdb.kmer.neighbors cannot calculate the {0}-mer neighbors of a {1}-mer".format(k, len(s)))
    else:
        left_neighbor_kmers = left_neighbors(kmer)
        right_neighbor_kmers = right_neighbors(kmer)
        """
        # TYPE 1: [first char removed ... + ... c  : ['A', "c", "g", "T"]]
        # TYPE 2: [["A", "C", "G", "T"] : c + last char removed  ]
        """
        right_neighbor_ids = list(map(kmer_to_id, right_neighbor_kmers))
        left_neighbor_ids = list(map(kmer_to_id, left_neighbor_kmers))
#         logger.debug(""" flower garden - joan G. Stark fir sur. fleicitaciones
#                                     wWWWw
#    vVVVv (___) wWWWw  wWWWw  (___)  vVVVv
#    (___)  ~Y~  (___)  vVVVv   ~Y~   (___)
#     ~Y~   \|    ~Y~   (___)    |/    ~Y~
#     \|   \ |/   \| /  \~Y~/   \|    \ |/
#    \\|// \\|// \\|/// \\|//  \\|// \\\|///
# jgs^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#             """,
#             """   computing the neighbor structure uwu ...
#                   _(_)_                          wWWWw   _
#       @@@@       (_)@(_)   vVVVv     _     @@@@  (___) _(_)_
#      @@()@@ wWWWw  (_)\    (___)   _(_)_  @@()@@   Y  (_)@(_)
#       @@@@  (___)     `|/    Y    (_)@(_)  @@@@   \|/   (_)\
#        /      Y       \|    \|/    /(_)    \|      |/      |
#     \ |     \ |/       | / \ | /  \|/       |/    \|      \|/
# jgs \\|//   \\|///  \\\|//\\\|/// \|///  \\\|//  \\|//  \\\|// 
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# """,
#               "",
#         )

        
        if quiet is not True:
            logger.debug("kmerdb.kmer.neighbors creating neighbors for kmer_id : {0}\nkmer : \"  {1}  \"\nneighbors : \n\n{2}\n{3}\nids: \n\n{4}\n{5}".format(kmer_id, kmer, left_neighbors, right_neighbors, left_neighbor_ids, right_neighbor_ids))
        return left_neighbor_ids + right_neighbor_ids


def validate_seqRecord_and_detect_IUPAC(seqRecord:SeqRecord, k:int, quiet_iupac_warning:bool=True):
    """
    Helper method for validating seqRecord and warnings for non-standard IUPAC residues for nucleic-acid and amino-acid sequences
    
    :param seqRecord: a BioPython SeqRecord object containing a nucleic acid string
    :type seqRecord: Bio.SeqRecord.SeqRecord
    :param k: The choice of k being used to shred k-mers
    :type int
    :param quiet_iupac_warning: verbosity parameter
    :type quiet_iupac_warning: bool
    :raises TypeError: if quiet_iupac_warning keyword argument is not a bool
    :raises TypeError: if the input (first argument) is not a Bio.SeqRecord
    :raises ValueError: if *non* IUPAC symbols are detected in the seqRecord object.
    :raises ValueError: if the input (first argument) is shorter than k
    :returns: a set of the letters detected, and length of the validated Bio.SeqRecord sequence
    :rtype: tuple


    """
    if type(quiet_iupac_warning) is not bool:
        raise TypeError("kmerdb.kmer.Kmers.validate_seqRecord_and_detect_IUPAC expects keyword argument 'quiet_iupac_warning' to be a bool")
    elif isinstance(seqRecord, SeqRecord) is False:
        raise TypeError("kmerdb.kmer.Kmers.validate_seqRecord_and_detect_IUPAC expects a Bio.SeqRecord.SeqRecord object as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.kmer.validate_seqRecord_and_detect_IUPAC expects an int as its second positional argument")
    letters = set(seqRecord.seq)
    seqlen = len(seqRecord.seq)

        
    if seqlen < k:
        logger.error("Offending sequence ID: {0}".format(seqRecord.id))
        raise ValueError("kmerdb.kmer.validate_seqRecord_and_detect_IUPAC expects that each input sequence is longer than k.")

    if is_sequence_na(seqRecord) is True and is_sequence_aa(seqRecord) is False:
        extended_iupac_symbols = list(letters.intersection(permitted_NA_characters) - standard_lettersNA)
        all_non_iupac_symbols = list(letters - permitted_NA_characters_with_N)
    elif is_sequence_aa(seqRecord) is True and is_sequence_na(seqRecord) is False:
        extended_iupac_symbols = list(letters.intersection(permitted_AA_characters) - standard_lettersAA)
        all_non_iupac_symbols = list(letters - permitted_AA_characters_extended)
        
    if len(all_non_iupac_symbols) > 0:
        logger.warning(all_non_iupac_symbols)
        raise ValueError("Non-IUPAC symbols detected in the sequence '{0}'".format(seqRecord.id))
    elif len(extended_iupac_symbols) > 0: # FIXME:  7/5/25 Unknown if fixed. The module is working correctly so far
        if quiet_iupac_warning is False:
            logger.warning("Will completely refuse to include k-mers with 'N' (or 'X' for AA sequences)")
            logger.warning("All counts for k-mers including N will be discarded")
            logger.warning("Other IUPAC symbols will be replaced with their respective pairs/triads")
            logger.warning("And a count will be given to each, rather than a half count")
        elif quiet_iupac_warning is True:
            logger.warning("Suppressing warning that non-standard IUPAC residues (including N/X) are detected.")
    return (letters, extended_iupac_symbols, all_non_iupac_symbols, seqlen)





def shred(seq:SeqRecord, k:int, replace_with_none:bool=False, quiet_iupac_warning:bool=True):
    """
    Take a seqRecord fasta/fastq object and slice according to the IUPAC charset.
    Doublets become replace with two counts, etc.
    If replace_with_none:bool is True, then kmer_id of None is appended to kmer_id for ambiguous IUPAC characters
    Otherwise, a random k-mer id from the expansion is selected and added to the kmer_ids array. Other ids are collected in the bonus_kmer_ids array

    :param seqRecord:
    :type Bio.SeqRecord.SeqRecord:
    :raise ValueError: Non IUPAC characters detected
    :raise ValueError: Internal Error: an invalid non-standard IUPAC character was not categoried properly
    :returns: a 3-tuple of k-mer ids, sequence ids, and positions
    :rtype: tuple
    """
    if type(k) is not int:
        raise TypeError("kmerdb.kmer.shred() expects an int as its second positional argument")
    elif not isinstance(seq, SeqRecord):
        raise TypeError("kmerdb.kmer.shred() expects a Bio.SeqRecord.SeqRecord as its first positional argument")
    seq_id = seq.id
    is_na = is_sequence_na(seq)
    is_aa = is_sequence_aa(seq)
    letters, iupac, non_iupac, seqlen = validate_seqRecord_and_detect_IUPAC(seq, k, quiet_iupac_warning=quiet_iupac_warning)
    #kmers = [seq.seq[i:(i+k)] for i in range(seqlen - k + 1)]
    seq_ids = []
    kmer_ids = []
    data = []
    for i in range(seqlen - k + 1):
        kmer = str(seq.seq[i:(i+k)])
        kmer_id = kmer_to_id(kmer) # NOTE: Could be None

        if is_na is True:
            nonstandard = set(kmer) - standard_lettersNA
            num_nonstandard = len(nonstandard)
            n_count = kmer.count("N")
            num_all_nonstandard = n_count + num_nonstandard

            if len(non_iupac) > 0:
                raise ValueError("Non-IUPAC symbols detected in '{0}'".format(seq_id))
            if num_all_nonstandard == 0:
                data.append({"seq_id": seq_id, "kmer": kmer, "kmer_id": kmer_id, "pos": i})
                continue
            elif num_all_nonstandard > 0:
                if replace_with_none is True:
                    kmer_ids.append(None)
                    continue
                """
                Here I expand all Nucleic Acid doublets, triplets, and N residues 
                """

                # Start by only creating substitutions for base kmer
                doublets_subd = _substitute_na_doublets(kmer)
                # Now substitute triplets
                triplets_subd = []
                for _kmer in doublets_subd:
                    for t in _substitute_na_triplets(_kmer):
                        triplets_subd.append(t)
                # Finally, substitute N residues
                all_subd = []
                for _kmer in triplets_subd:
                    for n in substitute_residue_with_chars(_kmer, "N", list(standard_lettersNA)):
                        all_subd.append(n)
                #all_subd = [n for n in substitute_residue_with_chars(_kmer, "N", standard_lettersNA) for _kmer in triplets_subd]
                d = []
                for _kmer in all_subd:
                    d.append({"seq_id": seq_id, "kmer": _kmer, "kmer_id": kmer_to_id(_kmer), "pos": i})
                data += d

        elif is_aa is True:
            raise RuntimeError("Internal error: Unimplemented for amino-acids")
        else:
            raise ValueError("Internal error: Unknown situation regarding kmer {0} and IUPAC expansion module")


    _kmer_ids = [_kmer["kmer_id"] for _kmer in data]
    _seq_ids = [_kmer["seq_id"] for _kmer in data]
    _pos = [_kmer["pos"] for _kmer in data]
    #print("Num kmer ids: {0}\nNum seq ids: {1}\nNum pos: {2}".format(len(_kmer_ids), len(_seq_ids), len(_pos)))
    return _kmer_ids, _seq_ids, _pos








def substitute_residue_with_chars(kmer:str, char:str, substitutions:list):
    """
    Substitute all position of residue 'char' with the provided 'substitutions' combinatorially.

    """
    if type(kmer) is not str:
        raise TypeError("kmerdb.kmer.substitute_residue_with_chars() expects a str as its first argument")
    elif type(char) is not str or len(char) > 1:
        raise TypeError("kmerdb.kmer.substitute_residue_with_chars() expects a single character as its second argument")
    elif type(substitutions) is not list or not all(type(c) is str and len(c) == 1 for c in substitutions):
        raise TypeError("kmerdb.kmer.substitute_residue_with_chars() expects a list of single characters as its third argument")

    kmers = []
    to_replace = ""
    num_replacements = kmer.count(char)
    if num_replacements == 0:
        return [kmer]
    if num_replacements == 1:
        return list(map(lambda c: kmer.replace("N", c), substitutions))
    else:
        for j in range(num_replacements):
            if to_replace == '':
                to_replace = str(substitutions)
            else:
                to_replace += ", " + str(substitutions)
        try:
            prod = list(eval("product({0})".format(to_replace)))
            for to_replace in prod: # List of characters to replace R with
                _nukmer = kmer
                for c in to_replace: # Characters substituted with R
                    _nukmer = _nukmer.replace(char, c, 1)
                assert _nukmer.count(char) == 0, "Internal error: kmer '{0}' has multiple ambiguous residues remaining after substitution".format(_nukmer)
                kmers.append(_nukmer)
        except Exception as e:
            raise e
    return kmers







    
def _substitute_na_doublets(kmer:str):
    to_sub = ''
    for c in "RYSWKM":
        if c in kmer:
            to_sub += c
    match to_sub:
        case "":
            return [kmer]
        case "R":
            return substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])
        case "Y":
            return substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])
        case "S":
            return substitute_residue_with_chars(kmer, "S", IUPAC_NA_DOUBLETS["S"])
        case "W":
            return substitute_residue_with_chars(kmer, "W", IUPAC_NA_DOUBLETS["W"])
        case "K":
            return substitute_residue_with_chars(kmer, "K", IUPAC_NA_DOUBLETS["K"])
        case "M":
            return substitute_residue_with_chars(kmer, "M", IUPAC_NA_DOUBLETS["M"])
        case "RY":
            return [ry for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RS":
            return [rs for rs in substitute_residue_with_chars(r, "S", IUPAC_NA_DOUBLETS["S"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RW":
            return [rw for rw in substitute_residue_with_chars(r, "W", IUPAC_NA_DOUBLETS["W"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RK":
            return [rk for rk in substitute_residue_with_chars(r, "K", IUPAC_NA_DOUBLETS["K"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RM":
            return [rm for rm in substitute_residue_with_chars(r, "M", IUPAC_NA_DOUBLETS["M"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "YS":
            return [ys for ys in substitute_residue_with_chars(y, "S", IUPAC_NA_DOUBLETS["S"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YW":
            return [yw for yw in substitute_residue_with_chars(y, "W", IUPAC_NA_DOUBLETS["W"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YK":
            return [yk for yk in substitute_residue_with_chars(y, "K", IUPAC_NA_DOUBLETS["K"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YM":
            return [ym for ym in substitute_residue_with_chars(y, "M", IUPAC_NA_DOUBLETS["M"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "SW":
            return [sw for sw in substitute_residue_with_chars(s, "W", IUPAC_NA_DOUBLETS["W"]) for s in substitute_residue_with_chars(kmer, "S", IUPAC_NA_DOUBLETS["S"])]
        case "SK":
            return [sk for sk in substitute_residue_with_chars(s, "K", IUPAC_NA_DOUBLETS["K"]) for s in substitute_residue_with_chars(kmer, "S", IUPAC_NA_DOUBLETS["S"])]
        case "SM":
            return [sm for sm in substitute_residue_with_chars(s, "M", IUPAC_NA_DOUBLETS["M"]) for s in substitute_residue_with_chars(kmer, "S", IUPAC_NA_DOUBLETS["S"])]
        case "WK":
            return [wk for wk in substitute_residue_with_chars(w, "K", IUPAC_NA_DOUBLETS["K"]) for w in substitute_residue_with_chars(kmer, "W", IUPAC_NA_DOUBLETS["W"])]
        case "WM":
            return [wm for wm in substitute_residue_with_chars(w, "M", IUPAC_NA_DOUBLETS["M"]) for w in substitute_residue_with_chars(kmer, "W", IUPAC_NA_DOUBLETS["W"])]
        case "KM":
            return [km for km in substitute_residue_with_chars(k, "M", IUPAC_NA_DOUBLETS["M"]) for k in substitute_residue_with_chars(kmer, "K", IUPAC_NA_DOUBLETS["K"])]



        
        case "RYS":
            return [rys for rys in substitute_residue_with_chars(ry, "S", IUPAC_NA_DOUBLETS["S"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RYW":
            return [ryw for ryw in substitute_residue_with_chars(ry, "W", IUPAC_NA_DOUBLETS["W"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RYK":
            return [ryk for ryk in substitute_residue_with_chars(ry, "K", IUPAC_NA_DOUBLETS["K"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RYM":
            return [rym for rym in substitute_residue_with_chars(ry, "M", IUPAC_NA_DOUBLETS["M"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RSW":
            return [rsw for rsw in substitute_residue_with_chars(rs, "W", IUPAC_NA_DOUBLETS["W"]) for rs in substitute_residue_with_chars(r, "S", IUPAC_NA_DOUBLETS["S"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RSK":
            return [rsk for rsk in substitute_residue_with_chars(rs, "K", IUPAC_NA_DOUBLETS["K"]) for rs in substitute_residue_with_chars(r, "S", IUPAC_NA_DOUBLETS["S"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RSM":
            return [rsm for rsm in substitute_residue_with_chars(rs, "M", IUPAC_NA_DOUBLETS["M"]) for rs in substitute_residue_with_chars(r, "S", IUPAC_NA_DOUBLETS["S"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RWK":
            return [rwk for rwk in substitute_residue_with_chars(rw, "K", IUPAC_NA_DOUBLETS["K"]) for rw in substitute_residue_with_chars(r, "W", IUPAC_NA_DOUBLETS["W"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RWM":
            return [rwm for rwm in substitute_residue_with_chars(rw, "K", IUPAC_NA_DOUBLETS["K"]) for rw in substitute_residue_with_chars(r, "W", IUPAC_NA_DOUBLETS["W"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RKM":
            return [rkm for rkm in substitute_residue_with_chars(rk, "M", IUPAC_NA_DOUBLETS["M"]) for rk in substitute_residue_with_chars(r, "K", IUPAC_NA_DOUBLETS["K"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "YSW":
            return [ysw for ysw in substitute_residue_with_chars(ys, "W", IUPAC_NA_DOUBLETS["W"]) for ys in substitute_residue_with_chars(y, "S", IUPAC_NA_DOUBLETS["S"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YSK":
            return [ysk for ysk in substitute_residue_with_chars(ys, "K", IUPAC_NA_DOUBLETS["K"]) for ys in substitute_residue_with_chars(y, "S", IUPAC_NA_DOUBLETS["S"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YSM":
            return [ysm for ysm in substitute_residue_with_chars(ys, "M", IUPAC_NA_DOUBLETS["M"]) for ys in substitute_residue_with_chars(y, "S", IUPAC_NA_DOUBLETS["S"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YWK":
            return [ywk for ywk in substitute_residue_with_chars(yw, "K", IUPAC_NA_DOUBLETS["K"]) for yw in substitute_residue_with_chars(y, "W", IUPAC_NA_DOUBLETS["W"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YWM":
            return [ywm for ywm in substitute_residue_with_chars(yw, "M", IUPAC_NA_DOUBLETS["M"]) for yw in substitute_residue_with_chars(y, "W", IUPAC_NA_DOUBLETS["W"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YKM":
            return [ykm for ykm in substitute_residue_with_chars(yk, "M", IUPAC_NA_DOUBLETS["M"]) for yk in substitute_residue_with_chars(y, "K", IUPAC_NA_DOUBLETS["K"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "SWK":
            return [swk for swk in substitute_residue_with_chars(sw, "K", IUPAC_NA_DOUBLETS["K"]) for sw in substitute_residue_with_chars(s, "W", IUPAC_NA_DOUBLETS["W"]) for s in substitute_residue_with_chars(kmer, "S", IUPAC_NA_DOUBLETS["S"])]
        case "SWM":
            return [swm for swm in substitute_residue_with_chars(sw, "M", IUPAC_NA_DOUBLETS["M"]) for sw in substitute_residue_with_chars(s, "W", IUPAC_NA_DOUBLETS["W"]) for s in substitute_residue_with_chars(kmer, "S", IUPAC_NA_DOUBLETS["S"])]
        case "SKM":
            return [skm for skm in substitute_residue_with_chars(sk, "M", IUPAC_NA_DOUBLETS["M"]) for sk in substitute_residue_with_chars(s, "K", IUPAC_NA_DOUBLETS["K"]) for s in substitute_residue_with_chars(kmer, "S", IUPAC_NA_DOUBLETS["S"])]
        case "WKM":
            return [wkm for wkm in substitute_residue_with_chars(wk, "M", IUPAC_NA_DOUBLETS["M"]) for wk in substitute_residue_with_chars(w, "K", IUPAC_NA_DOUBLETS["K"]) for w in substitute_residue_with_chars(kmer, "W", IUPAC_NA_DOUBLETS["W"])]


        
        case "RYSW":
            return [rysw for rysw in substitute_residue_with_chars(rys, "W", IUPAC_NA_DOUBLETS["W"]) for rys in substitute_residue_with_chars(ry, "S", IUPAC_NA_DOUBLETS["S"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RYSK":
            return [rysk for rysk in substitute_residue_with_chars(rys, "K", IUPAC_NA_DOUBLETS["K"]) for rys in substitute_residue_with_chars(ry, "S", IUPAC_NA_DOUBLETS["S"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RYSM":
            return [rysm for rysm in substitute_residue_with_chars(rys, "M", IUPAC_NA_DOUBLETS["M"]) for rys in substitute_residue_with_chars(ry, "S", IUPAC_NA_DOUBLETS["S"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RYWK":
            return [rywk for rywk in substitute_residue_with_chars(ryw, "K", IUPAC_NA_DOUBLETS["K"]) for ryw in substitute_residue_with_chars(ry, "W", IUPAC_NA_DOUBLETS["W"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RYWM":
            return [rywm for rywm in substitute_residue_with_chars(ryw, "M", IUPAC_NA_DOUBLETS["M"]) for ryw in substitute_residue_with_chars(ry, "W", IUPAC_NA_DOUBLETS["W"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RYKM":
            return [rykm for rykm in substitute_residue_with_chars(ryk, "M", IUPAC_NA_DOUBLETS["M"]) for ryk in substitute_residue_with_chars(ry, "W", IUPAC_NA_DOUBLETS["W"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RSWK":
            return [rswk for rswk in substitute_residue_with_chars(rsw, "K", IUPAC_NA_DOUBLETS["K"]) for rsw in substitute_residue_with_chars(rs, "W", IUPAC_NA_DOUBLETS["W"]) for rs in substitute_residue_with_chars(r, "S", IUPAC_NA_DOUBLETS["S"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RSWM":
            return [rswm for rswm in substitute_residue_with_chars(rsw, "M", IUPAC_NA_DOUBLETS["M"]) for rsw in substitute_residue_with_chars(rs, "W", IUPAC_NA_DOUBLETS["W"]) for rs in substitute_residue_with_chars(r, "S", IUPAC_NA_DOUBLETS["S"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RSKM":
            return [rskm for rskm in substitute_residue_with_chars(rsk, "M", IUPAC_NA_DOUBLETS["M"]) for rsk in substitute_residue_with_chars(rs, "W", IUPAC_NA_DOUBLETS["W"]) for rs in substitute_residue_with_chars(r, "S", IUPAC_NA_DOUBLETS["S"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "RWKM":
            return [rwkm for rwkm in substitute_residue_with_chars(rwk, "M", IUPAC_NA_DOUBLETS["M"]) for rwk in substitute_residue_with_chars(rw, "K", IUPAC_NA_DOUBLETS["K"]) for rw in substitute_residue_with_chars(r, "W", IUPAC_NA_DOUBLETS["W"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case "YSWK":
            return [yswk for yswk in substitute_residue_with_chars(ysw, "K", IUPAC_NA_DOUBLETS["K"]) for ysw in substitute_residue_with_chars(ys, "W", IUPAC_NA_DOUBLETS["W"]) for ys in substitute_residue_with_chars(y, "S", IUPAC_NA_DOUBLETS["S"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YSWM":
            return [yswm for yswm in substitute_residue_with_chars(ysw, "M", IUPAC_NA_DOUBLETS["M"]) for ysw in substitute_residue_with_chars(ys, "W", IUPAC_NA_DOUBLETS["W"]) for ys in substitute_residue_with_chars(y, "S", IUPAC_NA_DOUBLETS["S"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YSKM":
            return [yskm for yskm in substitute_residue_with_chars(ysk, "M", IUPAC_NA_DOUBLETS["M"]) for ysk in substitute_residue_with_chars(ys, "K", IUPAC_NA_DOUBLETS["K"]) for ys in substitute_residue_with_chars(y, "S", IUPAC_NA_DOUBLETS["S"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "YWKM":
            return [ywkm for ywkm in substitute_residue_with_chars(ywk, "M", IUPAC_NA_DOUBLETS["M"]) for ywk in substitute_residue_with_chars(yw, "K", IUPAC_NA_DOUBLETS["K"]) for yw in substitute_residue_with_chars(y, "W", IUPAC_NA_DOUBLETS["W"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]
        case "SWKM":
            return [swkm for swkm in substitute_residue_with_chars(swk, "M", IUPAC_NA_DOUBLETS["M"]) for swk in substitute_residue_with_chars(sw, "K", IUPAC_NA_DOUBLETS["K"]) for sw in substitute_residue_with_chars(s, "W", IUPAC_NA_DOUBLETS["W"]) for s in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]

        case 'RYSWK':
            return [ryswk for ryswk in substitute_residue_with_chars(rysw, "K", IUPAC_NA_DOUBLETS["K"]) for rysw in substitute_residue_with_chars(rys, "W", IUPAC_NA_DOUBLETS["W"]) for rys in substitute_residue_with_chars(ry, "S", IUPAC_NA_DOUBLETS["S"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]


        
        case 'RYSWM':
            return [ryswm for ryswm in substitute_residue_with_chars(rysw, "M", IUPAC_NA_DOUBLETS["M"]) for rysw in substitute_residue_with_chars(rys, "W", IUPAC_NA_DOUBLETS["W"]) for rys in substitute_residue_with_chars(ry, "S", IUPAC_NA_DOUBLETS["S"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case 'RYSKM':
            return [ryskm for ryskm in substitute_residue_with_chars(rysk, "M", IUPAC_NA_DOUBLETS["M"]) for rysk in substitute_residue_with_chars(rys, "K", IUPAC_NA_DOUBLETS["K"]) for rys in substitute_residue_with_chars(ry, "S", IUPAC_NA_DOUBLETS["S"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case 'RYWKM':
            return [rywkm for rywkm in substitute_residue_with_chars(rywk, "M", IUPAC_NA_DOUBLETS["M"]) for rywk in substitute_residue_with_chars(ryw, "K", IUPAC_NA_DOUBLETS["K"]) for ryw in substitute_residue_with_chars(ry, "W", IUPAC_NA_DOUBLETS["W"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case 'RSWKM':
            return [rswkm for rswkm in substitute_residue_with_chars(rswk, "M", IUPAC_NA_DOUBLETS["M"]) for rswk in substitute_residue_with_chars(rsw, "K", IUPAC_NA_DOUBLETS["K"]) for rsw in substitute_residue_with_chars(rs, "K", IUPAC_NA_DOUBLETS["K"]) for rs in substitute_residue_with_chars(r, "S", IUPAC_NA_DOUBLETS["S"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case 'YSWKM':
            return [yswkm for yswkm in substitute_residue_with_chars(yswk, "M", IUPAC_NA_DOUBLETS["M"]) for yswk in substitute_residue_with_chars(ysw, "k", IUPAC_NA_DOUBLETS["K"]) for ysw in substitute_residue_with_chars(ys, "W", IUPAC_NA_DOUBLETS["W"]) for ys in substitute_residue_with_chars(y, "S", IUPAC_NA_DOUBLETS["S"]) for y in substitute_residue_with_chars(kmer, "Y", IUPAC_NA_DOUBLETS["Y"])]


        case "RYSWKM":
            return [ryswkm for ryswkm in substitute_residue_with_chars(ryswk, "M", IUPAC_NA_DOUBLETS["M"]) for ryswk in substitute_residue_with_chars(rysw, "K", IUPAC_NA_DOUBLETS["K"]) for rysw in substitute_residue_with_chars(rys, "W", IUPAC_NA_DOUBLETS["W"]) for rys in substitute_residue_with_chars(ry, "S", IUPAC_NA_DOUBLETS["S"]) for ry in substitute_residue_with_chars(r, "Y", IUPAC_NA_DOUBLETS["Y"]) for r in substitute_residue_with_chars(kmer, "R", IUPAC_NA_DOUBLETS["R"])]
        case _:
            raise ValueError("Internal error: substitution of all positions in kmer '{0}' with residues '{1}' failed...".format(kmer, to_sub))














def _substitute_aa_doublets(kmer:str):
    to_sub = ''
    for c in "JBZ":
        if c in kmer:
            to_sub += c
    match to_sub:
        case "":
            return [kmer]
        case "JB":
            return [jb for jb in substitute_residue_with_chars(j, "B", IUPAC_AA_DOUBLETS["B"]) for j in substitute_residue_with_chars(kmer, "J", IUPAC_AA_DOUBLETS["J"])]
        case "JZ":
            return [jz for jz in substitute_residue_with_chars(j, "Z", IUPAC_AA_DOUBLETS["Z"]) for j in substitute_residue_with_chars(kmer, "J", IUPAC_AA_DOUBLETS["J"])]
        case "BZ":
            return [bz for bz in substitute_residue_with_chars(b, "Z", IUPAC_AA_DOUBLETS["Z"]) for b in substitute_residue_with_chars(kmer, "B", IUPAC_AA_DOUBLETS["B"])]
        case "JBZ":
            return [jbz for jbz in substitute_residue_with_chars(jb, "Z", IUPAC_AA_DOUBLETS["Z"]) for jb in substitute_residue_with_chars(j, "B", IUPAC_AA_DOUBLETS["B"]) for j in substitute_residue_with_chars(kmer, "J", IUPAC_AA_DOUBLETS["J"])]
        case _:
            raise ValueError("Internal Error: substitution of all positions in kmer '{0}' with residues '{1}' failed...".format(kmer, to_sub))
        
def _substitute_na_triplets(kmer:str):
    to_sub = ''
    for c in "BDHV":
        if c in kmer:
            to_sub += c
    match to_sub:
        case "":
            return [kmer]
        case "B":
            return substitute_residue_with_chars(kmer, "B", IUPAC_NA_TRIPLETS["B"])
        case "D":
            return substitute_residue_with_chars(kmer, "D", IUPAC_NA_TRIPLETS["D"])
        case "H":
            return substitute_residue_with_chars(kmer, "H", IUPAC_NA_TRIPLETS["H"])
        case "V":
            return substitute_residue_with_chars(kmer, "V", IUPAC_NA_TRIPLETS["V"])
        case "BD":
            return [bd for bd in substitute_residue_with_chars(b, "D", IUPAC_NA_TRIPLETS["D"]) for b in substitute_residue_with_chars(kmer, "B", IUPAC_NA_TRIPLETS["B"])]
        case "BH":
            return [bh for bh in substitute_residue_with_chars(b, "H", IUPAC_NA_TRIPLETS["H"]) for b in substitute_residue_with_chars(kmer, "B", IUPAC_NA_TRIPLETS["B"])]
        case "BV":
            return [bv for bv in substitute_residue_with_chars(b, "V", IUPAC_NA_TRIPLETS["V"]) for b in substitute_residue_with_chars(kmer, "B", IUPAC_NA_TRIPLETS["B"])]
        case "DH":
            return [dh for dh in substitute_residue_with_chars(d, "H", IUPAC_NA_TRIPLETS["H"]) for d in substitute_residue_with_chars(kmer, "D", IUPAC_NA_TRIPLETS["D"])]
        case "DV":
            return [dv for dv in substitute_residue_with_chars(d, "V", IUPAC_NA_TRIPLETS["V"]) for d in substitute_residue_with_chars(kmer, "D", IUPAC_NA_TRIPLETS["D"])]
        case "HV":
            return [hv for hv in substitute_residue_with_chars(h, "V", IUPAC_NA_TRIPLETS["V"]) for h in substitute_residue_with_chars(kmer, "H", IUPAC_NA_TRIPLETS["H"])]
        case "BDH":
            return [bdh for bdh in substitute_residue_with_chars(bd, "H", IUPAC_NA_TRIPLETS["H"]) for bd in substitute_residue_with_chars(b, "D", IUPAC_NA_TRIPLETS["D"]) for b in substitute_residue_with_chars(kmer, "B", IUPAC_NA_TRIPLETS["B"])]
        case "BDV":
            return [bdh for bdh in substitute_residue_with_chars(bd, "H", IUPAC_NA_TRIPLETS["H"]) for bd in substitute_residue_with_chars(b, "D", IUPAC_NA_TRIPLETS["D"]) for b in substitute_residue_with_chars(kmer, "B", IUPAC_NA_TRIPLETS["B"])]
        case "BHV":
            return [bhv for bhv in substitute_residue_with_chars(bh, "V", IUPAC_NA_TRIPLETS["V"]) for bh in substitute_residue_with_chars(b, "H", IUPAC_NA_TRIPLETS["H"]) for b in substitute_residue_with_chars(kmer, "B", IUPAC_NA_TRIPLETS["B"])]
        case "DHV":
            return [dvh for dvh in substitute_residue_with_chars(dv, "H", IUPAC_NA_TRIPLETS["H"]) for dv in substitute_residue_with_chars(d, "V", IUPAC_NA_TRIPLETS["V"]) for d in substitute_residue_with_chars(kmer, "D", IUPAC_NA_TRIPLETS["D"])]
        case "BDHV":
            return [bdvh for bdvh in substitute_residue_with_chars(bdv, "H", IUPAC_NA_TRIPLETS["H"]) for bdv in substitute_residue_with_chars(bd, "V", IUPAC_NA_TRIPLETS["V"]) for bd in substitute_residue_with_chars(b, "D", IUPAC_NA_TRIPLETS["D"]) for b in substitute_residue_with_chars(kmer, "B", IUPAC_NA_TRIPLETS["B"])]
        case _:
            raise ValueError("Internal Error: substitution of all positions in kmer '{0}' with residues '{1}' failed...".format(kmer, to_sub))

