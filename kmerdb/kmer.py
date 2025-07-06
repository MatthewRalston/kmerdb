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
from itertools import chain, repeat

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


permitted_AA_characters = standard_lettersAA
permitted_AA_characters_extended = standard_lettersAA.union(IUPAC_AA_DOUBLET_CHARS)



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
        all_non_iupac_symbols = list(letters - permitted_NA_characters)
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
    :returns: a 2-tuple of k-mer ids and bonus k-mer ids, depending on replace_with_none parameterization
    :rtype: tuple
    """
    seq_id = seq.id
    is_na = is_sequence_na(seq)
    is_aa = is_sequence_aa(seq)
    letters, iupac, non_iupac, seqlen = validate_seqRecord_and_detect_IUPAC(seq, k, quiet_iupac_warning=quiet_iupac_warning)
    kmers = [seq.seq[i:(i+k)] for i in range(seqlen - k + 1)]

    kmer_ids = []
    bonus_kmer_ids = []
    for i in range(seqlen - k + 1):
        kmer = seq.seq[i:(i+k)]
        if len(non_iupac) > 0:
            raise ValueError("Non-IUPAC symbols detected in '{0}'".format(seq_id))
        elif is_na is True:
            nonstandard_chars = set(kmer) - permitted_NA_characters_with_N
        elif is_aa is True:
            nonstandard_chars = set(kmer) - permitted_AA_characters_with_N

            
        if is_na is True and "N" in iupac:
            logger.warning("Omitting k-mer with ambiguous N residue from kmer '{0}' in sequence '{1}' at position {2} from consideration".format(kmer, seq_id, i))
            kmer_ids.append(None)
            continue
        elif is_aa is True and "X" in iupac:
            logger.warning("Omitting k-mer with ambiguous X residue from kmer '{0}' in sequence '{1}' at position {2} from consideration".format(kmer, seq_id, i))
            kmer_ids.append(None)
            continue

        if len(nonstandard_chars) == 0:
            kmer_ids.append(kmer_to_id(kmer))
        elif len(nonstandard_chars) > 0 and replace_with_none is True:
            kmer_ids.append(None) # Use None for a non-standard, valid IUPAC kmer id under certain conditions
        elif len(nonstandard_chars) == 1 and replace_with_none is False:
            # Otherwise, subsitute doublets or triplets
            logger.warning("Non-standard IUPAC character {0} detected in sequence '{1}' at position {2}: {3}".format(c, seq_id, i, kmer))
            """
            Doublet/triplet expansion and random selection of kmer_id with ambiguous valid IUPAC residue
            """
            if is_na is True:
                if replace_with_none is True:
                    kmer_ids.append(None)
                elif c in IUPAC_NA_DOUBLET_CHARS:
                    d1, d2 = (kmer.replace(c, x) for x in IUPAC_NA_DOUBLETS[c])
                    random_kmer = random.choice( (d1, d2) )
                    other_kmer = d1 if random_kmer == d2 else d2
                    random_kmer_id = kmer_to_id(random_kmer)
                    other_kmer_id = kmer_to_id(other_kmer)
                    
                    logger.warning("Subsituting a non-standard IUPAC doublet residue in k-mer {0}  : {2} and {2}".format(kmer, random_kmer, other_kmer))
                    kmer_ids.append(random_kmer_id)
                    bonus_kmer_ids.append(other_kmer_id)
                elif c in IUPAC_NA_TRIPLET_CHARS:
                    t1, t2, t3 = (kmer.replace(c, x) for x in IUPAC_NA_TRIPLETS[c])
                    random_kmer = random.choice( (t1, t2, t3) )
                    random_kmer_id = kmer_to_id(random_kmer)
                    kmer_ids.append(random_kmer)
                    adds = []
                    add_ids = []
                    for t in (t1, t2, t3):
                        if t != random_kmer:
                            t_kmer_id = kmer_to_id(t)
                            adds.append(t)
                            adds_ids.append(t_kmer_id)
                    adds_ids += add_ids
                    logger.warning("Subsituting a non-standard IUPAC triplet residue in k-mer {0}  : {2} and {2} and {3}".format(kmer, random_kmer, adds[0], adds[1]))
            elif is_aa is True:
                if replace_with_none is True:
                    kmer_ids.append(None)
                elif c in IUPAC_AA_DOUBLET_CHARS:
                    d1, d2 = (kmer.replace(c, x) for x in IUPAC_AA_DOUBLETS[c])
                    random_kmer = random.choice( (d1, d2 ) )
                    other_kmer = d1 if kmer == d2 else d2
                    random_kmer_id = kmer_to_id(random_kmer)
                    other_kmer_id = kmer_to_id(other_kmer)
                    logger.warning("Subsituting a non-standard IUPAC double residue in k-mer {0}  : {2} and {2}".format(kmer, random_kmer, other_kmer))
                    kmer_ids.append(random_kmer_id)
                    bonus_kmer_ids.append(other_kmer_id)
            else:
                raise ValueError("Internal Error: Somehow character {0} from kmer {1} in sequence '{2}' at position {3} wasn't categorized properly as a non-standard IUPAC symbol expansion into a doublet or triplet properly. ".format(c, kmer, seq_id, i))
        elif len(nonstandard_chars) > 1:
            logger.warning("Multiple characters are ambiguous in kmer '0' and prohibit simple replacement strategies for doublet or triplet expansion of this kmer. Too much ambiguity, using k-mer id 'None' and moving on...")
            continue
        else:
            raise ValueError("Internal error: Unknown situation regarding kmer {0} and IUPAC expansion module")
    return kmer_ids, bonus_kmer_ids
