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
from itertools import chain, repeat
import logging
logger = logging.getLogger(__file__)
from Bio import SeqIO, Seq
import Bio



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

binaryToLetterAA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y"] + ["!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ","] # original 20 elements + 12 dummy characters

standard_lettersAA = set(binaryToLetterAA)

IUPAC_AA_DOUBLETS = {
    "B": ["R", "N"], # Aspartate (R) or Aspartamine (N)
    "Z": ["E", "Q"], # Glutamate (E) or Glutamine (Q)
    "J": ["L", "I"]  # Leucine (L) or Isoleucine (I)
}
IUPAC_AA_DOUBLET_CHARS=set(IUPAC_AA_DOUBLETS)

get_doublet_chars_AA = lambda x: IUPAC_AA_DOUBLETS[x]
get_triplet_chars_AA = lambda x: IUPAC_AA_DOUBLETS[x]
permitted_AA_characters = standard_lettersAA
permitted_AA_characters_extended = standard_lettersAA.union(IUPAC_AA_DOUBLET_CHARS)



#############################
#
# Key quick lambdas
#
#############################
replace_char = lambda seq, x, y: seq.replace(x, y) # str

get_doublet_chars_AA = lambda x: IUPAC_AA_DOUBLETS[x]
get_triplet_chars_AA = lambda x: IUPAC_AA_DOUBLETS[x]
get_doublet_chars_NA = lambda x: IUPAC_NA_DOUBLETS[x] # list
get_triplet_chars_NA = lambda x: IUPAC_NA_TRIPLETS[x]

def is_sequence_na(s:str):
    if type(s) is str:
        letters = set(str(s))
    elif isinstance(s, Bio.SeqRecord.SeqRecord):
        letters = set(str(s.seq))
    else:
        raise TypeError("kmer.is_sequence_na expects a str/Bio.SeqRecord as its only positional argument")

    if letters.difference(permitted_NA_characters_with_N) == set():
        return True
    elif letters.difference(permitted_AA_characters_extended) == set():
        return False
    else:
        raise ValueError("Unable to determine IUPAC nature of the following sequence\n{0}".format(s))


def is_sequence_aa(s:str):
    if type(s) is str:
        letters = set(str(s))
    elif isinstance(s, Bio.SeqRecord.SeqRecord):
        letters = set(str(s.seq))
    else:
        raise TypeError("kmer.is_sequence_na expects a str/Bio.SeqRecord as its only positional argument")

    if letters.difference(permitted_AA_characters_extended) == set():
        return True
    elif letters.difference(permitted_NA_characters_with_N) == set():
        return False
    else:
        raise ValueError("Unable to determine IUPAC nature of the following sequence\n{0}".format(s))
    

#############################
#
# Class definition
#
#############################

class Kmers:
    """A wrapper class to pass variables through the multiprocessing pool. Methods in this class are (c) Matthew Ralston as above/throughout. That said, the k-mer bit-shifting trick
    
    :ivar k: The choice of k to shred with
    :ivar strand_specific: Include k-mers from forward strand only (TO BE DEPRECATED)
    :ivar verbose: print extra logging in the case of fasta files
    :ivar all_metadata: return extra metadata about sequence locations
    """
    def __init__(self, k, strand_specific:bool=True, amino_acid:bool=False, verbose:bool=False, all_metadata:bool=False):
        """

        :param k: The choice of k to shred with
        :type k: int
        :param strand_specific: Include k-mers from forward strand only (TO BE DEPRECATED)
        :type strand_specific: bool
        :param verbose: Whether or not to print extra logging in the case of fasta and not in the case of fastq
        :type verbose: bool
        :param all_metadata: Whether or not to pass back extra metadata (TO BE DEPRECATED)
        :type all_metadata: bool


        """
        if type(k) is not int:
            raise TypeError("kmerdb.kmer.Kmers expects an int as its first positional argument")
        elif type(strand_specific) is not bool:
            raise TypeError("kmerdb.kmer.Kmers expects a bool as its first keyword argument 'strand_specific'")
        elif type(amino_acid) is not bool:
            raise TypeError("kmerdb.kmer.Kmers expects a bool as its second keyword argument 'amino_acid'")
        elif type(verbose) is not bool:
            raise TypeError("kmerdb.kmer.Kmers expects a bool as its third keyword argument 'verbose'")
        elif type(all_metadata) is not bool:
            raise TypeError("kmerdb.kmer.Kmers expects a bool as its fourth keyword argument 'all_metadata' [DEPRECATED]")
        self.k = k 
        self.strand_specific = strand_specific
        self.amino_acid = amino_acid
        self.verbose = verbose
        self.all_metadata = all_metadata

        #self.__all_permitted_NA_characters = __permitted_NA_characters_with_N
        if self.amino_acid is False:
            self.n = "N"
            self.IUPAC_DOUBLET = IUPAC_NA_DOUBLETS
            self.IUPAC_TRIPLET = IUPAC_NA_TRIPLETS
            self.get_doublet_chars = get_doublet_chars_NA
            self.get_triplet_chars = get_triplet_chars_NA
            
            self._permitted_characters = permitted_NA_characters_with_N
            self._permitted_characters_with_N = permitted_NA_characters_with_N
            self._permitted_characters_extended = permitted_NA_characters
            self.standard_letters = standard_lettersNA
        else:
            self.n = "X"
            self.IUPAC_DOUBLE = IUPAC_AA_DOUBLETS
            self.IUPAC_TRIPLET = {}
            self.get_doublet_chars = get_doublet_chars_AA
            self.get_triplet_chars = get_doublet_chars_AA
            
            self._permitted_characters = permitted_AA_characters
            self._permitted_characters_with_N = permitted_AA_characters
            self._permitted_characters_extended = permitted_AA_characters_extended
            self.standard_letters = standard_lettersAA

    def validate_seqRecord_and_detect_IUPAC(self, seqRecord:Bio.SeqRecord.SeqRecord, is_fasta:bool=True, is_aa:bool=False, quiet_iupac_warning:bool=True):
        """
        Helper method for validating seqRecord and warnings for non-standard IUPAC residues for Nucleic Acid sequences

        :param seqRecord: a BioPython SeqRecord object containing a nucleic acid string
        :type seqRecord: Bio.SeqRecord.SeqRecord 
        :param is_fasta: are the inputs ALL .fasta?
        :type is_fasta: bool
        :param quiet_iupac_warning: verbosity parameter
        :type quiet_iupac_warning: bool
        :raise ValueError: if *non* IUPAC symbols are detected in the seqRecord object.
        :returns: a set of the letters detected, and length of the validated Bio.SeqRecord sequence
        :rtype: tuple


        """
        if type(quiet_iupac_warning) is not bool:
            raise TypeError("kmerdb.kmer.Kmers.validate_seqRecord_and_detect_IUPAC expects keyword argument 'quiet_iupac_warning' to be a bool")
        elif type(is_aa) is not bool:
            raise TypeError("kmerdb.kmer.Kmers.validate_seqRecord_and_detect_IUPAC expects keyword argument 'is_aa' to be a bool")
        elif not isinstance(seqRecord, Bio.SeqRecord.SeqRecord):
            raise TypeError("kmerdb.kmer.Kmers.validate_seqRecord_and_detect_IUPAC expects a Bio.SeqRecord.SeqRecord object as its first positional argument")
        letters = set(str(seqRecord.seq)) 
        seqlen = len(str(seqRecord.seq))  

        
        if seqlen < self.k:
            logger.error("Offending sequence ID: {0}".format(seqRecord.id))
            raise ValueError("kmerdb.kmer.validate_seqRecord_and_detect_IUPAC_NA expects that each input sequence is longer than k.")
        
        extended_iupac_symbols = list(letters.intersection(self._permitted_characters_extended) - self.standard_letters)
        #has_N = self._permitted_NA_characters_with_N 
        all_non_iupac_symbols = list(letters - self._permitted_characters_extended)
        
        if len(all_non_iupac_symbols) > 0:
            logger.warning("One or more unexpected characters in the {0} sequence".format("fasta" if self.is_fasta else "fastq"))
            logger.warning(all_non_iupac_symbols)
            raise ValueError("Non-IUPAC symbols detected in the sequence '{0}'".format(seqRecord.id))
        elif len(extended_iupac_symbols) > 0: # FIXME: 
            if quiet_iupac_warning is False:
                logger.warning("Will completely refuse to include k-mers with 'N' (or 'X' for AA sequences)")
                logger.warning("All counts for k-mers including N will be discarded")
                logger.warning("Other IUPAC symbols will be replaced with their respective pairs/triads")
                logger.warning("And a count will be given to each, rather than a half count")
            elif quiet_iupac_warning is True:
                logger.warning("Suppressing warning that non-standard IUPAC residues (including N/X) are detected.")
        return (letters, seqlen)


    def _shred_for_graph(self, seqRecord:Bio.SeqRecord.SeqRecord):
        """
        Introduced in 0.7.7, required for valid assembly. Simply omits the rejection of IUPAC non-standard residues.

        :param seqRecord: a BioPython sequence object to shred into k-mers
        :type seqRecord: Bio.SeqRecord.SeqRecord
        :raise RuntimeError: Non-IUPAC and non-standard IUPAC (other than ATGC+N) residues detected
        :return: k-mer ids
        :rtype: list
        """
        kmers = []
        
        letters, seqlen = self.validate_seqRecord_and_detect_IUPAC(seqRecord, quiet_iupac_warning=True)


        # Iterate over *each* k-mer in the sequence by index
        for i in range(seqlen - self.k + 1):
            kmer_seq = seqRecord.seq[i:(i+self.k)] # Creates the k-mer as a slice of a seqRecord.seq
            # No non-standard IUPAC residues allowed for _shred for use in graph.py
            nonstandard_iupac_symbols = list(set(kmer_seq) - self.standard_letters)
            non_iupac_symbols = list(set(kmer_seq) - self._permitted_characters_extended)

            if len(non_iupac_symbols) > 1:
                logger.error("Non-IUPAC symbols:")
                logger.error(non_iupac_symbols)
                raise RuntimeError("Non-IUPAC symbol(s) detected...")
            elif len(nonstandard_iupac_symbols) == 0:
                #logger.debug("Perfect sequence content (ATCG) detected...")
                kmers.append(kmer_to_id(s))
            elif len(nonstandard_iupac_symbols) == 1 and iupac_symbols[0] == self.n:
                #logger.debug("Only N-abnormal sequence content detected (aside from ATCG)...")
                kmers.append(kmer_to_id(s))
            elif len(nonstandard_iupac_symbols) > 1:
                logger.error("Assembly cannot continue with this k-mer, and it will be discarded, possibly affecting assembly process")
                raise RuntimeError("Non-standard IUPAC symbols detected in .fa/.fq file...")
        return kmers
            
    def shred(self, seqRecord, is_aa:bool=False):
        """
        Take a seqRecord fasta/fastq object and slice according to the IUPAC charset.
        Doublets become replace with two counts, etc.

        :param seqRecord: 
        :type seqRecord: Bio.SeqRecord.SeqRecord
        :raise RuntimeError: No IUPAC characters detected or Non-IUPAC characters detected
        :raise ValueError: Non-IUPAC characters detected
        :returns: a dictionary of information about the Bio.SeqRecord shredded and the k-mers produced (including the IUPAC doublet/triplet expansion
        :rtype: dict

        """

        letters, seqlen = self.validate_seqRecord_and_detect_IUPAC(seqRecord, is_aa=is_aa, quiet_iupac_warning=True)

            
        kmers = []

        kmers = [seqRecord.seq[i:(i+self.k)] for i in range(seqlen - self.k + 1)]

        for s in kmers:
            iupac_symbols = list(set(list(s)).intersection(self._permitted_characters_extended))
            non_iupac_symbols = list(set(list(s)).difference(self._permitted_characters_extended))

            ######################################################################33
            #       I U P A C        m o d u l e
            ######################################################################33

            
            seqs = []

            if self.verbose is not False:
                sys.stderr.write("running IUPAC module. generating sequences. L252\n") # 10/1/24
                sys.stderr.write("system message: k-mer sequence -        [| {0} |]\n".format(s))
                sys.stderr.write("Iterating over k-mer for doublet/triplet iupac expansion\n")
            nonstandard_chars = (set(list(s)) - set(self.n)) - self.standard_letters
            for c in list(nonstandard_chars):

                sys.stderr.write("Non-permitted character found: {0}".format(c))
                raise RuntimeError("can we get in here actually?")
                sys.stderr.write("character: {0}\n".format(c))
                sys.stderr.write("K-mer: {0}\n".format(s))
                sys.stderr.write("IUPAC symbols in the whole sequence: {0}\n".format(iupac_symbols))
                sys.stderr.write("Non-IUPAC symbols in the whole sequence: {0}\n".format(non_iupac_symbols))
                if len(iupac_symbols) > 0:
                    sys.stderr.write("expanding iupac doublets/triplets...\n")
                #print("seq0: {0}".format(seqs))
                if c in self.standard_letters:
                    raise ValueError("Cannot have standard letters to this point")


                # print(c)
                # print("doublet expansions:")
                if c in self.IUPAC_DOUBLET.keys():
                    for x in self.get_doublet_chars(c):
                        doublet1, doublet2 = replace_char(s, c, x)
                        seqs.append(doublet1)
                        seqs.append(doublet2)

                if c in self.IUPAC_TRIPLET.keys():
                    for x in self.get_triplet_chars(c):
                        trip1, trip2, trip3 = replace_char(s, c, x)
                        seqs.append(trip1)
                        seqs.append(trip2)
                        seqs.append(trip3)
            
            if seqs is None and len(non_iupac_symbols) > 0:
                raise ValueError("Non-IUPAC symbols detected in {0}".format(s))
            elif seqs is None and len(iupac_symbols) == 4 and len(set(iupac_symbols) - self.standard_letters) == 0:
                e = ValueError("standard iupac symbols likely. ATCG or a subset for this k-mer. Likely not generated from sequencing data")
                #raise e
                continue # 9/30/24 whooopsshuaazzzaaah

            elif seqs is None and len(iupac_symbols) > 0 and len(iupac_symbols) < 4:
                e = ValueError("just... stnadard iupac symbols likely. ATCG or a subset for this k-mer. Likely not generated from sequencing data") # 10/3/24
                #raise e
                continue
            # elif seqs is None and len(iupac_symbols) >0 and len(iupac_symbols) > 4:
            #     e = ValueError("extra IUPAC characters detected in {0}".format(s))
            #     #raise e

            #raise RuntimeError("Still debug. bona suarte")
            """
            You shall not pass! ... gas!
            Gandalyf the greatest
            """
            # if seqs is None:
            #     sys.stderr.write("K-mer Sequence: {0}\n".format(s))
            #     sys.stderr.write("IUPAC symbols found in sequence: {0}\n".format(iupac_symbols))
            #     sys.stderr.write("Non-IUPAC symbols found in sequence: {0}\n".format(non_iupac_symbols))
            #     raise RuntimeError("seqs array was not populated. Needs debugging...")

            # if seqs is None:
            #     raise RuntimeError("Cannot pass here without have seqs generated from shredding process")
            # else:
            #     print("Squeakuences")
            #     print("Seeeeqquueeencees:")
            #     print(seqs)


                
            for x in seqs:
                logger.debug(x)
                kmers.append(kmer_to_id(x))
                if not self.strand_specific: # Reverse complement by default
                    kmers.append(kmer_to_id(Seq.Seq(x).reverse_complement()))

                    # OKAY you are supposed to remember that
                    # the number of k-mer ids will no-longer match the appended metadata,
                    # So even though the lambda is correct, you have to do separate metadata for each id
            # elif len(non_iupac_symbols) > 0:
            #     sys.stderr.write("TEST CONDITIONS:\n")
            #     sys.stderr.write("Non-permitted characters should be 0: {0}".format(len(set(str(s)) - self._permitted_NA_characters)) + "\n")
            #     sys.stderr.write("IUPAC symbols should be > 0: {0}".format(len(iupac_symbols)) + "\n")
            #     sys.stderr.write("Non-IUPAC characters should be 0: {0}".format(len(non_iupac_symbols)) + "\n")
            #     sys.stderr.write("Non-IUPAC symbols: {0}".format(non_iupac_symbols) + "\n")
            #     sys.stderr.write("IUPAC symbols: {0}".format(iupac_symbols) + "\n")
            #     sys.stderr.write("K-mer: '{0}'".format(str(s)) + "\n")
            #     raise ValueError("uuuuuuuuuupppppppsssssssss")
            # else:
            #     try:
            #         kmers.append(kmer_to_id(s))
            #         if not self.strand_specific: # Reverse complement by default
            #             kmers.append(kmer_to_id(s.reverse_complement()))
            #     except KeyError as e:
            #         sys.stderr.write("One or more non-IUPAC letters found their way into a part of the code they're not supposed to go\n")
            #         sys.stderr.write("We officially support IUPAC but the statement challenging the sequence content failed, causing a genuine runtime error\n")
            #         sys.stderr.write("This caused the following KeyError\n")
            #         sys.stderr.write(e.__str__() + "\n")
            #         sys.stderr.write("Letters in the sequence: {0}".format(letters) + "\n")
            #         sys.stderr.write("Permitted letters: {0}".format(self._permitted_NA_characters) + "\n")
            #         raise RuntimeError("IUPAC standard extra base pairs (R, B, etc.) or non-IUPAC characters detected in the sequence")
            del s
        if self.verbose:
            sys.stderr.write("            --- ~~~ --- ~~~  shredded ~~~ --- ~~~ ---\n")
            sys.stderr.write("a {0}bp long sequence was shredded into L-k+1 {1} total and {2} unique k-mers\n\n{3} were discarded due to sequence content\n\n".format(len(seqRecord.seq), len(str(seqRecord.seq))-self.k+1, len(list(set(kmers))), len([x for x in kmers if x is None])) + "\n")




        # def kolmogorov_complexity(self, kmers:list):
        #     """

        #     """
        #     kmer_strings = list(map(lambda s: id_to_kmer(s), kmers))
        #     for i, s in enumerate(kmer_strings):


        #         letters = set(map(lambda c: s.count(c), standard_letters))

        #         match len(letters):
        #             case 0:
        #                 raise ValueError("Error during kmerdb.kmer.Kmer.shred.kolmogorov_complexity. Cannot calculate complexity of sequence with no letters")
        #             case 1:
        #                 return 1
        #             case _:
        #                 return _kolmogorov_complexity(s)
                              
        return {'id': seqRecord.id, 'kmers': kmers, 'read_length': seqlen}
            
        #return {'id': seqRecord.id, 'kmers': kmers, "seqids": repeat(seqRecord.id, len(starts)), "starts": starts, 'reverses': reverses}



#############################
#
# Functions
#
#############################

    

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
    :raise TypeError: str argument required, non-str type detected
    :raise KeyError: Non-standard (ATCG) character detected
    :returns: The kPal-inspired binary encoding (thanks!)
    :rtype: int

    """
        # and not isinstance(s, Bio.Seq.Seq)  # Shortened the checker as this is a common function.
    if type(s) is not str and not isinstance(s, Bio.SeqRecord.SeqRecord) and not isinstance(s, Bio.Seq.Seq): # Typecheck the input k-mer
        logger.error("Seq: {0}, type: {1}".format(s, type(s)))
        raise TypeError("kmerdb.kmer.kmer_to_id expects a str as its positional argument")
    """
    Determine if sequence is NA or AA
    Offloaded to a bool
    
    Omitting the is_aa bool typecheck as this is a frequently used function
    """
    idx = 0
    if is_aa is True and s.find("X") != -1:
        idx = None
    elif is_aa is False and s.find('N') != -1: # k-mer with 'N' do not have a binary encoding
        #logger.debug(TypeError("kdb.kmer.kmer_to_id expects the letters to contain only nucleotide symbols ATCG"))
        idx = None

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

    if isinstance(s, Bio.Seq.Seq) or isinstance(s, Bio.SeqRecord.SeqRecord):
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
    :raise TypeError: int argument required, non-int type detected
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
    elif type(is_aa) is not bool:
        raise TypeError("kmerdb.id_to_kmer expects the keyword argument 'is_aa' to be a bool")
    
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

        if kmer_len != 0:
            raise ValueError("kmer.id_to_kmer returned an empty k-mer. The kmer-id was '{0}', the choice of k is {1}".format(id, k))
        elif kmer_len != k:
            raise ValueError("kmer.id_to_kmer encountered an inconsistency error. Expected lenght of k-mer sequence was {0}".format(k, kmer_len))


        #just_reversed = kmer.reverse()
        kmer.reverse()
        return ''.join(kmer)


def neighbors(kmer:str, kmer_id:int,  k:int, quiet:bool=True):
    """

    3/11/24 revived. given a k-mer of length k, give its neighbors.

    This is so ugly.

    But rock on :
    

    :param kmer: The sequence as a string that will be sliced for k-mers
    :type kmer: str
    :param kmer_id: The k-mer id for neighbor generation
    :type kmer_id: int
    :param k: The int k 
    :type k: int
    :raise TypeError: Requires either a Bio.SeqRecord or str
    :raise TypeError: Requires k to be an int
    :raise ValueError: k must match the length of the input k-mer
    :returns: A list of 'neighboring' or 'adjacent' k-mers from derived from the input k-mer
    :rtype: list
    """
    if not isinstance(kmer, str):
        raise TypeError("kmerdb.kmer.neighbors expects a Biopython Seq object as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.kmer.neighbors expects an int as its second positional argument")
    elif len(kmer) != k:
        raise ValueError("kmerdb.kmer.neighbors cannot calculate the {0}-mer neighbors of a {1}-mer".format(k, len(s)))
    else:
        import copy
        letters1 = copy.deepcopy(binaryToLetterNA)
        letters2 = copy.deepcopy(letters1)
    
        firstCharRemoved = kmer[1:]
        lastCharRemoved  = kmer[:-1]

        new_type1 = list(map(lambda c: firstCharRemoved + c, letters1)) # Form the up neighbors

        new_type2 = list(map(lambda c: c + lastCharRemoved, letters2))  # Form the down neighbors
        
        """
        # TYPE 1: [first char removed ... + ... c  : ['A', "c", "g", "T"]]

        # TYPE 2: [["A", "C", "G", "T"] : c + last char removed  ]
        """
        new_type1_ids = list(map(kmer_to_id, new_type1))
        new_type2_ids = list(map(kmer_to_id, new_type2))

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
            sys.stderr.write("kmerdb.kmer.neighbors creating neighbors for k-mer {0} : '{1}'...".format(kmer_id, kmer) + "\n")
            sys.stderr.write(" ========================\n\n")
            sys.stderr.write(", ".join(list(map(str, new_type1_ids))) + "\n")
            sys.stderr.write(", ".join(list(map(str, new_type2_ids))) + "\n\n")
            sys.stderr.write(" ------------------------\n\n")

            sys.stderr.write("""
        k-id : {0}
        kmer : \"    {1}        \"

        'neighbors'

        {2}
        {3}

        'ids':

        {4}
        {5}
\n""".format(kmer_id, kmer, new_type1, new_type2, new_type1_ids, new_type2_ids))
        #return {"appended_first_char_all_ommitted": new_type1_ids, "prepended_last_char_all_omitted": new_type2_ids}


        
        return new_type1_ids + new_type2_ids
