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

letterToBinary={ # Unicode UTF-8 byte codes for ACGT
    65: 0,
    67: 1,
    71: 2,
    84: 3
}
binaryToLetter=['A', 'C', 'G', 'T']
standard_letters=set("ACTG")


"""
IUPAC support mappings for the k-mer counter (NOT USED IN THE ASSEMBLER)
"""
IUPAC_TO_DOUBLET = {
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"]
}
IUPAC_DOUBLET_CHARS=IUPAC_TO_DOUBLET.keys()
IUPAC_TO_TRIPLET = {
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"]
}
IUPAC_TRIPLET_CHARS=IUPAC_TO_TRIPLET.keys()
# *extremely* necessary variable
n = "N"
#############################
#
# Lambdas
#
#############################


get_doublet_chars = lambda x: IUPAC_TO_DOUBLET[x] # list
get_triplet_chars = lambda x: IUPAC_TO_TRIPLET[x]
replace_char = lambda seq, x, y: seq.replace(x, y) # str

#############################
#
# Class definition
#
#############################

class Kmers:
    """A wrapper class to pass variables through the multiprocessing pool
    
    :ivar k: The choice of k to shred with
    :ivar strand_specific: Include k-mers from forward strand only (TO BE DEPRECATED)
    :ivar verbose: print extra logging in the case of fasta files
    :ivar all_metadata: return extra metadata about sequence locations
    """
    def __init__(self, k, strand_specific:bool=True, verbose:bool=False, all_metadata:bool=False):
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
            raise TypeError("kmerdb.kmer.Kmers expects a bool as its second positional argument")
        elif type(verbose) is not bool:
            raise TypeError("kmerdb.kmer.Kmers expects a bool as its third positional argument")
        elif type(all_metadata) is not bool:
            raise TypeError("kmerdb.kmer.Kmers expects a bool as its fourth positional argument")
        self.k = k 
        self.strand_specific = strand_specific
        self.verbose = verbose
        self.all_metadata = all_metadata
        self.__permitted_chars = set("ACTGRYSWKMBDHVN")

    def validate_seqRecord_and_detect_IUPAC(self, seqRecord:Bio.SeqRecord.SeqRecord, is_fasta:bool=True, quiet_iupac_warning:bool=True):
        """
        Helper method for validating seqRecord and warnings for non-standard IUPAC residues.

        :param seqRecord: a BioPython SeqRecord object 
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
            raise TypeError("kmerdb.kmer.validate_seqRecord_and_detect_IUPAC expects keyword argument 'quiet_iupac_warning' to be a bool")
        elif not isinstance(seqRecord, Bio.SeqRecord.SeqRecord):
            raise TypeError("kmerdb.kmer.validate_seqRecord_and_detect_IUPAC expects a Bio.SeqRecord.SeqRecord object as its first positional argument")
        letters = set(seqRecord.seq) # This is ugly. Should really be explicitly cast to str
        seqlen = len(seqRecord.seq)  # `` Ditto. Should be explicitly cast to str

        
        if seqlen < self.k:
            logger.error("Offending sequence ID: {0}".format(seqRecord.id))
            raise ValueError("kmerdb.kmer.validate_seqRecord_and_detect_IUPAC expects that each input sequence is longer than k.")

        all_iupac_symbols = list(letters.intersection(self.__permitted_chars) - standard_letters)
        all_non_iupac_symbols = list(letters - self.__permitted_chars)
        
        if len(all_non_iupac_symbols) > 0:
            logger.warning("One or more unexpected characters in the {0} sequence".format("fasta" if self.is_fasta else "fastq"))
            logger.warning(list(letters - self.__permitted_chars))
            raise ValueError("Non-IUPAC symbols detected in the fasta file")
        elif len(all_iupac_symbols) > 0:
            if quiet_iupac_warning is False:
                logger.warning("Will completely refuse to include k-mers with 'N'")
                logger.warning("All counts for k-mers including N will be discarded")
                logger.warning("Other IUPAC symbols will be replaced with their respective pairs/triads")
                logger.warning("And a count will be given to each, rather than a half count")
            elif quiet_iupac_warning is True:
                logger.warning("Suppressing warning that non-standard IUPAC residues (including N) are detected.")
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
            s = seqRecord.seq[i:(i+self.k)] # Creates the k-mer as a slice of a seqRecord.seq
            # No non-standard IUPAC residues allowed for _shred for use in graph.py
            nonstandard_iupac_symbols = list(set(s) - standard_letters)
            non_iupac_symbols = list(set(s) - self.__permitted_chars)

            if len(non_iupac_symbols) > 1:
                logger.error("Non-IUPAC symbols:")
                logger.error(non_iupac_symbols)
                raise RuntimeError("Non-IUPAC symbol(s) detected...")
            elif len(nonstandard_iupac_symbols) == 0:
                #logger.debug("Perfect sequence content (ATCG) detected...")
                kmers.append(kmer_to_id(s))
            elif len(nonstandard_iupac_symbols) == 1 and iupac_symbols[0] == n:
                #logger.debug("Only N-abnormal sequence content detected (aside from ATCG)...")
                kmers.append(kmer_to_id(s))
            elif len(nonstandard_iupac_symbols) > 1:
                logger.error("Assembly cannot continue with this k-mer, and it will be discarded, possibly affecting assembly process")
                raise RuntimeError("Non-standard IUPAC symbols detected in .fa/.fq file...")
        return kmers
            
    def shred(self, seqRecord):
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
        letters, seqlen = self.validate_seqRecord_and_detect_IUPAC(seqRecord, quiet_iupac_warning=True)

            
        kmers = []
        starts = []
        reverses = []
        # Each of the n-k+1 string slices become the k-mers
        for i in range(seqlen - self.k + 1):
            s = seqRecord.seq[i:(i+self.k)]
            #logger.debug(letters - self.__permitted_chars)
            iupac_symbols = list(set(s).intersection(self.__permitted_chars) - standard_letters)
            non_iupac_symbols = list(set(s) - self.__permitted_chars)

            ######################################################################33
            #       I U P A C        m o d u l e
            ######################################################################33

            if len(set(str(s)) - self.__permitted_chars) == 0 and len(iupac_symbols) > 0 and len(non_iupac_symbols) == 0:
                seqs = None
                for c in list(set(s) - standard_letters):
                    # logger.debug("character: {0}".format(c))
                    # logger.debug("K-mer: {0}".format(s))
                    # logger.debug("IUPAC symbols in the whole sequence: {0}".format(iupac_symbols))
                    # logger.debug("Non-IUPAC symbols in the whole sequence: {0}".format(non_iupac_symbols))
                    if c in IUPAC_DOUBLET_CHARS:
                        #logger.info("Doublets being replaced")
                        if seqs is None:
                            seqs = list(map(lambda x: replace_char(str(s), c, x), get_doublet_chars(c)))
                        else:
                            seqs = list(chain.from_iterable(map(lambda s: list(map(lambda x: replace_char(s, c, x), get_doublet_chars(c))), seqs)))
                    elif c in IUPAC_TRIPLET_CHARS:
                        #logger.info("Triplets being replaced")
                        if seqs is None:
                            seqs = list(map(lambda x: replace_char(str(s), c, x), get_triplet_chars(c)))
                        else:
                            seqs = list(chain.from_iterable(map(lambda s: list(map(lambda x: replace_char(s, c, x), get_triplet_chars(c))), seqs)))
                    elif c in standard_letters:
                        raise RuntimeError("Standard DNA residue detected in IUPAC module")
                    elif c == n:
                        logger.warning("N content detected")
                        continue
                    else:
                        sys.stderr.write(str(seqRecord) + "\n")
                        sys.stderr.write("Full sequence above" + "\n")
                        sys.stderr.write("K-mer: {0}".format(s) + "\n")
                        sys.stderr.write(s + "\n")
                        raise RuntimeError("kmerdb.kmer.shred: Non-IUPAC character '{0}' made it into sequences generated from IUPAC doublet/triplet counting of the k-mer '{1}'".format(c, s))
                if seqs is None and "N" not in s:
                    sys.stderr.write("Permitted IUPAC symbols: {0}".format(iupac_symbols))
                    sys.stderr.write("Non-IUPAC symbols detected in sequence '{0}': {1}".format(s, non_iupac_symbols))
                    sys.stderr.write("The following k-mer did not have any IUPAC letters: {0}".format(s))
                    raise RuntimeError("kmerdb.kmer.shred: A sequence was rejected from non-IUPAC symbols")
                elif seqs is None and "N" in s:
                    continue
                one_word = "".join(seqs)
                if len(set(one_word) - self.__permitted_chars) > 0:
                    sys.stderr.write(one_word + "\n")
                    sys.stderr.write("Doublets/triplets produced from k-mer '{0}': \n{1}".format(s, seqs) + "\n")
                    raise ValueError("kmerdb.kmer.shred: at least one non-IUPAC symbol was found during doublets/triplet expansion of the k-mer '{0}'".format(s))
                else:
                    for x in seqs:
                        logger.debug(x)
                        kmers.append(kmer_to_id(x))
                        if self.all_metadata:
                            reverses.append(False)
                            starts.append(i)
                        if not self.strand_specific: # Reverse complement by default
                            kmers.append(kmer_to_id(Seq.Seq(x).reverse_complement()))
                            if self.all_metadata:
                                reverses.append(True)
                                starts.append(i)

                    # OKAY you are supposed to remember that
                    # the number of k-mer ids will no-longer match the appended metadata,
                    # So even though the lambda is correct, you have to do separate metadata for each id
            elif len(non_iupac_symbols) > 0:
                sys.stderr.write("TEST CONDITIONS:\n")
                sys.stderr.write("Non-permitted characters should be 0: {0}".format(len(set(str(s)) - self.__permitted_chars)) + "\n")
                sys.stderr.write("IUPAC symbols should be > 0: {0}".format(len(iupac_symbols)) + "\n")
                sys.stderr.write("Non-IUPAC characters should be 0: {0}".format(len(non_iupac_symbols)) + "\n")
                sys.stderr.write("Non-IUPAC symbols: {0}".format(non_iupac_symbols) + "\n")
                sys.stderr.write("IUPAC symbols: {0}".format(iupac_symbols) + "\n")
                sys.stderr.write("K-mer: '{0}'".format(str(s)) + "\n")
                raise ValueError("kmerdb.kmer.shred: Non-IUPAC symbol(s) detected")
            else:
                try:
                    kmers.append(kmer_to_id(s))
                    if self.all_metadata:
                        reverses.append(False)
                        starts.append(i)
                    if not self.strand_specific: # Reverse complement by default
                        kmers.append(kmer_to_id(s.reverse_complement()))
                        if self.all_metadata:
                            reverses.append(True)
                            starts.append(i)
                except KeyError as e:
                    sys.stderr.write("One or more non-IUPAC letters found their way into a part of the code they're not supposed to go\n")
                    sys.stderr.write("We officially support IUPAC but the statement challenging the sequence content failed, causing a genuine runtime error\n")
                    sys.stderr.write("This caused the following KeyError\n")
                    sys.stderr.write(e.__str__() + "\n")
                    sys.stderr.write("Letters in the sequence: {0}".format(letters) + "\n")
                    sys.stderr.write("Permitted letters: {0}".format(self.__permitted_chars) + "\n")
                    raise RuntimeError("IUPAC standard extra base pairs (R, B, etc.) or non-IUPAC characters detected in the sequence")
            del s
        if self.verbose:
            sys.stderr.write("            --- ~~~ --- ~~~  shredded ~~~ --- ~~~ ---\n")
            sys.stderr.write("a {0}bp long sequence was shredded into L-k+1 {1} total and {2} unique k-mers\n\n{3} were discarded due to sequence content\n\n".format(len(seqRecord.seq), len(str(seqRecord.seq))-self.k+1, len(list(set(kmers))), len([x for x in kmers if x is None])) + "\n")
        return {'id': seqRecord.id, 'kmers': kmers, "seqids": repeat(seqRecord.id, len(starts)), "starts": starts, 'reverses': reverses}



#############################
#
# Functions
#
#############################

    

def kmer_to_id(s):
    """Convert a fixed length k-mer string to the binary encoding parameterized upon that same k

    Note that the conversion of a k-mer string to an id integer
    is consistent regardless of k, 
    because the k is implicit in the k-mer string's size.

    Therefore, this method does not need to be wrapped in the k-mer class



    :param s: The input k-mer as string
    :type s: str
    :raise TypeError: str argument required, non-str type detected
    :raise KeyError: Non-standard (ATCG) character detected
    :returns: The kPal-inspired binary encoding (thanks!)
    :rtype: int

    """

    if not isinstance(s, str) and type(s) is not str and not isinstance(s, Bio.Seq.Seq) and not isinstance(s, Bio.SeqRecord.SeqRecord): # Typecheck the input k-mer
        logger.error("Seq: {0}, type: {1}".format(s, type(s)))
        raise TypeError("kmerdb.kmer.kmer_to_id expects a str as its argument")
    elif s.find('N') != -1: # k-mer with 'N' do not have a binary encoding
        #logger.debug(TypeError("kdb.kmer.kmer_to_id expects the letters to contain only nucleotide symbols ATCG"))
        return None
    else: 
        idx = 0
        if isinstance(s, Bio.Seq.Seq) or isinstance(s, Bio.SeqRecord.SeqRecord):
            s = str(s)
        for c in bytes(s, "UTF-8"): # Use byteshifting for fast conversion to binary encoding
            idx = idx << 2
            try:
                idx = idx | letterToBinary[c]
            except KeyError as e:
                sys.stderr.write("Entire sequence: {0}".format(s) + "\n")
                sys.stderr.write("Problematic character: {0}".format(c) + "\n")
                raise e
        return idx


def id_to_kmer(id, k):
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
    else:
        kmer = ""
        for i in range(k):
            kmer += binaryToLetter[id & 0x03]
            id = id >> 2
        kmer = list(kmer)
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
        letters1 = copy.deepcopy(binaryToLetter)
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

#         logger.debug(""" flower garden - joan G. Stark

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

        sys.stderr.write("kmerdb.kmer.neighbors creating neighbors for k-mer {0} : '{1}'...".format(kmer_id, kmer) + "\n")
        sys.stderr.write(" ========================\n\n")
        sys.stderr.write(", ".join(list(map(str, new_type1_ids))) + "\n")
        sys.stderr.write(", ".join(list(map(str, new_type2_ids))) + "\n\n")
        sys.stderr.write(" ------------------------\n\n")
        
        if quiet is not True:
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
