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

letterToBinary={ # Unicode UTF-8 byte codes for ACGT
    65: 0,
    67: 1,
    71: 2,
    84: 3
}
binaryToLetter=['A', 'C', 'G', 'T']
standard_letters=set("ACTG")


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
    :ivar strand_specific: Include k-mers from forward strand only
    """
    def __init__(self, k, strand_specific=True, fasta=False, all_metadata=False):
        """

        :param k: The choice of k to shred with
        :type k: int
        :param strand_specific: Include k-mers from forward strand only
        :type strand_specific: bool
        :param fasta: Whether or not to print extra logging in the case of fasta and not in the case of fastq
        :type fasta: bool
        :param all_metadata: Whether or not to pass back extra metadata
        :type all_metadata: bool


        """
        if type(k) is not int:
            raise TypeError("kmerdb.kmer.Kmers.__init__() expects an int as its first positional argument")
        elif type(strand_specific) is not bool:
            raise TypeError("kmerdb.kmer.Kmers.__init__() expects a bool as its second positional argument")
        elif type(fasta) is not bool:
            raise TypeError("kmerdb.kmer.Kmers.__init__() expects a bool as its third positional argument")
        elif type(all_metadata) is not bool:
            raise TypeError("kmerdb.kmer.Kmers.__init__() expects a bool as its fourth positional argument")
        self.k = k 
        self.strand_specific = strand_specific
        self.fasta = fasta
        self.all_metadata = all_metadata
        self.__permitted_chars = set("ACTGRYSWKMBDHVN")

    def shred(self, seqRecord):
        """

        :param seqRecord: 
        :type seqRecord: Bio.SeqRecord.SeqRecord
        :returns: 
        :rtype: 

        """
        
        
        if not isinstance(seqRecord, Bio.SeqRecord.SeqRecord):
            raise TypeError("kmerdb.kmer.Kmers expects a Bio.SeqRecord.SeqRecord object as its first positional argument")
        letters = set(seqRecord.seq)
        seqlen = len(seqRecord.seq)
        if seqlen < self.k:
            logger.error("Offending sequence ID: {0}".format(seqRecord.id))
            raise ValueError("kmerdb expects that each input sequence is longer than k.")
        all_iupac_symbols = list(letters.intersection(self.__permitted_chars) - standard_letters)
        all_non_iupac_symbols = list(letters - self.__permitted_chars)
        if len(all_non_iupac_symbols) > 0:
            logger.warning("One or more unexpected characters in the {0} sequence".format("fasta" if self.fasta else "fastq"))
            logger.warning(list(letters - self.__permitted_chars))
            raise ValueError("Non-IUPAC symbols detected in the fasta file")
        elif len(all_iupac_symbols) > 0:
            logger.warning("Will completely refuse to include k-mers with 'N'")
            logger.warning("All counts for k-mers including N will be discarded")
            logger.warning("Other IUPAC symbols will be replaced with their respective pairs/triads")
            logger.warning("And a count will be given to each, rather than a half count")
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
                    logger.debug("character: {0}".format(c))
                    logger.debug("K-mer: {0}".format(s))
                    logger.debug("IUPAC symbols in the whole sequence: {0}".format(iupac_symbols))
                    logger.debug("Non-IUPAC symbols in the whole sequence: {0}".format(non_iupac_symbols))
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
                        logger.error(str(seqRecord))
                        logger.error("Full sequence above")
                        logger.error("K-mer: {0}".format(s))
                        logger.error(s)
                        raise RuntimeError("Non-IUPAC character '{0}' made it into the IUPAC extension of kmerdb.kmer.Kmer.shred".format(c))
                if seqs is None and "N" not in s:
                    logger.debug(iupac_symbols)
                    logger.debug(non_iupac_symbols)
                    logger.error(s)
                    logger.error("The k-mer above did not have any IUPAC letters")
                    raise RuntimeError("No sequences generated in IUPAC module of kmerdb.kmer.Kmer.shred")
                elif seqs is None and "N" in s:
                    continue
                one_word = "".join(seqs)
                if len(set(one_word) - self.__permitted_chars) > 0:
                    logger.error(one_word)
                    logger.error(seqs)
                    raise ValueError("Still IUPAC characters in the strings")
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
                logger.debug("TEST CONDITIONS")
                logger.debug("Non-permitted characters should be 0: {0}".format(len(set(str(s)) - self.__permitted_chars)))
                logger.debug("IUPAC symbols should be >= 0: {0}".format(len(iupac_symbols)))
                logger.debug("Non-IUPAC characters should be 0: {0}".format(len(non_iupac_symbols)))
                logger.debug("Non-IUPAC symbols: {0}".format(non_iupac_symbols))
                logger.debug("IUPAC symbols: {0}".format(iupac_symbols))
                logger.debug("K-mer: {0}".format(str(s)))
                raise ValueError("Non-IUPAC symbol detected")
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
                    logger.debug("One or more non-IUPAC letters found their way into a part of the code they're not supposed to go")
                    logger.debug("We officially support IUPAC but the statement challenging the sequence content failed, causing a genuine runtime error")
                    logger.debug("This caused the following KeyError")
                    logger.debug(e)
                    logger.error("Letters in the sequence: {0}".format(letters))
                    logger.errror("Permitted letters: {0}".format(self.__permitted_chars))
                    logger.error("Unexpected behavior, rerun with -vv (DEBUG verbosity) to see more information")
                    raise RuntimeError("IUPAC standard extra base pairs (R, B, etc.) or non-IUPAC characters detected in the sequence")
            del s
        if self.fasta:
            sys.stderr.write("            --- ~~~ --- ~~~  shredded ~~~ --- ~~~ ---\n")
            sys.stderr.write("a {0}bp long sequence was shredded into L-k+1 {1} total and {2} unique k-mers\n\n{3} were discarded due to sequence content\n\n".format(len(seqRecord.seq), len(str(seqRecord.seq))-self.k+1, len(list(set(kmers))), len([x for x in kmers if x is None])))
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

    :param s: The input k-mer
    :type s: str
    :returns: The kPal-inspired binary encoding
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
                logger.error("Entire sequence: {0}".format(s))
                logger.error("Problematic character: {0}".format(c))
                raise e
        return idx


def id_to_kmer(id, k):
    if type(id) is not int:
        raise TypeError("kmerdb.id_to_kmer expects an int as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.id_to_kmer expects an int as its second positional argument")
    else:
        kmer = ""
        for i in range(k):
            kmer += binaryToLetter[id & 0x03]
            id = id >> 2
        kmer = list(kmer)
        kmer.reverse()
        return ''.join(kmer)


def neighbors(s, k):
    if not isinstance(s, str):
        raise TypeError("kmerdb.kmer.neighbors expects a Biopython Seq object as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.kmer.neighbors expects an int as its second positional argument")
    elif len(s) != k:
        raise TypeError("kmerdb.kmer.neighbors cannot calculate the {0}-mer neighbors of a {1}-mer".format(k, len(s)))
    else:
        import copy
        letters1 = copy.deepcopy(binaryToLetter)
        letters2 = copy.deepcopy(letters1)
    
        firstChar = s[0]
        lastChar  = s[-1]
        suffix    = s[1:]
        prefix    = s[:-1]
        rightNeighborIds = dict((c, kmer_to_id(suffix+c)) for c in letters1)
        leftNeighborIds = dict((c, kmer_to_id(c+prefix)) for c in letters2)
        return {"suffixes": rightNeighborIds, "prefixes": leftNeighborIds}
