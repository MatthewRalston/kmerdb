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
from itertools import repeat
import logging
logger = logging.getLogger(__file__)
from Bio import SeqIO, Seq
import Bio


letterToBinary={ # Unicode UTF-8 byte codes for ACGT
    65: 0,
    67: 1,
    71: 2,
    84: 3
}
binaryToLetter=['A', 'C', 'G', 'T']



class Kmers:
    """A wrapper class to pass variables through the multiprocessing pool
    
    :ivar k: The choice of k to shred with
    :ivar strand_specific: Include k-mers from forward strand only
    """
    def __init__(self, k, strand_specific=True):
        """

        :param k: The choice of k to shred with
        :type k: int
        :param strand_specific: Include k-mers from forward strand only
        :type strand_specific: bool


        """
        if type(k) is not int:
            raise TypeError("kmerdb.kmer.Kmers.__init__() expects an int as its first positional argument")
        elif type(strand_specific) is not bool:
            raise TypeError("kmerdb.kmer.Kmers.__init__() expects a bool as its second positional argument")
        self.k = k 
        self.strand_specific = strand_specific

    def shred(self, seqRecord):
        """

        :param seqRecord: 
        :type seqRecord: Bio.SeqRecord.SeqRecord
        :returns: 
        :rtype: 

        """
        if not isinstance(seqRecord, Bio.SeqRecord.SeqRecord):
            raise TypeError("kmerdb.kmer.Kmers expects a Bio.SeqRecord.SeqRecord object as its first positional argument")
        kmers = []
        starts = []
        reverses = []
        # Each of the n-k+1 string slices become the k-mers
        for c in range(len(seqRecord.seq) - self.k + 1):
            s = seqRecord.seq[c:(c+self.k)]
            kmers.append(str(s))
            reverses.append(False)
            starts.append(c)
            if not self.strand_specific: # Reverse complement by default
                kmers.append(str(s.reverse_complement()))
                reverses.append(True)
                starts.append(c)
        
        sys.stderr.write("            --- ~~~ --- ~~~  shredded ~~~ --- ~~~ ---\n")
        sys.stderr.write("a {0}bp long sequence was shredded into L-k+1 {1} total and {1} unique k-mers\n\n".format(len(seqRecord.seq), len(seqRecord.seq)-self.k+1, len(kmers)))
        output = {'id': seqRecord.id, 'kmers': list(filter(lambda x: x is not None, map(kmer_to_id, kmers))), "seqids": repeat(seqRecord.id, len(starts)), "starts": starts, 'reverses': reverses}
        #logger.debug(output['seqids'])
        return output
        #return list(filter(lambda x: x is not None, map(kmer_to_id, kmers)))


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
