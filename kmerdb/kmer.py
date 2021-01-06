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
        # Each of the n-k+1 string slices become the k-mers
        for c in range(len(seqRecord.seq) - self.k + 1):
            s = seqRecord.seq[c:(c+self.k)]
            kmers.append(str(s))
            if self.strand_specific: # Reverse complement by default
                kmers.append(str(s.reverse_complement()))
        return {'id': seqRecord.id, 'kmers': list(filter(lambda x: x is not None, map(kmer_to_id, kmers)))}
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

    if not isinstance(s, str): # Typecheck the input k-mer
        raise TypeError("kmerdb.kmer.kmer_to_id expects a Biopython Seq object as its argument")
    elif s.find('N') != -1: # k-mer with 'N' do not have a binary encoding
        #logger.debug(TypeError("kdb.kmer.kmer_to_id expects the letters to contain only nucleotide symbols ATCG"))
        return None
    else: 
        idx = 0
        for c in bytes(s, "UTF-8"): # Use byteshifting for fast conversion to binary encoding
            idx = idx << 2
            idx = idx | letterToBinary[c]
        return idx



