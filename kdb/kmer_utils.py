import logging
logger = logging.getLogger(__file__)
import re


class Utils:
    
    def __init__(self, k: int):
        """
        A utility class for converting sequences to binary representation.
        A class was needed because the declaration of k for decoding purposes is dynamic.
        The class also contains a get_neighbors function which was homeless.
        """
        if type(k) is not int: # Typecheck the value of k
            raise TypeError("kdb.kmer_utils.Utils expects an int as its first positional argument")
        
        self.k = k
        self.letterToBinary={ # Unicode UTF-8 byte codes for ACGT
            65: 0,
            67: 1,
            71: 2,
            84: 3
        }
        self.binaryToLetter=['A', 'C', 'G', 'T']

    def sequence_to_binary(self, s):
        """ Convert a fixed length k-mer string to the binary encoding parameterized upon that same k

        :param s: The input k-mer
        :type s: str
        :returns: The kPal-inspired binary encoding
        :rtype: int

        """
        if type(s) is not str: # Typecheck the input k-mer
            raise TypeError("kdb.kmer_utils.sequence_to_binary expects a str as its argument")
        elif re.sub("[ACTG]", '', s) != '': # Check if non ACGT characters exist with regex
            if "N" in s: # k-mer with 'N' do not have a binary encoding
                return None
            else:
                raise TypeError("kdb.kmer_utils.sequence_to_binary expects the letters to contain only nucleotide symbols ATCG")
        idx = 0
        for c in bytes(s, "UTF-8"): # Use byteshifting for fast conversion to binary encoding
            idx = idx << 2
            idx = idx | self.letterToBinary[c]

        return idx

    def binary_to_sequence(self, b):
        """ Decode a k-mer from its encoding parameterized upon k

        :param b: The bytes encoding the string, the unique id of the k-mer in the 2-bit space. A 2bit byte space is the most compression DNA information can be stored in.
        :type b: int
        :returns: The decoded k-mer value
        :rtype: str

        """
        if type(b) is not int: # Type check the integer to decode
            raise TypeError("kdb.kmer_utils.binary_to_sequence expects an int as its first positional argument")
        result=''
        for i in range(self.k): # Use byte-shifting to convert bit pairs to letters
            result += self.binaryToLetter[b & 0x03]
            b = b >> 2
        result = list(result)
        result.reverse()
        return ''.join(result)

    def get_neighbors(self, kid):
        """ Acquire the neighboring k-mer ids. Doesn't guarantee the k-mer has been observed or exists in the biological space we're trying to constrain by scope, it only presents the id of possible k-mer ids, not the ones spanning the data observed under consideration by the users. That would be handled after the fileutil class during graph operations. Commands involving simpler querying of the counts in the file need to make use of ids presented by get_neighbors, instead of that being handled by the fileutil class. The size of the graph information is polynomial in memory limitations defined by the way the profile is used in memory and accessed randomly on disk via the index. Basically, any use of this function is contingent on the definition of questions asked of the software with specific biological examples.

        :param kid: A starting k-mer id
        :returns: A list of neighboring k-mer ids
        :rtype: list

        """
        if type(kid) is not int:
            raise TypeError("kdb.kmer_utils.get_neighbors expects an int as its first positional argument")
        kmer = self.binary_to_sequence(kid)

        first_chars = self.binaryToLetter.copy()
        last_chars = self.binaryToLetter.copy()
        first_chars.remove(kmer[0])
        last_chars.remove(kmer[-1])
        return list(map(self.sequence_to_binary, list(map(lambda x: x + kmer[1:], first_chars)) + list(map(lambda x: kmer[:-1] + x, last_chars))))


