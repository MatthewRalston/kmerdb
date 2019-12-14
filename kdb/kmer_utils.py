import logging
logger = logging.getLogger(__file__)
import re

letterToBinary={ # Unicode UTF-8 byte codes for ACGT
    65: 0,
    67: 1,
    71: 2,
    84: 3
    }
binaryToLetter=['A', 'C', 'G', 'T']


def sequence_to_binary(s):
    """ Convert a fixed length k-mer string to the binary encoding parameterized upon that same k

    :param s: The input k-mer
    :type s: str
    :returns: The kPal-inspired binary encoding
    :rtype: int

    """
    if type(s) is not str:
     raise TypeError("kdb.kmer_utils.sequence_to_binary expects a str as its argument")
    elif re.sub("[ACTG]", '', s) != '':
        if "N" in s:
            return None
        else:
            raise TypeError("kdb.kmer_utils.sequence_to_binary expects the letters to contain only nucleotide symbols ATCG")
    idx = 0
    for c in bytes(s, "UTF-8"):
        idx = idx << 2
        idx = idx | letterToBinary[c]

    return idx

def binary_to_sequence(b, k):
    """ Decode a k-mer from its encoding parameterized upon k

    :param b: The bytes encoding the string, the unique id of the k-mer in the 2-bit space. A 2bit byte space is the most compression DNA information can be stored in.
    :type b: int
    :param k: The value to decode the k-mer string with from the bytes b.
    :type k: int
    :returns: The decoded k-mer value
    :rtype: str

    """
    if type(b) is not int:
        raise TypeError("kdb.kmer_utils.binary_to_sequence expects an int as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kdb.kmer_utils.binary_to_sequence expects an int as its second positional argument")
    result=''
    for i in range(k):
        result += binaryToLetter[b & 0x03]
        b = b >> 2
    result = list(result)
    result.reverse()
    return ''.join(result)

def get_neighbors(kid, k):
    """ Acquire the neighboring k-mer ids. Doesn't guarantee the k-mer has been observed or exists in the biological space we're trying to constrain by scope, it only presents the id of possible k-mer ids, not the ones spanning the data observed under consideration by the users. That would be handled after the fileutil class during graph operations. Commands involving simpler querying of the counts in the file need to make use of ids presented by get_neighbors, instead of that being handled by the fileutil class. The size of the graph information is polynomial in memory limitations defined by the way the profile is used in memory and accessed randomly on disk via the index. Basically, any use of this function is contingent on the definition of questions asked of the software with specific biological examples.

    :param kid: 
    :param k: 
    :returns: 
    :rtype: 

    """
    if type(kid) is not int:
        raise TypeError("kdb.kmer_utils.get_neighbors expects an int as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kdb.kmer_utils.get_neightbors expects an int as its second positional argument")
    kmer = binary_to_sequence(kid, k)

    first_chars = binaryToLetter.copy()
    last_chars = binaryToLetter.copy()
    first_chars.remove(kmer[0])
    last_chars.remove(kmer[-1])
    return list(map(sequence_to_binary, list(map(lambda x: x + kmer[1:], first_chars)) + list(map(lambda x: kmer[:-1] + x, last_chars))))

