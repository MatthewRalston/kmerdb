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
    if type(s) is not str:
     raise TypeError("kdb.kmer_converter.sequence_to_binary expects a str as its argument")
    elif re.sub("[ACTG]", '', s) != '':
     raise TypeError("kdb.kmer_converter.sequence_to_binary expects the letters to contain only nucleotide symbols ATCG")
    idx = 0
    for c in bytes(s, "UTF-8"):
        idx = idx << 2
        idx = idx | letterToBinary[c]

    return idx

def binary_to_sequence(b, k):
    if type(b) is not int:
        raise TypeError("kdb.kmer_converter.binary_to_sequence expects an int as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kdb.kmer_converter.binary_to_sequence expects an int as its second positional argument")
    result=''
    for i in range(k):
        result += binaryToLetter[b & 0x03]
        b = b >> 2
    result = list(result)
    result.reverse()
    return ''.join(result)

def get_neighbors(kid, k):
    if type(kid) is not int:
        raise TypeError("kdb.kmer_converter.get_neighbors expects an int as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kdb.kmer_converter.get_neightbors expects an int as its second positional argument")
    kmer = binary_to_sequence(kid, k)

    first_chars = binaryToLetter.copy()
    last_chars = binaryToLetter.copy()
    first_chars.remove(kmer[0])
    last_chars.remove(kmer[-1])
    return list(map(sequence_to_binary, list(map(lambda x: x + kmer[1:], first_chars)) + list(map(lambda x: kmer[:-1] + x, last_chars))))

