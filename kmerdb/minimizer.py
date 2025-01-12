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
import os

#import yaml
#from collections import OrderedDict

import numpy as np

from kmerdb import kmer, config, util

def select_lexicographical_minimizers(seq, k, window_size, kmer_ids):
    """
    Select k-mers based on lexicographical minimizers and return a binary array indicating selected k-mers.
    
    :param seq: Input DNA sequence
    :type seq: str
    :param k: Choice of k
    :type k: int
    :param window_size: Size of the sliding window
    :type window_size: int
    :param kmer_ids: Array of kmer IDs (1 to 4^k)
    :type kmer_ids: list

    :returns: Binary array where 1 indicates selected k-mer, 0 otherwise
    :rtype: np.array
    """

    N = len(kmer_ids)
    
    if len(seq) < k:
        raise ValueError("Sequence length was less than k")

    minimizers = np.zeros(N, dtype="int16")

    for i in range(N - k + 1):
        if i % window_size == 0:
            subseq = seq[i:i+k]
            kmer_id = kmer.kmer_to_id(subseq)
            minimizers[kmer_id] = 1


    return minimizers

