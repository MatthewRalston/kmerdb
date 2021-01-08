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
import numpy as np
import yaml
import math
import functools
import time


import logging
logger = logging.getLogger(__file__)



from kmerdb import fileutil, index, kmer
import Bio

def markov_probability(seq:Bio.SeqRecord.SeqRecord, kdbrdr:fileutil.KDBReader, kdbidx:index.IndexReader):
    """
    :param seq:
    :type SeqRecord: Bio.SeqRecord.SeqRecord
    :returns
    :rtype:
    """
    if not isinstance(seq, Bio.SeqRecord.SeqRecord):
        raise TypeError("kmerdb.probability.MarkovProbability expects a Bio.SeqRecord.SeqRecord as its first positional argument")

    elif not isinstance(kdbrdr, fileutil.KDBReader):
        raise TypeError("kmerdb.probability.MarkovProbability expects a kmerdb.fileutil.KDBReader as its second positionl argument")
    elif not isinstance(kdbidx, index.IndexReader):
        raise TypeError("kmerdb.probability.MarkovProbability expects a kmerdb.index.IndexReader as its third positional argument")

        
    k = kdbrdr.header['k']
    N = 4**k
    total_kmers = sum([f['total_kmers'] for f in kdbrdr.header['files']])
    if len(seq) < k:
        raise ValueError("kmerdb.probability.MarkovProbability expects the sequence to be at least k={0} in length".format(k))
    x = seq.seq[:k]
    kmer_id_x = kmer.kmer_to_id(x)
    if kmer_id_x > N:
        raise ValueError("kmerdb.probability.MarkovProbability expected the k-mer id {0} of the subsequence x='{1}' to have an id less than N = 4^k = {2}".format(kmer_id_x, x, N))
    if not isinstance(kdbidx.index, np.ndarray):
        raise ValueError("kmerdb.probability.MarkovProbability expects the kdbidx to be a Numpy array")
    elif kdbidx.index.size != N:
        raise ValueError("kmerdb.probability.MarkovProbability expects the kdbidx to have exactly {0} items, found {1}".format(N, kdbidx.size))

    
    kmer_id, count, neighbors = index.read_line(kdbrdr, kdbidx, kmer_id_x)
    if kmer_id is None or count is None or neighbors is None:
        logger.error("K-mer id: {0}\n".format(kmer_id_x, ))
        logger.error(kmer.id_to_kmer(kmer_id_x, k))
        time.sleep(1)
        raise RuntimeError("Index encountered an invalid initial k-mer id index value.")
    px =  count / float(total_kmers)
    #kmer_seq = kmer.id_to_kmer(km
    # Find the correct neighbor and calculate the transition probability
    prefix = x
    total = 0
    product = 1
    for i in range(k, len(seq.seq)):
        sum_qses = 0
        for s in neighbors["suffixes"]:
            kmer_id, count, _ = index.read_line(kdbrdr, kdbidx, s)
            if kmer_id is None or count is None or _ is None:
                raise ValueError("k-mer id '{0}' had an offset of zero in the index, it was not observed in the genome. Effective probability of sequence is 0 and Log-odds ratio of the sequence being generated from the k-mer profile is effectively 0 as well.")
            if s != kmer_id:
                raise ValueError("k-mer id read from file did not match the provided neighbor suffix id")
            sum_qses += count/float(total_kmers)
        # Prefix is the next prefix sequence to find a suffix for
        prefix = prefix[1:] + seq.seq[i]
        # We take the prefixes count and new neighbors
        kmer_id, count, neighbors = index.read_line(kdbrdr, kdbidx, kmer.kmer_to_id(prefix))
        # Now we calculate qt similar to px, the frequency (count / total number of k-mers, read from the header)
        qt = count/float(total_kmers)
        # Now we calculate the transition probability of sequence s to sequence t
        # as qt divided by the sum of each qsc
        ast = qt / sum_qses
        # For the final log odds ratio, we sum the logs of the transition probabilities
        total += math.log(ast)
        # For the final probability, we multiply the transition probabilities
        product *= ast
    # Fianlly qxis represents our null model in a way.
    # Suppose that the null model consisted of the same number of k-mer as observed,
    # but uniformly distributed amongst all the k-mers
    qxis = [total_kmers/N for x in range(k, len(seq.seq))]
    # We put this on the bottom of the Log-odds Ratio because the null model contrasts the transition probabilities
    lor = log(px) + total / sum([math.log(qx) for qx in qxis])
    # 
    pseq = px*product#/functools.reduce(lambda x,y: x*y, qxis)

    return {"seq": seq, "log_odds_ratio": lor, "p_of_seq": pseq}
