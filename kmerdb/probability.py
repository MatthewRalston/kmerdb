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
        raise TypeError("kmerdb.probability.markov_probability expects a Bio.SeqRecord.SeqRecord as its first positional argument")

    elif not isinstance(kdbrdr, fileutil.KDBReader):
        raise TypeError("kmerdb.probability.markov_probability expects a kmerdb.fileutil.KDBReader as its second positionl argument")
    elif not isinstance(kdbidx, index.IndexReader):
        raise TypeError("kmerdb.probability.markov_probability expects a kmerdb.index.IndexReader as its third positional argument")

        
    k = kdbrdr.metadata['k']
    mononucleotides = kdbrdr.metadata['files'][0]['mononucleotides']
    N = 4**k
    total_kmers = sum([f['total_kmers'] for f in kdbrdr.metadata['files']])
    um = 1/N # This is the uniform frequency for the space
    if len(seq) < k:
        raise ValueError("kmerdb.probability.markov_probability expects the sequence to be at least k={0} in length".format(k))
    x = seq.seq[:k]
    kmer_id_x = kmer.kmer_to_id(x)
    if kmer_id_x > N:
        raise ValueError("kmerdb.probability.markov_probability expected the k-mer id {0} of the subsequence x='{1}' to have an id less than N = 4^k = {2}".format(kmer_id_x, x, N))
    if not isinstance(kdbidx.index, np.ndarray):
        raise ValueError("kmerdb.probability.markov_probability expects the kdbidx to be a Numpy array")
    elif kdbidx.index.size != N:
        raise ValueError("kmerdb.probability.markov_probability expects the kdbidx to have exactly {0} items, found {1}".format(N, kdbidx.size))

    
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
    total2 = 0
    product = 1
    product2 = 1
    total_nucleotides = sum(mononucleotides.values())
    logger.debug(seq.seq)
    logger.debug(len(seq.seq))
    for i in range(k, len(seq.seq)):

        sum_qses = 0
        sum_aij = 0

        for char, idx in neighbors["suffixes"].items():
            kmer_id, count, _ = index.read_line(kdbrdr, kdbidx, idx)
            if kmer_id is None or count is None or _ is None:
                kmer_id = idx
                count = 0
                #raise ValueError("k-mer id '{0}' had an offset of zero in the index, it was not observed in the genome. Effective probability of sequence is 0 and Log-odds ratio of the sequence being generated from the k-mer profile is effectively 0 as well.".format(s))
            aij = mononucleotides[char] / total_nucleotides
            sum_aij += mononucleotides[char] / total_nucleotides
            sum_qses += count/float(total_kmers)
        # Prefix is the next prefix sequence to find a suffix for
        if i >= len(seq.seq):
            break
            #raise ValueError("kmerdb.probability.markov_probability encountered a value of i={0} larger than the length of the sequence {1}".format(i, len(seq.seq)))
        elif seq.seq[i] == "":
            logger.warning("Reached the end of the sequence, i={0} is producing no letters".format(i))
            break
        else:
            logger.debug("{0} : prefix: '{1}' + '{2}'".format(i, prefix[1:], seq.seq[i]))
            prefix = prefix[1:] + seq.seq[i]

        if len(prefix) != k:
            break
        else:
            logger.debug("New k-mer: '{0}'".format(prefix))
            logger.debug("Length of prefix: {0}".format(len(prefix)))
            try:
                kmerid = kmer.kmer_to_id(prefix)
            except KeyError as e:
                logger.debug(i)
                logger.debug("this should be the new letter: '{0}'".format(seq.seq[i]))
                logger.debug("this should be the existing prefix: '{0}'".format(prefix))
                raise e
            logger.debug("Iteration number {0}, calculating probabilities for k-mer id={1} seq={2}".format(i, kmerid, prefix))
            logger.debug("Reading counts/neighbors from file according to k-mer id")
            # We take the prefixes count and new neighbors
            kmer_id, count, neighbors = index.read_line(kdbrdr, kdbidx, kmerid)

            # Now we calculate qt similar to px, the frequency (count / total number of k-mers, read from the metadata)
        
            qt = count/float(total_kmers)
            # Now we calculate the transition probability of sequence s to sequence t
            # as qt divided by the sum of each qsc
            ast = qt / sum_qses

            # For the final log odds ratio, we sum the logs of the transition probabilities
            total += math.log10(ast)
            total2 += math.log10(sum_aij) # The ratio of the uniform frequency to the first order Markov background
            # For the final probability, we multiply the transition probabilities
            product *= ast
            product2 *= aij
    # Fianlly qxis represents our null model in a way.
    # Suppose that the null model consisted of the same number of k-mer as observed,
    # but uniformly distributed amongst all the k-mers
    qxis = [um for x in range(k, len(seq.seq))]
    # We put this on the bottom of the Log-odds Ratio because the null model contrasts the transition probabilities
    null_model = math.log10(um) + total2
    lor = (math.log10(px) + total) / null_model
    # 
    pseq = px*product#/functools.reduce(lambda x,y: x*y, qxis)
    prand = product2/N
    return {"seq": seq, "log_odds_ratio": lor, "p_of_seq": pseq}
