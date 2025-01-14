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


def minimizers_from_fasta_and_kdb(fasta_file, kdb_file):
    from Bio import SeqIO
    import numpy as np
    from kmerdb import fileutil, config, util, minimizer
    import json
    metadata = None
    N = None
    kmer_ids = None
    fasta_seqs = []
    mins = None

    fa_sfx = os.path.splitext(fasta_file)[-1]
    kdb_sfx = os.path.splitext(arguments.kdb)[-1]

    
    if kdb_sfx != ".kdb" and kdb_sfx != ".kdbg": # A filepath with invalid suffix
        raise IOError("Input .kdb filepath '{0}' does not end in '.kdb'".format(kdb_file))
    elif not os.path.exists(kdb_file):
        raise IOError("Input .kdb filepath '{0}' does not exist on the filesystem".format(kdb_file))
    elif fa_sfx != ".fa" and fa_sfx != ".fna" and fa_sfx != ".fasta":
        raise IOError("Input .fasta filepath '{0}' does not exist on the filesystem".format(arguments.fasta))

    with open(fasta_file, mode="r") as ifile:
        for record in SeqIO.parse(ifile, "fasta"):
            fasta_seqs.append(str(record.seq))

    with fileutil.open(arguments.kdb, mode='r', slurp=True) as kdb_in:
        metadata = kdb_in.metadata

            
        kmer_ids_dtype = metadata["kmer_ids_dtype"]
        N = 4**metadata["k"]
        if metadata["version"] != config.VERSION:
            logger.log_it("KDB version is out of date, may be incompatible with current KDBReader class", "WARNING")
        kmer_ids = kdb_in.kmer_ids

        for seq in fasta_seqs:

            minimizers = minimizer.select_lexicographical_minimizers(seq, metadata['k'], arguments.window_size, kmer_ids)

            # TODO: FIXME
            
            mins_as_list = list(map(int, minimizers))

            # 
            if mins is None:
                mins=mins_as_list
            else:
                for kmer_is_selected, i in enumerate(mins_as_list):
                    if mins_as_list[i] == 1:
                        pass
                    elif kmer_is_selected == 1:
                        mins[i] = 1
                
    return minimizers

def make_alignment(reference_fasta_seqs, query_fasta_seqs, k):


    reference_minimizers = select_lexicographical_minimizers()


def read_minimizer_kdbi_index_file(kdbi1):
    """
    Read a minimizer file into memory
    """
    minimizers = []
    with open(kdbi1, 'r') as minimizer_index_file:
        for l in minimizer_index_file:

            line = l.strip("\n").split("\t")

            assert len(line) == 2, "Index file corrupted. Use 'kmerdb minimizer' to regenerate the index file."
            kmer_id, is_minimizer = line
            
            if is_minimizer == 1:
                minimizers.append(kmer_id)
    return minimizers
