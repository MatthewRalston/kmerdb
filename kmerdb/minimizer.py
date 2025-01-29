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

import math


#import yaml
#from collections import OrderedDict




import numpy as np

from kmerdb import kmer, config, util

#from kmerdb import logger as kdbLogger
#global logger = 




def select_lexicographical_minimizers(seq, seq_id, k, window_size, kmer_ids):
    """
    Select k-mers based on lexicographical minimizers and return a binary array indicating selected k-mers.
    
    :param seq: Input DNA sequence
    :type seq: str
    :param seq_id: Sequence ID as fasta header string
    :type seq_id: str
    :param k: Choice of k
    :type k: int
    
    :param window_size: Size of the sliding window
    :type window_size: int
    :param kmer_ids: Array of kmer IDs (1 to 4^k)
    :type kmer_ids: list

    :returns: Binary array where 1 indicates selected k-mer, 0 otherwise
    :rtype: int, np.array
    """

    
    
    N = len(kmer_ids)
    if not type(seq) is str:
        raise TypeError("kmerdb.minimizer.select_lexicographical_minimizers expects a str as its first positional argument")        
    if len(seq) < k:
        raise ValueError("Sequence length was less than k")
    
    sequence_length = len(seq)
    
    is_minimizer = np.zeros(sequence_length, dtype="uint32")

    L = N - k + 1


    coords = []

    for i in range(L): # Number of minimizers = floor(sequence_length/window_size)
        if i % window_size == 0:
            is_min = 1
            is_minimizer[i] = 1
            is_min = 1
        else:
            is_min = 0
            is_minimizer[i] = 0
            
                
        subseq = seq[i:i+k]
        kmer_id = kmer.kmer_to_id(subseq)

        coords.append((seq_id, i, kmer_id, is_min))
    
    # coords is [(seq_id, i:L, kmer_id, is_minimizer), ...]            
    return math.floor(L / window_size), is_minimizer, coords # coords is a list of kmer start positions in the input fasta sequence, referencing the kmer_id of the k-mer at the minimizer position, is_minimizer has dimension L = N-k+1, and the number of minimizers selected as first positional return, is_minimizer array as second positional return, and coordinate array with kmer_id values and 1:L/window_size number of elements pointing to the kmer_ids associated with the threetuple of (fasta_id, fasta_coord, kmer_id)



    # len(coord.subtract.keys)

def read_sequences_from_fastafile(fasta_file):
    """
    Parse sequence identifiers and sequences into two arrays.
    """
    fasta_seqs = []
    fasta_ids = []
    from Bio import SeqIO
    with open(fasta_file, mode="r") as ifile:
        for record in SeqIO.parse(ifile, "fasta"):
            fasta_seqs.append(str(record.seq))
            fasta_ids.append(str(record.id))
    return fasta_seqs, fasta_ids


def get_minimizers_from_kdbfile_and_input_fastafile(kdbfile:str, fasta_file:str, window_size, seqs:list=None, ids:list=None):


    from kmerdb import fileutil
    if seqs is None or ids is None:
        fasta_seqs, fasta_ids = read_sequences_from_fastafile(fasta_file)
    elif not all(type(fasta_id) is str for fasta_id, i in ids):
        raise TypeError("fasta_ids list is not a list of fasta sequence ids")
    elif not all(type(fasta_seq) is str for fasta_seq, j in seqs):
        raise ValueError("error no sequences")
    else:
        # Unimplemented
        fasta_seqs = seqs
        fasta_ids = ids

    coord_map = []
    coord_obj_array = []
    is_min = {}
        
    if fasta_file is None or type(fasta_file) is not str:
        raise TypeError("")
    elif type(kdbfile) is not str or not os.path.exists(kdbfile):
        raise ValueError("kmerdb.minimizer.get_minimizers_from_kdbfile_and_input_fastafile .kdb path does not exist on filesystem")
    else:
        with fileutil.open(kdbfile, mode='r', slurp=True) as kdb_in:
            metadata = kdb_in.metadata

            k = metadata["k"]
            kmer_ids_dtype = metadata["kmer_ids_dtype"]
            N = 4**metadata["k"]
            if metadata["version"] != config.VERSION:
                sys.stderr.write("KDB version is out of date, may be incompatible with current KDBReader class\n")
            kmer_ids = kdb_in.kmer_ids


            
            for i, seq_id in enumerate(fasta_ids):
                seq = fasta_seqs[i] 
                # print(seq)
                # print(seq_id)
                # print(i)

                # sys.exit(1)
                
                # is minimizer is a boolean array 0 or 1
                # coords is a numpy 32 array of kmer_ids
                sed_id = fasta_ids[i]
                seqlen, is_minimizer, coords = select_lexicographical_minimizers(seq, seq_id, metadata['k'], window_size, kmer_ids)

                is_min[seq_id] = is_minimizer
                """
                Fixme! coordmap should be primarily indexed on coordinate
                """
                coord_map = [[] for _ in range(seqlen)]
                #coord_map[kmer_ids[i]] = [] # coordmap becomes hashmap structure with kmer_id, fasta_seq_id1, start_coord1 , fasta_seq_id2, start_coord2

                
                for j in range(seqlen):
                    subseq = seq[j:j+k]
                    kmer_id = kmer.kmer_to_id(subseq)
                    coord_map.append({"fasta_id": fasta_ids[i], "sequence_length": seqlen, "coord": j, "kmer_id": kmer_id, "is_minimizer": is_minimizer[j]}) # Is_min should be the jth sequence element as a k-mer, and whether or not the jth k-mer is used as a minimizer


                    # fasta_id[i] where i is number of queries, len, the coord, and the is_minimizer bool


                    coord_obj_array.append((fasta_ids[i], seqlen, j, kmer_id, is_minimizer[j]))
                    # As hashmapped structure
                    #coord_obj_array[fasta_ids[i]].push((fasta_ids[i], seqlen, j, is_minimizer[j]))
                
                    
        return k, kmer_ids, fasta_ids, is_min, coord_map, coord_obj_array # kmer_ids are kmer actg ids, as read from the kdb file, fasta_ids is an array of inputs sequence ids, coordmap is a list of hashmaps


def print_coordmap(coordinate_mapping):
    """
    A hash of hashes, 2D associative array, not implemented further.


    
    
    """


    for fasta_id, sequence_length, basepair, is_min in coordinate_mapping:
        assert type(fasta_id) is str, "FASTA ID was not a sequence id string"
        assert type(sequence_length) is int, "Sequence length input was not a number"
        assert type(basepair) is int, " - error no coordinate \_( |>\">| )_/ ... ( ...wtf?)"
        assert type(is_min) is int, "should be typed numeric 1s and 0s"
        
        print("\t".join((fasta_id, sequence, basepair, is_min)))
    return






def minimizers_from_fasta_and_kdb(fasta_file, kdb_file, window_size):
    from Bio import SeqIO
    import numpy as np
    from kmerdb import fileutil, config, util, minimizer
    import json
    metadata = None
    N = None
    kmer_ids = None
    fasta_seqs = []
    fasta_ids = []
    mins = {}

    fa_sfx = os.path.splitext(fasta_file)[-1]
    kdb_sfx = os.path.splitext(kdb_file)[-1]

    
    if kdb_sfx != ".kdb" and kdb_sfx != ".kdbg": # A filepath with invalid suffix
        raise IOError("Input .kdb filepath '{0}' does not end in '.kdb'".format(kdb_file))
    elif not os.path.exists(kdb_file):
        raise IOError("Input .kdb filepath '{0}' does not exist on the filesystem".format(kdb_file))
    elif fa_sfx != ".fa" and fa_sfx != ".fna" and fa_sfx != ".fasta":
        raise IOError("Input .fasta filepath '{0}' does not exist on the filesystem".format(arguments.fasta))

    fasta_ids, fasta_seqs = read_sequences_from_fastafile(fasta_file)

    
            

    """
    Fixme
    """ 

            
    k, kmer_ids, _, is_minimizer, coordmap, coord_obj_array = get_minimizers_from_kdbfile_and_input_fastafile(kdb_file, fasta_file, window_size)

    return k, kmer_ids, fasta_ids, is_minimizer, coordmap, coord_obj_array # kmer_ids are kmer actg ids, as read from the kdb file, fasta_ids are a hashmap of inputs sequence ids, coordmap is a list of coordinate-as-index based minimizer locations. coordmap = [{"fasta_id": ..., "seqlen": n, }, ...] 



    
def read_minimizer_kdbi_index_file(kdbi1):
    """
    Read a minimizer file into memory
    """

    if kdbi_sfx != ".kdbi.1": # A filepath with invalid suffix
        raise IOError("Input .kdb index filepath '{0}' does not end in '.kdb'".format(arguments.kdbi1))

    minimizers = []
    with open(kdbi1, 'r') as minimizer_index_file:
        for l in minimizer_index_file:

            line = l.strip("\n").split("\t")

            assert len(line) == 2, "Index file corrupted. Use 'kmerdb minimizer' to regenerate the index file."
            kmer_id, is_minimizer = line
            
            if is_minimizer == 1:
                minimizers.append(kmer_id)
    return minimizers

def print_minimizer_kdbi_index_file(mins, kmer_ids, filename):
    """
    Print a minimizer array to a .kdbi.1 index file
    """

    if type(mins) is not list:
        raise TypeError("kmerdb.minimizer.print_minmimizer_kdbi_index_file epects its first positional argument to be a list")
    elif type(kmer_ids) is not list:
        raise TypeError("kmerdb.minimizer.print_minmimizer_kdbi_index_file epects its second positional argument to be a list")
    elif len(mins) != len(kmer_ids):
        raise ValueError("kmerdb.minimizer.print_minmimizer_kdbi_index_file expects the number of minimizers should match the number of potential kmers")
    elif type(filename) is not str:
        raise TypeError("kmerdb.minimizer.print_minmimizer_kdbi_index_file expects the third positional argument to be a str")
    
    with open(filename, 'w') as ofile:
        for kmer, i in enumerate(kmer_ids):
            ofile.write("{0}\t{1}\n".format(kmer, mins[i]))

    return filename
            
