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
import gzip
import math
import logging

logger = logging.getLogger(__file__)

import numpy as np

from kmerdb import kmer, config, util



def select_lexicographical_minimizers(seq, seq_id, k, window_size):
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
    :rtype: int, np.array, np.array
    """
    N = 4**k
    if type(seq) is not str:
        raise TypeError("kmerdb.minimizer.select_lexicographical_minimizers() expects a str as its first positional argument")
    elif type(seq_id) is not str:
        raise TypeError("kmerdb.minimizer.select_lexicographical_minimizers() expects a str as its second positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.minimizer.select_lexicographical_minimizers() expects an int as its third positional argument")
    elif type(window_size) is not int:
        raise TypeError("kmerdb.minimizer.select_lexicographical_minimizers() expects an int as its fourth positional argument")
    
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
        else:
            is_min = 0
            is_minimizer[i] = 0
        # 
        subseq = seq[i:i+k]
        if len(subseq) != k:
            raise ValueError("WHOOPS")
        kmer_id = kmer.kmer_to_id(subseq)

        coords.append((seq_id, i, kmer_id, is_min))
    """
    # coords is [(seq_id, i:L, kmer_id, is_minimizer), ...]
    """
    return coords # 4-tuple of (seq_id:str, i:int, kmer_id:int, is_min:bool/int)
    #return math.floor(L / window_size), is_minimizer, coords # coords is a list of kmer start positions in the input fasta sequence, referencing the kmer_id of the k-mer at the minimizer position, is_minimizer has dimension L = N-k+1, and the number of minimizers selected as first positional return, is_minimizer array as second positional return, and coordinate array with kmer_id values and 1:L/window_size number of elements pointing to the kmer_ids associated with the threetuple of (fasta_id, fasta_coord, kmer_id)


def read_fasta_sequences(fasta_fh):
    """
    Parse sequence identifiers and sequences into two arrays.

    
    """
    fastas = []
    from Bio import SeqIO
    try:
        for record in SeqIO.parse(fasta_fh, "fasta"):
            fastas.append({
                "seq": str(record.seq),
                "id": str(record.id)
            })
    except Exception as e:
        logger.error("kmerdb.minimizer.read_sequences_from_fastafile() encountered an error while reading this file and needs to exit")
        raise e
    finally:
        fasta_fh.close()
    return fastas


#def get_minimizers_from_kdbfile_and_input_fastafile(kdbfile:str, fasta_file:str, window_size, seqs:list=None, ids:list=None):
def read_minimizers(kdbfile:str, fasta_file:str, window_size):

    from kmerdb import fileutil

    if type(kdbfile) is not str:
        raise TypeError("kmerdb.minimizer.read_minimizers() expects a str as its first positional argument")
    elif type(fasta_file) is not str:
        raise TypeError("kmerdb.minimizer.read_minimizers() expects a str as its second positional argument")
    elif type(window_size) is not int:
        raise TypeError("kmerdb.minimizer.read_minimizers() expects an int as its third positional argument")
    
    if not os.path.exists(kdbfile) or not os.access(kdbfile, os.R_OK):
        raise IOError("kmerdb.minimizer.read_minimizers() could not read '{0}' from the filesystem".format(kdbfile))
    elif not os.path.exists(fasta_file) or not os.access(fasta_file, os.R_OK):
        raise IOError("kmerdb.minimizer.read_minimizers() could not read '{0}' from the filesystem".format(fasta_file))
    
    if seqs is None or ids is None:
        if fasta_file.endswith(".gz"):
            with gzip.open(fasta_file, mode='r') as ifile:
                fastas = read_fasta_sequences(fasta_file)
        else:
            with open(fasta_file, mode='r') as ifile:
                fastas = read_fasta_sequences(fasta_file)
                
    coord_map = []
    coord_obj_array = []
    is_min = {}

    raise RuntimeError("WOOPS")
    
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
            # is minimizer is a boolean array 0 or 1
            # coords is a numpy 32 array of kmer_ids
            sed_id = fasta_ids[i]
            coords = select_lexicographical_minimizers(seq, seq_id, metadata['k'], window_size, kmer_ids)

            raise RuntimeError("This needs fixing too...")
            """
            # coords is [(seq_id, i:L, kmer_id, is_minimizer), ...]
            """

            #is_min = [c[3] for c in coords]

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
    for fasta_id, i, kmer_id, is_min in coordinate_mapping:
        assert type(fasta_id) is str, "FASTA ID was not a sequence id string"
        assert type(i) is int, "Sequence coordinate was not a number"
        assert type(kmer_id) is int, " - error no coordinate \_( |>\">| )_/ ... ( ...wtf?)"
        assert type(is_min) is int, "should be typed numeric 1s and 0s"
        print("\t".join((fasta_id, sequence, basepair, is_min)))
    return






def make_minimizers(fasta_file:str, kdb_file:str, window_size:int):
    """
    :param fasta_file:
    :type str:
    :param kdb_file:
    :type int:
    :param window_size:
    :type int:
    :returns: (coordmap, k) # coordmap is a dictionary keyed on sequence id, whose value is list of 4-tuples of (seq_id, coordinate, kmer_id, is_min:bool)
    """

    from kmerdb import fileutil, config, util

    metadata = None
    N = None

    coordmap = {}
    mins = {}

    fasta_sfx = os.path.splitext(fasta_file)[-1]
    kdb_sfx = os.path.splitext(kdb_file)[-1]

    if type(fasta_file) is not str:
        raise TypeError("kmerdb.minimizer.make_minimizers() expects a str as its first positional argument")
    elif type(kdb_file) is not str:
        raise TypeError("kmerdb.minimizer.make_minimizers() expects a str as its second positional argument")

    if fasta_sfx == ".gz" and (fasta_file.endswith(".fasta.gz") or fasta_file.endswith(".fna.gz") or fasta_file.endswith(".fa.gz")):
        pass
    elif fasta_sfx != ".fasta" and fasta_sfx != ".fna" and fasta_sfx != ".fa": 
        raise IOError("Reference .fasta filepath does not end in '.fasta'/'.fna.'.fa'")
    elif not os.path.exists(fasta_file) or not os.access(fasta_file, os.R_OK):
        raise IOError("Reference .fasta filepath '{0}' does not exist on the filesystem".format(fasta_file))
    elif kdb_sfx != ".kdb":
        raise IOError("Reference .kdb filepath does not end in '.kdb'")
    elif not os.path.exists(kdb_file) or not os.access(kdb_file, os.R_OK):
        raise IOError("Reference .kdb filepath '{0}' does not exist on the filesystem".format(kdb_file))

    if fasta_file.endswith(".gz"):
        logger.info("Reading .gz compressed fasta file '{0}'...".format(fasta_file))
        with gzip.open(fasta_file, mode='r') as ifile:
            fastas = read_fasta_sequences(ifile)
    else:
        logger.info("Reading fasta file '{0}'...".format(fasta_file))
        with open(fasta_file, mode='r') as ifile:
            fastas = read_fasta_sequences(ifile)
    """
    Fixme
    """
    logger.info("Reading .kdb file '{0}'...".format(kdb_file))
    with fileutil.open(kdb_file, mode='r', slurp=True) as kdb_in:
        metadata = kdb_in.metadata

        k = metadata["k"]
        kmer_ids_dtype = metadata["kmer_ids_dtype"]
        N = 4**metadata["k"]
        if metadata["version"] != config.VERSION:
            sys.stderr.write("KDB version is out of date, may be incompatible with current KDBReader class\n")
        kmer_ids = kdb_in.kmer_ids
        for s in fastas:
            seq = s["seq"]
            seq_id = s["id"]
            logger.debug("Collecting lexicographical minimizers for sequence '{0}'...".format(seq_id))
            coords = select_lexicographical_minimizers(seq, seq_id, metadata['k'], window_size)
            coordmap[seq_id] = coords
    logger.info("Completed loading minimizers from .kdb file '{0}'...".format(kdb_file))

    return coordmap, k




    
def read_minimizer_kdbi(kdbi_file):
    """
    Read a minimizer file into memory
    """
    kdbi_sfx = os.path.splitext(kdbi_file)[-1]
    
    if kdbi_sfx != ".kdbi.1": # A filepath with invalid suffix
        raise IOError("Input .kdb index filepath '{0}' does not end in '.kdbi.1'".format(kdbi_file))

    minimizers = []
    with open(kdbi1, 'r') as minimizer_index_file:
        for line in minimizer_index_file.readlines():
            line = lin.split("\t")

            assert len(line) == 4, "Index file corrupted. Use 'kmerdb minimizers' to regenerate the index file."
            seq_id, i, kmer_id, is_min = line
            minimizers.append( (seq_id, i, kmer_id, is_min) )
    return minimizers

def print_minimizer_kdbi(coordmap, filename, dump_all:bool=False):
    """
    Print a minimizer array to a .kdbi.1 index file
    """
    if type(coordmap) is not dict:
        raise TypeError("kmerdb.minimizer.print_minimizer_kdbi() expects its first positional argument to be a list")
    elif not all((type(c) is tuple and type(c[0]) is str and type(c[1]) is int and type(c[2]) is int and type(c[3]) is int) for seq_id, coords in coordmap.items() for c in coords):
        raise ValueError("kmerdb.minimizer.print_minimizer_kdbi() expects its first positional argument to be a coordinate-mapping of minimizers")
    elif type(filename) is not str:
        raise TypeError("kmerdb.minimizer.print_minimizer_kdbi() expects its second positional argument to be a str")
    elif type(dump_all) is not bool:
        raise TypeError("kmerdb.minimizer.print_minimizer_kdbi() expects the keyword argument 'dump_all' to be a bool")

    
    with open(filename, 'w') as ofile:
        for seq_id, coords in coordmap.items():
            for c in coords:
                if c[-1] == 0 and dump_all is False:
                    continue
                ofile.write("\t".join(list(map(str, c))) + "\n")

    return filename
            
