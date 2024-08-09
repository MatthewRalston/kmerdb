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
import io
import gzip

import psutil
import tempfile

import yaml, json
from collections import deque, OrderedDict
import math
import re
from itertools import chain


from builtins import open as _open

import jsonschema
from Bio import SeqIO, Seq, bgzf


import numpy as np
#import networkx as nx

from kmerdb import fileutil, parse, kmer, config, util

# Logging configuration
global logger
logger = None



def open(filepath, mode="r", metadata=None, slurp:bool=False, logger=None):
    """
    Opens a file for reading or writing. Valid modes are 'xrwbt'. 
    Returns a lazy-loading KDBGReader object or a KDBGWriter object.
    The data may be force loaded with 'slurp=True'


    Opens a .kdbg file, still only view goals.

    Opens for reading or writing, legacy function.


    fn sig: (

    filepath
    mode(r/w param)
    metadata:None
    slurp: False

    )


    Stands in front of the kdb classes.

     also available in more general fileutil.py                                        

    metadata data vector (jsonschema)         described in config.py

    squash target unassigned on todo (slurp)


    nodes
    weights

    dtype defs, avail in json

    graph description avail in json

    traversal not implemented.


    


    :param filepath:
    :type filepath: str
    :param mode:
    :type mode: str
    :param metadata: The file header/metadata dictionary to write to the file.
    :type metadata: dict
    :param slurp: Immediately load all data into KDBGReader
    :type slurp: bool
    :raise ValueError: mode parameter related errors
    
    :returns: kmerdb.fileutil.KDBGReader/kmerdb.fileutil.KDBGWriter
    :rtype: kmerdb.fileutil.KDBGReader
    """
    if type(filepath) is not str:
        raise TypeError("kmerdb.graph.open expects a str as its first positional argument")
    # elif not os.access(filepath, os.R_OK):
    #     raise ValueError("kmerdb.graph.open expects an existing filepath as its first positional argument")
    elif type(mode) is not str:
        raise TypeError("kmerdb.graph.open expects the keyword argument 'mode' to be a str")
    elif (mode == "w" or mode == "x") and (metadata is not None and (isinstance(metadata, OrderedDict) or type(metadata) is dict)):
        pass
    elif (mode == "w" or mode == "x") and metadata is not None:
        raise TypeError("kmerdb.graph.open expects an additional metadata dictionary")
    elif type(slurp) is not bool:
        raise TypeError("kmerdb.graph.open expects a boolean for the keyword argument 'slurp'")
    modes = set(mode)
    if modes - set("xrwbt") or len(mode) > len(modes):
        raise ValueError("invalid mode: {}".format(mode))

    # oooo....


    # .......oooOO.ooo.....



    # ..........o.o.o.....O..


    # 
    creating = "x" in modes
    reading  = "r" in modes
    writing  = "w" in modes
    binary   = "b" in modes
    text     = "t" in modes

    if text and binary:
        raise ValueError("can't have text and binary mode at once")
    elif not (creating or reading or writing):
        raise ValueError("must have exactly one or read/write")

    if "r" in mode.lower():
        return KDBGReader(filename=filepath, mode=mode, slurp=slurp)
    elif "w" in mode.lower() or "x" in mode.lower():
        return KDBGWriter(filename=filepath, mode=mode, metadata=metadata)
    else:
        raise ValueError("Bad mode %r" % mode)


    self.logger = logger
    self._loggable = logger is not None


def bytesize_of_metadata(metadata):
    """
    defers to util.get_bytesize_of_metadata
    """
    return util.bytessize_of_metadata

    
def parse_kdbg_table_line(kdbg_table_line:str, sort_arr:bool=False, row_dtype:str="uint64"):
    """
    Parses a line according to the expected .kdbg syntax, and returns the python data types expected as a tuple.

    .kdbg


    table def
    =========================
    vector 'e'

    i          uint64
    node1_id         uint64
    node2_id         uint64
    w          uint64

    :param line:
    :type kdbg_table_line: str
    :raise TypeError: .kdbg table line was not str
    :raise TypeError: sort_arr was not bool
    :raise ValueError: .kdbg table has 4 columns - index, node1_id, node2_id, edge_weight
    :returns: node1_id, node2_id, weight
    :rtype: tuple
    """

    if type(kdbg_table_line) is not str:
        raise TypeError("kmerdb.graph.parse_kdbg_table_line expects a str as its first positional argument")
    elif type(kdbg_table_line) is str and kdbg_table_line == "":
        raise StopIteration("empty table line")
    
    if type(sort_arr) is not bool:
        raise TypeError("kmerdb.graph.parse_kdbg_table_line needs 'sort_arr' to be bool")
    try:
        np.dtype(row_dtype)
    except TypeError as e:
        sys.stderr.write("Invalid numpy dtype\n")
        raise e
    finally:



        linesplit = kdbg_table_line.rstrip().split("\t")
        
        if len(linesplit) != 4:
            sys.stderr.write("Full line:\n{0}".format(line) + "\n")
            raise ValueError("kmerdb.graph.parse_kdbg_table_line encountered a .kdbg line without 4 columns. Invalid format for edge list")
        else:
            i, node1_id, node2_id, weight = linesplit
            i, node1_id, node2_id, weight = int(i), int(node1_id), int(node2_id), int(weight)
        
            e = np.array([node1_id, node2_id, weight], dtype=row_dtype)
            
            # Still needs the sort_array default functionality as switchable.
            # Still needs the solver (ss-sp, or bfs)
            
            # # TODO: bfs would necessitate the objective
            # vanilla cpu implementation
            
            # # TODO:
            # Needs the node delexer.
            # is this a pre or post sort?
            # should be a post
            # Needs the signature and the implementation
            
            # Objective function def
            return [i] + list(e)
    

            

def parsefile(filepath:str, k:int, quiet:bool=True, b:int=50000, both_strands:bool=False, logger=None):
    """
    Parse a single sequence file.
    
    :param filepath: Path to a fasta or fastq file
    :type filepath: str
    :param k: Choice of k to shred k-mers with
    :type k: int
    :param b: Number of reads (per block) to process in parallel
    :type b: int
    :param both_strands: Strand specificity argument for k-mer shredding process
    :type both_strands: bool
    :raise TypeError: first argumentshould be a valid filepath
    :raise TypeError: second argument k should be an int
    :raise TypeError: quiet keyword arg should be a bool
    :raise TypeError: *unimplemented* - both_strands should be a bool
    :raise ValueError: Invalid length-of-0 kmer_id array produced during k-mer shredding
    :raise AssertionError: relationship 4^k = unique_kmers + (unique)nullomers found invalid
    :returns: (edge_list, header, counts, nullomer_array) header_dictionary is the file's metadata for the header block
    :rtype: (numpy.ndarray, dict, list, list)

    """

    """
    Type Checking

    """

    from kmerdb import kmer
    if type(filepath) is not str:
        raise TypeError("kmerdb.graph.parsefile expects a str as its first positional argument")
    elif not os.path.exists(filepath):
        raise OSError("kmerdb.graph.parsefile could not find the file '{0}' on the filesystem".format(filepath))
    elif type(k) is not int:
        raise TypeError("kmerdb.graph.parsefile expects an int as its second positional argument")
    elif type(b) is not int:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'b' to be an int")
    elif type(both_strands) is not bool:
        raise TypeError("kmerdb.graph.parsefile expects the keyword argument 'both_strands' to be a bool")

    _loggable = logger is not None
    

    N = 4**k
    counts = np.zeros(N, dtype="uint64")


    if _loggable:
        logger.log_it("Initializing edge list fileparser...", "DEBUG")
    seqprsr = parse.SeqParser(filepath, b, k)
    fasta = not seqprsr.fastq # Look inside the seqprsr object for the type of file
    
    # Initialize the kmer class for sequence shredding activities
    Kmer = kmer.Kmers(k) # A wrapper class to shred k-mers with
    
    # Actually load in records 'recs'
    recs = [r for r in seqprsr]


    if _loggable:
        logger.log_it("Read {0} sequences from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"), "INFO")
        
        
    while len(recs):
        num_recs = len(recs)
        #logger.debug("\n\nAcquiring list of all k-mer ids from {0} sequence records...\n\n".format(num_recs))


        list_of_dicts = list(map(Kmer._shred_for_graph, recs))

        if _loggable:
            logger.log_it("k-mer shredding a block of {0} reads/sequences".format(num_recs), "INFO")
        kmer_ids = list(chain.from_iterable(list_of_dicts))

        
        num_kmers = len(kmer_ids)

        if num_kmers == 0:
            raise ValueError("No k-mers to add. Something likely went wrong. Please report to the issue tracker")
        else:
            sys.stderr.write("\nAccumulating all k-mers from this set of records...\n")
            for kmer in kmer_ids:
                counts[kmer] += 1
        

        # END WHILE routine
        # load_more_records, according to 'block size'. Legacy syntax
        recs = [r for r in seqprsr] # The next block of exactly 'b' reads
        # This will be logged redundantly with the sys.stderr.write method calls at line 141 and 166 of parse.py (in the _next_fasta() and _next_fastq() methods)
        #sys.stderr("\n")
        # JUST LOGGING
        if _loggable:
            logger.log_it("Read {0} more records from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"), "DEBUG")




    if _loggable:
        logger.log_it("Read {0} k-mers from the input file".format(len(kmer_ids)), "INFO")

    sys.stderr.write("\n\n\n")
    if _loggable:
        logger.log_it("K-mer counts read from input(s)", "INFO")
    sys.stderr.write("="*40 + "\n")

    if _loggable:
        logger.log_it("Constructing weighted edge list...", "INFO")
    edge_list = make_graph(kmer_ids, k, quiet=quiet, logger=logger)

    if _loggable:
        logger.log_it("Number of edges: {0}".format(len(edge_list)), "INFO")
        logger.log_it("Edge list generation completed...", "INFO")

    unique_kmers = int(np.count_nonzero(counts))
    # FIXME
    num_nullomers = int(len(list(np.where(counts == 0)[0])))


    all_kmer_ids = list(range(N))
    is_nullomer = np.where(counts == 0)
    kmer_ids_np = np.array(all_kmer_ids, dtype="uint64")
    
    nullomer_array = list(kmer_ids_np[is_nullomer])
    
    assert num_nullomers == N - unique_kmers, "kmerdb.graph.parsefile found inconsistencies between two ways of counting nullomers. Internal error."

    
    seqprsr.total_kmers = int(np.sum(counts))
    seqprsr.unique_kmers = unique_kmers
    seqprsr.nullomers = num_nullomers
    seqprsr.nullomer_array = nullomer_array

    
    return (edge_list, seqprsr.header_dict(), counts, seqprsr.nullomer_array)

def make_graph(kmer_ids:list, k:int=None, quiet:bool=True, logger=None):
    """
    Make a neighbor graph from a k-mer id list.

    More specifically, this is a edge list

    Even more so, an edge list of technically adjacent k-mers in a .fa or .fq file.

    Taks from that what you will. It means its a 1D list with implicit, i.e. testable order, and no metadata.

    This particular adjacency list will have redundant relationships:

    i.e. they have been observed again elsewhere in the .fa sequence (1 or multiple) or the .fq file (millions of reads)

    So they are simply tallied.


    :param kmer_ids: list of k-mer ids
    :type kmer_ids: list
    :param k: Choice of k to shred k-mers with
    :type k: int
    :param quiet: Verbosity parameter
    :type quiet: bool
    :raise TypeError: first argument was not a list
    :raise TypeError: keyword argument k was not an inteer
    :raise TypeError: quiet keyword arg should be a bool
    :raise AssertionError: A *suppressed* error indicating whether two subsequent ids in the kmer_ids array were found to be unrelated by their sequences. Suppressed to support multi-sequence .fa/.fq files.
    :returns: An edge list, keyed on a k-mer id 2-tuple e.g. edges[(kmer_id1, kmer_id2)] = weight
    :rtype: dict


    """


    """
    DEBUG

    bad_pair

    A debugging parameter for assessing specific k-mer ids. The debugging process for self-assessing the k-mer pairs and the nature of the non-adjacency is not developed.
    
    3/18/24

    More specifically, i am still encountering the bug and haven't developed the algorithm beyond this point yet. 
    I'm still working on developing logging and processing throughout based on this debug variable set at this point in the variable namespace, and lexical order of the code.
    """

    bad_pair = (3050519, 3908357) # A debugging constant used throughout


    if kmer_ids is None or type(kmer_ids) is not list:
        raise TypeError("kmerdb.graph.make_graph expects a list as its only positional argument")
    elif k is None or type(k) is not int:
        raise TypeError("kmerdb.graph.make_graph expects the keyword argument k to be an int")
    elif quiet is None or type(quiet) is not bool:
        raise TypeError("kmerdb.graph.make_graph expects the keyword argument quiet to be a bool")

    _loggable = logger is not None
    # Initialize total_list of possible k-mer neighbor pairs from the edge list.
    # Step 1. Complete

    # # NIICE.
    # neighbors = None
    # #neighbors = 



    # arrangements[e] 


    # # this is an accumulator for lesser edges (edges where both count/frequency is exceedingly exceedingly rare.) air.air.
    
    # for e in edges_ids:
    #     if (arrangements[e] == 1):
    #         singletons[e] +=1
    #     elif (arrangements[e] == 0):
    #         nullomers[e] += 1
    #     elif (arrangements[e] == max(arrangements.values())):
    #         maxes = max(arrangements.values())

    #     all_neighbors[e] = 0

    # all_all_all_neighbors = None    
    # """
    # Iterate over all k-mer ids, creating lists of neighbors derived from the kmer.
    # Omit the self relationship, as the kmer is not a neighbor of itself.

    # 3/18/24

    # BACKLOG
    # TODO item
    # this isn't complete, as prepopulating the possible edges in k-space, is still in process among other things.
    # Convert to k-mer ids and add the adjacency tuple to the output

    # A minimal neighbor graph is a dictionary, keyed by the pair
    # """

    all_edges_in_kspace = {}


    """
    First we calc the possible neighbor space for the subspace of k-space observed in the file, and all of its 8 neighbors 
    obtained by taking off the leading character, and appending characters to the right of the k-1-mer.

    And vice versa.

    8*num_kmers_observed

    The first loop to populate that subspace p is simple enough that simplification or better than constant time would not be necessary.

    """



    """
    ["Actggtca"] = ["
    """
    neighbors = {}



    for i in range(len(kmer_ids)):

        k_id = kmer_ids[i]

               
        current_kmer = kmer.id_to_kmer(kmer_ids[i], k)
        common_seq = current_kmer[1:-1]

        local_neighbors = kmer.neighbors(current_kmer, k_id, k, quiet=quiet)
        
        neighbors[k_id] = local_neighbors


        # logger.error("             ----------------- 15633431 + 12202077 --------------- END ")


        # logger.error(" |||| =========     + = || = - = || = + = || = - = || + =          =========")

        # 3/18/24
        # LEGACY


        # marked as legacy
        
        # if k_id == 15633431 or k_id == 12202077:
        #     logger.debug("-"*11)
            
        #     logger.debug("K-mer {0} and all (non-self) neighbors:".format(k_id))
        #     logger.debug("{0}  -  {1}".format(current_kmer, "k1 + c"))
        #     logger.debug("-"*11)
        #     logger.debug("<:" + ",".join(char_first) + " |||| " + common_seq + "  ||||  " + ",".join(char_last) + ">:")
        #     logger.debug(",".join(list(map(str, char_first_ids))) + " |||| " + ",".join(list(map(str, char_last_ids))))


        #     logger.debug("Okay, so this is just the k-mer + either char first (left) or char last (right)")
        # for pair_ids in just_eight_edges:

        #     if 15633431 in pair_ids or 12202077 in pair_ids: # ACGT
        #         print(pair_ids)

            
        #     all_edges_in_kspace[pair_ids] = 0

        """
        Old win condition. I wanted to see that a neighbor pair that was originally not found:
        The 12-mers with the following ids, which had this subsequence in common 'GTGGATAACCT', was not found in the keys of the edge list dictionary.
        """

        # missing_key = (15633431, 12202077)
        # if missing_key in all_edges_in_kspace.keys():
        #     print("ITS THERE! LIAR")
        #     print("EDGE OF KSPACE")
        #     print("{0} : {1}".format(missing_key, all_edges_in_kspace[missing_key]))
        #     logger.error(" |||| =========     + = || = - = || = + = || = - = || + =          =========")
        
        # if (15633431, 12202077) in just_eight_edges:
        #     print("Sequence:")
        #     print(current_kmer)
        #     print(just_eight_edges)
            
        #     print("SUCCESS. This k-mer, current k-mer : {0}|{1} has 8 neighbors. {2} . None should be self".format(k_id, current_kmer, just_eight_edges))
        #     #print(all_edges_in_kspace)
        #     sys.exit(1)
        if _loggable:
            logger.log_it("adjacency space completed...", "INFO")
        pass

        
    """
    At this point, the adjacency structure is a simple list of tuples. Redundant tuples may exist in disparate parts of the graph.

    The structure of the neighbor_graph is [(from_node_id, to_node_id), ...]

    Next, the redundant k-mer id pairs are tallied, so that the edge weights are calculated
    """

    error_count = 0

    for i in range(len(kmer_ids)):
        if i+1 == len(kmer_ids):
            break
        else:
            """
            The following constitutes a 'neighboring' pair/edge of two k-mers (nodes):
            the k-mers are sequential in the original k-mer id list (implicitly sorted by virtue of the orientation of the read in the .fasta/.fastq file) 
            this list of k-mers actual neighbors can be validated by inspecting debug level output below (-vv).
            """
            # 3/18/24
            # premature optimization found. not gonna discuss.
            # 3/18/24
            # actually not? the point of the 'adjacency' in the simple list form from the KDBGReader invocation of the kmers structure from the /.fa/ format...
            # so its capturing the raw kmer array (of ids) from .fa and/or .fq data by virtue of some existing features
            # and it's collecting the proper ordering of the k-mer relationship integrity across the lines of the .fq
            # pair = (the current k-mer and its nearest neighbor in the sequence file i.e. whether or not it was adjacenct in the 'sequence' proper or whether the k-mer is in a list with multiple sequences, as in given from a .fq file, such that adjacent *reads* in the file produce a k-mer input list that has k-mer ids in a list that represents residue positions in both a .fa and a .fq file. Implicit order is understood by the programmer using this, but the computer itself still has no knowledge of whether a k-mer (or residue) in the file is adjacent in the genome coordinates to its literal neihbor in the list, vs on a subsequent chromosome, or Illumina read in the case of .fq input.
            """
            Basically at this point we haven't guaranteed anything about the input k-mer list. It's just been parsed and passed in from either a .fa or .fq file.
            Okay so now we're creating an artificial pair of k-mers, a relationship between two ids, that also constitutes a k-1-mer in common, regardless of the orientation of the 'neighbor' relationship. Additional assert statements below

            Essentially, this guarantees that, by having the k-2-mer in common, that the current_kmer and any of its 8 possible neighbors have at least a k-2mer in common among all types of neighbor relaionships.
            In addition, the 'neighbor' relationship between the current_kmer and a single neighbor, also has a k-1-mer in common.
            Initialize the k-mer pair (duple of ids) in the all_edges_in_kspace var, capturing the pair_id space in the keys iterable from the .keys()
            """

            current_kmer_id = kmer_ids[i]
            current_kmer = kmer.id_to_kmer(current_kmer_id, k)
            pair = (current_kmer, kmer.id_to_kmer(kmer_ids[i+1], k)) 
            pair_ids = tuple(map(kmer.kmer_to_id, pair))
            
            common_seq = current_kmer[1:]


            if quiet is False:
                sys.stderr.write("Kmer_id: = kmer_ids[i]     =        {0}\n".format(current_kmer_id))
                sys.stderr.write(" -----pair: = '{0}'  <-> '{1}'\n".format(pair[0], pair[1]))
            all_edges_in_kspace[pair_ids] = 0
            # Create the neighbor lists from the list of chars, prepended or appended to the k-mer suffix/prefix
            """
            BOGUS ASSERTIONS AND NONGUARANTEES...
            """
            assert len(current_kmer) == k
            assert len(common_seq) == k-1
            
            #logger.debug("Validations completed for k-mer {0}.".format(current_kmer_id))
            k1 = current_kmer[1:] # 11 characters (k - 1) remove the left most char, so new char goes on end to make the neighbor
            # not_the_king_but_the_rook = None
            k2 = current_kmer[:-1] # 11 characters ( k - 1) remove the final character so new char is prepended
            # Do we even assert these?

            # Comments:
            # The char_last requires its first char deleted.
            # Then a character is prepended to the beginning to the reverse of the sequence, or appending to the forward sequence.
            # The char_first requires its terminal char deleted, and then 4 residue characters are prepended to the forward sequence, or appended to the childless reverse
            # we call char_first 'childless' because it is assigned responsibility for the 4 interchangables on its prepend, while it's oldest parent is removed.
            """
            current_kmer
            pair : the pair of adjacent k-mers with the k-1mer as identical sequence
            pair_ids  : ids
            k1 : k-1mer   current_kmer[1:] left-most char ommitted, (neighbors *will* get appended to right side)
            k2 : k-1mer   current_kmer[:-1] right-most char ommitted, neighbors prepended

            # e.g. ACTGACTG
            # pair0 = CTGACTG (k-1 mer via first char removed. k-1 mer that receives appends)
            # pair1 = ACTGACT (k-1 mer via last char removed. k-1 mer that gets at its prepend.)

            DECLARATIONS FINISHED
            """
            # Commented 3/11/24 ish thanks nic
            # Working on some bills, then back to the drawing board.
            # I mean, aside from a lab job, if i could wfh or work localish then taking a month or so off to go to the mountains for soil sampling.
            # I'd really want to go with a guide or expert.
            # That way I could sample the most diverse ecosystems for sampling
            # But prior to that I'd need the go-kit prototype.
            # Maybe a soil sampling protocol, expert backpack, sampling fluids/solvents, portable lab cooler.
            # Then run the samples off of some grant.
            """
            LOCAL

            B A D         P A I R S
            debuging problematic pair of k-mer ids
            """
            if pair_ids == bad_pair:
                sys.stderr.write("Problematic pair of k-mer ids identified...\n")
                sys.stderr.write("Pair: ({0})\n".format(", ".join(pair_ids)))
                raise ValueError("DEBUGGING PURPOSES ONLY: Identified a pair of k-mers noted as problematic. If you see this error, file an issue at https://github.com/MatthewRalston/kmerdb")

            try:
                """
                This is real concise. performance focused, not too much style. 
                but this is the best way in python to concisely assert that the string slice is equivalent from k-mer to k-mer neighbor. Ensuring this throughout is damn routine.
                essentially, this code block is added as indication that these assertions hold on a single sequence .fa file
                by ensuring that sequential k-mer ids retrieved (via sliding window) from a single sequence .fa are, in fact, related by a common sequence.
                when the assertion error is triggered, it is most likely that this is due to subsequent reads in a .fastq file, resulting in an erroneous 'pair' of k-mers from the algorithm above, targeting the n'th and n+1'th k-mers in the kmer_id array.
                """
                assert pair[0][1:] == pair[1][:-1], "kmerdb.graph expects neighboring k-mers to have at least k-1 residues in common. Must be from a fastq file."
                assert (common_seq in pair[0]) and (common_seq in pair[1]), "kmerdb.graph expects the common_seq be in both sequences of the pair."
            except AssertionError as e:
                error_count += 1
                # logger.error(pair)
                # logger.warning("This k-mer pair should have this subsequence in common: {0}".format(common_seq))
                # logger.warning("k1: '{0}'           k2: '{1}'".format(k1, k2))
                # logger.warning("It is explicitly assumed this is from .fastq input")
                # #logger.debug("Types: pair[0] : {0} pair[1] : {1}".format(type(pair[0]), type(pair[1])))
                # logger.debug("'{0}' == '{1}'".format(pair[0][1:], pair[1][:-1]))
                # logger.debug("Is pair0[1:] == pair1[:-1]? : {0}".format(pair[0][1:] == pair[1][:-1]))
                # logger.debug(e)
                pass
            try:
                # Accumulate for each edge pair
                #if (15633431, 12202077) == pair_ids:
                if bad_pair == pair_ids:
                    sys.stderr.write("\n\nin accumulator... cannot identify the problem for this pairing again. It is in the wrong order in the key submission, and can't be recovered for some reason\n\n")
                    sys.stderr.write(bad_pair)

                    raise ValueError("Internal Error. Related to debugging problematic pairs during .fastq inter-read pair removal from edge list.")
                """
                ACCUMULATOR           --------- DONE
                """
                all_edges_in_kspace[pair_ids] += 1
                """

                You might be thinking "HOLY CRAP THATS A LOT OF MESSAGES" without the --quiet flag.

                That's right. Also:


                "Well, why can't you just fix the parser so it *only* returns valid edges from the .fq file?"

                i.e. check whether each k-mer... from each read... is at the final index of the read... right?


                So instead, I iterate over the k-mer array one time, and remove the edges. The performance difference between the approaches
                would be negligible. 

                """

                
                if quiet is False:
                    sys.stderr.write("Edge: {0} => {1}".format(pair[0], pair[1]))
#                     sys.stderr.write("""
#  ╱|、
# (˚ˎ 。7  
#  |、˜〵          
# じしˍ,)ノ
# \n""")                        
            except KeyError as e:
                sys.stderr.write("Invalid key(pair) used to access edge list\n")
                sys.stderr.write("PAIR: {0} => {1}\n".format(pair[0], pair[1]))
                sys.stderr.write("pair ids: {0}\n".format(pair_ids))
                sys.stderr.write(30*"=" + "\n")
                # for i, p in enumerate(all_edges_in_kspace.keys()):
                #     if bad_pair == pair_ids:
                #         # Double checking
                #         logger.error("Debugging the 'bad pair', problematic k-mer id pair...")
                #         logger.error(pair_ids)
                raise e
            # print("heeeeeey, wow thens doo the prin funshen")
            # print(pair_ids, pair, "wow thans prin funshun ( ---------)_______________________prin fun-shun_______________________________________(    --------------)")
            # print("print funshun you're so smart and tall")
            # print("print funshun why did you sit next to me")
            # print("function")


            maxlen = 15 # Placeholder
            p1str = str(pair_ids[0])
            p2str = str(pair_ids[1])
            p1str = p1str + (maxlen - len(p1str))*" "
            p2str = p2str + (maxlen - len(p2str))*" "

            if quiet is False:
                sys.stderr.write("k-mer 'pair' ({0}, {1}) adjacency from input file(s) verified".format(p1str, p2str))
                sys.stderr.write("\r")
            pass


    if error_count > 0:
        sys.stderr.write("\n\n\nNOTE: ADJACENCY ERRORS DETECTED: Found {0} 'improper' k-mer pairs/adjacencies from the input file(s),\n where subsequent k-mers in the k-mer id array (produced from the sliding window method over input seqs/reads) did not share k-1 residues in common.\n These *may* be introduced in the array from the last k-mer of one seq/read (in the .fa/.fq) and the first k-mer of the next seq/read.\n".format(error_count))
    sys.stderr.write("\n\n\n")
    sys.stderr.write("All edges populated *and* accumulated across k-mer id arrays from inputs.\n")
    sys.stderr.write("\n\n\n'Graph' generation complete\n\n")

            
    # This variable is unrestricted. A larger k may inflate this var beyond the memory limits of the computer where the observed k-space can be represented with enough resolution through enough profile diversity.
    """
    node_pairs and the k-space diversity 


    Depending on inputs, this may become large.
    """
    return all_edges_in_kspace



def create_graph(nodes:list, edge_tuples:list, gpu:bool=False):

    import networkx as nx

    if nodes is None or type(nodes) is not list or not all(type(n) is not int for n in nodes):
        raise TypeError("kmerdb.graph.create_graph expects the first argument to be a list of ints")
    elif edge_tuples is None or type(edge_tuples) is not list or not all(len(e) != 3 for e in edge_tuples) or not all((type(e[0]) is int and type(e[1]) is int) for e in edge_tuples):
        raise TypeError("kmerdb.graph.create_graph expects the second argument to be a list of tuples of length 2")
    elif gpu is None or type(gpu) is not bool:
        raise TypeError("kmerdb.graph.create_graph expects the keyword argument gpu to be a bool")


    """
    Now we make the networkx graph
    """
    G = nx.Graph()

    G.add_nodes_from(nodes)
    
        
    G.add_edges_from(edge_tuples)


    

    

def w_lexer():
    pass

            
def v_order_lexer(e:tuple, asc:bool=False):
    """
    Uh.
    """
    if type(asc) is not bool:
        raise TypeError("kmerdb.graph requires valid bool keyword arg 'asc'")
    
    if type(e) is not tuple:
        logger.error("kmerd.graph.threetuple_v_order_lexer needs a valid tuple")
        
        sys.exit(1)
        exit(1)
    elif len(e) != 2:
        raise IndexError("kmerdb.graph.v_order_lexer requires length 2")
    elif asc is not bool:
        logger.error("kmerdb.graph.v_order_lexer accepts the sort order as an input")
        raise TypeError("kmerdb.graph requires a lexical order")
    elif asc is False:
        raise RuntimeError("requires the order")

    
    if len(e) != 3:
        logger.error("kmerdb.graph invalid input")
        raise TypeError("kmerdb.graph requires a 3-tuple")

    n1, n2, w = e

    # Maybe...
    n1 = int(n1)
    n2 = int(n2)

    # This ... smh
    if asc is True:
        if n1 > n2:
            return (n2, n1, w)
        elif n2 > n1:
            return (n1, n2, w)
        else:
            raise RuntimeError('x')
    elif asc is False:
        if n1 > n2:
            return (n1, n2, w)
        elif n2 > n1:
            return (n2, n1, w)
        else:
            raise RuntimeError('x')
    else:
        
        logger.error("idgaf")
        return (n1, n2, w)


def randomsort_lexer(g):
    return (random.randint(0, 1428945), random.randint(0, 1428945), random.randint(0, 1428945))

def wrupsort_lexer(g):
    return sorted([g[0], g[1]])


    
def wrongsort_lexer(g, asc:bool=False):
    if asc is True:
        if n1 > n2:
            return (n1, n2, w)
        elif n2 > n1:
            return (n2, n1, w)
        else:
            raise RuntimeError('x')
    elif asc is False:
        if n1 > n2:
            return (n2, n1, w)
        elif n2 > n1:
            return (n1, n2, w)
        else:
            raise RuntimeError('x')




class KDBGReader(bgzf.BgzfReader):
    """
    A class to read .kdbg files.




    1 args

    8 keyword args

    fileobj
    mode
    n1_dtype
    n2_dtype
    weights_dtype
    sort
    slurp

    :ivar filename: The choice of k to shred with
    :ivar fileobj: An existing fileobject from io.IOBase
    :ivar mode: read/write mode
    :ivar n1_dtype: NumPy dtype
    :ivar n2_dtype: NumPy dtype
    :ivar weights_dtype: NumPy dtype
    :ivar sort: *unimplemented* - uh.. sort the .kdbg edges somehow
    :ivar slurp: Autoload the .kdbg file

    
    """
    def __init__(self, filename:str, fileobj:io.IOBase=None, mode:str="r", n1_dtype:str="uint64", n2_dtype:str="uint64", weights_dtype:str="uint64", sort:bool=False, slurp:bool=False, logger=None):
        """
        A wrapper around Bio.bgzf.BgzfReader

        :param filename: A valid filepath
        :type filename: str
        :param fileobj: An existing fileobject from io.IOBase
        :type fileobj: io.IOBase
        :param mode: read/write mode
        :type mode: str
        :param n1_dtype: NumPy dtype
        :type n1_dtype: str
        :param n2_dtype: NumPy dtype
        :type n2_dtype: str
        :param weights_dtype: NumPy dtype
        :type weights_dtype: str
        :param sort: *unimplemented* - uh.. sort the .kdbg edges somehow
        :type sort: bool
        :param slurp: Autoload the .kdbg file
        :type slurp: bool
        :raise TypeError: fileobj was not a file-type interface, inheriting from io.IOBase
        :raise TypeError: filename was not a str
        :

        """
        """

         fileobj type checking
        """


        if fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'fileobj' to be a file object")
        if filename is not None and type(filename) is not str:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'filename' to be a str")
        elif not os.access(filename, os.R_OK):
            raise ValueError("kmerdb.graph.open expects an existing filepath as its first positional argument")

        elif mode is not None and type(mode) is not str:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'mode' to be a str")
        elif n1_dtype is not None and type(n1_dtype) is not str:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'n1_dtype' to be a str")
        elif n2_dtype is not None and type(n2_dtype) is not str:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'n2_dtype' to be a str")
        elif weights_dtype is not None and type(weights_dtype) is not str:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'weights_dtype' to be a str")

        elif sort is not None and type(sort) is not bool:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'sort' to be a bool")
        elif slurp is not None and type(slurp) is not bool:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'slurp' to be a bool")

        # if logger is None:
        #     raise ValueError("graph.KDBGReader expects the keyword argument 'logger' to be valid")
        
        """

        Handle fileobj or 
        """
        
        if fileobj:
            assert filename is None
            handle = fileobj
            assert "b" in handle.mode.lower()
        else:
            if "w" in mode.lower() or "a" in mode.lower():
                raise ValueError("Must use read mode (default), not write or append mode")
            handle = _open(filename, "rb")

        self._text = "b" not in mode.lower()
        if self._text:
            self._newline = "\n"
        else:
            self._newline = b"\n"
        self._handle       = handle
        self._filepath     = self._handle.name

        #       Placeholder
        self.max_cache     = 100
        
        self._buffers      = {}
        self._block_start_offset = None
        self._block_raw_length = None

        self.logger = logger
        self._loggable = logger is not None
        """

        np.array defs

        n1
        n2
        weights

        """
        
        self.n1         = None
        self.n1_dtype   = None
        self.n2         = None
        self.n2_dtype   = None
        self.weights       = None
        self.weights_dtype = None

        self.read_and_validate_kdbg_header()

        if slurp is True:
            if self._loggable:
                self.logger.log_it("Reading .kdbg contents into memory", "INFO")
            self.slurp(n1_dtype=n1_dtype, n2_dtype=n2_dtype, weights_dtype=weights_dtype)

        self.is_int = True
        if handle is not None:
            handle.close()
        if self._handle is not None:
            self._handle.close()
        self._handle = None
        handle = None
        fileobj=None
        return

    def read_and_validate_kdbg_header(self):

        """

        KDBGReader PyYAML                metadata ingestion

        """

        '''
        Here we want to load the metadata blocks. We want to load the first two lines of the file: the first line is the version, followed by the number of metadata blocks
        '''
        # 0th block
        if self._loggable:
            self.logger.log_it("Loading the 0th block from '{0}'...".format(self._filepath), "INFO")
        self._load_block(self._handle.tell())

        self._buffer = self._buffer.rstrip(config.header_delimiter)
        
        initial_header_data = OrderedDict(yaml.safe_load(self._buffer))

        # Placeholder
        num_header_blocks = None

        if type(initial_header_data) is str:
            if self._loggable:
                self.logger.log_it("Inappropriate type for the header data.", "ERROR")
            #logger.info("Um, what the heckin' is this in my metadata block?")
            raise TypeError("kmerdb.graph.KDBGReader could not parse the YAML formatted metadata in the first blocks of the file")
        elif type(initial_header_data) is OrderedDict and self._loggable:
            self.logger.log_it("Successfully parsed the 0th block of the file, which is expected to be the first block of YAML formatted metadata", "INFO")
            self.logger.log_it("Assuming YAML blocks until delimiter reached.", "INFO")

        """

        KDBGReader
        
        YAML metadata spec validation
        """
            
        if "version" not in initial_header_data.keys():
            raise TypeError("kmerdb.graph.KDBGReader couldn't validate the header YAML")
        elif "metadata_blocks" not in initial_header_data.keys():
            raise TypeError("kmerdb.graph.KDBGReader couldn't validate the header YAML")

        if self._loggable:
            self.logger.log_it(initial_header_data, "INFO")
        
        if initial_header_data["metadata_blocks"] != 1:
            if self._loggable:
                self.logger.log_it("More than 1 metadata block: uhhhhh are we loading any additional blocks", "ERROR")
            raise IOError("Cannot read more than 1 metadata block yet")

        for i in range(initial_header_data["metadata_blocks"] - 1):
            if self._loggable:
                self.logger.log_it("Multiple metadata blocks read, most likely from a composite edge-graph...", "WARNING")
            self._load_block(self._handle.tell())
            addtl_header_data = yaml.safe_load(self._buffer.rstrip(config.header_delimiter))
            if type(addtl_header_data) is str:
                if self._loggable:
                    self.logger.log_it(str(addtl_header_data), "ERROR")
                raise TypeError("kmerdb.graph.KDBGReader determined the data in the {0} block of the header data from '{1}' was not YAML formatted".format(i, self._filepath))
            elif type(addtl_header_data) is dict:
                sys.stderr.write("\r")
                sys.stderr.write("Successfully parsed {0} blocks of YAML formatted metadata".format(i))
                initial_header_data.update(addtl_header_data)
                num_header_blocks = i
            else:
                if self._loggable:
                    self.logger.log_it(str(addtl_header_data), "ERROR")
                raise RuntimeError("kmerdb.graph.KDBGReader encountered a addtl_header_data type that wasn't expected when parsing the {0} block from the .kdb file '{1}'.".format(i, self._filepath))

        #raise RuntimeError("kmerdb.graph.KDBGReader encountered an unexpected type for the header_dict read from the .kdb header blocks")

        if self._loggable:
            self.logger.log_it("Validating the header YAML...", "INFO")

        try:
            jsonschema.validate(instance=initial_header_data, schema=config.graph_schema)
            self.metadata = dict(initial_header_data)
            
            self.k = self.metadata['k']
            self.n1_dtype = self.metadata["n1_dtype"]
            self.n2_dtype = self.metadata["n2_dtype"]
            self.weights_dtype = self.metadata["weights_dtype"]
            self.sorted = self.metadata["sorted"]

        except jsonschema.ValidationError as e:
            if self._loggable:
                self.logger.log_it("kmerdb.graph.KDBGReader couldn't validate the header/metadata YAML from {0} header blocks".format(num_header_blocks), "ERROR")

                self.logger.log_it("""


Failed to validate YAML header.

-------------------------------



 You can store unique k-mer counts, total nullomer counts, and other metadata alongside the weighted edge list with the 'kmerdb graph' command


 ...then try again to view the header with 'kmerdb header'



""", "ERROR")

                
                self.logger.log_it(e.__str__, "ERROR")
            raise ValueError("Requires kmerdb v{0} format YAML header. Body is .tsv format table, .bgzf file.       - weighted edge list        (idx, node1, node2, weight)")
        self.metadata["header_offset"] = self._handle.tell()
        #logger.debug("Handle set to {0} after reading header, saving as handle offset".format(self.metadata["header_offset"]))
        #self._reader = gzip.open(self._filepath, 'r')
        self._offsets = deque()
        for values in bgzf.BgzfBlocks(self._handle):
            #logger.debug("Raw start %i, raw length %i, data start %i, data length %i" % values)
            self._offsets.appendleft(values) # raw start, raw length, data start, data length
        if len(self._offsets) == 0:
            raise IOError("kmerdb.graph.KDBGReader opened an empty file")
        # Skip the zeroth block
        self._load_block()
        # print(str(self._buffer)) # 1
        # print(self.readline())
        # self._load_block()
        # print(self._buffer) # 2
        # print(self.readline())


    def read_line(self):
        """
        Read and parse a single line from the .kdbg file

        :returns: node1_id, node2_id, weight
        :rtype: tuple
        
        """
        
        line = self.readline()


        if self.n1_dtype == "uint64" and self.n2_dtype == "uint64" and self.weights_dtype == "uint64":
            return parse_kdbg_table_line(line, row_dtype="uint64")
        else:
            raise IOError("kmerdb.graph.KDBGReader.read_line Cannot determine file dtype")



    def slurp(self, n1_dtype:str="uint64", n2_dtype:str="uint64", weights_dtype:str="uint64"):
        """
        Autoload the .kdbg file into memory

        :param n1_dtype: NumPy dtype
        :type str:
        :param n2_dtype: NumPy dtype
        :type str:
        :param weights_dtype: NumPy dtype
        :type str:

        """
        
        if type(n1_dtype) is not str:
            raise TypeError("kmerdb.graph.KDBGReader.slurp expects n1_dtype to be a str")
        if type(n2_dtype) is not str:
            raise TypeError("kmerdb.graph.KDBGReader.slurp expects n2_dtype to be a str")
        if type(weights_dtype) is not str:
            raise TypeError("kmerdb.graph.KDBGReader.slurp expects weights_dtype to be a str")

        try:
            np.dtype(n1_dtype)
            np.dtype(n2_dtype)
            np.dtype(weights_dtype)
        except TypeError as e:
            raise TypeError("kmerdb.graph.KDBGReader.slurp expects each dtype keyword argument to be a valid numpy data type")

        vmem = psutil.virtual_memory()

        i = 0

        idx = []
        arrs = []
        node1 = []
        node2 = []
        weights = []
        reading=True
        while reading is True:

            try:
                line = next(self)

            except StopIteration as e:
                if self._loggable:
                    self.logger.log_it("Finished loading .kdbg through slurp (on init)", "ERROR")
                #raise e
                pass
            if line is None:
                if self._loggable:
                    self.logger.log_it("'next' returned None. Panic", "WARNING")

                raise RuntimeError("kmerdb.graph.KDBGReader.slurp Panicked on new line")


                  
            try:
                
                arr = self.read_line()
                arrs.append(arr)

                
                i, n1, n2, weight = arr

                idx.append(i)
                node1.append(n1)
                node2.append(n2)
                weights.append(weight)
            except StopIteration as e:
                reading = False

        
        self.n1 = np.array(node1, dtype=n1_dtype)
        self.n2 = np.array(node2, dtype=n2_dtype)
        self.weights = np.array(weights, dtype=weights_dtype)
        return



    
class KDBGWriter(bgzf.BgzfWriter):
    """
    A wrapper class around Bio.bgzf.BgzfWriter to write a .kdbg file to disk.


    :ivar filename: valid filepath
    :ivar mode: read/write mode
    :ivar metadata: header information
    :ivar both_strands: *unimplemented*
    :ivar fileobj: io.IOBase
    :ivar compresslevel: compression parameter
    """
    
    def __init__(self, filename:str=None, mode:str="w", metadata:OrderedDict=None, both_strands:bool=False, fileobj:io.IOBase=None, compresslevel:int=6, logger=None):
        """
        A wrapper around Bio.bgzf.BgzfWriter

        :param filename: A valid filepath
        :type filename: str
        :param mode: read/write mode
        :type mode: str
        :param metadata: A metadata header for the .kdbg file
        :type metadata: collections.OrderedDict
        :param both_strands: *unimplemented*
        :type both_strands: bool
        :param compresslevel: compression parameter
        :type compresslevel: int
        :raise TypeError: filename was not a str
        :raise TypeError: mode was invalid
        :raise TypeError: both_strands was invalid
        :raise TypeError: fileobj was invalid
        :raise TypeError: compresslevel was invalid
        :raise ValueError: mode was invalid
        :raise NotImplementedError: append mode was invalid
        """
        
        if filename is None or type(filename) is not str:
            raise TypeError("kmerdb.graph.KDBGWriter expects the filename to be a str")
        elif mode is None or type(mode) is not str:
            raise TypeError("kmerdb.graph.KDBGWriter expects the mode to be a str")
        elif metadata is None or (type(metadata) is not OrderedDict and type(metadata) is not dict):
            raise TypeError("kmerdb.graph.KDBGWriter - invalid metadata argument")
        elif both_strands is None or type(both_strands) is not bool:
            raise TypeError("kmerdb.graph.KDBGWriter expects the keyword argument 'both_strands' to be a bool")
        elif fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kmerdb.graph.KDBGWriter expects the keyword argument 'fileobj' to be an instance of io.IOBase")
        elif compresslevel is None or type(compresslevel) is not int:
            raise TypeError("kmerdb.graph.KDBGWriter expects the keyword argument 'compresslevel' to be an int")

        self.logger = logger
        self._loggable = logger is not None

        
        if fileobj:
            assert filename is None, "kmerdb.graph expects filename to be None is fileobj handle is provided"
            handle = fileobj
        else:
            if "w" not in mode.lower() and "a" not in mode.lower():
                raise ValueError("Must use write or append mode, not %r" % mode)
            elif "wb" == mode:
                pass
            elif mode == "w":
                pass
            elif "a" in mode.lower():
                raise NotImplementedError("Append mode is not implemented yet")
                # handle = _open(filename, "ab")
            else:
                raise RuntimeError("Unknown mode for .kdbg file writing class kmerdb.graph.KDBGWriter: {0}".format(mode))
        self._text = "b" not in mode.lower()
        self._handle = _open(filename, "wb")
        self._buffer = b"" if "b" in mode.lower() else ""
        self.compresslevel = compresslevel


        self.logger = logger
        self._loggable = logger is not None
        
        """
        Write the header to the file
        """

        if self._loggable:
            self.logger.log_it("Constructing a new .kdbg file '{0}'...".format(self._handle.name), "INFO")


        # 3-04-2024
        yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)

        self.metadata = metadata

        #self._write_block(metadata_slice)
        if "b" in mode.lower():
            metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
            metadata_plus_delimiter_in_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')
            self.metadata["metadata_blocks"] = math.ceil( sys.getsizeof(metadata_plus_delimiter_in_bytes) / ( 2**16 ) ) # First estimate
            metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
            metadata_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')
            self.metadata["metadata_blocks"] = math.ceil( sys.getsizeof(metadata_bytes) / ( 2**16 ) ) # Second estimate
            metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
            metadata_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')
            if self._loggable:
                self.logger.log_it("Writing the {0} metadata blocks to the new file".format(self.metadata["metadata_blocks"]), "INFO")
                self.logger.log_it("Header is being written as follows:\n{0}".format(yaml.dump(self.metadata, sort_keys=False)), "DEBUG")

            # 01-01-2022 This is still not a completely functional method to write data to bgzf through the Bio.bgzf.BgzfWriter class included in BioPython
            # I've needed to implement a basic block_writer, maintaining compatibility with the Biopython bgzf submodule.
            #self.write(bytes(yaml.dump(metadata, sort_keys=False), 'utf-8'))
        
            for i in range(self.metadata["metadata_blocks"]):
                metadata_slice = metadata_bytes[:65536]
                metadata_bytes = metadata_bytes[65536:]
                self._write_block(metadata_slice)

                #self._write_block
                self._buffer = b""
                self._handle.flush()
        elif "w" == mode.lower() or "x" == mode.lower():
            metadata_yaml_str = yaml.dump(metadata, sort_keys=False)
            metadata_bytes = bytes(metadata_yaml_str, 'utf-8')
            mode = "b"
            self._buffer = b""
            self.write(metadata_bytes)

            self._handle.flush()
        else:
            raise RuntimeError("Could not determine proper encoding for write operations to .kdb file")
        

    

        

class Parseable:
    def __init__(self, arguments, logger):
        self.arguments = arguments
        self.logger = logger
        
    def parsefile(self, filename):
        """Wrapper function for graph.parsefile to keep arguments succinct for deployment through multiprocessing.Pool
            
        :param filename: the filepath of the fasta(.gz)/fastq(.gz) to process with kmerdb.graph.parsefile
        :type filename: str
        """
        return parsefile(filename, self.arguments.k, quiet=self.arguments.quiet, logger=self.logger)



