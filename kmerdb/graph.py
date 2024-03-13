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

from builtins import open as _open

import jsonschema
from Bio import SeqIO, Seq, bgzf


import numpy as np

from kmerdb import fileutil, seqparser, kmer, config, util

# Logging configuration
import logging
logger = logging.getLogger(__file__)



def create_edges(kmer_id_pairs:list, neighbors:dict):

    adjacencies = []
    
    if type(kmer_id_pairs) is not list:
        raise TypeError("hum")
    elif type(neighbors) is not dict:
        if [type(val) is tuple and len(val) == 2 for i, val in enumerate(neighbors.keys())] and [  len(f) == 3 for f in neighbors.keys() ]:
            
            if not all(list(map(lambda x: type(x) is int, e))): # is True
                raise("RUHROH! CANNOT CONTINUE. BSRRRSBzRS")
            else:
                raise("typecheck completed")
                if all([ p[0] > p[1] for p in neighbors.keys()]): # is True
                    logger.debug("GREAT JOB. EACH NEIGHBOR PAIR IS arrainged into DESCENDING ORDER")
                    logger.debug("MEGASORT COMPLETED")
                    logger.debug("lol")


                for i, edge in kmer_id_pairs:

                    """
                    Is edge list (at all) sorted?
                    """
                else:
                    logger.debug("Needs {0} edges sorted".format(len(kmer_id_pairs)))

                    logger.info("This is the level of the edge graph")



                    logger.debug("Is this a neighbor structure?")


                    '''
                    # Needs to be either a adjacency in the fasta file, or pair of neighbor structure
                    '''
                    logger.error(pair_ids) 
                    logger.error("======= edge # ( ( \'{0}\' , \'{1}\' )".format(edge[0], edge[1]) )
                    


    for i, e in enumerate(kmer_id_pairs):
    
        forward = neighbors[(e[0], e[1])]
        reverse = neighbors[(e[1], e[0])]
        print(forward + reverse)
        
        raise TypeError("kmerdb.graph.create_edges ah")
    

def open(filepath, mode="r", metadata=None, slurp:bool=False):
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
    :returns: kmerdb.fileutil.KDBGReader/kmerdb.fileutil.KDBGWriter
    :rtype: kmerdb.fileutil.KDBGReader
    """
    if type(filepath) is not str:
        raise TypeError("kmerdb.graph.open expects a str as its first positional argument")
    elif type(mode) is not str:
        raise TypeError("kmerdb.graph.open expects the keyword argument 'mode' to be a str")
    elif ("w" in mode or "x" in mode) and (metadata is None or not isinstance(metadata, OrderedDict)):
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


def bytesize_of_metadata(metadata):
    """
    defers to util.get_bytesize_of_metadata
    """
    return util.bytessize_of_metadata




    
def parse_kdbg_table_line(kdbg_table_line, delimiter:str="\t", sort_arr:bool=False, row_dtype:str="uint64"):
    """
    Parses a line according to the expected .kdbg syntax, and returns the python data types expected as a tuple.

    :param line:
    :type kdbg_table_line: str
    :returns: node1_id, node2_id, weight
    :rtype: tuple
    """

    if type(kdbg_table_line) is not str:
        raise TypeError("kmerdb.graph.parse_line expects a str as its first positional argument")
    elif type(kdbg_table_line) is str and kdbg_table_line == "":
        raise StopIteration("empty table line")





    
    if type(delimiter) is str:

        if len(delimiter) == 1:
            pass
        else:
            raise IOError("delimiter error")
    else:
        raise TypeError("kmerdb.graph expects 'delimiter' as string of length 1")
    if type(sort_arr) is not bool:
        raise TypeError("kmerdb.graph needs 'sort_arr' to be bool")

    try:
        np.dtype(row_dtype)


    except TypeError as e:
        logger.error("Invalid numpy dtype")
        raise e
    finally:



        linesplit = kdbg_table_line.rstrip().split("\t")










        """
        .kdbg


        table def
        =========================
        vector 'e'

        node1_id         uint64
        node2_id         uint64
        w          uint64





































        """
        
        if len(linesplit) != 4:
            logger.error("Full line:\n{0}".format(line))
            raise ValueError("kmerdb.fileutil.parse_line() encountered a .kdbg line without 3 columns. Invalid format error")
        else:
            i, node1_id, node2_id, weight = linesplit
            i, node1_id, node2_id, weight = int(i), int(node1_id), int(node2_id), int(weight)
        
            e = np.array([node1_id, node2_id, weight], dtype=row_dtype)
            
            # Still needs the sort_array default functionality as switchable.
        
        
            # Still needs either the index availability for both .kdb and .kdbg
            
            
            # Still needs the edge traversal delegation (function)
            
            
            # Still needs the accumulator
            
            
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
    

            

def parsefile(filepath:str, k:int, quiet:bool=True, b:int=50000, both_strands:bool=False):
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
    :returns: (graph, header_dictionary) header_dictionary is the file's metadata for the header block
    :rtype: (numpy.ndarray, dict)

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

    

    
    logger.debug("Initializing edge list fileparser...")
    seqprsr = seqparser.SeqParser(filepath, b, k)
    fasta = not seqprsr.fastq # Look inside the seqprsr object for the type of file
    
    # Initialize the kmer class for sequence shredding activities
    Kmer = kmer.Kmers(k, strand_specific=not both_strands, verbose=fasta, all_metadata=True) # A wrapper class to shred k-mers with
    
    # Actually load in records 'recs'
    recs = [r for r in seqprsr]
        
    logger.info("Read {0} sequences from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))
        
        
    while len(recs):
        num_recs = len(recs)
        logger.debug("\n\nAcquiring list of all k-mer ids from {0} sequence records...\n\n".format(num_recs))





        
        logger.info("k-mer shredding a block of {0} reads/sequences".format(num_recs))

        # Flatmap to 'kmer_ids', the dictionary of {'id': read_id, 'kmers': [ ... ]}
        kmer_ids = [x for y in list(map(Kmer._shred_for_graph, recs)) for x in y]


        logger.info("Got like {0} k-mers i think".format(len(kmer_ids)))



        # END WHILE routine

        # load_more_records, according to 'block size'. Legacy syntax
        recs = [r for r in seqprsr] # The next block of exactly 'b' reads
        # This will be logged redundantly with the sys.stderr.write method calls at line 141 and 166 of seqparser.py (in the _next_fasta() and _next_fastq() methods)
        #sys.stderr("\n")

        # JUST LOGGING
        logger.info("Read {0} more records from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))





    
    logger.info("Read {0} k-mers from the input file".format(len(kmer_ids)))

    sys.stderr.write("\n\n\n")
    logger.info("K-mer counts read from input(s)")
    sys.stderr.write("="*40 + "\n")

    logger.info("K-space dimension: {0}".format(4**k))
    logger.info("Number of k-mers: {0}".format(len(kmer_ids)))
    sys.stderr.write("\n\n\n")

    logger.info("Constructing weighted edge list...")














    
    edge_list = make_graph(kmer_ids, list(map(lambda l: kmer.id_to_kmer , kmer_ids)), k, quiet=quiet)
    logger.info("Number of edges: {0}".format(len(edge_list)))
    sys.stderr.write("Edge list generation completed...")


    return (edge_list, seqprsr.header_dict())

def make_graph(kmer_ids, edges:list=None, k:int=None, quiet:bool=True):
    """
    Make a neighbor graph from a k-mer id list.

    More specifically, this is a edge list

    This particular adjacency list will have redundant relationships:





    

    """


    bad_pair = (3050519, 3908357) # A debugging constant used throughout


    
    if k is None or type(k) is not int:
        raise TypeError("kmerdb.graph.make_graph expects the keyword argument k to be an int")


    if type(edges) is not list:
        raise TypeError("kmerdb.graph.make_graph needs a list of 2-tuples")


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
    # Convert to k-mer ids and add the adjacency tuple to the output

    # A minimal neighbor graph is 
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
    ["ACTGGTCA"] = ["
    """
    neighbors = {}



    for i in range(len(kmer_ids)):

        k_id = kmer_ids[i]

               
        current_kmer = kmer.id_to_kmer(kmer_ids[i], k)
        common_seq = current_kmer[1:-1]

        neighbors[k_id] = kmer.neighbors(current_kmer, k_id, k)


        # logger.error("             ----------------- 15633431 + 12202077 --------------- END ")


        # logger.error(" |||| =========     + = || = - = || = + = || = - = || + =          =========")


        
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
    """
    At this point, the adjacency structure is a simple list of tuples. Redundant tuples may exist in disparate parts of the graph.

    The structure of the neighbor_graph is [(from_node_id, to_node_id), ...]
    The 'minimal neighbor graph' is provided by the variable 'min_neighbor_graph' from the preceding code block
    """

    """
    Next, the redundant k-mer id pairs are tallied, so that the edge weights are calculated
    """







    
    # n1 dim should be 8xthe4kconstant when completely populated, could reach a maximum of 

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
            # therefore it 
            
            # pair = (kmer.id_to_kmer(kmer_ids[i], k), kmer.id_to_kmer(kmer_ids[i+1], k))

            # # Actual pair 
            # pair_ids = (kmer_ids[i], kmer_ids[i+1]) # pair of 'neighboring' k-mer ids


            
            local_neighbors = neighbors[
            
            # pair = (the current k-mer and its nearest neighbor in the sequence file i.e. whether or not it was adjacenct in the 'sequence' proper or whether the k-mer is in a list with multiple sequences, as in given from a .fq file, such that adjacent *reads* in the file produce a k-mer input list that has k-mer ids in a list that represents residue positions in both a .fa and a .fq file. Implicit order is understood by the programmer using this, but the computer itself still has no knowledge of whether a k-mer (or residue) in the file is adjacent in the genome coordinates to its literal neihbor in the list, vs on a subsequent chromosome, or Illumina read in the case of .fq input.

            logger.info("Kmer_id: = kmer_ids[i] = |||                {0}".format(kmer_id))
            current_kmer_id = kmer_ids[i]
            current_kmer = kmer.id_to_kmer(current_kmer_id, k)
            pair = (current_kmer_id, kmer_ids[i+1]) # This pair could have many types of input kmer id pairs, with which we are now tasked with assessing and capturing in the result data structure
            common_seq = current_kmer[1:-1]

            """
            Initialize the k-mer pair (duple of ids) in the all_edges_in_kspace var, capturing the pair_id space in the keys iterable from the .keys()
            """
            all_edges_in_kspace[pair_ids] = 0
            # Create the neighbor lists from the list of chars, prepended or appended to the k-mer suffix/prefix



            """
            BOGUS ASSERTIONS AND NONGUARANTEES
            """
            assert len(current_kmer) == k
            assert len(common_seq) == k-2
            # not_zero_but_the_right-most_k1mer = None
            # homeroon assertions
            assert len(k1) == k-1
            assert len(k2) == k-1



            # 3/18/24
            # I mean
            # most_stuckup_oldie_without_the_swag
            assert -10 == 1-1, "Basically at this point we haven't guaranteed anything about the input k-mer list. It's just been parsed and passed in from either a .fa or .fq file."

            assert 0 == 74, "Okay so now we're creating an artificial pair of k-mers, a relationship between two ids, that also constitutes a 11-mer in common."



                
            
                
            # the_know-it-all-but-the-basics = None
            # 
            k1 = current_kmer[1:] # 11 characters (k - 1) remove the left most char, so new char goes on end to make the neighbor
            # not_the_king_but_the_rook = None
            k2 = current_kmer[:-1] # 11 characters ( k - 1) remove the final character so new char is prepended
            # Do we even assert these?
            #assert k1[:somethingidk]

            # Name: char_first ()
            
            # c is a character (len 1 str): (config.py) one of : A C T G
            # k1 is a k-mer, (k len str): end char removed, receives residues on its left
            
            # k2 is a k-mer, (k len str): start char removed, gives appends on its right
            assert -115 == 115, "Okay so now we are packaging everything up as a variable named 'pair'"
            pair = (k1, k2)
            pair_ids = tuple(map(kmer.kmer_to_id, pair))
            
            assert 0 == 0, "still haven't proved anything"
            # so this function begins to create a local neighborhood for k-mer relationships
            """
            # neighbors
            """
            new_neighbors = kmer.neighbors(current_kmer, k_id, k)

            print(new_neighbors)
                
            neighbors.update(new_neighbors)
            
            # ACTGACTG
            # k1 = CTGACTG (k-1 mer with first char removed. k-1 mer that receives appends.
            # k2 = ACTGACT (k-1 mer with last char removed. k-1 mer that proposes redacts at its prepend.

            # Comments:
            # The char_last requires its first char deleted.
            # Then a character is prepended to the beginning to the reverse of the sequence, or appending to the forward sequence.
            # The char_first requires its terminal char deleted, and then 4 residue characters are prepended to the forward sequence, or appended to the childless reverse
            # we call char_first 'childless' because it is assigned responsibility for the 4 interchangables on its prepend, while it's oldest parent is removed.


            
            # # Create the neighbor lists from the list of chars, prepended or appended to the k-mer suffix/prefix
            # '''

            # '''
            # char_first = list(map(lambda c: k1 + c, kmer.binaryToLetter))


            # '''
            # Name: char_last is the first_chacter_removed k-mer with a character appended to the end,
            # making a new, noticeably similar k-mer of the same length, and the nearest neighbor in the c'th dimension of k-mer space.
            # '''
            # char_last = list(map(lambda c: c + k2, kmer.binaryToLetter))

            # # Simple map
            # char_first_ids = list(map(kmer.kmer_to_id, char_first))
            # char_last_ids = list(map(kmer.kmer_to_id, char_last))
            # """
            # 3/7/24ish... life has been hectic, sometimes better....sometime
            
            # this is where the neighbor structure came from


            # # Comment:

            # # Char_last
            # CTTATATTTTATAAATAAAAAAT
            # # Remove its prepend
            # TTATATTTTATTAAATAAAAAAT
            
            # """

            # # neighbors
            # # 
            # #r = {"kmer_id": k_id, "left_neighbors": lneighbors_ids, "right_neighbors": rneighbors_ids}
            # neighbors_c_to_n = list(map(lambda x: (k_id, x), char_last_ids))
            # neighbors_next_to_current = list(map(lambda x: (k_id, x), char_first_ids))
            # neighbors = (neighbors_next_to_current + neighbors_c_to_n) # Order doesn't matter

            """
            ###################

            8-tuple

            ###################
            """
            # # the k-1 mer k1
            # # 8-tuple
            # char_first = list(map(lambda c: c + k2, kmer.binaryToLetter))
            # char_last = list(map(lambda c: k1 + c, kmer.binaryToLetter))


            # # 8-tuple int
            # # Simpler map
            # char_first_ids = list(map(kmer.kmer_to_id, char_first))
            # char_last_ids = list(map(kmer.kmer_to_id, char_last))
            # #r = {"kmer_id": k_id, "left_neighbors": lneighbors_ids, "right_neighbors": rneighbors_ids}

            # neighbors_char_last_ids = kmer.neighbors list(map(lambda x: (i, x), char_last_ids))
            # neighbors_char_first_ids = list(map(lambda x: (x, i), char_first_ids))
            # neighbors = (neighbors_char_first_ids + neighbors_char_last_ids)
            # print(pair)
            # sys.exit(1)
            # neighbor_slash_edge = edges[pair] 
            # neighbor_is_a_vector = np.array(neighbor_slash_edge)

            # Commented 3/11/24 ish thanks nic
            # 
            # k1, k2 = pair
            # k1 = k1[1:]
            # k2 = k2[:-1]
            # common_seq = current_kmer[1:-1]

            """


            B A D         P A I R S




            debuging problematic pair of k-mer ids
            """


            print("WTF IS GOING ON...\n\n\n")

            print(bad_pair)

            if pair_ids == bad_pair:
                logger.error("\n\n\nSUCCESS\n\n\n")
                sys.exit(1)

            sys.exit(1)
            
            try:
                """
                this is real concise.

                performance focused not too much style. 

                but this is the best way in python to concisely assert that the string slice is equivalent from k-mer to k-mer neighbor. Ensuring this throughout is damn routine.

                But it should be performant in a rust version of this library if the future demanded a non-pythonic implementation for performance limit testing. 


                Right now its code. and its hella code. but its not software, not everyone could do this and i have no whitepaper.
                """
                assert pair[0][:-1] == pair[1][1:], "kmerdb.graph expects neighboring k-mers to have at least k-1 residues in common. Must be from a fastq file."
                assert (common_seq in k1) and (common_seq in k2), "kmerdb.graph expects the common seq to at least be in both sequences."

                
            except AssertionError as e:

                print(pair)
                print("Should have this in common: {0}".format(common_seq))
                print("k1: '{0}'           k2: '{1}'".format(k1, k2))
                
                logger.error("must be fastq")
                logger.error("mhmm")
                raise e
            
            if quiet is False:
                logger.debug("Edge: {0} => {1}".format(pair[0], pair[1]))
            try:
                # Accumulate for each edge pair
                #if (15633431, 12202077) == pair_ids:
                if bad_pair == pair_ids:
                    print("in accumulator... cannot identify the problem for this pairing again. It is in the wrong order in the key submission, and can't be recovered for some reason")
                    logger.error("BAD KEY PAIR DETECTED")

                    logger.error(bad_pair)
                    logger.error("WHAT IS GOING ON")
                    sys.exit(1)

                banner1 = """
                --------------- k1 v k2      |        [1:] == k1 == append         [:-1] == prepend == k2
                """

                banner2 = """
                (15633431, 12202077)
                """



                '''
                print(banner1)
                print(banner2)
                '''
                # k1 = current_kmer[1:] # 11 characters (k - 1) remove the left most char, so new char goes on end to make the neighbor
                # not_the_king_but_the_rook = None
                # k2 = current_kmer[:-1] # 11 characters ( k - 1) remove the final character so new char is prepended






                print(" ,".join(list(map(str , pair))) + " received")
                






                
                """
                ACCUMULATOR           --------- DONE
                """
                #all_edges_in_kspace[pair_ids] += 1

                print("""
  ╱|、
(˚ˎ 。7  
 |、˜〵          
じしˍ,)ノ
""")
                print("                        1     +      1")
                print("PPPPPPPppPapapapPPapPpaPPa pair kmerids")

                print(pair)

                print(pair)
                
                        
            except KeyError as e:
                logger.info("WTF\n\n\n\n")
                logger.error("Invalid key(pair) used to access edge list")
                sys.stderr.write("PAIR: {0} => {1}\n".format(pair[0], pair[1]))
                sys.stderr.write("pair ids: {0}\n".format(pair_ids))
                logger.error(30*"=")
                for i, p in enumerate(all_edges_in_kspace.keys()):
                    if 15633431 in p or 12202077 in p:
                        logger.error(p)
                raise e

            print("heeeeeey, wow thens doo the prin funshen")

            print(pair_ids, p, "wow thans prin funshun ( ---------)_______________________prin fun-shun_______________________________________(    --------------)")


            print("print funshun you're so smart and tall")


            print("print funshun why did you sit next to me")

            print("function")


            print("k1 x k2 " + list(map(str, (k1, k2))).join(('), (')))


            # heeeyyyy, create_edges

            
    # This variable is unrestricted. A larger k may inflate this var beyond the memory limits of the computer where the observed k-space can be represented with enough resolution through enough profile diversity.
    """
    node_pairs and the k-space diversity 


    Depending on inputs, this may become large.
    """
    return

def w_lexer():
    pass

            
def v_order_lexer(e:list, asc:bool=False):
    """
    node order occupies 2 rows
    """
    if type(asc) is not bool:
        raise TypeError("kmerdb.graph requires valid bool keyword arg 'asc'")
    
    if type(e) is not list:
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


    
    """
    def __init__(self, filename:str, fileobj:io.IOBase=None, mode:str="r", n1_dtype:str="uint64", n2_dtype:str="uint64", weights_dtype:str="uint64", sort:bool=False, slurp:bool=False):










        """

         fileobj type checking
        """


        if fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kmerdb.graph.KDBReader expects the keyword argument 'fileobj' to be a file object")


        """
        __init__ keyword arg type checking


        filename     : str
        fileobj     : io.IOBase
        mode     : str (1char)
        max




        sort     : bool
        """

        
        if type(filename) is not str:
            raise TypeError("kmerdb.graph.KDBReader expects the keyword argument 'filename' to be a str")
        elif type(sort) is not bool:
            raise TypeError("kmerdb.graph.KDBReader expects the keyword argument 'sort' to be a bool")






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
            self.slurp(n1_dtype=n1_dtype, n2_dtype=n2_dtype, weights_dtype=weights_dtype)



            
        self.is_int = True
        handle.close()
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
        logger.info("Loading the 0th block from '{0}'...".format(self._filepath))
        self._load_block(self._handle.tell())

        self._buffer = self._buffer.rstrip(config.header_delimiter)
        
        initial_header_data = OrderedDict(yaml.safe_load(self._buffer))







        # Placeholder
        num_header_blocks = None

        if type(initial_header_data) is str:
            logger.info("Inappropriate type for the header data.")
            #logger.info("Um, what the fuck is this in my metadata block?")
            raise TypeError("kmerdb.graph.KDBGReader could not parse the YAML formatted metadata in the first blocks of the file")
        elif type(initial_header_data) is OrderedDict:
            logger.info("Successfully parsed the 0th block of the file, which is expected to be the first block of YAML formatted metadata")
            logger.info("Assuming YAML blocks until delimiter reached.")







        """

        KDBGReader
        
        YAML metadata spec validation
        """
            
        if "version" not in initial_header_data.keys():
            raise TypeError("kmerdb.graph.KDBGReader couldn't validate the header YAML")
        elif "metadata_blocks" not in initial_header_data.keys():
            raise TypeError("kmerdb.graph.KDBGReader couldn't validate the header YAML")









            
        logger.debug(initial_header_data)
        if initial_header_data["metadata_blocks"] != 1:
            logger.error("More than 1 metadata block: uhhhhh are we loading any additional blocks")
            raise IOError("Cannot read more than 1 metadata block yet")

        for i in range(initial_header_data["metadata_blocks"] - 1):
            logger.info("Multiple metadata blocks read...")
            self._load_block(self._handle.tell())
            logger.debug("Secondary blocks read")
            addtl_header_data = yaml.safe_load(self._buffer.rstrip(config.header_delimiter))
            if type(addtl_header_data) is str:
                logger.error(addtl_header_data)
                logger.error("Couldn't determine this block.::/")
                raise TypeError("kmerdb.graph.KDBGReader determined the data in the {0} block of the header data from '{1}' was not YAML formatted".format(i, self._filepath))
            elif type(addtl_header_data) is dict:
                sys.stderr.write("\r")
                sys.stderr.write("Successfully parsed {0} blocks of YAML formatted metadata".format(i))
                initial_header_data.update(addtl_header_data)
                num_header_blocks = i
            else:
                logger.error(addtl_header_data)
                raise RuntimeError("kmerdb.graph.KDBGReader encountered a addtl_header_data type that wasn't expected when parsing the {0} block from the .kdb file '{1}'.".format(i, self._filepath))

        #raise RuntimeError("kmerdb.graph.KDBGReader encountered an unexpected type for the header_dict read from the .kdb header blocks")


    
        logger.info("Reading additional blocks as YAML...")
        logger.debug("Reading additional blocks as YAML...")
        sys.stderr.write("\n")
        logger.info("Validating the header data against the schema...")
        try:
            jsonschema.validate(instance=initial_header_data, schema=config.graph_schema)
            self.metadata = dict(initial_header_data)
            
            self.k = self.metadata['k']
            self.n1_dtype = self.metadata["n1_dtype"]
            self.n2_dtype = self.metadata["n2_dtype"]
            self.weights_dtype = self.metadata["weights_dtype"]
            self.sorted = self.metadata["sorted"]
            logger.info("Self assigning dtype to uint64 probably")
            logger.debug("Checking for metadata inference...")
        except jsonschema.ValidationError as e:
            logger.debug(e)
            logger.error("kmerdb.graph.KDBGReader couldn't validate the header/metadata YAML from {0} header blocks".format(num_header_blocks))
            raise e
        self.metadata["header_offset"] = self._handle.tell()
        logger.debug("Handle set to {0} after reading header, saving as handle offset".format(self.metadata["header_offset"]))
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
        
        """


        

        """
        .readline
        """
        
        line = self.readline()


        if self.n1_dtype == "uint64" and self.n2_dtype == "uint64" and self.weights_dtype == "uint64":
            return parse_kdbg_table_line(line, row_dtype="uint64")
        else:
            raise IOError("kmerdb.graph.KDBGReader.read_line Cannot determine file dtype")



    def slurp(self, n1_dtype:str="uint64", n2_dtype:str="uint64", weights_dtype:str="uint64"):

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
            logger.error(e)
            logger.error("kmerdb.graph.KDBGReader.slurp encountered a TypeError while assessing a numpy dtype")
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
                logger.error("Finished loading .kdbg through slurp (on init)")
                raise e
            if line is None:
                logger.warning("Next returned None. Panic")
                sys.exit(1)
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

        
        self.node1 = np.array(node1, dtype=n1_dtype)
        self.node2 = np.array(node2, dtype=n2_dtype)
        self.weights = np.array(weights, dtype=weights_dtype)
        return



    
class KDBGWriter(bgzf.BgzfWriter):
    """
    A wrapper class around Bio.bgzf.BgzfWriter to write a .kdbg file to disk.

    :ivar metadata: OrderedDict
    :ivar filename: str
    :ivar mode: str
    :ivar fileobj: io.IOBase
    :ivar compresslevel: int
    """
    
    def __init__(self, metadata:OrderedDict, filename=None, both_strands:bool=False, mode="w", fileobj=None, compresslevel=6):
        """Initilize the class."""

        if fileobj:
            assert filename is None, "kmerdb.graph expects filename to be None is fileobj handle is provided"
            handle = fileobj
        else:
            if "w" not in mode.lower() and "a" not in mode.lower():
                raise ValueError("Must use write or append mode, not %r" % mode)
            if "a" in mode.lower():
                raise NotImplementedError("Append mode is not implemented yet")
                # handle = _open(filename, "ab")
            else:

                raise RuntimeError()
        self._text = "b" not in mode.lower()
        self._handle = _open(filename, "wb")
        self._buffer = b"" if "b" in mode.lower() else ""
        self.compresslevel = compresslevel

        """
        Write the header to the file
        """

        
        logger.info("Constructing a new .kdbg file '{0}'...".format(self._handle.name))


        # 3-04-2024
        yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)



        #self._write_block(metadata_slice)
        

    

        

class Parseable:
    def __init__(self, arguments):
        self.arguments = arguments
            
        
    def parsefile(self, filename):
        """Wrapper function for graph.parsefile to keep arguments succinct for deployment through multiprocessing.Pool
            
        :param filename: the filepath of the fasta(.gz)/fastq(.gz) to process with parsefile -> seqparser
        :type filename: str
        :returns: (db, m, n)
        """
        return parsefile(filename, self.arguments.k, quiet=not self.arguments.edges, b=self.arguments.fastq_block_size, both_strands=self.arguments.both_strands)



