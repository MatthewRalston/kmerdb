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
import numpy as np
import networkx as nx

from kmerdb import fileutil, parse, kmer, config, util



def graph_convert(n1:np.ndarray, n2:np.ndarray, weights:np.ndarray, k:int=None):
    """
    :param n1: A numpy.ndarray of node 1 kmer ids from an edge list
    :type numpy.ndarray:
    :param n2: A numpy.ndarray of node 2 kmer ids from an edge list
    :type numpy.ndarray:
    :param weights: A numpy.ndarray of weights corresponding to the edge
    :type numpy.ndarray:
    """

    if type(k) is not int:
        raise TypeError("kmerdb.graphlib.graph_convert() expects the keyword argument 'k' to be an int")
    elif isinstance(n1, np.ndarray) is False:
        raise TypeError("kmerdb.graphlib.graph_convert() expects a numpy.ndarray as its first positional argument")
    elif isinstance(n2, np.ndarray) is False:
        raise TypeError("kmerdb.graphlib.graph_convert() expects a numpy.ndarray as its second positional argument")
    elif isinstance(weights, np.ndarray) is False:
        raise TypeError("kmerdb.graphlib.graph_convert() expects a numpy.ndarray as its third positional argument")
    elif n1.shape == () or n2.shape == () or weights.shape == ():
        raise ValueError("kmerdb.graphlib.graph_convert() found that numpy.ndarrays were empty")
    elif n1.shape != n2.shape:
        raise ValueError("kmerdb.graphlib.graph_convert() found that the node1 and node2 array lengths did not match |  node1: {0} node2: {1}".format(n1.shape, n2.shape))
    elif n2.shape != weights.shape:
        raise ValueError("kmerdb.graphlib.graph_convert() found that the node2 and weights array lengths did not match |  node2: {0} weights: {1}".format(n2.shape, weights.shape))

    N = 4**k
    # edge_list = list(zip(n1.tolist(), n2.tolist(), weights.tolist()))

    
    edge_list_length = n1.shape[0] #
    edge_list = [ (int(n1[i]), int(n2[i]), {'weight': int(weights[i])}) for i in range(edge_list_length)]
    G = nx.Graph()


    G.add_nodes_from(range(N)) # Adds a list of k-mer ids as nodes to the graph
    G.add_edges_from(edge_list)

    print("Number of nodes: {0}".format(G.number_of_nodes()))
    print("Number of edges: {0}".format(G.number_of_edges()))


    
def make_coo_graph(seq_filepath:str, k:int, quiet:bool=True):
    if type(seq_filepath) is not str:
        raise TypeError("kmerdb.graph.make_coo_graph() expects a fasta/fastq sequence filepath as a str as its first positional_regument")
    elif type(k) is not int:
        raise TypeError("kmerdb.graph.make_coo_graph() expects an int for k as the second positional argument")
    elif type(quiet) is not bool:
        raise TypeError("kmerdb.graph.make_coo_graph() expects the keyword argument 'slurp' to be a bool")

    if os.path.exists(filename) is False or os.access(filename, os.R_OK) is False:
        raise ValueError("kmerdb.graph.make_coo_graph() expects the filepath to be be readable on the filesystem")
    N = 4**k
    dtype = "uint32" if k < 17 else "uint64"
    """
    # {'i': tkmer_id, 'kmer_id': kmer_id, 'seq_id': seq_id, 'pos': pos}
    # (i, kmer_id, seq_id, pos)
    # We need a list of tuples with the following attributes:
    #      - the index of the kmer in the TOTAL k-mer sequence space
    #      - the kmer_id for the k-mer starting that position
    #      - the seq_id for the sequence that was shredded
    #      - and an ALMOST identical index of the position of the k-mer id
    """
    
    col1 =[] # IS THE INDEX not the k-mer id of the k-mer
    col2 =[] # 
    vals =[]
    # col1 = np.zeros(range(N), dtype=dtype)
    # col2 = np.zeros(range(N), dtype=dtype)
    # vals = np.zeros(range(N), dtype=dtype)
    for seq_id, seq in parse_sequence_file(seq_filepath, return_tuple=True):
        """
        Shred sequences and get a list of kmer ids in order. We do not count kmers at the moment.
        Each k-mer needs to be converted to pairs of the id and all neighbors, 8 records total, corresponding to the neighboring kmer_ids by taking a single base from the left/right side of the kmer_id of interest.

        The maximum number (maximum k-mer id) that can be stored in a uint32 numpy.ndarray is 429967295 exactly
        therefore unit32 can only be used when k <= 16
        Store the results in numpy uint32 (~4B is the max size of unit32) if k<16 EXACTLY or uint64 otherwise
        
        """
        
    
