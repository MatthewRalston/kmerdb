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



import logging
import argparse
import os
import sys
import yaml
import json
import time
import signal


import pdb

from multiprocessing import cpu_count

from collections import OrderedDict


from Bio import bgzf

#import concurrent.futures

from kmerdb import logger as kdbLogger

global logger
logger = None


global step
step = 0

global feature
feature = 0

global exit_code
exit_code = -1


global exit_summary
exit_summary = None


def version(arguments):
    from kmerdb import config

    print(config.VERSION)


def print_argv():
    argv = sys.argv
    sys.stderr.write(" ".join(argv[0:4]) + " ...\n")

def citation(arguments):

    MODULE_ROOT = os.path.dirname(__file__)
    citation_file = os.path.join(MODULE_ROOT,  'CITATION.txt')
    if os.access(citation_file, os.R_OK):

        sys.stderr.write("Removing '{0}'\n".format(citation_file))
        os.remove(citation_file)

    sys.stderr.write("On the real, gotta eat.\n")
    sys.stderr.write("Consider a +1 on Github to keep it real...\n\n")

def index_file(arguments):
    from kmerdb import fileutil, index
    from kmerdb.config import DONE


    idx = []
    with fileutil.open(arguments.kdb, mode='r', slurp=True, with_index=True) as kdb:
        k = kdb.metadata['k']
        metadata = kdb.metadata
            
        line_index = index.open(filepath=arguments.kdb, indexfilepath=arguments.kdb + "i", mode="u", idx=kdb.index, force=arguments.force, DEBUG=arguments.debug)
    sys.stderr.write(DONE)


def shuf(arguments):
    import random
    
    from kmerdb import fileutil, index
    from kmerdb.config import DONE

    if os.path.splitext(arguments.kdb_in)[-1] != ".kdb":
        raise ValueError("kdb shuf requires an indexed .kdb file. Received '{0}'".format(arguments.kdb))
    else:
        arguments.kdb_index = arguments.kdb_in + "i"
    
        with fileutil.open(arguments.kdb_in, mode='r') as kdb:
            with index.open(arguments.kdb_index, mode='r') as kdbi:
                shuffled = list(range(len(kdbi.index)))
                random.shuffle(shuffled)
                random.shuffle(shuffled)
                with fileutil.open(arguments.kdb_out, mode='w', metadata=kdb.metadata) as kdb_out:
                    x = 1
                    for i in shuffled:
                        kmer_id, count, kmer_metadata = index.read_line(kdb, kdbi, i)
                        kdb_out.write("{0}\t{1}\t{2}".format(x, count, kmer_metadata))
                        x += 1
    sys.stderr.write(DONE)
    
# def markov_probability(arguments):
#     """
#     A very old function that is probably broken in the indexing function.
    
#     :param arguments: argparse Namespace
#     :type arguments:
#     """
#     import pandas as pd
#     import numpy as np
#     from kmerdb import fileutil, index, probability, parse
#     from kmerdb.config import DONE

#     if os.path.splitext(arguments.kdb)[-1] != ".kdb":
#         raise IOError("Model .kdb filepath does not end in '.kdb'")


#     if index.has_index(arguments.kdb):
#         arguments.kdbi = arguments.kdb + "i"
#         #df = pd.DataFrame([], columns=["SequenceID", "Log_Odds_ratio", "p_of_seq"])
#         profiles = np.array([], dtype="int64")
#         with fileutil.open(arguments.kdb, 'r', slurp=True) as kdb:
#             k = kdb.metadata['k']
#             with index.open(arguments.kdbi, 'r') as kdbi:
#                 with parse.SeqParser(arguments.seqfile, arguments.fastq_block_size, k) as seqprsr:
#                     recs = [r for r in seqprsr]
#                     if seqprsr.fastq:
#                         logger.debug("Read exactly b=={0}=={1} records from the {2} seqparser object".format(b, len(recs), s))
#                         assert len(recs) == b, "The seqparser should return exactly {0} records at a time".format(b)
#                     else:
#                         logger.debug("Read {0} sequences from the {1} seqparser object".format(len(recs), seqprsr))
#                         logger.debug("Skipping the block size assertion for fasta files")

#                     while len(recs): # While the seqprsr continues to produce blocks of reads
#                         # Do something here
                    
#                         markov_probs = list(map(lambda p: [p["seq"].name, p["log_odds_ratio"], p["p_of_seq"]], [probability.markov_probability(seq, kdb, kdbi) for seq in recs]))

#                         sys.stderr.write(json.dumps(markov_probs))
#                         if profiles.shape == (0,):
#                             profiles = np.array(markov_probs)
#                         else:
#                             np.append(profiles, markov_probs, axis=0)

#                         recs = [r for r in seqprsr] # Essentially accomplishes an iteration in the file, wrapped by the parse.SeqParser class
#         df = pd.DataFrame(profiles, columns=["SequenceID", "Log_Odds_ratio", "p_of_seq"])
#         df.to_csv(sys.stdout, sep=arguments.delimiter, index=False)

#     else:
#         raise IndexError(".kdb file '{0}' has no corresponding index file. Please run 'kdb index -h' for details on index generation".format(arguments.kdb))

    
#     sys.stderr.write(DONE)






def expanded_help(arguments):

    import sys
    
    argv = sys.argv
    
    from kmerdb import config, appmap

    kmerdb_appmap = appmap.kmerdb_appmap( argv )

    #kmerdb_appmap.print_program_header()
    
    if arguments.method not in config.subcommands:
        raise ValueError("unsupported method")
    elif arguments.method == "profile":
        kmerdb_appmap.print_profile_header()
    elif arguments.method == "graph":
        kmerdb_appmap.print_graph_header()
    elif arguments.method == "index":
        #kmerdb_appmap.print_index_header()
        raise ValueError("Unsupported method")
    elif arguments.method == "shuf":
        #kmerdb_appmap.print_shuf_header()
        raise ValueError("Unsupported method")
    elif arguments.method == "matrix":
        kmerdb_appmap.print_matrix_header()
    elif arguments.method == "distance":
        kmerdb_appmap.print_distance_header()
    elif arguments.method == "minimizers":
        kmerdb_appmap.print_minimizers_header()

        
    sys.stderr.write("\n\nUse --help for expanded usage\n")

def get_minimizers(arguments):

    from Bio import SeqIO
    import numpy as np
    from kmerdb import fileutil, config, util, minimizer
    import json
    metadata = None
    N = None
    kmer_ids = None
    fasta_seqs = []
    coordmap = None
    
    assert type(arguments.kdb) is str, "kdb_in must be a str"
    kdb_sfx = os.path.splitext(arguments.kdb)[-1]
    fa_sfx = os.path.splitext(arguments.fasta)[-1]
    
    if kdb_sfx != ".kdb": # A filepath with invalid suffix
        raise IOError("Input .kdb filepath '{0}' does not end in '.kdb'".format(arguments.kdb))
    elif not os.path.exists(arguments.kdb):
        raise IOError("Input .kdb filepath '{0}' does not exist on the filesystem".format(arguments.kdb))
    elif fa_sfx != ".fa" and fa_sfx != ".fna" and fa_sfx != ".fasta":
        raise IOError("Input .fasta filepath '{0}' does not exist on the filesystem".format(arguments.fasta))
    

    # coordmap = [(fasta_id_str, i:L, kmer_id, is_minimizer), ...]
    # kmer_ids length : input
    # is_minimizer length : m = L/window_size hashmap of sequence coords and is_minimizer bool
    # coordmap 4-tuple (fasta_id_str, i:floor(m), kmer_id)
    k, kmer_ids, fasta_ids, is_minimizer, coordmap, coord_obj_array = minimizer.minimizers_from_fasta_and_kdb(arguments.fasta, arguments.kdb, arguments.window_size)
    
    output_index_file = "{0}.{1}.kdbi.1".format(str(arguments.kdb).strip(".kdb"), k)
    logger.log_it("\n\nPrinting results of minimizers to '{0}'\n\n".format(output_index_file), "INFO")

    m = num_minimizers = len(is_minimizer)



    """
    Writes the fourtuple to a plain-text index-file (.kdbi.1)

    (fasta_id_str, 1:num_possible_minimizers, kmer_id, is_minimizer_bool)
    
    if --include-unselected-minimizers (-a) is provided, will write all L=N-k+1 possible minimizer selections from input sequences
    """
    with open(output_index_file, "w") as ofile:
        # if arguments.include_unselected_minimizers:
        #     for seq in is_minimizer.keys():

        #         """
        #         FIXME
        #         """
        #         for i in range(len(is_minimizer[seq])):
        #             """
        #             fasta_id_str, coord, kmer_id_at_coord, is_minimizer
        #             """
        #             ofile.write("\t".join((fasta_ids[i], i, coordmap[i][2], is_minimizer[i], )))
        # else:
        
        # for minimizer_sequence_coordinate, i in enumerate(range(len(coordmap))):

        #     ofile.write("\t".join((fasta_ids[i], minimizer_sequence_coordinate, coordmap[i][2], True)) + "\n")


        for i, ktuple in enumerate(coord_obj_array):
            ofile.write("\t".join(list(map(str, ktuple))) + "\n")
            if not arguments.quiet:
                print(ktuple)

        
        
    logger.log_it("Wrote minimizers to .kdbi.1 index file '{0}'...".format(output_index_file), "INFO")

            

def get_alignments(arguments):
    """
    Produces Smith-Waterman alignment by loading minimizers

    """
    import tempfile
    import copy
    import gzip
    import shutil

    
    from Bio import SeqIO

    from kmerdb import minimizer, fileutil, parse, kmer, util

    
    reference_minimizers = []
    query_minimizers = []

    reference_fasta_seqs = []
    query_fasta_seqs = []

    reference_fasta_ids = []
    query_fasta_ids = []

    kmer_ids = None
    k = None
    
    query_fasta_sfx = os.path.splitext(arguments.query)[-1]
    #reference_kdb_sfx = os.path.splitext(arguments.kdb)[-1]
    #reference_kdbi_sfx = ".".join(os.path.splitext(arguments.kdbi)[-2:])
    reference_fasta_sfx = os.path.splitext(arguments.reference)[-1]
    reference_kdb_sfx = os.path.splitext(arguments.reference_kdb)[-1]

    # Create as temporary file, hold minimizer array in memory



            
    if reference_fasta_sfx != ".fasta" and reference_fasta_sfx != ".fna" and reference_fasta_sfx != ".fa": # A filepath with invalid suffix
        raise IOError("Viewable fasta filepath does not end in '.fasta'/'.fna.'.fa'")
    elif not os.path.exists(arguments.reference):
        raise IOError("Viewable .fasta filepath '{0}' does not exist on the filesystem".format(arguments.reference))
    elif query_fasta_sfx != ".fasta" and query_fasta_sfx != ".fna" and query_fasta_sfx != ".fa": 
        raise IOError("Viewable .fasta filepath does not end in '.fasta'/'.fna.'.fa'")
    elif not os.path.exists(arguments.query):
        raise IOError("Viewable .fasta filepath '{0}' does not exist on the filesystem".format(arguments.query))
    elif reference_kdb_sfx != ".kdb":
        raise IOError("Viewable .kdb filepath does not end in '.kdb'")
    elif not os.path.exists(arguments.reference_kdb):
        raise IOError("Viewable .kdb filepath '{0}' does not exist on the filesystem".format(arguments.reference_kdb_sfx))


    if str(arguments.reference).endswith(".fasta.gz") or str(arguments.reference).endswith(".fna.gz") or str(arguments.reference).endswith(".fa.gz"):
        with gzip.open(arguments.reference, mode="r") as ifile:
            for record in SeqIO.parse(ifile, "fasta"):
                reference_fasta_seqs.append(str(record.seq))
                reference_fasta_ids.append(str(record.id))
    else:
        with open(arguments.reference, mode="r") as ifile:
            for record in SeqIO.parse(ifile, "fasta"):
                reference_fasta_seqs.append(str(record.seq))
                reference_fasta_ids.append(str(record.id))
    if str(arguments.query).endswith(".fasta.gz") or str(arguments.query).endswith(".fna.gz") or str(arguments.query).endswith(".fa.gz"):
        with gzip.open(arguments.query, mode="r") as ifile:
            for record in SeqIO.parse(ifile, "fasta"):
                query_fasta_seqs.append(str(record.seq))
                query_fasta_ids.append(str(record.id))
    else:
        with open(arguments.query, mode="r") as ifile:
            for record in SeqIO.parse(ifile, "fasta"):
                query_fasta_seqs.append(str(record.seq))
                query_fasta_ids.append(str(record.id))

                



    # Make Query kdb


    # Calculate minimizers

     
    if reference_kdb_sfx == ".kdb":
        with fileutil.open(arguments.reference_kdb, mode='r', slurp=True) as kdb_in:
            metadata = kdb_in.metadata

            
            kmer_ids_dtype = metadata["kmer_ids_dtype"]
            N = 4**metadata["k"]
            k = metadata["k"]
            if metadata["version"] != config.VERSION:
                logger.log_it("KDB version is out of date, may be incompatible with current KDBReader class", "WARNING")
            logger.log_it("Reading from file...", "INFO")

            kmer_ids = kdb_in.kmer_ids

    
    newargs = copy.deepcopy(arguments)

    newargs.k = k
    newargs.minK = None
    newargs.maxK = None

    
    #newargs.output_name = ".".join(os.path.splitext(str(arguments.query))[:-1])
    
    newargs.output_name = "temporary_kdb_file"
    
    newargs.parallel = 1
    newargs.input = [str(newargs.query)]

    newargs.show_hist = False
    newargs.sorted = False
    newargs.quiet = True
    #print(newargs)
    
    _profile(newargs)


    query_kdb_filepath = "{0}.{1}.kdb".format(newargs.output_name, k)

    if not os.path.exists(query_kdb_filepath):
        OSError("Could not locate temporary kdb filepath '{0}' on filesystem".format(query_kdb_filepath))
    else:
        logger.log_it("Completed generating .kdb file '{0}' from query fasta file '{1}'.".format(query_kdb_filepath, newargs.query), "INFO")

    kmer_ids, reference_fasta_ids, reference_fasta_coords, reference_minimizers = minimizer.minimizers_from_fasta_and_kdb(arguments.reference, arguments.reference_kdb, arguments.window_size)    
    kmer_ids, query_fasta_ids, query_fasta_coords, query_minimizers = minimizer.minimizers_from_fasta_and_kdb(arguments.query, query_kdb_filepath, arguments.window_size)


    print("unimplemented")
    #alignment.smith_waterman_with_minimizers(query_fasta_ids, reference_fasta_ids, query_fasta_seqs, reference_fasta_seqs, arguments.reference_kdb, query_kdb_filepath, query_coords, reference_coords, reference_minimizers, query_minimizers)
        
    

    print(config.DONE)
    
    
def distances(arguments):
    """
    An end-user function to provide CLI access to certain distances
    """
    from multiprocessing import Pool
    import multiprocessing as mp

    import pandas as pd
    import numpy as np
    from scipy.spatial.distance import pdist, squareform


    from kmerdb import fileutil, config, util, distances




    global feature
    global step


    

    has_cython = False
    try:


        from kmerdb.distance import pearson_correlation
        from kmerdb import python_distances as distance
        
        logger.log_it("Importing distance module...", "INFO")
        has_cython = True
    except ImportError as e:
        logger.log_it(e.__str__(), "ERROR")
        logger.log_it("Could not import Cython distance submodule", "ERROR")
        raise RuntimeError("Cannot import pyx library")

    n = len(arguments.input)

    if len(arguments.input) > 1:

        feature = 1
        step = 1

        
        file_reader = fileutil.FileReader()
        files = list(map(lambda f: fileutil.open(f, 'r', slurp=False), arguments.input))
        
        logger.log_it("Files: {0}".format(files), "DEBUG")

        
        if not all(os.path.splitext(kdb)[-1] == ".kdb" for kdb in arguments.input):
            raise IOError("One or more parseable .kdb filepaths did not end in '.kdb'")
        dtypes = [x.metadata["kmer_ids_dtype"] for x in files]
        suggested_dtype = dtypes[0]


        logger.log_it(", ".join(dtypes), "ERROR")
        logger.log_it("Suggested dtype: {0}".format(suggested_dtype), "ERROR")

        logger.log_it("\n\n\n\n{0}\n\n\n\n".format(suggested_dtype), "INFO")

        
        ks = [kdbrdr.k for kdbrdr in files]
        suggested_k = ks[0]
        if not all(kdbrdr.k == suggested_k for kdbrdr in files):
            logger.log_it("Files: {0}".format(files), "ERROR")

            logger.log_it("Choices of k: {0}".format(ks), "ERROR")
            logger.log_it("By default the inferred value of k is used when k is not specified at the command line, which was {0}".format(suggested_k), "ERROR")

            
            raise TypeError("One or more files did not have k set to be equal to {0}".format(arguments.k))


        
        logger.log_it("Files: {0}".format(files), "DEBUG")

        
        if not all(kdbrdr.dtype == suggested_dtype for kdbrdr in files):
            raise TypeError("One or more files did not have dtype = {0}".format(suggested_dtype))
        if arguments.parallel > 1:
            with Pool(processes=arguments.parallel) as pool:
                files = pool.map(file_reader.load_file, arguments.input)
        else:
            files = list(map(lambda f: fileutil.open(f, 'r', slurp=True), arguments.input))


        data = [kdbrdr.slurp() for kdbrdr in files]



        
        logger.log_it("Suggested dtype {0}".format(suggested_dtype), "ERROR")
        
        profiles = np.array(data, dtype=suggested_dtype)

        #profiles = np.transpose(profiles)

        # The following does *not* transpose a matrix defined as n x N=4**k
        # n is the number of independent files/samples, (fasta/fastq=>kdb) being assessed
        # N=4**k is the dimensionality of the vector, sparse or not, that makes up the perceived profile, the descriptors is the column dimension.
        # The names for the columns are taken as the basenamed filepath, with all extensions (substrings beginning with '.') stripped.

        if arguments.column_names is None:
            columns = list(map(lambda kdbrdr: os.path.basename(kdbrdr._filepath).split(".")[0], files))
        else:
            with open(arguments.column_names, 'r') as column_names:
                columns = [line.rstrip() for line in column_names]
        if len(columns) != len(files):
            raise RuntimeError("Number of column names {0} does not match number of input files {1}...".format(len(columns, len(files))))

        
        logger.log_it("Shape: {0}".format(profiles.shape), "DEBUG")
        logger.log_it("Converting arrays of k-mer counts into a pandas DataFrame...", "INFO")
        
        #columns = list(df.columns) # I hate this language, it's hateful.
        n = len(columns)
    elif len(arguments.input) == 1 and (arguments.input[0] == "STDIN" or arguments.input[0] == "/dev/stdin"):

        feature = 2
        step = 2


        
        logger.log_it("Reading input as tsv/csv from STDIN", "INFO")
        try:
            df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.log_it(e.__str__(), "ERROR")
            logger.log_it("Pandas error on DataFrame reading. Perhaps a null dataset being read?", "ERROR")
            
            raise e
        profiles = np.array(df)
        #profiles = np.transpose(np.array(df))
        columns = list(df.columns)



        # profiles.shape = [row number, column number]
        #
        
        assert profiles.shape[1] == len(columns), "distance | Number of columns of dataframe did not match the number of rows of the profile"
        n = len(columns)

    elif len(arguments.input) == 1 and (os.path.splitext(arguments.input[0])[-1] == ".tsv" or os.path.splitext(arguments.input[0])[-1] == ".csv"):

        feature = 2
        step = 2
        
        try:
            df = pd.read_csv(arguments.input[0], sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.log_it("Pandas error on DataFrame reading. Perhaps a null dataset being read?", "ERROR")
            raise e
        except FileNotFoundError as e:
            logger.log_it(e.__str__(), "ERROR")
            logger.log_it("Input file not found", "ERROR")
            raise e
        #profiles = np.transpose(np.array(df))
        profiles = np.array(df)
        columns = list(df.columns) # I'm sorry I ever made this line. Please forgive me.
        # This is just gratuitous code and language. I'm really really not sure what I want to express here.
        #assert profiles.shape[0] == len(columns), "distance | Number of columns of dataframe did not match the number of rows of the profile"
        n = len(columns)
    elif len(arguments.input) == 1 and os.path.splitext(arguments.input[0])[-1] == ".kdb":
        logger.log_it("Not sure why you'd want a singular distance.", "ERROR")
        logger.log_it("kdb distance requires more than one .kdb file as positional inputs", "ERROR")
        sys.exit(1)
    else:
        logger.log_it("bin/kdb.distances() received {0} arguments as input, which were not supported.".format(len(arguments.input)), "ERROR")
        sys.exit(1)



    feature = 3
    step = 3
    logger.log_it("Custom calculating a {0}x{0} '{1}' distance matrix...".format(n, arguments.metric), "INFO")    

    #print(profiles.shape)
    #raise RuntimeError("phooey")
    
    if arguments.metric in ["pearson", "spearman"]:

        profiles = np.transpose(np.array(df))
        
        #files = list(map(lambda f: fileutil.open(f, 'r', slurp=True), arguments.input))
        data = [['' for x in range(n)] for y in range(n)]
        for i in range(n):
            for j in range(n):
                logger.log_it("Calculating the {0}x{1} cell".format(i, j), "DEBUG")
                logger.log_it("Calculating {0} distance between {1}:({2}) and {3}:({4}) columns...".format(arguments.metric, columns[i], i, columns[j], j), "INFO")

                if i == j:
                    logger.log_it("Was self, reverting to identity.", "INFO")
                    data[i][j] = distance.identity[arguments.metric]
                elif i > j:
                    data[i][j] = None
                    logger.log_it("Refusing to recompute custom distance metric for symmetric matrix", "INFO")

                elif i < j:
                    #print("Info: ixj {0}x{1}")
                    logger.log_it("Info: ixj {0}x{1}".format(i, j), "DEBUG")
                    
                    if arguments.metric == "pearson":
                        logger.log_it("Computing custom Pearson correlation coefficient...", "INFO")

                        
                        #data[i][j] = distance.correlation(profiles[i], profiles[j])
                        if has_cython is True:
                            idstr = "{0}x{1}".format(i, j)
                            r = mp.Value('d', 0.0)
                            p = mp.Process(target=pearson_correlation, args=(profiles[i], profiles[j], profiles[i].size, r))
                            p.start()
                            #data[i][j] = correlation(profiles[i], profiles[j], profiles[i].size)
                            p.join()
                            data[i][j] = r.value
                        else:
                            raise RuntimeError("Cannot calculate pearson correlation without NumPy and Cython")
                    elif arguments.metric == "spearman":
                        cor, pval = distance.spearman(profiles[i], profiles[j])
                        logger.log_it("Computing SciPy Spearman distance", "INFO")
                        
                        # FIXME! also print pval matrices
                        data[i][j] = cor
                    elif arguments.metric == "EMD":
                        logger.log_it("Custom EMD distance is deprecated", "ERROR")
                        raise RuntimeError("Genuine bug, please report the usage of deprecated custom distance functions")
                        data[i][j] = distance.EMD(profiles[i], profiles[j])
                    elif arguments.metric == "d2s":
                        logger.log_it("Custom d2s distance is deprecated", "ERROR")
                        raise RuntimeError("Genuine bug, please report the usage of deprecated custom distance functions")
                        data[i][j] = distance.d2s(profiles[i], profiles[j])
                    else:
                        logger.log_it("Other distances are not implemented yet", "ERROR")
                        sys.exit(1)
                # This double loop quickly identifies empty cells and sets the data correctly from the permutation above
        logger.log_it("Joining processes for parallel correlations", "INFO")
        logger.log_it("\n\n\nDone joining processes...", "INFO")


        

        logger.log_it("Completing matrix from shared process values", "INFO")
        for i in range(n):
            for j in range(n):
                if not i < j:
                    data[i][j] = None
        logger.log_it("Filling in symmetrical matrix...", "INFO")
        for i in range(n):
            for j in range(n):
                if data[i][j] is None:
                    data[i][j] = data[j][i]
                    #logger.info("\n\n\nSwerve swerve\n\n\n")

                    
        dist = np.array(data)
    else:

        
        dist = pdist(np.transpose(profiles), metric=arguments.metric)
        dist = squareform(dist)


    if dist.shape == (2,2):
        print(dist[0][1])
    else:

        df = pd.DataFrame(dist, columns=columns)
        
        ## FIXME: CUSTOM sorting code, not commiting to git repo
        #suffixes = [(int(x.split("_")[1]), i) for i, x in enumerate(column_names)] # A list of a 2-tuple of the correct sort order and the index
        #suffixes.sort(key=lambda x: x[0])
        #sorted_column_names = [column_names[s[1]] for s in suffixes]
        #df.sort_values(sorted_column_names)
        

        df.to_csv(sys.stdout, sep=arguments.output_delimiter, index=False)

def get_matrix(arguments):


    
    import logging as pylog
    pylog.getLogger('matplotlib.font_manager').disabled = True
    pylog.getLogger('matplotlib').setLevel(logging.WARNING)

    from multiprocessing import Pool
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    
    from kmerdb import fileutil, config



    global logger
    global feature
    global step


    step = 0
    feature = 0

    
    if len(arguments.input) > 1:

        
        files = list(map(lambda f: fileutil.open(f, 'r', slurp=False, logger=logger), arguments.input))

        if arguments.k is None:
            arguments.k = files[0].k
        if not all(os.path.splitext(kdb)[-1] == ".kdb" for kdb in arguments.input):
            raise IOError("One or more parseable .kdb filepaths did not end in '.kdb'")
        ks = [kdbrdr.k for kdbrdr in files]
        suggested_k = ks[0]
        if not all(k == suggested_k for k in ks):
            logger.log_it("Files: {0}".format(files), "ERROR")
            
            logger.log_it("Choices of k: {0}".format(ks), "ERROR")

            logger.log_it("By default the inferred value of k is used when k is not specified at the command line, which was {0}".format(arguments.k), "ERROR")

            logger.log_it("Default choice of k, since not specified at the CLI is {0}".format(arguments.k), "ERROR")

            logger.log_it("Proceeding with k set to {0}".format(suggested_k), "ERROR")

            
            raise TypeError("One or more files did not have k set to be equal to {0}".format(arguments.k))
        dtypes = [kdbrdr.metadata["count_dtype"] for kdbrdr in files]
        suggested_dtype = dtypes[0]
        if not all(kdbrdr.count_dtype == suggested_dtype for kdbrdr in files):
            raise TypeError("One of more files did not have dtype = {0}".format(suggested_dtype))

        if arguments.parallel > 1:
            infile = fileutil.FileReader(arguments, logger=logger)
            
            with Pool(processes=arguments.parallel) as pool:
                try:


                    
                    files = pool.map(infile.parsefile, arguments.input)
                    
                except StopIteration as e:
                    logger.log_it("Files: {)}".format(", ".join(arguments.input)), "ERROR")
                    raise e
        else:
            files = []

            for f in arguments.input:
                logger.log_it("Reading input from '{0}'...Please be patient...".format(f), "INFO")

                fh = fileutil.open(f, 'r', slurp=True, logger=logger) # no metadata on matrix command
                files.append(fh)


        data = [kdbrdr.counts for kdbrdr in files]
                     
        pure_data = np.array(data, dtype=suggested_dtype)
        profiles = np.transpose(pure_data)
        # The following does *not* transpose a matrix defined as n x N=4**k
        # n is the number of independent files/samples, (fasta/fastq=>kdb) being assessed
        # N=4**k is the dimensionality of the vector, sparse or not, that makes up the perceived profile, the descriptors is the column dimension.

        # The names for the columns are taken as the basenamed filepath, with all extensions (substrings beginning with '.') stripped.
        logger.log_it("============================================================", "INFO")
        logger.log_it("Converting arrays of k-mer counts into a pandas DataFrame...", "INFO")
        logger.log_it("============================================================", "INFO")
        


        if arguments.column_names is None:
            columns = list(map(lambda kdbrdr: os.path.basename(kdbrdr._filepath).split(".")[0], files))
        else:
            with open(arguments.column_names, 'r') as column_names:
                rows = [line.rstrip() for line in column_names]
                if len(rows) == 1 and type(rows[0]) == str:
                    columns = rows[0].split("\t")
                elif len(rows) > 1:
                    columns = list(map(lambda r: r.rstrip().rstrip(".kdb"),rows))
                else:
                    raise ValueError("Could not properly parse .tsv column_names accessory file")
        if len(columns) != len(files):
            raise RuntimeError("Number of column names {0} does not match number of input files {1}...".format(len(columns), len(files)))
        suggested_metadata = files[0].metadata
        expected, num_columns = 4**suggested_k, len(columns)
        if profiles.shape != (expected, num_columns):
            #logger.info("Time is money.")
            logger.log_it("Check shape self...", "ERROR")

            logger.log_it("Expected shape: {0} x {1}".format(expected, num_columns), "ERROR")

            logger.log_it("Actual shape: {0}".format(profiles.shape), "ERROR")

            
            raise RuntimeError("Raw profile shape (a NumPy array) doesn't match expected dimensions")
        df = pd.DataFrame(profiles, columns=columns)
        is_uint32 = False
        is_uint64 = False
        is_float32 = False
        is_float64 = False
        is_uint32 = all([x.dtype.name == 'uint32' for x in profiles])
        is_uint64 = all([x.dtype.name == 'uint64' for x in profiles])
        is_float32 = all([x.dtype.name == 'float32' for x in profiles])
        is_float64 = all([x.dtype.name == 'float64' for x in profiles])

    elif len(arguments.input) == 0 or (len(arguments.input) == 1 and (arguments.input[0] == "STDIN" or arguments.input[0] == "/dev/stdin")):

        feature = 1
        step = 1
        
        logger.log_it("Reading input as tsv/csv from STDIN", "INFO")
        try:
            df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.log_it(e.__str__(), "INFO")

            logger.log_it("Pandas error on DataFrame reading. Perhaps a null dataset being read?", "ERROR")
            raise e
        columns = list(df.columns)
    elif len(arguments.input) == 1 and (os.path.splitext(arguments.input[0])[-1] == ".tsv" or os.path.splitext(arguments.input[0])[-1] == ".csv"):

        feature = 1
        step = 1
        
        
        logger.log_it("Reading input file as tsv/csv...", "INFO")
        try:
            df = pd.read_csv(arguments.input[0], sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.log_it(e.__str__(), "ERROR")


            logger.log_it("Pandas error on DataFrame reading. Perhaps a null dataset being read?", "ERROR")

            
            raise e
        columns = list(df.columns)
    else:
        logger.log_it(arguments, "ERROR")

        logger.log_it(str(arguments.input), "ERROR")


        
        raise ArgumentError("kmerdb matrix received {0} arguments as input, and this is not supported.".format(len(arguments.input)))

    final_df = None
    if arguments.method == "DESeq2":
        
        feature = 2
        step = 4

        
        # if arguments.normalize_with == "ecopy":
        #     import ecopy as ep
        #     logger.info("Normalizing the DataFrame for sample size with ecopy...")
        #     normalized = ep.transform(df, method="normalize", axis=0)
        # # Normalizing across the 0th axis will normalize across the row meaning between samples. To arrive at this I checked the column sums of both to be sure I understood what was being normalized.
        # # We don't necessarily need to standardize, and not standardizing will make rationalizations or sensitivity to the real number range more clear.
        # if arguments.normalize_with == "DESeq2":

        logger.log_it("Normalizing the DataFrame with DESeq2 NB-compatible count normalization via rpy2...", "INFO")

        '''
        The following was shown to me by the following medium post.
        https://w1ndy.medium.com/calling-r-libraries-from-python-5ffbf3c3e5a8
        '''
        from rpy2.robjects.packages import importr, PackageNotInstalledError
        from rpy2.robjects import FloatVector, DataFrame, r, pandas2ri, Formula
        pandas2ri.activate()
        try:
            #fitdistrplus = importr('fitdistrplus')
            deseq = importr('DESeq2')
            count_matrix = pandas2ri.py2rpy(df)
            ncol = r['ncol']
            seq = r['seq']
            colData = pd.DataFrame(np.transpose(columns), columns=["Species"])
            
            logger.log_it("colData:\n{0}".format(colData), "INFO")

                
            dds = r['DESeqDataSetFromMatrix'](count_matrix, colData, Formula('~ Species'))
            logger.log_it("DESeq dataset:\n{0}".format(str(dds)), "INFO")
            
            dds = r['estimateSizeFactors'](dds)
                
            normalized = r['counts'](dds, normalized = True)
            if not arguments.no_normalized_ints:
                normalized = np.array(np.rint(normalized), dtype="int64")
            normalized = pd.DataFrame(normalized, columns=columns)
            
            logger.log_it("Normalized matrix shape: {0}".format(normalized.shape), "INFO")


            #     if arguments.method == "Unnormalized":
            #         cf = r['descdist'](r_dataframe, discrete=True)
            #     else:
            #         cf = r['descdist'](r_dataframe, discrete=False)
        except PackageNotInstalledError as e:
            logger.log_it(e.__str__(), "ERROR")

            logger.log_it("One or more R packages were not installed. Please see the install documentation about installing the associated R packages", "ERROR")

            logger.log_it("Use -vv for suggested installation instructions.", "ERROR")

            
            raise e

        #normalized.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)

        # logger.error("False ending of Normalized")
        # logger.debug("final_df should be set as normalized")
        # sys.exit(1)
        final_df = normalized
    elif arguments.method == "from":

        feature = 2
        step = 4
        
        #df.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)
        final_df = df
    elif arguments.method == "Frequency":

        feature = 2
        step = 4
        
        
        final_df = df # 4/5/24 - frequencies are given in the standard format
        #k = suggested_metadata['k']
        #total_kmers = suggested_metadata["total_kmers"]
        #final_df = df.div(total_kmers)
    elif arguments.method == "PCA":



        feature = 2
        step = 2


        
        # This method is actually dimensionality reduction via SVD, and this process is used during principal components analysis.
        # We generate the elbow graph in this step if the required dimensionality parameter '-n' is not supplied.
        # In this case we assume they have not provided the parameter because they'd like to use the so-called 'auto-detection'
        
        # I've specified a feature called 'auto-detect' such that the elbow graph is produced
        # they supply the number of dimensions they'd like after viewing the graph,
        # And in this way by viewing the graph, the 'correct' number of dimensions is attained.
        #
        #################################
        # PCA dimensionality reduction
        #################################
    
        logger.log_it("Performing preliminary dimensionality reduction to the MLE", "INFO")


        
        pca = PCA()#n_components="mle", svd_solver='auto')
        pca.fit(np.transpose(df)) # PCA of the normalized matrix or its transpose?
        # We want the number of k-mers, the number of features reduced, so we transpose the original matrix


        plt.plot(range(1, len(pca.explained_variance_ratio_)+1), pca.explained_variance_ratio_.cumsum(), marker='o', linestyle="--")
        plt.title("Explained variance by components")
        plt.xlabel("Number of components")
        plt.ylabel("Cumulative explained variance")
        plt.savefig(config.pca_variance_fig_filepath)

        if arguments.n is not None:
            logger.log_it("Using selected PCA dimensionality to reduce the transpose matrix/DataFrame again for use in 'kdb kmeans'", "INFO")
            pca = PCA(n_components=arguments.n)
            pca.fit(np.transpose(df))
            sys.stderr.write("\n\n\n")
            sys.stderr.write("-"*30+ "\n")
            
            sys.stderr.write("Explained variances: {0}\n".format(pca.explained_variance_ratio_))
            sys.stderr.write("Log-likelihoods: {0}\n".format(pca.score_samples(np.transpose(df))))
            sys.stderr.write("Log-likelihood of all samples: {0}\n".format(pca.score(np.transpose(df))))
            sys.stderr.write("MLE estimate of components for dimensionality reduction produced this shape: {0}\n".format(pca.components_.shape))

            score_matrix = pca.transform(np.transpose(df))
            score_df = pd.DataFrame(np.transpose(score_matrix), columns=columns)

            #score_df.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)
            final_df = score_df
        else:
            logger.log_it("You must look at '{0}' to decide on the choice of '-n' for specify when generating the dimensionally reduced matrix for clustering function".format(config.pca_variance_fig_filepath), "WARNING")


            
        #logger.debug("Components:\n{0}".format(pca.components_))
        # The shape for my first run was 10x11 and I have 11 samples. So I decided to not transpose this matrix to simply add it to a DataFrame and rename the samples.
        #df = pd.DataFrame(pca.components_, columns=list(map(lambda kdbrdr: os.path.splitext(os.path.basename(kdbrdr._filepath))[0], files)))
    elif arguments.method == "tSNE":


        feature = 2
        step = 3




        
        '''
        In t-SNE, perplexity or k is defined as 2^S, where S is the Shannon entropy of the conditional probability distribution.
        '''
        if arguments.n is None:
            raise TypeError("'kdb matrix tSNE' requires a keyword argument '-n' equal to the number of components of the subspace for tSNE to project the data into. A choice of 2 is recommended")
        tsne = TSNE(n_components=arguments.n, perplexity=arguments.perplexity).fit_transform(np.transpose(df))
        tsne_df = pd.DataFrame(np.transpose(tsne), columns=columns)
        #tsne_df.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)
        final_df = tsne_df

    ## FIXME: CUSTOM SORTING. WILL NOT COMMIT TO GIT REPO
    #suffixes = [(int(x.split("_")[1]), i) for i, x in enumerate(column_names)] # A list of a 2-tuple of the correct sort order and the index
    #suffixes.sort(key=lambda x: x[0])
    #sorted_column_names = [column_names[s[1]] for s in suffixes]
    #final_df.sort_values(sorted_column_names)

    final_df.to_csv(sys.stdout, sep=arguments.output_delimiter, index=arguments.with_index)

    
    logger.log_it("Done printing {0} matrix to STDOUT".format(arguments.method), "INFO")

    #logger.info("Beginning distribution analysis in R...")
    #logger.warn("Not implemented in Python...")
    
    sys.stderr.write(config.DONE)



def kmeans(arguments):
    import logging
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    from io import TextIOWrapper
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn import metrics
    from sklearn.cluster import KMeans
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    from kmerdb import config

    
    if isinstance(arguments.input, argparse.FileType) or isinstance(arguments.input, TextIOWrapper):
        df = pd.read_csv(arguments.input, sep=arguments.delimiter)
        column_names = df.columns
        df = df.transpose()
        #df.set_index('index')
    elif arguments.input is None or arguments.input == "/dev/stdin" or arguments.input == "STDIN":
        df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        column_names = df.columns
        df = df.transpose()
        #df.set_index('index')
    else:
        logger.log_it("An unknown IO type was detected in kmeans", "ERROR")
        raise IOError("Unable to process input argument")

        
    num_samples, num_features = df.shape
    logger.log_it("Input DataFrame shape: {0}".format(df.shape), "INFO")
    
    t0 = time.time()

    
    wcss = []
    for i in range(1, num_samples+1):
        kmeans_pca = KMeans(n_clusters = i, init="k-means++", random_state=42)
        kmeans_pca.fit(df)
        wcss.append(kmeans_pca.inertia_)

    plt.plot(range(1, num_samples+1), wcss, marker='o', linestyle='--')
    plt.xlabel("Number of Clusters")
    plt.ylabel("WCSS")
    plt.title("K-means with PCA")
    plt.savefig(config.kmeans_elbow_graph_fig_filepath)
    sys.stderr.write("="*42 + "\n\n")
    sys.stderr.write(" " * 10 + "The range 1:{0} was swept to find an elbow for the clustering\n".format(num_samples))
    sys.stderr.write("The figure was written to '{0}'\n\n".format(config.kmeans_elbow_graph_fig_filepath))
    sys.stderr.write("="*42 + "\n\n")

    
    if arguments.k is not None and arguments.method == "sklearn":
        final_kmeans = KMeans(n_clusters = arguments.k, init="k-means++", random_state=42)
        final_kmeans.fit(df)
        labels = final_kmeans.labels_
        df.insert(0, "Cluster", labels)
        #labels = None
        # final_df.replace({"Cluster": {
        #     1: "One"
        # }})
        df.to_csv(arguments.output, sep=arguments.output_delimiter)
        colnames = list(df.columns)
        fig, ax = plt.subplots()
        scatter = ax.scatter(df.iloc[:, 1], df.iloc[:, 2], c=df.iloc[:, 0])
        legend = ax.legend(*scatter.legend_elements(), loc="upper right", title="Cluster")

        print("K-means coordinates")
        for x,y,z in zip(np.array(df.iloc[:, 1]), np.array(df.iloc[:, 2]), column_names):
            
            print("{0}\t{1}\t{2}".format(x,y,z))
            #sys.stderr.write("{0}\t{1}\t{2}\n".format(x,y,z))
            ax.annotate(z, xy=(x, y), xycoords='data', xytext=(0, -45), textcoords='offset points', arrowprops=dict(facecolor='black', shrink=0.05), horizontalalignment='left', verticalalignment='bottom')
        ax.set_xlabel("Dim1")
        ax.set_ylabel("Dim2")
        ax.set_title("K-means clustering")
        ax.add_artist(legend)
        ax.grid(True)


        plt.savefig(config.kmeans_clustering_fig_filepath)
        
        #centroids = kmeans.cluster_centers_
        logger.log_it(str(list(column_names)), "INFO")
        
        logger.log_it(", ".join(list(map(str, labels))), "INFO")
    
        '''
        Clustering score quality estimation from 
        https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_digits.html#define-our-evaluation-benchmark
    
        '''
        name = "K-means metrics"
        fit_time = time.time() - t0
        results = [name, fit_time, final_kmeans.inertia_]

        clustering_metrics = [
            metrics.homogeneity_score,
            metrics.completeness_score,
            metrics.v_measure_score,
            metrics.adjusted_rand_score,
            metrics.adjusted_mutual_info_score,
        ]
        
        results += [m(labels, final_kmeans.labels_) for m in clustering_metrics]

        results += [metrics.silhouette_score(df, final_kmeans.labels_, metric="euclidean", sample_size=300,)]

        formatter_result = ("\n{:9s}\t{:.3f}s\t{:.0f}\t{:.3f}\t{:.3f}"
                            "\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n\n\n")
        sys.stderr.write("\n" + "\t".join(["Name", "Total_time", "Inertia", "Homogeneity_score", "Completeness_score", "V_measure_score", "Adj_rand_score", "Adj_mutual_info_score", "Silhouette"]) + "\n")
        sys.stderr.write(formatter_result.format(*results))

        #kmer_pca = pca(df, scale=False) # It's already been normalized (scaled) in get_matrix()
        #print(kmer.summary_rot())


    elif arguments.k is not None and arguments.method == "Biopython":
        from Bio import Cluster

        clust, error, nfound = Cluster.kcluster(df, nclusters=arguments.k, dist=arguments.distance)

        labels = list(clust)

        
        df.insert(0, "Cluster", labels)
        #labels = None
        # final_df.replace({"Cluster": {
        #     1: "One"
        # }})

        df.to_csv(arguments.output, sep=arguments.output_delimiter)
        fig, ax = plt.subplots()
        scatter = ax.scatter(df.iloc[:, 1], df.iloc[:, 2], c=df.iloc[:, 0])
        legend = ax.legend(*scatter.legend_elements(), loc="upper right", title="Cluster")
        for x,y,z in zip(np.array(df.iloc[:, 1]), np.array(df.iloc[:, 2]), column_names):
            sys.stderr.write("{0}\t{1}\t{2}".format(x,y,z))
            ax.annotate(z, xy=(x, y), xycoords='data', xytext=(0, -45), textcoords='offset points', arrowprops=dict(facecolor='black', shrink=0.05), horizontalalignment='left', verticalalignment='bottom')

        
        ax.set_xlabel("Dim1")
        ax.set_ylabel("Dim2")
        ax.set_title("K-means clustering")
        ax.add_artist(legend)
        ax.grid(True)
        plt.savefig(config.kmeans_clustering_fig_filepath)
        
        #centroids = kmeans.cluster_centers_
        logger.log_it(", ".join(list(column_names)), "INFO")

        logger.log_it(", ".join(labels), "INFO")

        
    sys.stderr.write("Generated {0} clusters projected onto reduced dimension 1 and 2 of the input dataset\n".format(arguments.k))
    sys.stderr.write("The figure was written to {0}\n\n".format(config.kmeans_clustering_fig_filepath))

    logger.log_it("The annotated sample matrix can be printed by specifying an output files with [ -o|--output OUTFILE ]", "INFO")
    
    sys.stderr.write(config.DONE)



def hierarchical(arguments):
    """
    Thanks to https://www.analyticsvidhya.com/blog/2019/05/beginners-guide-hierarchical-clustering/
    for the inspiration
    """
    import logging
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    from io import TextIOWrapper
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import scipy.cluster.hierarchy as shc

    from Bio.Cluster import treecluster
    
    from kmerdb import config


    
    if isinstance(arguments.input, argparse.FileType) or isinstance(arguments.input, TextIOWrapper):
        df = pd.read_csv(arguments.input, sep=arguments.delimiter)
        column_names = list(df.columns)
        df = df.transpose()
    elif arguments.input is None or arguments.input == "STDIN" or arguments.input == "/dev/stdin":
        df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        column_names = list(df.columns)
        df = df.transpose()
    else:
        logger.log_it("An unknown IO type was detected in hierarchical", "INFO")
        raise IOError("Unable to process input")


    
    num_samples, num_features = df.shape
    logger.log_it("Input DataFrame shape: {0}".format(df.shape), "INFO")

    #tree = treecluster(distancematrix=np.array(df)) # Can this be exported to Newick?
    plt.figure(figsize=(16, 9)) # figsize=(10, 7)
    plt.title("Dendrogram")
    dend = shc.dendrogram(shc.linkage(df, method=arguments.method), labels=column_names, leaf_rotation=90)
    #plt.show()
    plt.savefig(config.hierarchical_clustering_dendrogram_fig_filepath)

    sys.stderr.write("Saving the matplotlib dendrogram to '{0}'...".format(config.hierarchical_clustering_dendrogram_fig_filepath))

    data = np.array(df).tolist()

    n = len(data)
    m = len(data[0])
    for i in range(n):
        for j in range(m):
            if i < j:
                data[i].pop(-1)
    #print(data)
    
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from Bio import Phylo
    dm = DistanceMatrix(column_names, data)
    constructor = DistanceTreeConstructor()
    upgmatree = constructor.upgma(dm)
    print(upgmatree)
    Phylo.draw(upgmatree, branch_labels=lambda c: c.branch_length)
    Phylo.write(upgmatree, config.spearman_upgma_tree_phy, "phyloxml")
    sys.stderr.write("Saving the phylip format tree to '{0}'".format(config.spearman_upgma_tree_phy))
    sys.stderr.write(config.DONE)

# def rarefy(arguments):
#     logging.getLogger('matplotlib.font_manager').disabled = True
#     logging.getLogger('matplotlib').setLevel(logging.WARNING)
#     from io import TextIOWrapper
#     import numpy as np
#     import pandas as pd
#     import matplotlib.pyplot as plt

    
#     from ecopy.diversity import rarefy
#     import ecopy as ep
#     from kmerdb import config


#     if isinstance(arguments.input, argparse.FileType) or isinstance(arguments.input, TextIOWrapper):
#         df = pd.read_csv(arguments.input, sep=arguments.delimiter)
#         column_names = df.columns
#         df = df.transpose()
#     elif arguments.input is None:
#         df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
#         column_names = df.columns
#         df = df.transpose()
#     else:
#         logger.error("An unknown IO type was detected in bin/kdb.cluster()")
#         sys.exit(1)


#     num_samples, num_features = df.shape
#     t0 = time.time()

#     logger.warning("Multiple files: Untested behavior with ecopy.diversity.rarefy")
#     test = rarefy(df, config.ecopy_rarefaction_fig_filepath, 'rarecurve')
#     #arguments.output.write(test)
#     if arguments.output is None:
#         print(test)
#     else:
#         arguments.output.write(test)

#     sys.stderr.write(config.DONE)
        
        

def header(arguments):
    """
    Another end-user function that takes an argparse Namespace object.
    This function just reads the metadata header, can print in json.
    """
    from kmerdb import fileutil, config, util, graph

    sfx = os.path.splitext(arguments.kdb)[-1]
    metadata = None

    if sfx != ".kdb" and sfx != ".kdbg": # A filepath with invalid suffix
        raise IOError("Viewable .kdb(g) filepath does not end in '.kdb' or '.kdbg'")
    elif not os.path.exists(arguments.kdb):
        raise IOError("Input .kdb(g) filepath '{0}' does not exist on the filesystem".format(arguments.kdb_in))

    if sfx == ".kdb":
        kdb = fileutil.open(arguments.kdb, mode='r', sort=False, slurp=False)
        metadata = kdb.metadata
            
        kmer_ids_dtype = metadata["kmer_ids_dtype"]
        N = 4**metadata["k"]
        if metadata["version"] != config.VERSION:
            logger.log_it(".kdb version is out of date, may be incompatible with current KDBReader class", "WARNING")

        assert kdb.kmer_ids.size == N, "view | read kmer_ids size did not match N from the header metadata"
        assert kdb.counts.size == N, "view | read counts size did not match N from the header metadata"
        assert kdb.frequencies.size == N, "view | read frequencies size did not match N from the header metadata"
        metadata = kdb.metadata
    elif sfx == ".kdbg":
        kdb = graph.open(arguments.kdb, mode='r', slurp=False)
        if kdb.metadata["version"] != config.VERSION:
            logger.log_it(".kdb file version is out of date, may be incompatible with current fileutil.KDBReader class", "WARNING")


            
        N = 4**kdb.metadata["k"]
        metadata = kdb.metadata

    else:
        raise ValueError("Input to 'kmerdb header' is a .kdb or .kdbg (gzipped .tsv) file of a count vector or a graph. Requires YAML metadata header.Try creating a k-mer count profile with 'kmerdb profile' or edge list with 'kmerdb graph'")

    if arguments.json:
        print(dict(kdb.metadata))
    else:
        yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
        print(yaml.dump(metadata))
        print(config.header_delimiter)
            
def view(arguments):
    """
    Another end-user function that takes an argparse Namespace object as its only input.
    This function facilitates the primary view of the dataset. It is a flat file after all, so we just need
    to display data properly to the user. This function is responsible for opening a KDBReader object,
    slurping the contents into memory, assert format sanity checks, and then iterate over those contents to reproduce the file to stdout.
    This function also has --un-sort and --re-sort features for sanity checking the sort order. Un-sort will un-sort a file onto just the k-mer id
    The un-sort function will produce what 'kmerdb profile -k $k input.fa output.kdb' would produce without the --sorted flag.
    The re-sort function will re-sort the list lexically according to count, and display the results on stdout.
    The default action is to display the .kdb file exactly as it is, assuming that the parsing KDBReader stands up.


    :param arguments: argparse namespace
    :type arguments: 
    """
    
    import numpy as np
    from kmerdb import fileutil, config, util, kmer, graph
    import json
    metadata = None
    N = None
    def get_header(line, header):
        """
        A little helper recurrence function for grabbing the additions to the header.
        I don't really know if this is germane... Untested  function.
        """
        

        if line.rstrip() != config.end_header_line:
            header += line
            return header
        else:
            header_dict = yaml.safe_load(header)
            if type(header_dict) is not dict:
                logger.log_it("Tried to parse:\n{0}\n".format(header), "DEBUG")

                
                raise ValueError("Could not parse YAML formatted header")
            else:
                return header_dict

    assert type(arguments.kdb_in) is str, "kdb_in must be a str"
    sfx = os.path.splitext(arguments.kdb_in)[-1]




    
    if sfx != ".kdb" and sfx != ".kdbg": # A filepath with invalid suffix
        raise IOError("Viewable .kdb(g) filepath does not end in '.kdb' or '.kdbg'")
    elif not os.path.exists(arguments.kdb_in):
        raise IOError("Viewable .kdb(g) filepath '{0}' does not exist on the filesystem".format(arguments.kdb_in))
    if sfx == ".kdb":
        with fileutil.open(arguments.kdb_in, mode='r', sort=arguments.re_sort, slurp=True) as kdb_in:
            metadata = kdb_in.metadata

            
            kmer_ids_dtype = metadata["kmer_ids_dtype"]
            N = 4**metadata["k"]
            if metadata["version"] != config.VERSION:
                logger.log_it("KDB version is out of date, may be incompatible with current KDBReader class", "WARNING")
            if arguments.kdb_out is None or (arguments.kdb_out == "/dev/stdout" or arguments.kdb_out == "STDOUT"): # Write to stdout, uncompressed
                if arguments.header:
                    yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
                    print(yaml.dump(metadata, sort_keys=False))
                    print(config.header_delimiter)
            logger.log_it("Reading from file...", "INFO")

            
            
            try:
                if not arguments.un_sort and arguments.re_sort and metadata["sorted"] is True:
                    kmer_ids_sorted_by_count = np.lexsort((kdb_in.counts, kdb_in.kmer_ids))
                    reverse_kmer_ids_sorted_by_count = np.flipud(kmer_ids_sorted_by_count)
                    for i, idx in enumerate(kmer_ids_sorted_by_count):
                        kmer_id = kdb_in.kmer_ids[i]
                        logger.log_it("The first is an implicit row-index. The second is a k-mer id, then the counts and frequencies.", "DEBUG")
                        
                        logger.log_it("{0}\t{1}\t{2}\t{3}".format(i, kmer_id, kdb_in.counts[kmer_id], kdb_in.frequencies[kmer_id]), "DEBUG")


                        print("{0}\t{1}\t{2}\t{3}".format(i, kmer_id, kdb_in.counts[kmer_id], kdb_in.frequencies[kmer_id]))
                else:
                    for i, idx in enumerate(kdb_in.kmer_ids):
                        kmer_id = kdb_in.kmer_ids[idx]
                        logger.log_it("The row in the file should follow this order:", "DEBUG")
                        
                        logger.log_it("The first is an implicit row-index. The second is a k-mer id, then the counts and frequencies.", "DEBUG")
                        
                        logger.log_it("{0}\t{1}\t{2}\t{3}".format(i, kmer_id, kdb_in.counts[kmer_id], kdb_in.frequencies[kmer_id]), "DEBUG")

                        
                        try:
                            if arguments.un_sort is True:
                                assert kmer_id == idx, "view | kmer_id {0} didn't match the expected k-mer id.".format(idx, kmer_id)
                                assert i == kmer_id, "view | kmer_id {0} didn't match the implicit index {1}".format(idx, i)
                            else:
                                #logger.debug("Not sorting, so skipping assertion about profile (col1, column 2)")
                                pass
                        except AssertionError as e:
                            logger.log_it(e.__str__(), "WARNING")
                            
                            logger.log_it("K-mer id {0} will be printed in the {1} row".format(idx, i), "WARNING")
                            
                        #raise e
                        logger.log_it("{0} line:".format(i), "DEBUG")
                        
                        logger.log_it("=== = = = ======= =  =  =  =  =  = |", "INFO")

                        if arguments.un_sort is True:
                            if arguments.no_singletons is True:
                                if kdb_in.counts[idx] > 2:
                                    
                                    print("{0}\t{1}\t{2}\t{3}".format(i, idx, kdb_in.counts[idx], kdb_in.frequencies[idx]))
                            else:
                                print("{0}\t{1}\t{2}\t{3}".format(i, idx, kdb_in.counts[idx], kdb_in.frequencies[idx]))                                
                        else:
                            if arguments.no_singletons is True:
                                if kdb_in.counts[idx] > 2:
                                    
                                    print("{0}\t{1}\t{2}\t{3}".format(i, idx, kdb_in.counts[idx], kdb_in.frequencies[idx]))
                            else:
                                print("{0}\t{1}\t{2}\t{3}".format(i, idx, kdb_in.counts[idx], kdb_in.frequencies[idx]))
                        #print("{0}\t{1}\t{2}\t{3}".format(i, idx, kdb_in.counts[kmer_id], kdb_in.frequencies[kmer_id]))
                # I don't think anyone cares about the graph representation.
                # I don't think this actually matters because I can't figure out what the next data structure is.
                # Is it a Cypher query and creation node set?
                # I need to demonstrate a capacity for graph based learning.
                # (:-|X) The dread pirate roberts got me.
                # :)
            except BrokenPipeError as e:
                logger.log_it(e.__str__(), "ERROR")
                raise e
        if arguments.kdb_out is not None and arguments.compress: # Can't yet write compressed to stdout
            logger.log_it("Can't write .kdb to stdout! We need to use a Bio.bgzf filehandle.", "ERROR")
            raise IOError("Can't write .kdb to stdout! We need to use a Bio.bgzf filehandle.")

        elif arguments.kdb_out is not None and type(arguments.kdb_out) is not str:
            raise ValueError("Cannot write a file to an argument that isn't a string")
        elif arguments.kdb_out is not None and (os.path.exists(arguments.kdb_out) or not os.path.exists(arguments.kdb_out)):
            if os.path.exists(arguments.kdb_out):
                logger.log_it("Overwriting '{0}'...".format(arguments.kdb_out), "WARNING")

            logger.log_it("Creating '{0}'...".format(arguments.kdb_out), "DEBUG")
            
            if arguments.kdb_out is not None:
                with fileutil.open(arguments.kdb_in, 'r', dtype=suggested_dtype, sort=arguments.sorted, slurp=True) as kdb_in:
                    assert kdb_in.kmer_ids.size == N, "view | read kmer_ids size did not match N from the header metadata"
                    assert kdb_in.counts.size == N, "view | read counts size did not match N from the header metadata"
                    assert kdb_in.frequencies.size == N, "view | read frequencies size did not match N from the header metadata"
                    with fileutil.open(arguments.kdb_out, metadata=metadata, mode='w') as kdb_out:
                        try:
                            for i, idx in enumerate(kdb_in.kmer_ids):
                                kmer_id = idx
                                seq = kmer.id_to_kmer(kmer_id, arguments.k)
                                kmer_metadata = kmer.neighbors(seq, arguments.k)
                                logger.log_it("The first is the actual row id. This is the recorded row-id in the file. This should always be sequential. Next is the k-mer id. ", "DEBUG")
                                kdb_out.write("{0}\t{1}\t{2}\t{3}\n".format(i, kmer_id, kdb_in.counts[kmer_id],  kdb_in.frequencies[kmer_id], kmer_metadata))
                    
                        except StopIteration as e:
                            logger.log_it(e.__str__(), "ERROR")
                            raise e
                        finally:
                            #kdb_out._write_block(kdb_out._buffer)
                            #kdb_out._handle.flush()
                            #kdb_out._handle.close()
                            sys.stderr.write(config.DONE)
    elif sfx == ".kdbg":
        kdbg_in = graph.open(arguments.kdb_in, mode='r', slurp=True)
        metadata = kdbg_in.metadata

        

        n1_dtype      = metadata["n1_dtype"]
        n2_dtype      = metadata["n2_dtype"]
        weights_dtype = metadata["weights_dtype"]

            
        if metadata["version"] != config.VERSION:
            logger.log_it("KDB version is out of date, may be incompatible with current KDBReader class", "WARNING")
            
        if arguments.kdb_out is None or (arguments.kdb_out == "/dev/stdout" or arguments.kdb_out == "STDOUT"): # Write to stdout, uncompressed
            if arguments.header:
                yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
                print(yaml.dump(metadata, sort_keys=False))
                print(config.header_delimiter)
        logger.log_it("Reading from file...", "INFO")
        
        logger.log_it("I cut off the json-formatted unstructured column for the main view.", "DEBUG")
        

        for i in range(len(kdbg_in.n1)):
            n1 = kdbg_in.n1[i]
            n2 = kdbg_in.n2[i]
            w  = kdbg_in.weights[i]
            logger.log_it("The row in the file should follow this order:", "DEBUG")
            
            logger.log_it("The first is an implicit row-index. The second and third are k-mer ids, then edge weight", "DEBUG")
            
            logger.log_it("{0}\t{1}\t{2}\t{3}".format(i, n1, n2, w), "DEBUG")
            
            logger.log_it("{0} line:".format(i), "DEBUG")

            logger.log_it("=== = = = ======= =  =  =  =  =  = |", "DEBUG")

            
            print("{0}\t{1}\t{2}\t{3}".format(i, n1, n2, w))
                # I don't think anyone cares about the graph representation.
                # I don't think this actually matters because I can't figure out what the next data structure is.
                # Is it a Cypher query and creation node set?
                # I need to demonstrate a capacity for graph based learning.
                # (:-|X) The dread pirate roberts got me.
                # :)
        if arguments.kdb_out is not None and arguments.compress: # Can't yet write compressed to stdout
            logger.log_it("Can't write kdb to stdout! We need to use a Bio.bgzf filehandle.", "ERROR")
            raise IOError("Can't write kdb to stdout! We need to use a Bio.bgzf filehandle.")

        elif arguments.kdb_out is not None and type(arguments.kdb_out) is not str:
            raise ValueError("Cannot write a file to an argument that isn't a string")
        elif arguments.kdb_out is not None and os.path.exists(arguments.kdb_out):
            
            logger.log_it("Overwriting '{0}'...".format(arguments.kdb_out), "WARNING")

            
        elif arguments.kdb_out is not None and not os.path.exists(arguments.kdb_out):
            
            logger.log_it("Creating '{0}'...".format(arguments.kdb_out), "DEBUG")
            
            if arguments.kdb_out is not None:
                with graph.open(arguments.kdb_out, metadata=metadata, mode='w') as kdb_out:
                    try:
                        for i in range(len(kdbg.n1)):
                            kdb_out.write("{0}\t{1}\t{2}\t{3}\n".format(i, kdbg_in.n1[i], kdbg_in.n2[i],  kdbg_in.w[i]))
                    
                    except StopIteration as e:
                        logger.log_it(e.__str__(), "ERROR")
                        raise e
                    finally:
                        #kdb_out._write_block(kdb_out._buffer)
                        #kdb_out._handle.flush()
                        #kdb_out._handle.close()
                        sys.stderr.write(config.DONE)
    else:
        raise ValueError("Input files in kmerdb are tab delimited .csv files, either a count vector or a edge list (block gzip compression). Requires a YAML metadata header. Try making a k-mer count profile/vector with 'kmerdb profile -k 12 <input_1.fa>' ")

            
def assembly(arguments):
    from kmerdb import graph








    assert type(arguments.kdbg) is str, "kdbg must be a str"
    sfx = os.path.splitext(arguments.kdbg)[-1]




    
    if sfx != ".kdb" and sfx != ".kdbg": # A filepath with invalid suffix
        raise IOError("Viewable .kdb(g) filepath does not end in '.kdb' or '.kdbg'")
    elif not os.path.exists(arguments.kdbg):
        raise IOError("Viewable .kdb(g) filepath '{0}' does not exist on the filesystem".format(arguments.kdbg))
    elif sfx == ".kdbg":
        with graph.open(arguments.kdbg, mode='r', slurp=True) as kdbg_in:
            metadata = kdbg_in.metadata

            N = len(kdbg_in.n1)            
            n1_dtype      = metadata["n1_dtype"]
            n2_dtype      = metadata["n2_dtype"]
            weights_dtype = metadata["weights_dtype"]



            
            if metadata["version"] != config.VERSION:
                logger.log_it(".kdb version is out of date, may be incompatible with current KDBReader class", "WARNING")

                
            if arguments.kdb_out is None or (arguments.kdb_out == "/dev/stdout" or arguments.kdb_out == "STDOUT"): # Write to stdout, uncompressed
                if arguments.header:
                    yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
                    print(yaml.dump(metadata, sort_keys=False))
                    print(config.header_delimiter)
            logger.log_it("Reading from file...", "INFO")

            
            try:



                nodes = list(set(kdbg_in.n1 + kdbg_in.n2))
                edges = list(zip(kdbg_in.n1, kdbg_in.n2, list(map(lambda w: {'weight': w} , kdbg_in.weights))))

                graph.create_graph(nodes, edges)
                
            except BrokenPipeError as e:
                logger.log_it(e.__str__(), "ERROR")
                raise e

    
                            
                                              
def make_graph(arguments):
    """
    Another ugly function that takes a argparse Namespace object as its only positional argument

    Basically, from a fasta file, I want to generate a file format that consists of unique row ids, k-mer ids, adjacency list, 

    and finally an int field (0 represents unused) representing the order of the row-kmer being used in a graph traversal.

    Note, that the .kdbg format may not be easily regenerated from a k-mer count vector alone, and thus a .fasta/.fastq is needed.

    The goal isn't to yet implement the traversal algorithm in this first commit. But instead I'd like to get the format specified.
    """
    from multiprocessing import Pool


    import jsonschema
    
    import numpy as np
    
    from kmerdb import kmer, util, graph
    from kmerdb.config import VERSION

    global logger
    global step
    global feature
    

    
    logger.log_it("Printing entire CLI argparse option Namespace...", "DEBUG")

    logger.log_it(str(arguments), "DEBUG")
    
    # The extension should be .kdb because I said so.
    logger.log_it("Checking extension of output file...", "INFO")

    
    if os.path.splitext(arguments.kdbg)[-1] != ".kdbg":
        raise IOError("Destination .kdbg filepath does not end in '.kdbg'")

    file_metadata = []
    theoretical_kmers_number = 4**arguments.k # Dimensionality of k-mer profile



    
    logger.log_it("Parsing {0} sequence files to generate a k-mer adjacency list...".format(len(list(arguments.seqfile))), "INFO")

    
    infile = graph.Parseable(arguments, logger=logger) # NOTE: Uses the kmerdb.graph module


    feature = 1
    step = 1

    
    logger.log_it("Processing {0} fasta/fastq files across {1} processors...".format(len(list(arguments.seqfile)), arguments.parallel), "INFO")
    
    logger.log_it("Parallel (if specified) mapping the kmerdb.parse.parsefile() method to the seqfile iterable", "DEBUG")
    
    logger.log_it("In other words, running the kmerdb.parse.parsefile() method on each file specified via the CLI", "DEBUG")
    
    if arguments.parallel > 1:
        with Pool(processes=arguments.parallel) as pool:
            # data files_metadata
            data = pool.map(infile.parsefile, arguments.seqfile) # Returns a list of k-mer ids
    else:
            # data files_metadata
            data = list(map(infile.parsefile, arguments.seqfile))



    """
    ######################################
    Step 1. Completed
    ######################################
    """



    step = 2


            
    """
    Summary statistics and metadata structure
    """

    nullomer_ids = []
    counts = []
    for d in data:
        nullomer_ids.extend(d[3])
        counts.extend(d[2])
        file_metadata.append(d[1])
    N = 4**arguments.k

    
    all_observed_kmers = int(np.sum(counts))
        
    unique_nullomers = len(set(nullomer_ids))
    
    unique_kmers = int(np.count_nonzero(counts))
    
    # Key assertion
    assert unique_kmers + unique_nullomers == theoretical_kmers_number, "kmerdb | internal error: unique nullomers ({0}) + unique kmers ({1}) should equal 4^k = {2} (was {3})".format(unique_nullomers, unique_kmers, theoretical_kmers_number, unique_kmers + unique_nullomers)
    #logger.info("created a k-mer composite in memory")

    
    metadata=OrderedDict({
        "version": VERSION,
        "metadata_blocks": 1,
        "k": arguments.k,
        "total_kmers": all_observed_kmers,
        "unique_kmers": unique_kmers,
        "unique_nullomers": unique_nullomers,
        "sorted": arguments.sorted,
        "n1_dtype": "uint64",
        "n2_dtype": "uint64",
        "weights_dtype": "uint64",
        "tags": [],
        "files": file_metadata
    })

    
    logger.log_it("Validating aggregated metadata structure for .kdbg header...", "INFO")

    try:
        np.dtype(metadata["n1_dtype"])
        np.dtype(metadata["n2_dtype"])
        np.dtype(metadata["weights_dtype"])
    except TypeError as e:
        logger.log_it(metadata, "ERROR")
        
        logger.log_it(e.__str__(), "ERROR")
        
        logger.log_it("kmerdb encountered a type error and needs to exit", "ERROR")
        
        raise TypeError("Incorrect dtype detected in metadata structure. Internal error.")
    
    try:
        jsonschema.validate(instance=metadata, schema=config.graph_schema)
    except jsonschema.ValidationError as e:
        logger.log_it(e.__str__(), "ERROR")

        
        logger.log_it("kmerdb.graph.KDBGReader couldn't validate the header/metadata YAML", "ERROR")

        
        raise e

    logger.log_it("Validation complete.", "DEBUG")
    
    logger.log_it("\n\n\tCompleted summation and metadata aggregation across all inputs...\n\n", "INFO")



    """
    ACCUMULATE ALL EDGE WEIGHTS ACROSS ALL FILES
    """
    """
    task 1: initialize the final edge datastructure, a hashmap, keyed on a pair of k-mer ids, and containing only an integer weight
    """





    step = 3



    
    # Initialize empty data structures
    all_edges_in_kspace = {}
    n1 = []
    n2 = []
    weights = []
    # Loop over the input files and initialize the dictionary on the valid pairs, setting them to 0.
    for d in data:
        edges, h, counts, nullomers = d
        pair_ids = edges.keys()
        for p in pair_ids:
            all_edges_in_kspace[p] = 0
    """
    task 2: Accumulate all edge weights across all files
    """
    for d in data:
        edges, h, counts, nullomers = d
        pair_ids = edges.keys()
        
        for p in pair_ids:
            try:
                all_edges_in_kspace[p] += edges[p]
            except KeyError as e:
                logger.log_it("Unknown edge detected, evaded prepopulation. Internal error.", "ERROR")
                raise e
    """
    task 3: Add edges (2 k-mer ids) and weight to lists for conversion to NumPy array.
    """
    for e in all_edges_in_kspace.keys():
        n1.append(e[0])
        n2.append(e[1])
        weights.append(all_edges_in_kspace[e])


    """
    #################################################
    Step 2. Completed
    #################################################
    """













        
    """
    N would be the number of edges, pairs of nodes, and weights.
    """
    N = len(n1)
    
    logger.log_it("Initializing Numpy arrays of {0} uint for the edges and corresponding weight...".format(N), "DEBUG")
    n1 = np.array(n1, dtype=metadata["n1_dtype"])
    n2 = np.array(n2, dtype=metadata["n2_dtype"])
    weights = np.array(weights, dtype=metadata["weights_dtype"])
    logger.log_it("Initialization of profile completed, using approximately {0} bytes per array".format(n1.nbytes), "INFO")

    """
    Write the YAML metadata header to the .kdbg file
    """
    yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
    sys.stderr.write(yaml.dump(metadata, sort_keys=False))

    sys.stderr.write("\n\n\n")
    
    logger.log_it("Wrote metadata header to .kdbg...", "INFO")
    

    
    """
    Cause pandas isn't edge-y enough.
    """
    # Convert the list of numpy arrays (the uint64 3-tuple of the edge) to pandas dataframe
    #df = pd.DataFrame(twoD_weighted_edge_list)
    #df.to_csv(sys.stdout, sep=arguments.output_delimiter, index=False)


    """
    Write the dataset (weighted edge list) to a file with '.kdbg' as its suffix.
    """


    feature = 2
    step = 4
    
    kdbg_out = graph.open(arguments.kdbg, mode='wb', metadata=metadata)
    
    try:
        sys.stderr.write("\n\n\nWriting edge list to {0}...\n\n\n".format(arguments.kdbg))
        for i, node1 in enumerate(n1):
            
            node2 = n2[i]
            w = weights[i]

            tupley = (i, node1, node2, w)
            tupley_dl = np.array(tupley, dtype="uint64")
            if arguments.quiet is False:
                print("\t".join(list(tupley_dl)))
                #print("{0}\t{1}\t{2}\t{3}".format(i, node1, node2, w))
            # i, node1, node2, weight
            kdbg_out.write("{0}\n".format("\t".join(list(tupley_dl))))
            #kdbg_out.write("{0}\t{1}\t{2}\t{3}\n".format(i, node1, node2, w))
    finally:
        kdbg_out._write_block(kdbg_out._buffer)
        kdbg_out._handle.flush()
        kdbg_out._handle.close()


        """
        Done around n birfday
        3/20/24        

        Final statistics to stderr
        """
        sys.stderr.write("\n\n\nFinal stats:\n\n\n")
        
        sys.stderr.write("Total k-mers processed: {0}\n".format(all_observed_kmers))
        sys.stderr.write("Unique nullomer count:   {0}\n".format(unique_nullomers))
        sys.stderr.write("Unique {0}-mer count:     {1}\n".format(arguments.k, unique_kmers))
        sys.stderr.write("Theoretical {0}-mer number (4^{0}):     {1}\n".format(arguments.k, theoretical_kmers_number))
        sys.stderr.write("="*30 + "\n")
        sys.stderr.write(".kdbg stats:\n")
        sys.stderr.write("-"*30 + "\n")
        sys.stderr.write("Edges in file:  {0}\n".format(N))
        sys.stderr.write("Non-zero weights: {0}\n".format(int(np.count_nonzero(weights))))
        sys.stderr.write("Sum of node degrees in file: {0}".format(2*N))
        sys.stderr.write("\nDone\n")

    logger.log_it("Done printing weighted edge list to .kdbg", "INFO")

    sys.stderr.write(config.DONE)

    
            
def profile(arguments):
    """
    A complex, near-end user function that handles a argparse Namespace as its only positional argument

    This function handles multiprocessing, NumPy type checking and array initialization, full metadata expansion if needed.

    It also manages count aggregation across multiple fasta/fastq files as default.
    
    The function produces a composite profile from multiple inputs, possibly in parallel,

    and provides functionality to produce a .kdb format file, while writing the same information (the full k-mer profile to stdout)

    This behavior can be suppressed with the --quiet CLI option.

    This function is also responsible for the sorting of output files.

    The help menu for this function is your friend.

    :param arguments: argparse Namespace
    :type arguments:
    """

    import copy
    
    logger.log_it("Printing entire CLI argparse option Namespace...", "DEBUG")
    logger.log_it(str(arguments), "DEBUG")

    # The number of specified processors should be available
    logger.log_it("Validating processor count for parallelization...", "DEBUG")
    if arguments.parallel <= 0 or arguments.parallel > cpu_count()+1:
        raise argparse.ArgumentError("-p, --parallel must be a valid processor count.")

    ## 6/11/24 removed because reasons
    # The extension should be .kdb because I said so.
    # logger.log_it("Checking extension of output file...", "INFO")
    # if os.path.splitext(arguments.kdb)[-1] != ".kdb":
    #     raise IOError("Destination .kdb filepath does not end in '.kdb'")
    samples = []


    
    if len(arguments.input) == 1:
        logger.log_it("Input filename is '{0}'".format(arguments.input[0]), "INFO")
        if ".fastq" in arguments.input[0] or ".fq" in arguments.input[0] or ".fastq.gz" in arguments.input[0] or ".fq.gz" in arguments.input[0]:
            logger.log_it("Input suffix is .fastq", "INFO")
        elif ".fasta" in arguments.input[0] or ".fa" in arguments.input[0] or ".fna" in arguments.input[0] or ".fasta.gz" in arguments.input[0] or ".fa.gz" in arguments.input[0]:
            logger.log_it("Input suffix is .fasta", "INFO")
        elif ".txt" in arguments.input[0] or ".tsv" in arguments.input[0]:
            logger.log_it("Input suffix is .txt, possibly samplesheet. will open as tsv", "INFO")
            # One sample per line
            samplesheet = arguments.input[0]
            with open(samplesheet, 'r') as ifile:
                for line in ifile:
                    sample = line.rstrip()
                    if os.access(sample, os.R_OK):
                        samples.append(sample)
                    else:
                        logger.log_it("Error while processing samplesheet '{0}'...".format(samplesheet), "ERROR")
                        raise ValueError("Couldn't open sample file '{0}' for reading".format(sample))
            arguments.input = samples
        else:

            raise ValueError("Could not determine file type of input")        
    else:
        raise ValueError("Could not determine POSIX access mode for one or more input files.")



    new_args = copy.deepcopy(arguments)


    
    if arguments.k is None:
        if arguments.minK is None or arguments.maxK is None:
            raise ValueError("In multi-k mode, arguments --minK and --maxK are required.")
        elif arguments.minK < 4 or arguments.maxK < 5 or arguments.minK > 25 or arguments.maxK > 30:
            raise ValueError("Valid --min|max between 3 and 30")
        else:
            logger.log_it("Doing things...", "INFO")
            for k in range(arguments.minK, arguments.maxK+1):

                new_args.k = k
                #raise RuntimeError("Multi-k mode")
                _profile(new_args)
                #_profile(arguments) # FIXME
                pass
    elif arguments.k is not None:
        logger.log_it("Running in single-k mode", "INFO")

        #raise RuntimeError("single-k mode")
        _profile(new_args)

    
    return 




def _profile(arguments):
    from multiprocessing import Pool
    import json
    import time
    import numpy as np
    from kmerdb import parse, fileutil, kmer, util
    from kmerdb.config import VERSION

    global logger
    global feature
    global step


    step = 0
    feature = 0
    
    file_metadata = []
    total_kmers = 4**arguments.k # Dimensionality of k-mer profile
    N = total_kmers
    theoretical_kmers_number = N

    logger.log_it("Parsing {0} sequence files to generate a composite k-mer profile...".format(len(list(arguments.input))), "INFO")
    nullomers = set()
    # Initialize arrays
    
    # kmerdb.parse.Parseable

    infile = parse.Parseable(arguments, logger=logger) #

    feature += 1
    step += 1

    
    logger.log_it("Processing {0} fasta/fastq files across {1} processors...".format(len(list(arguments.input)), arguments.parallel), "INFO")
    
    logger.log_it("Parallel (if specified) mapping the kmerdb.parse.parsefile() method to the seqfile iterable", "DEBUG")
    
    logger.log_it("In other words, running the kmerdb.parse.parsefile() method on each file specified via the CLI", "DEBUG")

    """
    This is .fa/.fq parsing. It is delegated to the parse.Parseable class for argument intake, which then refers to the module-level method
    'parsefile', not to be confused with parse.Parseable.parsefile, which simply executes the parsefile function on behalf of the wrapper class
    'Parseable' such that the function may be used with 'multiprocessing.Pool'
    """

    
    
    if arguments.parallel > 1:
        with Pool(processes=arguments.parallel) as pool:
            data = pool.map(infile.parsefile, arguments.input)
    else:
        data = list(map(infile.parsefile, arguments.input))



                       
        
    # the actual 'data' is now a list of 4-tuples
    # Each 4-tuple represents a single file
    # (counts, header_dictionary<dict>, nullomers<list>, list)

    # Construct a final_counts array for the composite profile across all inputs



    

    
    step += 1
    counts = np.zeros(N, dtype="uint64")
    # Complete collating of counts across files
    # This technically uses 1 more arrray than necessary 'final_counts' but its okay
    logger.log_it("Summing counts from individual fasta/fastq files into a composite profile...", "INFO")
    

    for d in data:
        counts = counts + d[0]
        file_metadata.append(d[1])

        
    if np.sum(counts) == 0:
        raise ValueError("Each element of the array of k-mer counts was 0. Sum of counts == 0. Likely an internal error. Exiting.")


    
    assert np.sum(counts) > 0, "cowardly refusing to print Array of all 0's. More likely an internal error."
        
        
        
    sys.stderr.write("\n\n\tCompleted summation and metadata aggregation across all inputs...\n\n")
    # unique_kmers = int(np.count_nonzero(counts))
    # #total_nullomers = total_kmers - unique_kmers

    # nullomer_ids = list(map(lambda n: n[2], data))
    # all_observed_kmers = int(np.sum(list(map(lambda h: h['total_kmers'], file_metadata))))



    
    """
    3/20/24
    *Massive* nullomer and count validation regression.
    """

    """
    Summary statistics and metadata structure
    """
    step +=1
    feature += 1
    
    # nullomer_ids = []
    # #counts = []
    # file_metadata = []
    # for d in data:
    #     nullomer_ids.extend(d[2])
    #     counts.extend(d[0])
    #     file_metadata.append(d[1])

    #     # Extract the *new* logs from parse.parsefile
    #     newlogs = d[3]
    #     for line in newlogs:
    #         try:
    #             logger.logs.index(line)
    #         except ValueError as e:
    #             logger.logs.append(line)

    
    all_observed_kmers = int(np.sum(counts))
    unique_kmers = int(np.count_nonzero(counts))
    
    unique_nullomers = theoretical_kmers_number - unique_kmers
    #unique_nullomers = len(set(nullomer_ids))

    from kmerdb import lexer


    
    no_singletons = []
    hist = util.get_histo(list(counts))

    #print(hist)

    
    if arguments.show_hist is True:
        sys.stderr.write("Full histogram:\n")
        sys.stderr.write("[{0}]".format(", ".join(list(map(lambda x: str(x), hist)))))


    i, cov = lexer.max(hist)
    kmer_coverage = cov
    
    logger.log_it("Theoretical k-mer number: {0} | {1}".format(N, theoretical_kmers_number), "DEBUG")
    logger.log_it("Length of count array: {0}".format(counts.size), "DEBUG")
    logger.log_it("Number of non-zeroes: {0}".format(unique_kmers), "DEBUG")
    logger.log_it("Number of nullomers: {0}".format(unique_nullomers), "DEBUG")
    
    # Key assertion

    """
    7/15/24 Okay so the assertions are working fine, but something else is getting garbled and the counts are coming out all zero. 
    """
    
    assert unique_kmers + unique_nullomers == theoretical_kmers_number, "kmerdb | internal error: unique nullomers ({0}) + unique kmers ({1}) should equal 4^k = {2} (was {3})".format(unique_nullomers, unique_kmers, theoretical_kmers_number, unique_kmers + unique_nullomers)
    #logger.info("created a k-mer composite in memory")
    

    # nullomer_ids = []

    # for ns in file_metadata:
    #     nullomer_ids.extend(ns["nullomer_array"])
    
    # unique_nullomers = len(set(nullomer_ids))
    # # 
    # unique_kmers = int(np.sum(list(map(lambda h: h["unique_kmers"], file_metadata))))

    # assert unique_kmers + unique_nullomers == N, "kmerdb | internal error: unique nullomers ({0}) + unique kmers ({1}) should equal 4^k = {2} (was {3})".format(unique_nullomers, unique_kmers, N, unique_kmers + unique_nullomers)
    
    # all_observed_kmers = int(np.sum(counts))

    logger.log_it("Initial counting process complete, creating BGZF format file (.kdb)...", "INFO")
    
    logger.log_it("Formatting master metadata dictionary...", "INFO")
    


    
    metadata=OrderedDict({
        "version": VERSION,
        "metadata_blocks": 1,
        "k": arguments.k,
        "total_kmers": all_observed_kmers,
        "unique_kmers": unique_kmers,
        "unique_nullomers": unique_nullomers,
        "metadata": False,
        "sorted": arguments.sorted,
        "kmer_coverage": kmer_coverage,
        "kmer_ids_dtype": "uint64",
        "profile_dtype": "uint64",
        "kmer_coverage_histogram_dtype": "uint64",
        "count_dtype": "uint64",
        "frequencies_dtype": "float64",
        "tags": [],
        "files": file_metadata
    })

    """
    Validate Python datatypes listed in metadata header by calling the dtype function of the module of the NumPy 1.21.2
    """
    try:
        np.dtype(metadata["kmer_ids_dtype"])
        np.dtype(metadata["profile_dtype"])
        np.dtype(metadata["count_dtype"])
        np.dtype(metadata["frequencies_dtype"])
    except TypeError as e:
        logger.log_it(e.__str__(), "ERROR")

        raise TypeError("Incorrect dtype.")

    logger.log_it("Initializing Numpy array of {0} uint zeroes for the final composite profile...".format(total_kmers), "DEBUG")
    
    kmer_ids = np.array(range(N), dtype=metadata["kmer_ids_dtype"])
    profile = np.array(range(N), dtype=metadata["profile_dtype"])
    counts = np.array(counts, dtype=metadata["count_dtype"])
    hist_ = list(hist)
    hist = np.array(hist, dtype=metadata["kmer_coverage_histogram_dtype"])

    frequencies = np.divide(counts, metadata["total_kmers"])


    


    logger.log_it("Initialization of profile completed, using approximately {0} bytes for profile".format(counts.nbytes), "INFO")
    
    yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
    sys.stderr.write(yaml.dump(metadata, sort_keys=False))

    
    step += 1

    """
    Print histogram as 0-100, with the 0's and 1's reinserted to specification of histogram
    """



    # Output takes the form --output-name + '.' + $k + '.kdb'
    output_filepath = "{0}.{1}.kdb".format(arguments.output_name, arguments.k)

    logger.log_it("Collapsing the k-mer counts across the various input files into the .kdb file '{0}'".format(output_filepath), "INFO")    
    kdb_out = fileutil.open(output_filepath, 'wb', metadata=metadata)
    
    try:
        sys.stderr.write("\n\nWriting outputs to {0}...\n\n".format(output_filepath))


        # Numpy indexing is trash
        #kmer_ids = np.zeros(total_kmers, dtype=metadata["kmer_ids_dtype"])
        #profile = np.zeros(total_kmers, dtype=metadata["profile_dtype"])
        #counts = np.zeros(total_kmers, dtype=metadata["count_dtype"])
        #frequencies = np.zeros(total_kmers, dtype=metadata["frequencies_dtype"])



        sys.stderr.write("Full histogram :\n{0}\n".format(list(map(lambda x: int(x), hist))))

        if arguments.quiet is not True and arguments.show_hist is True:
            print(hist_[0:100])
            print("Here's the (mini-) histogram you asked for!")



        
        i, m = lexer.max(hist)

        if arguments.quiet is not True:
            print("idx i: {0}, current maximum: {1}".format(i, m))
        sys.stderr.write("idx i: {0}, current maximum: {1}".format(i, m))


        #print(hist_)
        #sys.stderr.write("mini_hist: {0}".format(hist_))
        #kdb_out.write(list(hist)[0:100])
        #time.sleep(240)

        if arguments.sorted:

            kmer_ids_sorted_by_count = np.lexsort(kmer_ids)
            reverse_kmer_ids_sorted_by_count = list(kmer_ids_sorted_by_count)
            reverse_kmer_ids_sorted_by_count.reverse()

            
            logger.log_it("K-mer id sort example: {0}".format(reverse_kmer_ids_sorted_by_count[:30]), "INFO")

            
            for i, idx in enumerate(reverse_kmer_ids_sorted_by_count):

                kmer_id = int(kmer_ids[idx])
                seq = kmer.id_to_kmer(kmer_id, arguments.k)
                
                logger.log_it("{0}\t{1}\t{2}\t{3}".format(i, kmer_ids[idx], counts[idx], frequencies[idx]), "INFO")

                
                c[idx] = counts[idx]
                f[idx] = frequencies[idx]
                if arguments.quiet is not True:
                    print("{0}\t{1}\t{2}\t{3}".format(i, kmer_id, c, f))
                kdb_out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(i, kmer_id, c, f))
        else:
            for i, idx in enumerate(kmer_ids):


                kmer_id = int(kmer_ids[idx])
                seq = kmer.id_to_kmer(kmer_id, arguments.k)
                c = counts[idx]
                f = frequencies[idx]

                logger.log_it("{0}\t{1}\t{2}\t{3}".format(i, kmer_id, c, f), "INFO")
                
                if arguments.quiet is not True:
                    print("{0}\t{1}\t{2}\t{3}".format(i, kmer_id, c, f))
                kdb_out.write("{0}\t{1}\t{2}\t{3}\n".format(i, kmer_id, c, f))
                
        logger.log_it("Wrote 4^k = {0} k-mer counts + neighbors to the .kdb file.".format(total_kmers), "INFO")

            
    finally:
        kdb_out._write_block(kdb_out._buffer)
        kdb_out._handle.flush()
        kdb_out._handle.close()

        """
        Done around n birfday
        3/20/24        

        Final statistics to stderr
        """
        sys.stderr.write("\n\n\nFinal stats:\n\n\n")
        sys.stderr.write("Total k-mers processed: {0}\n".format(all_observed_kmers))
        sys.stderr.write("Unique nullomer count:   {0}\n".format(unique_nullomers))
        sys.stderr.write("Unique {0}-mer count:     {1}\n".format(arguments.k, unique_kmers))
        sys.stderr.write("Estimated genome size: N/A\n".format(all_observed_kmers/kmer_coverage))
        sys.stderr.write("K-mer coverage:  {0}\n".format(kmer_coverage))
        sys.stderr.write("Theoretical {0}-mer number (4^{0}):     {1}\n".format(arguments.k, theoretical_kmers_number))
        sys.stderr.write("="*30 + "\n")
        sys.stderr.write("\nDone\n")


    
        


def citation_info():
    citation = None

    MODULE_ROOT = os.path.dirname(__file__)
    citation_file = os.path.join(MODULE_ROOT,  'CITATION.txt')
    if os.access(citation_file, os.R_OK):
        with open(citation_file, 'r') as citation_f:
            citation = citation_f.read().rstrip()

        if citation == "":
            return
        else:
            sys.stderr.write("Printing citation notice to stderr. This will not interfere with the execution of the program in any way. Please see CITATION_FAQ.md for any questions.\n")
            sys.stderr.write(citation + "\n\n\n")
            sys.stderr.write("Run 'kmerdb citation' to silence.\n")
    else:
        sys.stderr.write("Thanks for using/citing kmerdb\n")
    # else:
    #     raise IOError("Cannot locate the extra package data file 'kmerdb/CITATION', which should have been distributed with the program")


def get_program_header(arguments):
    import appmap
    import sys


    argv = sys.argv

    kmerdb_appmap = appmap.kmerdb_appmap( argv )


    kmerdb_appmap.print_program_header()


    return kmerdb_appmap

def cli():

    import sys

    from kmerdb import config, appmap

    
    argv = sys.argv

    
    sys.stderr.write("Running kdb script from '{0}'\n".format(__file__))
    sys.stderr.write("Checking installed environment...\n")
    primary_path = sys.path[0]
    ver_info = sys.version_info
    py_version = "python{0}.{1}".format(ver_info.major, ver_info.minor)
    #sys.path.append(os.path.join(primary_path, "../lib/{0}/site-packages".format(py_version)))
    #sys.path.append(os.path.join(os.path.dirname(__file__), "../lib"))
    #site_packages = [p for p in sys.path if "kdb" in p]
    
    #sys.stderr.write("PYTHONPATH={0}".format(sys.path))
    #sys.path.remove(os.path.dirname(os.path.abspath(__file__)))
    citation_info()


    
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(help="Use --help with sub-commands")


    profile_parser = subparsers.add_parser("profile", help=appmap.command_1_description)
    profile_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")

    profile_parser.add_argument("-k",  type=int, help="Choose k-mer size")
    profile_parser.add_argument("--minK", type=int, help="Minimum k for output k-mer profiles")
    profile_parser.add_argument("--maxK", type=int, help="Maximum k for output k-mer profiles")
    profile_parser.add_argument("-o", "--output-name", type=str, required=True, help="File name pattern for outputs. e.g.: example with underscores and dashes input_samplename_1 => input_samplename_1.11.kdb, input_samplename_1.12.kdb")


    profile_parser.add_argument("-p", "--parallel", type=int, default=1, help="Shred k-mers from reads in parallel")    
    # profile_parser.add_argument("--batch-size", type=int, default=100000, help="Number of updates to issue per batch to PostgreSQL while counting")
    # profile_parser.add_argument("-b", "--fastq-block-size", type=int, default=100000, help="Number of reads to load in memory at once for processing")
    #profile_parser.add_argument("-n", type=int, default=1000, help="Number of k-mer metadata records to keep in memory at once before transactions are submitted, this is a space limitation parameter after the initial block of reads is parsed. And during on-disk database generation")
    profile_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help="Number of logged lines to print to stderr. Default: 50")
    profile_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    #profile_parser.add_argument("--keep-db", action="store_true", help=argparse.SUPPRESS)
    #profile_parser.add_argument("--both-strands", action="store_true", default=False, help="Retain k-mers from the forward strand of the fast(a|q) file only")
    
    #profile_parser.add_argument("--all-metadata", action="store_true", default=False, help="Include read-level k-mer metadata in the .kdb")
    profile_parser.add_argument("--show-hist", action="store_true", default=False, help="Print the histogram of counts to stderr.")
    profile_parser.add_argument("--sorted", action="store_true", default=False, help="Sort the output kdb file by count")
    profile_parser.add_argument("--quiet", action="store_true", default=False, help="Do not log the entire .kdb file to stdout")
    profile_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    #profile_parser.add_argument("--sparse", action="store_true", default=False, help="Whether or not to store the profile as sparse")
    # Seqfile todo

                          
    
    profile_parser.add_argument("input", nargs="+", type=str, help="A plain text samplesheet: <samplesheet.txt|.fa|.fq.gz> - one filepath per line. OR one or more input .fasta or .fastq (.gz supported) files.")
    
    # profile_parser.add_argument("seqfile", nargs="+", type=str, metavar="<.fasta|.fastq>", help="Fasta or fastq files")
    # profile_parser.add_argument("kdb", type=str, help="Kdb file")
    profile_parser.set_defaults(func=profile)

    graph_parser = subparsers.add_parser("graph", help=appmap.command_2_description)
    graph_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")

    graph_parser.add_argument("-k", default=12, type=int, help="Choose k-mer size (Default: 12)", required=True)

    graph_parser.add_argument("-p", "--parallel", type=int, default=1, help="Parallel file reading only.")
    graph_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help="Number of logged lines to print to stderr. Default: 50")
    graph_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    #graph_parser.add_argument("-b", "--fastq-block-size", type=int, default=100000, help="Number of reads to load in memory at once for processing")
    #graph_parser.add_argument("--both-strands", action="store_true", default=False, help="Retain k-mers from the forward strand of the fast(a|q) file only")
    graph_parser.add_argument("--quiet", action="store_true", default=False, help="Do not list all edges and neighboring relationships to stderr")
    graph_parser.add_argument("--sorted", action="store_true", default=False, help=argparse.SUPPRESS)
    graph_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    #graph_parser.add_argument("--sparse", action="store_true", default=False, help="Whether or not to store the profile as sparse")

    graph_parser.add_argument("seqfile", nargs="+", type=str, metavar="<.fasta|.fastq>", help="Fasta or fastq files")
    graph_parser.add_argument("kdbg", type=str, help=".kdbg file")
    graph_parser.set_defaults(func=make_graph)

    # assembly_parser = subparsers.add_parser("assemble", help="Use NetworkX (and/or cugraph) to perform deBruijn graphs")
    # assembly_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    # assembly_parser.add_argument("-g", "--gpu", action="store_true", default=False, help="Utilize GPU resources (requires CUDA library cugraph)")
    # assembly_parser.add_argument("kdbg", type=str, help=".kdbg file")
    # assembly_parser.set_defaults(func=assembly)

    minimizer_parser = subparsers.add_parser("minimizers", help="Index an array of minimizers (selected-or-not) to associate with a position in reference")

    minimizer_parser.add_argument("fasta", type=str, help="A fasta file (.fa) to generate minimizers from")
    minimizer_parser.add_argument("kdb", type=str, help="A .kdb file to generate minimizers")
    #minimizer_parser.add_argument("-k", type=int, help="Choice of k for minimizer-selection", required=True)
    minimizer_parser.add_argument("-w", "--window-size",  type=int, help="Choice of window size (in base-pairs) to select minimizers with", required=True)
    
    minimizer_parser.add_argument("-a", "--include-unselected-minimizers", action="store_true", default=False, help="Include minimizer positions that were not selected as minimizers (length: (N-k+1)/window_size)")
    minimizer_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    minimizer_parser.add_argument("--quiet", action="store_true", default=False, help="Do not produce additional stdout")
    minimizer_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    minimizer_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    minimizer_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help=argparse.SUPPRESS)
    minimizer_parser.set_defaults(func=get_minimizers)

    alignment_parser = subparsers.add_parser("alignment", help="Create a Smith-Waterman-like alignment using a reference fasta, its minimizers, and query sequence, and its minimizers.")
    alignment_parser.add_argument("reference", type=str, help="Reference sequences in .fa/.fna/.fasta format")
    alignment_parser.add_argument("reference_kdb", type=str, help="Reference requires kmer count vector in .kdb format for minimizer selection to seed alignment")
    alignment_parser.add_argument("query", type=str, help="Query sequences in .fa/.fna/.fasta format. Not expected to be .kdb or indexed")
    alignment_parser.add_argument("-n", "--num-mins", type=int, help="Number of minimizers to seed an alignment")
    alignment_parser.add_argument("--match-score", type=int, default=3, help="Score of extending an SW alignment from minimizer seed match by 1bp (int| default: 3) ")
    alignment_parser.add_argument("--mismatch-score", type=int, default=-1, help="Mismatch penalty for SW alignment (int: -1)")

    
    #alignment_parser.add_argument("--gap-opening", type=int, default=-10, help="Affine gap opening penalty (int| default: -10)")
    #alignment_parser.add_argument("--gap-extend", type=int, default=-1, help="gap mismatch extension penalty (int| default: -1)")
    
    alignment_parser.add_argument("-w", "--window-size", type=int, help="Window size for sliding window in minimizer selection.", required=True)
    
    alignment_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    alignment_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    alignment_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    alignment_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help=argparse.SUPPRESS)
    alignment_parser.set_defaults(func=get_alignments)
    
    usage_parser = subparsers.add_parser("usage", help="provide expanded usage information on parameters and functions provided")
    usage_parser.add_argument("method", type=str, choices=("usage", "help", "profile", "graph", "index", "shuf", "matrix", "distance"), help="Print expanded usage statement")
    usage_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    usage_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    usage_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    usage_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help=argparse.SUPPRESS)
    usage_parser.set_defaults(func=expanded_help)

    help_parser = subparsers.add_parser("help", help="provide expanded help section on parameters and functions provided")
    help_parser.add_argument("-m", "--method", type=str, choices=("usage", "help", "profile", "graph", "index", "shuf", "matrix", "distance"), required=True, help="Print expanded usage statement")
    help_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    help_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    help_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    help_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help=argparse.SUPPRESS)
    help_parser.set_defaults(func=expanded_help)
    
    
    view_parser = subparsers.add_parser("view", help="View the contents of the .kdb or .kdbg file")
    view_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    view_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    view_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    view_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    view_parser.add_argument("-H", "--header", action="store_true", help="Include header in the output")
    view_parser.add_argument("--re-sort", action="store_true", help="Sort the k-mer profile before displaying") # FIXME
    view_parser.add_argument("--un-sort", action="store_true", help="Unsort the k-mer profile before displaying") # FIXME 
    view_parser.add_argument("--dtype", type=str, default="uint64", help="Read in the profiles as unsigned integer 64bit NumPy arrays")
    view_parser.add_argument("-d", "--decompress", action="store_true", help="Decompress input? DEFAULT: ")
    view_parser.add_argument("-c", "--compress", action="store_true", help="Print compressed output")
    view_parser.add_argument("--no-singletons", action="store_true", help="Do not print nullomers or singletons")
    
    #view_parser.add_argument("-d", "--decompress", action="store_true", help="Decompress input? DEFAULT: ")
    #view_parser.add_argument("-c", "--compress", action="store_true", help="Print compressed output")
    view_parser.add_argument("kdb_in", type=str, nargs="?", default=None, help="A k-mer database file (.kdb) to be read (Optional)")
    view_parser.add_argument("kdb_out", type=str, nargs="?", default=None, help="A k-mer database file (.kdb) to be written to (Optional)")
    view_parser.set_defaults(func=view)

    header_parser = subparsers.add_parser("header", help="Print the YAML header of the .kdb or .kdbg file and exit")
    header_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    header_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    header_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    header_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    header_parser.add_argument("-j", "--json", help="Print as JSON. DEFAULT: YAML")
    header_parser.add_argument("kdb", type=str, help="A k-mer database file (.kdb)")
    header_parser.set_defaults(func=header)


    matrix_parser = subparsers.add_parser("matrix", help=appmap.command_3_description)
    matrix_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    matrix_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    matrix_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    matrix_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help="Number of logged lines to print to stderr. Default: 50")

    matrix_parser.add_argument("-p", "--parallel", type=int, default=1, help="Read files in parallel")
    matrix_parser.add_argument("--with-index", default=False, action="store_true", help="Print the row indices as well")
    matrix_parser.add_argument("--column-names", default=None, type=str, help="A filepath to a plaintext flat file of column names.")
    matrix_parser.add_argument("--delimiter", default="\t", type=str, help="The choice of delimiter to parse the input .tsv with. DEFAULT: '\t'")
    matrix_parser.add_argument("--output-delimiter", type=str, default="\t", help="The output delimiter of the final csv/tsv to write. DEFAULT: '\t'")


    #matrix_parser.add_argument("--normalize-with", type=str, choices=["ecopy", "DESeq2"], default="DESeq2", help="Normalize with which method? DEFAULT: DESeq2")
    matrix_parser.add_argument("--no-normalized-ints", action="store_true", default=False, help="Don't round normalized counts to the nearest integer")
    matrix_parser.add_argument("-k", default=None, type=int, help="The k-dimension that the files have in common")
    matrix_parser.add_argument("-n", default=None, type=int, help="The number of dimensions to reduce with PCA or t-SNE. DEFAULT: an elbow graph will be generated if -n is not provided to help the user choose -n")

    matrix_parser.add_argument("--perplexity", default=5, type=int, help="The choice of the perplexity for t-SNE based dimensionality reduction")
    matrix_parser.add_argument("method", choices=["PCA", "tSNE", "DESeq2", "from"], default=None, help="Choice of dimensionality reduction, normalization method (DESeq2), or matrix-from (collate data only)")
    matrix_parser.add_argument("input", nargs="*", default=[], metavar="<kdbfile1 kdbfile2 ...|input.tsv|STDIN>", help="Two or more .kdb files, or another count matrix in tsv/csv")
    matrix_parser.set_defaults(func=get_matrix)
    
    # rarefy_parser = subparsers.add_parser("rarefy", help="Generate rarefaction information using ecopy.diversity.rarefy for the supplied .kdb files")
    # rarefy_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    # rarefy_parser.add_argument("-d", "--delimiter", default="\t", type=str, help="The choice of delimiter to parse the DataFrame with")
    # rarefy_parser.add_argument("--output-delimiter", type=str, default="\t", help="The output delimiter of the final csv/tsv to write.")
    # rarefy_parser.add_argument("-i", "--input", type=argparse.FileType("r"), default=None, help="The input reduced dimension or simply normalized matrix to use with K-means clustering")
    # rarefy_parser.add_argument("-o", "--output", type=argparse.FileType("w"), default=None, help="THe output csv/tsv of rarefied data")
    # rarefy_parser.set_defaults(func=rarefy)

    kmeans_parser = subparsers.add_parser("kmeans", help=appmap.command_7_description)
    kmeans_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    kmeans_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    kmeans_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help="Number of logged lines to print to stderr. Default: 50")
    kmeans_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    kmeans_parser.add_argument("-d", "--delimiter", type=str, default="\t", help="The delimiter of the input csv/tsv to parse, presumably produced by 'kdb matrix'.")
    kmeans_parser.add_argument("--output-delimiter", type=str, default="\t", help="The output delimiter of the final csv/tsv to write.")
    kmeans_parser.add_argument("--distance", type=str, default='e', choices=['e', 'b', 'c', 'a', 'u', 'x', 's', 'k'], help="Different distance metrics offered by kcluster. It is highly advised that you check both this source and their documentation to see how this is implemented.")

    
    kmeans_parser.add_argument("-k", default=None, type=int, help="The choice of k for clustering", required=True)
    kmeans_parser.add_argument("-i", "--input", type=argparse.FileType("r"), default=None, help="The input reduced dimension or mereley normalized matrix to use with K-means clustering")
    kmeans_parser.add_argument("-o", "--output", type=argparse.FileType("w"), default=None, help="The output csv/tsv with added 'label' to use for graphing in R, if the matplotlib graphs are not sufficient.")
    kmeans_parser.add_argument("method", choices=["sklearn", "Biopython"], default="Biopython", help="The Python implementation of k-means to use. The Biopython method is selected for access to alternative distance metrics")
    kmeans_parser.set_defaults(func=kmeans)

    hierarchical_parser = subparsers.add_parser("hierarchical", help=appmap.command_8_description)
    hierarchical_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    hierarchical_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    hierarchical_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help="Number of logged lines to print to stderr. Default: 50")
    hierarchical_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")

    hierarchical_parser.add_argument("-d", "--delimiter", type=str, default="\t", help="The delimiter to use when reading the csv.")
    hierarchical_parser.add_argument("-i", "--input", type=argparse.FileType("r"), default=None, help="The input distance matrix for hierarchical clustering")
    hierarchical_parser.add_argument("-m", "--method", type=str, choices=["single", "complete", "average", "weighted", "centroid", "median", "ward"], default="ward", help="The method for linkage fitting to use")
    hierarchical_parser.set_defaults(func=hierarchical)
    
    dist_parser = subparsers.add_parser("distance", help=appmap.command_4_description)
    dist_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    dist_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    dist_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")        
    dist_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    dist_parser.add_argument("--output-delimiter", type=str, default="\t", help="The output delimiter of the final csv/tsv to write.")
    dist_parser.add_argument("-p", "--parallel", type=int, default=1, help="Read files in parallel")
    dist_parser.add_argument("--column-names", type=str, default=None, help="A filepath to a plaintext flat file of column names.")
    dist_parser.add_argument("--delimiter", type=str, default="\t", help="The delimiter to use when printing the csv.")
    dist_parser.add_argument("-k", default=None, type=int, help="The k-dimension that the files have in common")
    
    dist_parser.add_argument("metric", choices=[
        "braycurtis",
        "canberra",
        "chebyshev",
        "cityblock",
        "correlation",
        "cosine",
        "dice",
        "euclidean",
        "hamming",
        "jaccard",
        "jensenshannon",
        "kulsinski",
        "mahalanobis",
        "matching",
        "minkowski",
        "pearson",
        "rogerstanimoto"
        "russelrao",
        "seuclidean",
        "sokalmichener",
        "sokalsneath",
        "spearman",
        "sqeuclidean",
        "yule"], default="correlation", help="Choice of distance metric between two profiles")
    dist_parser.add_argument("input", nargs="*", default=[], metavar="<kdbfile1 kdbfile2 ...|input.tsv|STDIN>", help="Two or more .kdb files, or another count matrix in tsv/csv")
    dist_parser.set_defaults(func=distances)


    index_parser = subparsers.add_parser("index", help=appmap.command_9_description)

    index_parser.add_argument("--force", action="store_true", help="Force index creation (if previous index exists")
    index_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    index_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    index_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    index_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)

    
    index_parser.add_argument("kdb", type=str, help="A k-mer database file (.kdb)")
    index_parser.set_defaults(func=index_file)

    shuf_parser = subparsers.add_parser("shuf", help=appmap.command_10_description)
    shuf_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    shuf_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    shuf_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    shuf_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    shuf_parser.add_argument("kdb_in", type=str, help="An indexed k-mer database file (.kdb)")
    shuf_parser.add_argument("kdb_out", type=str, help="The output shuffled k-mer database file (.kdb)")
    shuf_parser.set_defaults(func=shuf)


    version_parser = subparsers.add_parser("version", help=appmap.command_13_name)
    version_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    version_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    version_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    version_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    version_parser.set_defaults(func=version)



    
    # markov_probability_parser = subparsers.add_parser("probability", help=u"""
# Calculate the log-odds ratio of the Markov probability of a given sequence from the product (pi) of the transition probabilities(aij) times the frequency of the first k-mer (P(X1)), given the entire k-mer profile of a species.

# See https://matthewralston.github.io/quickstart#kmerdb-probability for more details.

# 1. Durbin, R., Eddy, S.R., Krogh, A. and Mitchison, G., 1998. Biological sequence analysis: probabilistic models of proteins and nucleic acids. Cambridge university press.
# """)

    # markov_probability_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    # markov_probability_parser.add_argument("-d", "--delimiter", type=str, default="\t", help="The delimiter to use when reading the csv.")
    # markov_probability_parser.add_argument("-b", "--fastq-block-size", type=int, default=100000, help="Number of reads to load in memory at once for processing")
    # markov_probability_parser.add_argument("seqfile", type=str, metavar="<.fasta|.fastq>", default=None, help="Sequences to calculate standard Markov-chain probabilities from, either .fasta or .fastq")
    # markov_probability_parser.add_argument("kdb", type=str, help="An indexed k-mer database file (.kdb)")
    # markov_probability_parser.set_defaults(func=markov_probability)

    citation_parser = subparsers.add_parser("citation", help="Silence the citation notice on further runs")
    citation_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    citation_parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Do not format errors and condense log")
    citation_parser.add_argument("-nl", "--num-log-lines", type=int, choices=config.default_logline_choices, default=50, help=argparse.SUPPRESS)
    citation_parser.add_argument("-l", "--log-file", type=str, default="kmerdb.log", help="Destination path to log file")
    citation_parser.set_defaults(func=citation)

    args=parser.parse_args()
    
    global logger
    global exit_code

    
    global step
    global feature

    
    global logs



    sys.stderr.write("Constructed a logger for the program...\n")
    #logger.debug(sys.path)

    # Print program header
    sys.stderr.write("Starting the program run-time timer...\n\n\n")
    start = time.time()


    """
    Print detailed debugging information prior to program log.
    """


    #signal.signal(signal.SIGINT, graceful_interrupt)
    #signal.signal(signal.SIGTERM, graceful_termination)


    # Extract the __init__ function invoked from argparsed arguments.
    function_name = vars(args)['func'].__name__

    # Now, map to the correct subcommand name via config

    if function_name in ["usage", "help", "citation"]:
        subcommand_name = function_name
    else:
        subcommand_name = config.subcommands[config.subcommand_functions.index(function_name)]



    from kmerdb import config
    logger = kdbLogger.Loggah(logfile=args.log_file or None, level=args.verbose)
        

    from kmerdb import appmap
    kmerdb_appmap = appmap.kmerdb_appmap( argv , logger )


    
    kmerdb_appmap.print_program_header()
    sys.stderr.write("Beginning program...\n")
    kmerdb_appmap.print_verbosity_header()
    
    if args.debug is True:
        args.func(args)
    else:
        logger.log_it("Running with error summary feature enabled, bypass with --debug for vague/invalid exceptions", "DEBUG")
        
        try:
            args.func(args)
            
            exit_code = 0

        except TypeError as e:
            
            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = 1
        
            raise e
        except ValueError as e:
        
            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = 2
        
            raise e
        except KeyError as e:
        
            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = 3
        
            raise e
        except RuntimeError as e:
        
            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = 4
        
            raise e
        except OSError as e:

            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = 5
        
            raise e
        except IOError as e:
            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = 6
            raise e
        
        except ArgumentError as e:

            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = 7
        
            raise e
        except AssertionError as e:

            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = 8
        
            raise e
        except FileNotFoundError as e:
            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            
            exit_code = 9
            raise e
        except Exception as e:
            exit_summary = kmerdb_appmap.exit_gracefully(e , subcommand=subcommand_name, step=step, feature=feature, logs=logger.logs, n_logs=args.num_log_lines or None)
            exit_code = -1
        
            raise e
        finally:
            sys.stderr.write("Program ran for {0} seconds...\n\n\n".format(time.time() - start))
            sys.stderr.write(config.thanks)
            sys.stderr.write(config.DONE)
            sys.exit(exit_code)
