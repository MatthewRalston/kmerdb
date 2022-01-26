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
from multiprocessing import cpu_count
from collections import OrderedDict


from Bio import bgzf

#import concurrent.futures

global logger
logger = None


def print_argv():
    argv = sys.argv
    sys.stderr.write(" ".join(argv[0:4]) + " ...\n")

def citation(arguments):
    import pkg_resources
    citation = None
    if pkg_resources.resource_exists('kmerdb', 'CITATION'):
        citation_fname = pkg_resources.resource_filename('kmerdb', 'CITATION')
        with open(citation_fname, 'w') as citation_file:
            citation_file.write("")


    sys.stderr.write("On the real, we gotta eat.")
    sys.stderr.write("Consider a +1 to keep it real...")

def index_file(arguments):
    from kmerdb import fileutil, index
    from kmerdb.config import DONE
    
    with fileutil.open(arguments.kdb, mode='r') as kdb:
        k = kdb.metadata['k']
    line_index = index.open(arguments.kdb, mode="xt", k=k)
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
    
def markov_probability(arguments):
    """
    A very old function that is probably broken in the indexing function.
    
    :param arguments: argparse Namespace
    :type arguments:
    """
    import pandas as pd
    import numpy as np
    from kmerdb import fileutil, index, probability, seqparser
    from kmerdb.config import DONE

    if os.path.splitext(arguments.kdb)[-1] != ".kdb":
        raise IOError("Model .kdb filepath does not end in '.kdb'")


    if index.has_index(arguments.kdb):
        arguments.kdbi = arguments.kdb + "i"
        #df = pd.DataFrame([], columns=["SequenceID", "Log_Odds_ratio", "p_of_seq"])
        profiles = np.array([], dtype="int64")
        with fileutil.open(arguments.kdb, 'r', slurp=True) as kdb:
            k = kdb.metadata['k']
            with index.open(arguments.kdbi, 'r') as kdbi:
                with seqparser.SeqParser(arguments.seqfile, arguments.fastq_block_size, k) as seqprsr:
                    recs = [r for r in seqprsr]
                    if seqprsr.fastq:
                        logger.debug("Read exactly b=={0}=={1} records from the {2} seqparser object".format(b, len(recs), s))
                        assert len(recs) == b, "The seqparser should return exactly {0} records at a time".format(b)
                    else:
                        logger.debug("Read {0} sequences from the {1} seqparser object".format(len(recs), seqprsr))
                        logger.debug("Skipping the block size assertion for fasta files")

                    while len(recs): # While the seqprsr continues to produce blocks of reads
                        # Do something here
                    
                        markov_probs = list(map(lambda p: [p["seq"].name, p["log_odds_ratio"], p["p_of_seq"]], [probability.markov_probability(seq, kdb, kdbi) for seq in recs]))

                        sys.stderr.write(json.dumps(markov_probs))
                        if profiles.shape == (0,):
                            profiles = np.array(markov_probs)
                        else:
                            np.append(profiles, markov_probs, axis=0)

                        recs = [r for r in seqprsr] # Essentially accomplishes an iteration in the file, wrapped by the seqparser.SeqParser class
        df = pd.DataFrame(profiles, columns=["SequenceID", "Log_Odds_ratio", "p_of_seq"])
        df.to_csv(sys.stdout, sep=arguments.delimiter, index=False)

    else:
        raise IndexError(".kdb file '{0}' has no corresponding index file. Please run 'kdb index -h' for details on index generation".format(arguments.kdb))

    
    sys.stderr.write(DONE)
                    
def distances(arguments):
    """
    An end-user function to provide CLI access to certain distances
    """
    from multiprocessing import Pool

    import pandas as pd
    import numpy as np
    from scipy.spatial.distance import pdist, squareform


    from kmerdb import fileutil, distance, config
    n = len(arguments.input)

    if len(arguments.input) > 1:
        file_reader = fileutil.FileReader()
        files = list(map(lambda f: fileutil.open(f, 'r', slurp=False), arguments.input))
        
        logger.debug("Files: {0}".format(files))
        if not all(os.path.splitext(kdb)[-1] == ".kdb" for kdb in arguments.input):
            raise IOError("One or more parseable .kdb filepaths did not end in '.kdb'")
        dtypes = [x.profile.dtype.name for x in files]
        suggested_dtype = dtypes[0]
        
        ks = [kdbrdr.k for kdbrdr in files]
        suggested_k = ks[0]
        if not all(kdbrdr.k == suggested_k for kdbrdr in files):
            logger.error("Files: {0}".format(files))
            logger.error("Choices of k: {0}".format(ks))
            logger.error("By default the inferred value of k is used when k is not specified at the command line, which was {0}".format(suggested_k))
            raise TypeError("One or more files did not have k set to be equal to {0}".format(arguments.k))
        logger.debug("Files: {0}".format(files))
        
        if not all(kdbrdr.dtype == suggested_dtype for kdbrdr in files):
            raise TypeError("One or more files did not have dtype = {0}".format(suggested_dtype))
        if arguments.parallel > 1:
            with Pool(processes=arguments.parallel) as pool:
                files = pool.map(file_reader.load_file, arguments.input)
        else:
            files = list(map(lambda f: fileutil.open(f, 'r', slurp=True), arguments.input))


        data = [kdbrdr.slurp() for kdbrdr in files]
        logger.info(data)
        pure_data = np.array(data, dtype=suggested_dtype)
        profiles = np.transpose(pure_data)

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
        logger.debug("Shape: {0}".format(profiles.shape))
        logger.info("Converting arrays of k-mer counts into a pandas DataFrame...")
        #columns = list(df.columns) # I hate this language, it's hateful.
        n = len(columns)
    elif len(arguments.input) == 1 and (arguments.input[0] == "STDIN" or arguments.input[0] == "/dev/stdin"):
        logger.info("Reading input as tsv/csv from STDIN")
        try:
            df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        profiles = np.transpose(np.array(df))
        columns = list(df.columns)
        assert profiles.shape[0] == len(columns), "distance | Number of columns of dataframe did not match the number of rows of the profile"
        n = len(columns)

    elif len(arguments.input) == 1 and (os.path.splitext(arguments.input[0])[-1] == ".tsv" or os.path.splitext(arguments.input[0])[-1] == ".csv"):
        try:
            df = pd.read_csv(arguments.input[0], sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        profiles = np.transpose(np.array(df))
        columns = list(df.columns) # I'm sorry I ever made this line. Please forgive me.
        # This is just gratuitous code and language. I'm really really not sure what I want to express here.
        assert profiles.shape[0] == len(columns), "distance | Number of columns of dataframe did not match the number of rows of the profile"
        n = len(columns)
    elif len(arguments.input) == 1 and os.path.splitext(arguments.input[0])[-1] == ".kdb":
        logger.error("Not sure why you'd want a singular distance.")
        logger.error("kdb distance requires more than one .kdb file as positional inputs")
        sys.exit(1)
    else:
        logger.error("bin/kdb.distances() received {0} arguments as input, which were not supported.".format(len(arguments.input)))
        sys.exit(1)

    # Masthead stuff
    sys.stderr.write(config.DEFAULT_MASTHEAD)
    if logger.level == logging.DEBUG or logger.level == logging.INFO:
        sys.stderr.write(config.DEBUG_MASTHEAD)
    sys.stderr.write(config.DISTANCE_MASTHEAD)

    logger.info("Custom calculating a {0}x{0} '{1}' distance matrix...".format(n, arguments.metric))
    
    if arguments.metric in ["pearson", "spearman"]:
        
        #files = list(map(lambda f: fileutil.open(f, 'r', slurp=True), arguments.input))
        data = [['' for x in range(n)] for y in range(n)]
        for i in range(n):
            for j in range(n):
                logger.debug("Calculating the {0}x{1} cell".format(i, j))
                logger.info("Calculating {0} distance between {1}:({2}) and {3}:({4}) columns...".format(arguments.metric, columns[i], i, columns[j], j))

                if i == j:
                    logger.info("Was self, reverting to identity.")
                    data[i][j] = distance.identity[arguments.metric]
                elif i > j:
                    data[i][j] = None
                    logger.info("Refusing to recompute custom distance metric for symmetric matrix")
                elif i < j:
                    #print("Info: ixj {0}x{1}")
                    logger.debug("Info: ixj {0}x{1}".format(i, j))
                    
                    if arguments.metric == "pearson":

                        logger.error("{0}".format(type(profiles[i])))
                        logger.info("Computing custom Pearson correlation coefficient...")
                        data[i][j] = distance.correlation(profiles[i], profiles[j])
                    elif arguments.metric == "euclidean":
                        logger.error("Custom Euclidean distance is DEPRECATED")
                        raise RuntimeError("Genuine bug, please report the usage of deprecated custom distance functions")
                        logger.info("Computing custom euclidean distance")
                        data[i][j] = distance.euclidean(profiles[i], profiles[j])
                    elif arguments.metric == "spearman":
                        cor, pval = distance.spearman(profiles[i], profiles[j])
                        logger.info("Computing SciPy Spearman distance")
                        # FIXME! also print pval matrices
                        data[i][j] = cor
                    elif arguments.metric == "EMD":
                        logger.error("Custom EMD distance is deprecated")
                        raise RuntimeError("Genuine bug, please report the usage of deprecated custom distance functions")
                        data[i][j] = distance.EMD(profiles[i], profiles[j])
                    elif arguments.metric == "d2s":
                        logger.error("Custom d2s distance is deprecated")
                        raise RuntimeError("Genuine bug, please report the usage of deprecated custom distance functions")
                        data[i][j] = distance.d2s(profiles[i], profiles[j])
                    else:
                        logger.error("Other distances are not implemented yet")
                        sys.exit(1)
                # This double loop quickly identifies empty cells and sets the data correctly from the permutation above
        logger.info("Filling in symmetrical matrix...")
        for i in range(n):
            for j in range(n):
                if data[i][j] is None:
                    data[i][j] = data[j][i]
                    #logger.info("\n\n\nSwerve swerve\n\n\n")
        dist = np.array(data)
    else:
        dist = pdist(np.transpose(profiles), metric=arguments.metric)
        dist = squareform(dist)
        # if arguments.metric == "correlation":
        #     ones = np.ones(dist.shape, dtype="int64")
        #     dist = np.subtract(ones, dist)
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
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    from multiprocessing import Pool
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    
    from kmerdb import fileutil, config


    
    if len(arguments.input) > 1:

        file_reader = fileutil.FileReader()
        files = list(map(lambda f: fileutil.open(f, 'r', slurp=False), arguments.input))

        if arguments.k is None:
            arguments.k = files[0].k
        if not all(os.path.splitext(kdb)[-1] == ".kdb" for kdb in arguments.input):
            raise IOError("One or more parseable .kdb filepaths did not end in '.kdb'")
        ks = [kdbrdr.k for kdbrdr in files]
        suggested_k = ks[0]
        if not all(k == suggested_k for k in ks):
            logger.error("Files: {0}".format(files))
            logger.error("Choices of k: {0}".format(ks))
            logger.error("By default the inferred value of k is used when k is not specified at the command line, which was {0}".format(arguments.k))
            logger.error("Default choice of k, since not specified at the CLI is {0}".format(arguments.k))
            logger.error("Proceeding with k set to {0}".format(suggested_k))
            raise TypeError("One or more files did not have k set to be equal to {0}".format(arguments.k))
        dtypes = [kdbrdr.count_dtype for kdbrdr in files]
        suggested_dtype = dtypes[0]
        if not all(kdbrdr.count_dtype == suggested_dtype for kdbrdr in files):
            raise TypeError("One of more files did not have dtype = {0}".format(suggested_dtype))
        if arguments.parallel > 1:
            with Pool(processes=arguments.parallel) as pool:
                files = pool.map(file_reader.load_file, arguments.input)
        else:
            files = list(map(lambda f: fileutil.open(f, 'r', slurp=True), arguments.input))

        data = [kdbrdr.slurp() for kdbrdr in files]
        logger.info(data)
        pure_data = np.array(data, dtype=suggested_dtype)
        profiles = np.transpose(pure_data)
        # The following does *not* transpose a matrix defined as n x N=4**k
        # n is the number of independent files/samples, (fasta/fastq=>kdb) being assessed
        # N=4**k is the dimensionality of the vector, sparse or not, that makes up the perceived profile, the descriptors is the column dimension.

        # The names for the columns are taken as the basenamed filepath, with all extensions (substrings beginning with '.') stripped.
        logger.info("============================================================")
        logger.info("Converting arrays of k-mer counts into a pandas DataFrame...")
        logger.info("============================================================")



        if arguments.column_names is None:
            columns = list(map(lambda kdbrdr: os.path.basename(kdbrdr._filepath).split(".")[0], files))
        else:
            with open(arguments.column_names, 'r') as column_names:
                columns = [line.rstrip() for line in column_names]
        if len(columns) != len(files):
            raise RuntimeError("Number of column names {0} does not match number of input files {1}...".format(len(columns), len(files)))
        suggested_metadata = files[0].metadata
        expected, num_columns = 4**suggested_k, len(columns)
        if profiles.shape != (expected, num_columns):
            #logger.info("Time is money.")
            logger.error("Check shape self...")
            logger.error("Expected shape: {0} x {1}".format(expected, num_columns))
            logger.error("Actual shape: {0}".format(profiles.shape))
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
        logger.info("Reading input as tsv/csv from STDIN")
        try:
            df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        columns = list(df.columns)
    elif len(arguments.input) == 1 and (os.path.splitext(arguments.input[0])[-1] == ".tsv" or os.path.splitext(arguments.input[0])[-1] == ".csv"):
        logger.info("Reading input file as tsv/csv")
        try:
            df = pd.read_csv(arguments.input[0], sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        columns = list(df.columns)
    else:
        logger.error(arguments)
        logger.error("bin/kdb.get_matrix() received {0} arguments as input, and this is not supported.".format(len(arguments.input)))
        logger.error(arguments.input)
        sys.exit(1)
    sys.stderr.write(config.DEFAULT_MASTHEAD)
    if logger.level == logging.DEBUG or logger.level == logging.INFO:
        sys.stderr.write(config.DEBUG_MASTHEAD)
    sys.stderr.write(config.MATRIX_MASTHEAD)


    final_df = None
    if arguments.method == "DESeq2":
        # if arguments.normalize_with == "ecopy":
        #     import ecopy as ep
        #     logger.info("Normalizing the DataFrame for sample size with ecopy...")
        #     normalized = ep.transform(df, method="normalize", axis=0)
        # # Normalizing across the 0th axis will normalize across the row meaning between samples. To arrive at this I checked the column sums of both to be sure I understood what was being normalized.
        # # We don't necessarily need to standardize, and not standardizing will make rationalizations or sensitivity to the real number range more clear.
        # if arguments.normalize_with == "DESeq2":

        logger.info("Normalizing the DataFrame with DESeq2 NB-compatible count normalization via rpy2...")
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
            logger.info("colData:\n{0}".format(colData))
                
            dds = r['DESeqDataSetFromMatrix'](count_matrix, colData, Formula('~ Species'))
            logger.info("DESeq dataset:\n{0}".format(str(dds)))
            dds = r['estimateSizeFactors'](dds)
                
            normalized = r['counts'](dds, normalized = True)
            if not arguments.no_normalized_ints:
                normalized = np.array(np.rint(normalized), dtype="int64")
            normalized = pd.DataFrame(normalized, columns=columns)
            logger.debug(normalized)
            logger.info("Normalized matrix shape: {0}".format(normalized.shape))

            #     if arguments.method == "Unnormalized":
            #         cf = r['descdist'](r_dataframe, discrete=True)
            #     else:
            #         cf = r['descdist'](r_dataframe, discrete=False)
        except PackageNotInstalledError as e:
            logger.error(e)
            logger.error("One or more R packages were not installed. Please see the install documentation about installing the associated R packages")
            logger.error("Use -vv for suggested installation instructions.")
            sys.exit(1)
        #normalized.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)

        # logger.error("False ending of Normalized")
        # logger.debug("final_df should be set as normalized")
        # sys.exit(1)
        final_df = normalized
    elif arguments.method == "pass":
        #df.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)
        final_df = df
    elif arguments.method == "Frequency":
        k = suggested_metadata['k']
        total_kmers = suggested_metadata["total_kmers"]
        final_df = df.div(total_kmers)
    elif arguments.method == "PCA":
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
    
        logger.info("Performing preliminary dimensionality reduction to the MLE")
        pca = PCA()#n_components="mle", svd_solver='auto')
        pca.fit(np.transpose(df)) # PCA of the normalized matrix or its transpose?
        # We want the number of k-mers, the number of features reduced, so we transpose the original matrix


        plt.plot(range(1, len(pca.explained_variance_ratio_)+1), pca.explained_variance_ratio_.cumsum(), marker='o', linestyle="--")
        plt.title("Explained variance by components")
        plt.xlabel("Number of components")
        plt.ylabel("Cumulative explained variance")
        plt.savefig(config.pca_variance_fig_filepath)

        if arguments.n is not None:
            logger.info("Using selected PCA dimensionality to reduce the transpose matrix/DataFrame again for use in 'kdb kmeans'")
            pca = PCA(n_components=arguments.n)
            pca.fit(np.transpose(df))
            #logger.debug("Explained variances: {0}".format(pca.explained_variance_ratio_))
            #logger.debug("Log-likelihoods: {0}".format(pca.score_samples(normalized)))
            #logger.debug("Overall log-likelihood of all samples: {0}".format(pca.score(normalized)))
            #logger.debug("MLE estimate of components for dimensionality reduction produced this shape: {0}".format(pca.components_.shape))

            score_matrix = pca.transform(np.transpose(df))
            score_df = pd.DataFrame(np.transpose(score_matrix), columns=columns)

            #score_df.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)
            final_df = score_df
        else:
            logger.warning("You must look at '{0}' to decide on the choice of '-n' for specify when generating the dimensionally reduced matrix for clustering function".format(config.pca_variance_fig_filepath))
            sys.exit(1)
        #logger.debug("Components:\n{0}".format(pca.components_))
        # The shape for my first run was 10x11 and I have 11 samples. So I decided to not transpose this matrix to simply add it to a DataFrame and rename the samples.
        #df = pd.DataFrame(pca.components_, columns=list(map(lambda kdbrdr: os.path.splitext(os.path.basename(kdbrdr._filepath))[0], files)))
    elif arguments.method == "tSNE":
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
    logger.info("Done printing {0} matrix to STDOUT".format(arguments.method))

    #logger.info("Beginning distribution analysis in R...")
    #logger.warn("Not implemented in Python...")
    
    sys.stderr.write(config.DONE)



def kmeans(arguments):
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
        logger.error("An unknown IO type was detected in bin/kdb.cluster()")
        sys.exit(1)

    sys.stderr.write(config.DEFAULT_MASTHEAD) # Not working
    if logger.level == logging.DEBUG or logger.level == logging.INFO: # Not working
        sys.stderr.write(config.DEBUG_MASTHEAD)
    sys.stderr.write(config.KMEANS_MASTHEAD) # Not working

        
    num_samples, num_features = df.shape
    logger.info("Input DataFrame shape: {0}".format(df.shape))
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
        logger.info(list(column_names))
        logger.info(labels)
    
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
        logger.info(list(column_names))
        logger.info(labels)
        
    sys.stderr.write("Generated {0} clusters projected onto reduced dimension 1 and 2 of the input dataset\n".format(arguments.k))
    sys.stderr.write("The figure was written to {0}\n\n".format(config.kmeans_clustering_fig_filepath))

    logger.info("The annotated sample matrix can be printed by specifying an output files with [ -o|--output OUTFILE ]")
    
    sys.stderr.write(config.DONE)



def hierarchical(arguments):
    """
    Thanks to https://www.analyticsvidhya.com/blog/2019/05/beginners-guide-hierarchical-clustering/
    for the inspiration
    """
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
        logger.error("An unknown IO type was detected in bin/kdb.cluster()")
        sys.exit(1)

    sys.stderr.write(config.DEFAULT_MASTHEAD)
    if logger.level == logging.DEBUG or logger.level == logging.INFO:
        sys.stderr.write(config.DEBUG_MASTHEAD)
    sys.stderr.write(config.HIERARCHICAL_MASTHEAD)
    
    num_samples, num_features = df.shape
    logger.info("Input DataFrame shape: {0}".format(df.shape))

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
#     from kdb import config


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

#     sys.stderr.write(config.DEFAULT_MASTHEAD)
#     if logger.level == logging.DEBUG:
#         sys.stderr.write(config.DEBUG_MASTHEAD)
#     sys.stderr.write(config.RAREFY_MASTHEAD)


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
    from kmerdb import fileutil, config, util

    if os.path.splitext(arguments.kdb)[-1] != ".kdb":
        raise IOError("Viewable .kdb filepath does not end in '.kdb'")

    with fileutil.open(arguments.kdb, mode='r') as kdb_in:
        if kdb_in.metadata["version"] != config.VERSION:
            logger.warning("KDB file version is out of date, may be incompatible with current fileutil.KDBReader class")
        N = 4**kdb_in.metadata["k"]

        assert kdb_in.profile.size == N, "view | read profile size did not match N from the header metadata"
        assert kdb_in.kmer_ids.size == N, "view | read kmer_ids size did not match N from the header metadata"
        assert kdb_in.counts.size == N, "view | read counts size did not match N from the header metadata"
        assert kdb_in.frequencies.size == N, "view | read frequencies size did not match N from the header metadata"

        if arguments.json:
            print(dict(kdb_in.metadata))
        else:
            yaml.add_representer(OrderedDict, util.represent_ordereddict)
            print(yaml.dump(kdb_in.metadata))
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
    from kmerdb import fileutil, config, util, kmer
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
                logger.debug("Tried to parse:\n{0}\n".format(header))
                raise ValueError("Could not parse YAML formatted header")
            else:
                return header_dict

    assert type(arguments.kdb_in) is str, "kdb_in must be a str"
    if os.path.splitext(arguments.kdb_in)[-1] != ".kdb": # A filepath with invalid suffix
        raise IOError("Viewable .kdb filepath does not end in '.kdb'")
    elif not os.path.exists(arguments.kdb_in):
        raise IOError("Viewable .kdb filepath '{0}' does not exist on the filesystem".format(arguments.kdb_in))
    with fileutil.open(arguments.kdb_in, mode='r', sort=arguments.re_sort, slurp=True) as kdb_in:
        metadata = kdb_in.metadata
        kmer_ids_dtype = metadata["kmer_ids_dtype"]
        N = 4**metadata["k"]
        if metadata["version"] != config.VERSION:
            logger.warning("KDB version is out of date, may be incompatible with current KDBReader class")
        if arguments.kdb_out is None or (arguments.kdb_out == "/dev/stdout" or arguments.kdb_out == "STDOUT"): # Write to stdout, uncompressed
            if arguments.header:
                yaml.add_representer(OrderedDict, util.represent_ordereddict)
                print(yaml.dump(metadata, sort_keys=False))
                print(config.header_delimiter)
        logger.info("Reading from file...")
        logger.debug("I cut off the json-formatted unstructured column for the main view.")
        #logger.debug(kdb_in.profile)
        try:
            if not arguments.un_sort and arguments.re_sort and metadata["sorted"] is True:
                kmer_ids_sorted_by_count = np.lexsort((kdb_in.counts, kdb_in.kmer_ids))
                reverse_kmer_ids_sorted_by_count = np.flipud(kmer_ids_sorted_by_count)
                for i, idx in enumerate(kmer_ids_sorted_by_count):
                    p = kdb_in.profile[i]
                    kmer_id = kdb_in.kmer_ids[kdb_in.profile[i]]
                    logger.debug("The first is an implicit row-index. The second is the corresponding profile. The third is a k-mer id, then the counts and frequencies.")
                    logger.debug("{0}\t{1}\t{2}\t{3}\t{4}".format(i, p, kmer_id, kdb_in.counts[kmer_id], kdb_in.frequencies[kmer_id]))

                    print("{0}\t{1}\t{2}\t{3}\t{4}".format(i, kmer_id, p, kdb_in.counts[kmer_id], kdb_in.frequencies[kmer_id]))
            else:
                for i, idx in enumerate(kdb_in.kmer_ids):
                    p = kdb_in.profile[i]
                    kmer_id = kdb_in.kmer_ids[kdb_in.profile[i]]
                    logger.debug("The row in the file should follow this order:")
                    logger.debug("The first is an implicit row-index. The second is the corresponding profile. The third is a k-mer id, then the counts and frequencies.")
                    logger.debug("{0}\t{1}\t{2}\t{3}\t{4}".format(i, p, kmer_id, kdb_in.counts[kmer_id], kdb_in.frequencies[kmer_id]))
                    try:
                        if arguments.un_sort is True:
                            assert p == idx, "view | kmer_id {0} didn't match the loaded profile {1}".format(idx, p)
                            assert kmer_id == idx, "view | kmer_id {0} didn't match the expected k-mer id.".format(idx, kmer_id)
                            assert i == kmer_id, "view | kmer_id {0} didn't match the implicit index {1}".format(idx, i)
                        else:
                            #logger.debug("Not sorting, so skipping assertion about profile (col1, column 2)")
                            pass
                    except AssertionError as e:
                        logger.warning(e)
                        logger.warning("K-mer id {0} will be printed in the {1} row".format(idx, i))
                        #raise e
                    logger.debug("{0} line:".format(i))
                    logger.debug("=== = = = ======= =  =  =  =  =  = |")
                    if arguments.un_sort is True:
                        print("{0}\t{1}\t{2}\t{3}\t{4}".format(i, i, idx, kdb_in.counts[idx], kdb_in.frequencies[idx]))
                    else:
                        print("{0}\t{1}\t{2}\t{3}\t{4}".format(i, p, idx, kdb_in.counts[kmer_id], kdb_in.frequencies[kmer_id]))
                # I don't think anyone cares about the graph representation.
                # I don't think this actually matters because I can't figure out what the next data structure is.
                # Is it a Cypher query and creation node set?
                # I need to demonstrate a capacity for graph based learning.
                # (:-|X) The dread pirate roberts got me.
                # :)
        except BrokenPipeError as e:
            logger.error(e)
            raise e
    if arguments.kdb_out is not None and arguments.compress: # Can't yet write compressed to stdout
        logger.error("Can't write kdb to stdout! We need to use a Bio.bgzf filehandle.")
        sys.exit(1)
    elif arguments.kdb_out is not None and type(arguments.kdb_out) is not str:
        raise ValueError("Cannot write a file to an argument that isn't a string")
    elif arguments.kdb_out is not None and os.path.exists(arguments.kdb_out):
        logger.warning("Overwriting '{0}'...".format(arguments.kdb_out))
    elif arguments.kdb_out is not None and not os.path.exists(arguments.kdb_out):
        logger.debug("Creating '{0}'...".format(arguments.kdb_out))
    if arguments.kdb_out is not None:
        with fileutil.open(arguments.kdb_in, 'r', dtype=suggested_dtype, sort=arguments.sorted, slurp=True) as kdb_in:
            assert kdb_in.profile.size == N, "view | read profile size did not match N from the header metadata"
            assert kdb_in.kmer_ids.size == N, "view | read kmer_ids size did not match N from the header metadata"
            assert kdb_in.counts.size == N, "view | read counts size did not match N from the header metadata"
            assert kdb_in.frequencies.size == N, "view | read frequencies size did not match N from the header metadata"
            with fileutil.open(arguments.kdb_out, metadata=metadata, mode='w') as kdb_out:
                try:
                    for i, idx in enumerate(kdb_in.kmer_ids):
                        kmer_id = idx
                        seq = kmer.id_to_kmer(kmer_id, arguments.k)
                        kmer_metadata = kmer.neighbors(seq, arguments.k)
                        kmer_id = kdb_in.kmer_ids[kdb_in.profile[i]]
                        p = kdb_in.profile[j]
                        logger.debug("The first is the actual row id. This is the recorded row-id in the file. This should always be sequential. Next is the k-mer id. ")
                        kdb_out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(i, p, kmer_id, kdb_in.counts[kmer_id],  kdb_in.frequencies[kmer_id], kmer_metadata))
                    
                except StopIteration as e:
                    logger.error(e)
                    raise e
                finally:
                    #kdb_out._write_block(kdb_out._buffer)
                    #kdb_out._handle.flush()
                    #kdb_out._handle.close()
                    sys.stderr.write(config.DONE)

            
def profile(arguments):
    """
    A complex, near-end user function that handles an arparse Namespace as its only positional argument

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
    from multiprocessing import Pool
    import json
    import time
    import numpy as np
    from kmerdb import parse, fileutil, kmer, util
    from kmerdb.config import VERSION

    logger.debug("Printing entire CLI argparse option Namespace...")
    logger.debug(arguments)

    # The number of specified processors should be available
    logger.debug("Validating processor count for parallelization...")
    if arguments.parallel <= 0 or arguments.parallel > cpu_count()+1:
        raise argparse.ArgumentError("-p, --parallel must be a valid processor count.")

    # The extension should be .kdb because I said so.
    logger.info("Checking extension of output file...")
    if os.path.splitext(arguments.kdb)[-1] != ".kdb":
        raise IOError("Destination .kdb filepath does not end in '.kdb'")




    file_metadata = []
    total_kmers = 4**arguments.k # Dimensionality of k-mer profile
    N = total_kmers

    logger.info("Parsing {0} sequence files to generate a composite k-mer profile...".format(len(list(arguments.seqfile))))
    nullomers = set()
    # Initialize arrays
    

    infile = parse.Parseable(arguments) #

    
    logger.info("Processing {0} fasta/fastq files across {1} processors...".format(len(list(arguments.seqfile)), arguments.parallel))
    logger.debug("Parallel (if specified) mapping the kmerdb.parse.parsefile() method to the seqfile iterable")
    logger.debug("In other words, running the kmerdb.parse.parsefile() method on each file specified via the CLI")
    if arguments.parallel > 1:
        with Pool(processes=arguments.parallel) as pool:
            data = pool.map(infile.parsefile, arguments.seqfile)
    else:
        data = list(map(infile.parsefile, arguments.seqfile))

    # 'data' is now a list of 4-tuples
    # Each 4-tuple represents a single file
    # (counts<numpy.array>, header_dictionary<dict>, nullomers<list>, all_kmer_metadata<list>)

    # Construct a final_counts array for the composite profile across all inputs
    logger.debug("Initializing large list for extended metadata")
    all_kmer_metadata = list([] for x in range(total_kmers)) if arguments.all_metadata else None
    logger.debug("Allocation completed")
    counts = np.zeros(N, dtype="uint64")
    # Complete collating of counts across files
    # This technically uses 1 more arrray than necessary 'final_counts' but its okay
    logger.info("Summing counts from individual fasta/fastq files into a composite profile...")
    for d in data:
        counts = counts + d[0]
        if arguments.all_metadata:
            logger.info("\n\nMerging metadata from all files...\n\n")
            all_kmer_metadata = util.merge_metadata_lists(arguments.k, all_kmer_metadata, d[3])

    sys.stderr.write("\n\n\tCompleted summation and metadata aggregation across all inputs...\n\n")
    unique_kmers = int(np.count_nonzero(counts))
    total_nullomers = total_kmers - unique_kmers
    all_observed_kmers = int(np.sum(counts))

    logger.info("Initial counting process complete, creating BGZF format file (.kdb)...")
    logger.info("Formatting master metadata dictionary...")

    metadata=OrderedDict({
        "version": VERSION,
        "metadata_blocks": 1,
        "k": arguments.k,
        "total_kmers": all_observed_kmers,
        "unique_kmers": unique_kmers,
        "unique_nullomers": total_nullomers,
        "metadata": arguments.all_metadata,
        "sorted": arguments.sorted,
        "kmer_ids_dtype": "uint64",
        "profile_dtype": "uint64",
        "count_dtype": "uint64",
        "frequencies_dtype": "float64",
        "tags": [],
        "files": [d[1] for d in data]
    })
        
    try:
        np.dtype(metadata["kmer_ids_dtype"])
        np.dtype(metadata["profile_dtype"])
        np.dtype(metadata["count_dtype"])
        np.dtype(metadata["frequencies_dtype"])
    except TypeError as e:
        logger.error(e)
        logger.error("kmerdb encountered a type error and needs to exit")
        raise TypeError("Incorrect dtype.")

    logger.debug("Initializing Numpy array of {0} uint zeroes for the final composite profile...".format(total_kmers))
    kmer_ids = np.array(range(N), dtype=metadata["kmer_ids_dtype"])
    profile = np.array(range(N), dtype=metadata["profile_dtype"])
    counts = np.array(counts, dtype=metadata["count_dtype"])
    frequencies = np.divide(counts, metadata["total_kmers"])
    logger.info("Initialization of profile completed, using approximately {0} bytes per profile".format(counts.nbytes))    
    sys.stderr.write(str(yaml.dump(metadata, sort_keys=False)))

    


    logger.info("Collapsing the k-mer counts across the various input files into the final kdb file '{0}'".format(arguments.kdb)) 
    kdb_out = fileutil.open(arguments.kdb, 'wb', metadata=metadata)
    try:
        sys.stderr.write("\n\nWriting outputs to {0}...\n\n".format(arguments.kdb))


        # Numpy indexing is trash
        #kmer_ids = np.zeros(total_kmers, dtype=metadata["kmer_ids_dtype"])
        #profile = np.zeros(total_kmers, dtype=metadata["profile_dtype"])
        #counts = np.zeros(total_kmers, dtype=metadata["count_dtype"])
        #frequencies = np.zeros(total_kmers, dtype=metadata["frequencies_dtype"])
        all_metadata = []
        if arguments.all_metadata:

            if arguments.sorted:
                j = 0
                list_of_duples = list(zip(kmer_ids, counts))
                logger.debug(list_of_duples[0:3])
                kmer_ids_sorted_by_count = np.lexsort(list_of_duples)
                reverse_kmer_ids_sorted_by_count = np.flipud(kmer_ids_sorted_by_count)
                for i, idx in enumerate(reverse_kmer_ids_sorted_by_count):
                    seq = kmer.id_to_kmer(idx, arguments.k)
                    kmer_metadata = kmer.neighbors(seq, arguments.k) # metadata is initialized by the neighbors
                    reads = []
                    starts = []
                    reverses = []
                    frequency = float(count)/N
                    for read, start, reverse in all_kmer_metadata[idx]:
                        reads.append(read)
                        starts.append(start)
                        reverses.append(reverse)
                        kmer_metadata["reads"] = reads
                        kmer_metadata["starts"] = starts
                        kmer_metadata["reverses"] = reverses
                    counts[idx] = counts[idx]
                    frequencies[idx] = frequencies[idx]
                    all_metadata.append(kmer_metadata)
                    assert j <= i + 1, "profile | row index was not sequential. Possibly read improperly."
                    assert idx == kmer_ids[i], "profile | row index did not match a kmer_id"
                    logger.info("First is the implicit row index, next is a k-mer id, next is the corresponding k-mer id to the row-index (may not match from sorting), next is the count and frequencies")
                    print("{0}\t{1}\t{2}\t{3}\t{4}".format(i, idx, kmer_ids[i], counts[idx], frequencies[idx]))
                    kdb_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i, idx, kmer_ids[i], counts[idx], frequencies[idx], json.dumps(kmer_metadata)))
                    j += 1
            else:
                j = 0
                for i, idx in enumerate(kmer_ids):
                    kmer_id = int(idx)
                    seq = kmer.id_to_kmer(kmer_id, arguments.k)
                    #logger.info("{0}\t{1}\t{2}\t{3}\t{4}".format(i, idx, kmer_ids[i], counts[kmer_id], frequencies[kmer_id]))
                    
                    kmer_metadata = kmer.neighbors(seq, arguments.k)
                    all_metadata.append(kmer_metadata)
                    counts[idx] = counts[idx]
                    frequencies[idx] = frequencies[idx]
                    logger.info("First is the implicit row index, next is a k-mer id, next is the corresponding k-mer id to the row-index (may not match from sorting), next is the count and frequencies")
                    if arguments.quiet is not True:
                        print("{0}\t{1}\t{2}\t{3}\t{4}".format(i, idx, kmer_ids[i], counts[idx], frequencies[idx]))
                    kdb_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i, idx, kmer_ids[i], counts[i], frequencies[i], json.dumps(kmer_metadata)))
                    j += 1
        else:
            if arguments.sorted:
                duple_of_arrays = (kmer_ids, counts)
                kmer_ids_sorted_by_count = np.lexsort(duple_of_arrays)
                reverse_kmer_ids_sorted_by_count = np.flipud(kmer_ids_sorted_by_count)
                logger.debug("K-mer id sort shape: {0}".format(len(kmer_ids_sorted_by_count)))
                for i, idx in enumerate(reverse_kmer_ids_sorted_by_count):

                    kmer_id = int(kmer_ids[idx])
                    logger.info("{0}\t{1}\t{2}\t{3}\t{4}".format(i, idx, kmer_ids[i], counts[idx], frequencies[idx]))
                    seq = kmer.id_to_kmer(kmer_id, arguments.k)
                    kmer_metadata = kmer.neighbors(seq, arguments.k)
                    all_metadata.append(kmer_metadata)
                    counts[idx] = counts[idx]
                    frequencies[idx] = frequencies[idx]
                    if arguments.quiet is not True:
                        print("{0}\t{1}\t{2}\t{3}\t{4}".format(i, idx, kmer_ids[idx], counts[idx], frequencies[idx]))
                    kdb_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i, idx, kmer_ids[i], counts[idx], frequencies[idx], json.dumps(kmer_metadata)))
            else:
                for i, idx in enumerate(kmer_ids):
                    kmer_id = int(idx)
                    logger.info("{0}\t{1}\t{2}\t{3}\t{4}".format(i, idx, kmer_ids[i], counts[idx], frequencies[idx]))
                    p = profile[i]
                    c = counts[i]
                    f = frequencies[i]
                    seq = kmer.id_to_kmer(int(kmer_id), arguments.k)
                    kmer_metadata = kmer.neighbors(seq, arguments.k) # metadata is initialized by the neighbors
                    all_metadata.append(kmer_metadata)
                    if arguments.quiet is not True:
                        print("{0}\t{1}\t{2}\t{3}\t{4}".format(i, idx, kmer_ids[i], counts[i], frequencies[i]))
                    kdb_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i, idx, kmer_ids[i], counts[i], frequencies[i], json.dumps(kmer_metadata)))
        logger.info("Wrote 4^k = {0} k-mer counts + neighbors to the .kdb file.".format(total_kmers))

        logger.info("Done")
            
    finally:
        kdb_out._write_block(kdb_out._buffer)
        kdb_out._handle.flush()
        kdb_out._handle.close()

        sys.stderr.write("Total k-mers processed: {0}\n".format(all_observed_kmers))
        sys.stderr.write("Final nullomer count:   {0}\n".format(total_nullomers))
        sys.stderr.write("Unique {0}-mer count:     {1}\n".format(arguments.k, unique_kmers))
        sys.stderr.write("Total {0}-mer count:     {1}\n".format(arguments.k, total_kmers))
    

        sys.stderr.write("\nDone\n")
    
        
def get_root_logger(level):
    levels=[logging.WARNING, logging.INFO, logging.DEBUG]
    if level < 0 or level > 2:
        raise TypeError("{0}.get_root_logger expects a verbosity between 0-2".format(__file__))
    logging.basicConfig(level=levels[level], format="%(levelname)s: %(asctime)s %(funcName)s L%(lineno)s| %(message)s", datefmt="%Y/%m/%d %I:%M:%S")
    root_logger = logging.getLogger()
    for name in logging.Logger.manager.loggerDict.keys():
        if ('boto' in name) or ('urllib3' in name) or ('s3' in name) or ('findfont' in name):
            logging.getLogger(name).setLevel(logging.WARNING)

    return root_logger


def citation_info():
    import pkg_resources
    citation = None
    if pkg_resources.resource_exists('kmerdb', 'CITATION'):
        citation = pkg_resources.resource_string('kmerdb', 'CITATION').decode('utf-8').rstrip()
        if citation == "":
            return
        else:
            sys.stderr.write("Printing citation notice to stderr. This will not interfere with the execution of the program in any way. Please see CITATION_FAQ.md for any questions.\n")
            sys.stderr.write(citation + "\n\n\n")
            sys.stderr.write("Run 'kmerdb citation' to silence.\n")
    else:
        raise IOError("Cannot locate the extra package data file 'kmerdb/CITATION', which should have been distributed with the program")


def cli():
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


    profile_parser = subparsers.add_parser("profile", help="Parse data into the database from one or more sequence files")
    profile_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    profile_parser.add_argument("-k", default=12, type=int, help="Choose k-mer size (Default: 12)")

    profile_parser.add_argument("-p", "--parallel", type=int, default=1, help="Shred k-mers from reads in parallel")

    profile_parser.add_argument("--batch-size", type=int, default=100000, help="Number of updates to issue per batch to PostgreSQL while counting")
    profile_parser.add_argument("-b", "--fastq-block-size", type=int, default=100000, help="Number of reads to load in memory at once for processing")
    profile_parser.add_argument("-n", type=int, default=1000, help="Number of k-mer metadata records to keep in memory at once before transactions are submitted, this is a space limitation parameter after the initial block of reads is parsed. And during on-disk database generation")
    #profile_parser.add_argument("--keep-S3-file", action="store_true", help="Download S3 file to the current working directory")
    profile_parser.add_argument("--keep-db", action="store_true", help=argparse.SUPPRESS)
    profile_parser.add_argument("--both-strands", action="store_true", default=False, help="Retain k-mers from the forward strand of the fast(a|q) file only")
    profile_parser.add_argument("--all-metadata", action="store_true", default=False, help="Include read-level k-mer metadata in the .kdb")
    profile_parser.add_argument("--sorted", action="store_true", default=False, help="Sort the output kdb file by count")
    profile_parser.add_argument("--quiet", action="store_true", default=False, help="Do not log the entire .kdb file to stdout")
    #profile_parser.add_argument("--sparse", action="store_true", default=False, help="Whether or not to store the profile as sparse")

    profile_parser.add_argument("seqfile", nargs="+", type=str, metavar="<.fasta|.fastq>", help="Fasta or fastq files")
    profile_parser.add_argument("kdb", type=str, help="Kdb file")
    profile_parser.set_defaults(func=profile)

    header_parser = subparsers.add_parser("header", help="Print the YAML header of the .kdb file and exit")
    header_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    header_parser.add_argument("-j", "--json", help="Print as JSON. DEFAULT: YAML")
    header_parser.add_argument("kdb", type=str, help="A k-mer database file (.kdb)")
    header_parser.set_defaults(func=header)
    
    view_parser = subparsers.add_parser("view", help="View the contents of the .kdb file")
    view_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    view_parser.add_argument("-H", "--header", action="store_true", help="Include header in the output")
    view_parser.add_argument("--re-sort", action="store_true", help="Sort the k-mer profile before displaying")
    view_parser.add_argument("--un-sort", action="store_true", help="Unsort the k-mer profile before displaying")
    view_parser.add_argument("--dtype", type=str, default="uint64", help="Read in the profiles as unsigned integer 64bit NumPy arrays")
    view_parser.add_argument("-d", "--decompress", action="store_true", help="Decompress input? DEFAULT: ")
    view_parser.add_argument("-c", "--compress", action="store_true", help="Print compressed output")
    view_parser.add_argument("kdb_in", type=str, nargs="?", default=None, help="A k-mer database file (.kdb) to be read (Optional)")
    view_parser.add_argument("kdb_out", type=str, nargs="?", default=None, help="A k-mer database file (.kdb) to be written to (Optional)")
    view_parser.set_defaults(func=view)


    matrix_parser = subparsers.add_parser("matrix", help="Generate a reduced-dimensionality matrix of the 4^k * n (k-mers x samples) data matrix.")
    matrix_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
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
    matrix_parser.add_argument("method", choices=["PCA", "tSNE", "DESeq2", "pass", "Frequency"], default=None, help="Choice of distance metric between two profiles")
    matrix_parser.add_argument("input", nargs="*", default=[], metavar="<kdbfile1 kdbfile2 ...|input.tsv|STDIN>", help="Two or more .kdb files, or another count matrix in tsv/csv")
    matrix_parser.set_defaults(func=get_matrix)
    
    # rarefy_parser = subparsers.add_parser("rarefy", help="Generate rarefaction information using ecopy.diversity.rarefy for the supplied .kdb files")
    # rarefy_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    # rarefy_parser.add_argument("-d", "--delimiter", default="\t", type=str, help="The choice of delimiter to parse the DataFrame with")
    # rarefy_parser.add_argument("--output-delimiter", type=str, default="\t", help="The output delimiter of the final csv/tsv to write.")
    # rarefy_parser.add_argument("-i", "--input", type=argparse.FileType("r"), default=None, help="The input reduced dimension or simply normalized matrix to use with K-means clustering")
    # rarefy_parser.add_argument("-o", "--output", type=argparse.FileType("w"), default=None, help="THe output csv/tsv of rarefied data")
    # rarefy_parser.set_defaults(func=rarefy)

    kmeans_parser = subparsers.add_parser("kmeans", help="Cluster the files according to their k-mer profile")
    kmeans_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    kmeans_parser.add_argument("-d", "--delimiter", type=str, default="\t", help="The delimiter of the input csv/tsv to parse, presumably produced by 'kdb matrix'.")
    kmeans_parser.add_argument("--output-delimiter", type=str, default="\t", help="The output delimiter of the final csv/tsv to write.")
    kmeans_parser.add_argument("--distance", type=str, default='e', choices=['e', 'b', 'c', 'a', 'u', 'x', 's', 'k'], help="Different distance metrics offered by kcluster. It is highly advised that you check both this source and their documentation to see how this is implemented.")

    
    kmeans_parser.add_argument("-k", default=None, type=int, help="The choice of k for clustering", required=True)
    kmeans_parser.add_argument("-i", "--input", type=argparse.FileType("r"), default=None, help="The input reduced dimension or mereley normalized matrix to use with K-means clustering")
    kmeans_parser.add_argument("-o", "--output", type=argparse.FileType("w"), default=None, help="The output csv/tsv with added 'label' to use for graphing in R, if the matplotlib graphs are not sufficient.")
    kmeans_parser.add_argument("method", choices=["sklearn", "Biopython"], default="Biopython", help="The Python implementation of k-means to use. The Biopython method is selected for access to alternative distance metrics")
    kmeans_parser.set_defaults(func=kmeans)

    hierarchical_parser = subparsers.add_parser("hierarchical", help="Use scipy.cluster.hierarchy to generate a dendrogram from a distance matrix")
    hierarchical_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    hierarchical_parser.add_argument("-d", "--delimiter", type=str, default="\t", help="The delimiter to use when reading the csv.")
    hierarchical_parser.add_argument("-i", "--input", type=argparse.FileType("r"), default=None, help="The input distance matrix for hierarchical clustering")
    hierarchical_parser.add_argument("-m", "--method", type=str, choices=["single", "complete", "average", "weighted", "centroid", "median", "ward"], default="ward", help="The method for linkage fitting in R to use")
    hierarchical_parser.set_defaults(func=hierarchical)
    
    dist_parser = subparsers.add_parser("distance", help="Calculate various distance metrics between profiles")
    dist_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
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


    index_parser = subparsers.add_parser("index", help="Create a index file that can be held in memory")
    index_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    index_parser.add_argument("kdb", type=str, help="A k-mer database file (.kdb)")
    index_parser.set_defaults(func=index_file)

    shuf_parser = subparsers.add_parser("shuf", help="Create a shuffled .kdb file")
    shuf_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    shuf_parser.add_argument("kdb_in", type=str, help="An indexed k-mer database file (.kdb)")
    shuf_parser.add_argument("kdb_out", type=str, help="The output shuffled k-mer database file (.kdb)")
    shuf_parser.set_defaults(func=shuf)

    
    markov_probability_parser = subparsers.add_parser("probability", help=u"""
Calculate the log-odds ratio of the Markov probability of a given sequence from the product (pi) of the transition probabilities(aij) times the frequency of the first k-mer (P(X1)), given the entire k-mer profile of a species.

See https://matthewralston.github.io/quickstart#kmerdb-probability for more details.

1. Durbin, R., Eddy, S.R., Krogh, A. and Mitchison, G., 1998. Biological sequence analysis: probabilistic models of proteins and nucleic acids. Cambridge university press.
""")

    markov_probability_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    markov_probability_parser.add_argument("-d", "--delimiter", type=str, default="\t", help="The delimiter to use when reading the csv.")
    markov_probability_parser.add_argument("-b", "--fastq-block-size", type=int, default=100000, help="Number of reads to load in memory at once for processing")
    markov_probability_parser.add_argument("seqfile", type=str, metavar="<.fasta|.fastq>", default=None, help="Sequences to calculate standard Markov-chain probabilities from, either .fasta or .fastq")
    markov_probability_parser.add_argument("kdb", type=str, help="An indexed k-mer database file (.kdb)")
    markov_probability_parser.set_defaults(func=markov_probability)

    citation_parser = subparsers.add_parser("citation", help="Silence the citation notice on further runs")
    citation_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    citation_parser.set_defaults(func=citation)
    
    args=parser.parse_args()
    global logger
    logger = get_root_logger(args.verbose)

    sys.stderr.write("Constructed a logger for the program...\n")
    #logger.debug(sys.path)
    args.func(args)
    
