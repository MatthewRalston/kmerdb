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
        with fileutil.open(arguments.kdb, 'r') as kdb:
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
    import pandas as pd
    import numpy as np
    from scipy.spatial.distance import pdist, squareform


    from kmerdb import fileutil, distance, config
    n = len(arguments.input)


    if len(arguments.input) > 1:
        files = list(map(lambda f: fileutil.open(f, 'r'), arguments.input))
        logger.debug("Files: {0}".format(files))
        if arguments.k is None:
            arguments.k = files[0].k
        if not all(os.path.splitext(kdb)[-1] == ".kdb" for kdb in arguments.input):
            raise IOError("One or more parseable .kdb filepaths did not end in '.kdb'")
        if not all(kdbrdr.k == arguments.k for kdbrdr in files):
            logger.error("Files: {0}".format(files))
            logger.error("Choices of k: {0}".format([kdbrdr.k for kdbrdr in files]))
            logger.error("By default the inferred value of k is used when k is not specified at the command line, which was {0}".format(arguments.k))
            raise TypeError("One or more files did not have k set to be equal to {0}".format(arguments.k))
        logger.debug("Files: {0}".format(files))
        data = [kdbrdr.slurp() for kdbrdr in files]
        logger.debug(data)
        profiles = np.transpose(np.array(data, dtype="int32"))
        # The following does *not* transpose a matrix defined as n x N=4**k
        # n is the number of independent files/samples, (fasta/fastq=>kdb) being assessed
        # N=4**k is the dimensionality of the vector, sparse or not, that makes up the perceived profile, the descriptors is the column dimension.

        # The names for the columns are taken as the basenamed filepath, with all extensions (substrings beginning with '.') stripped.
        logger.info("Converting arrays of k-mer counts into a pandas DataFrame...")

        column_names = list(map(lambda kdbrdr: os.path.basename(kdbrdr._filepath).split(".")[0], files))
        if len(column_names) != len(files):
            raise RuntimeError("Number of column names {0} does not match number of input files {1}...".format(len(column_names), len(files)))
        logger.debug("Shape: {0}".format(profiles.shape))
        expected = (4**arguments.k, len(column_names))
        if profiles.shape != expected:
            logger.error("Expected shape: {0}".format(expected))
            logger.error("Actual shape: {0}".format(profiles.shape))
            raise RuntimeError("Raw profile shape (a Numpy array) doesn't match expected dimensions")

        df = pd.DataFrame(profiles, columns=column_names)
    elif len(arguments.input) == 0 or (len(arguments.input) == 1 and (arguments.input[0] == "STDIN" or arguments.input[0] == "/dev/stdin")):
        logger.info("Reading input as tsv/csv from STDIN")
        try:
            df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        profiles = np.array(df)
        column_names = list(df.columns)
        n = len(column_names)
    elif len(arguments.input) == 1 and (os.path.splitext(arguments.input[0])[-1] == ".tsv" or os.path.splitext(arguments.input[0])[-1] == ".csv"):
        logger.info("Hidden: 1 argument. Reading input as tsv")
        try:
            df = pd.read_csv(arguments.input[0], sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        profiles = np.array(df)
        column_names = list(df.columns)
        n = len(column_names)
    elif len(arguments.input) == 1 and os.path.splitext(arguments.input)[-1] == ".kdb":
        logger.error("kdb distance requires more than one .kdb file as positional inputs")
        sys.exit(1)
    else:
        logger.error("bin/kdb.distances() received {0} arguments as input, which were not supported.".format(len(arguments.input)))
        sys.exit(1)

    # Masthead stuff
    sys.stderr.write(config.DEFAULT_MASTHEAD)
    if logger.level == logging.DEBUG:
        sys.stderr.write(config.DEBUG_MASTHEAD)
    sys.stderr.write(config.DISTANCE_MASTHEAD)

    logger.info("Calculating a {0}x{0} '{1}' distance matrix...".format(n, arguments.metric))
    
    if arguments.metric in ["spearman", "EMD", "d2s"]:
        data = [['' for x in range(n)] for y in range(n)]
        for i in range(n):
            for j in range(n):
                logger.info("Calculating {0} distance between {1} and {2}...".format(arguments.metric, column_names[i], column_names[j]))
                if i == j:
                    data[i][j] = distance.identity[arguments.metric]
                elif i > j:
                    data[i][j] = None
                elif i < j:
                    if arguments.metric == "correlation":
                        data[i][j] = distance.correlation(arguments.input[i], arguments.input[j])
                        # data[i][j] = distance.correlation(arguments.kdb[i], arguments.kdb[j])
                    elif arguments.metric == "euclidean":
                        data[i][j] = distance.euclidean(arguments.input[i], arguments.input[j])
                    elif arguments.metric == "spearman":
                        cor, pval = distance.spearman(profiles[i], profiles[j])
                        # FIXME! also print pval matrices
                        data[i][j] = cor
                    elif arguments.metric == "EMD":
                        data[i][j] = distance.EMD(profiles[i], profiles[j])
                    elif arguments.metric == "d2s":
                        data[i][j] = distance.d2s(profiles[i], profiles[j])
                    else:
                        logger.error("Other distances are not implemented yet")
                        sys.exit(1)
                # This double loop quickly identifies empty cells and sets the data correctly from the permutation above
        for i in range(n):
            for j in range(n):
                if data[i][j] is None:
                    data[i][j] = data[j][i]
        logger.info("Printing distance matrix...")
        logger.info(data)
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
        df = pd.DataFrame(dist, columns=column_names)


        ## FIXME: CUSTOM sorting code, not commiting to git repo
        #suffixes = [(int(x.split("_")[1]), i) for i, x in enumerate(column_names)] # A list of a 2-tuple of the correct sort order and the index
        #suffixes.sort(key=lambda x: x[0])
        #sorted_column_names = [column_names[s[1]] for s in suffixes]
        #df.sort_values(sorted_column_names)
        

        df.to_csv(sys.stdout, sep=arguments.output_delimiter, index=False)

def get_matrix(arguments):
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    
    from kmerdb import fileutil, config

    
    if len(arguments.input) > 1:
        files = list(map(lambda f: fileutil.open(f, 'r'), arguments.input))
        if arguments.k is None:
            arguments.k = files[0].k
        if not all(os.path.splitext(kdb)[-1] == ".kdb" for kdb in arguments.input):
            raise IOError("One or more parseable .kdb filepaths did not end in '.kdb'")
        if not all(kdbrdr.k == arguments.k for kdbrdr in files):
            logger.error("Files: {0}".format(files))
            logger.error("Choices of k: {0}".format([kdbrdr.k for kdbrdr in files]))
            logger.error("By default the inferred value of k is used when k is not specified at the command line, which was {0}".format(arguments.k))
            raise TypeError("One or more files did not have k set to be equal to {0}".format(arguments.k))
        profiles = np.transpose(np.array(list(map(lambda kdbrdr: kdbrdr.slurp(), files)), dtype="int32"))
        # The following does *not* transpose a matrix defined as n x N=4**k
        # n is the number of independent files/samples, (fasta/fastq=>kdb) being assessed
        # N=4**k is the dimensionality of the vector, sparse or not, that makes up the perceived profile, the descriptors is the column dimension.

        # The names for the columns are taken as the basenamed filepath, with all extensions (substrings beginning with '.') stripped.
        logger.info("============================================================")
        logger.info("Converting arrays of k-mer counts into a pandas DataFrame...")
        logger.info("============================================================")

        column_names = list(map(lambda kdbrdr: os.path.basename(kdbrdr._filepath).split(".")[0], files))
        df = pd.DataFrame(profiles, columns=column_names)
    elif len(arguments.input) == 0 or (len(arguments.input) == 1 and (arguments.input[0] == "STDIN" or arguments.input[0] == "/dev/stdin")):
        logger.info("Reading input as tsv/csv from STDIN")
        try:
            df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        column_names = list(df.columns)
    elif len(arguments.input) == 1 and (os.path.splitext(arguments.input[0])[-1] == ".tsv" or os.path.splitext(arguments.input[0])[-1] == ".csv"):
        logger.info("Reading input file as tsv/csv")
        try:
            df = pd.read_csv(arguments.input[0], sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        column_names = list(df.columns)
    else:
        logger.error(arguments)
        logger.error("bin/kdb.get_matrix() received {0} arguments as input, and this is not supported.".format(len(arguments.input)))
        logger.error(arguments.input)
        sys.exit(1)
    sys.stderr.write(config.DEFAULT_MASTHEAD)
    if logger.level == logging.DEBUG:
        sys.stderr.write(config.DEBUG_MASTHEAD)
    sys.stderr.write(config.MATRIX_MASTHEAD)


    final_df = None
    if arguments.method == "Normalized":
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
            colData = pd.DataFrame(np.transpose(column_names), columns=["Species"])
            logger.info("colData:\n{0}".format(colData))
                
            dds = r['DESeqDataSetFromMatrix'](count_matrix, colData, Formula('~ Species'))
            logger.info("DESeq dataset:\n{0}".format(str(dds)))
            dds = r['estimateSizeFactors'](dds)
                
            normalized = r['counts'](dds, normalized = True)
            if not arguments.no_normalized_ints:
                normalized = np.array(np.rint(normalized), dtype="int64")
            normalized = pd.DataFrame(normalized, columns=column_names)
            logger.debug(normalized)
            logger.info("Normalized matrix shape: {0}".format(normalized.shape))

            #     if arguments.method == "Unnormalized":
            #         cf = r['descdist'](r_dataframe, discrete=True)
            #     else:
            #         cf = r['descdist'](r_dataframe, discrete=False)
        except PackageNotInstalledError as e:
            logger.error(e)
            logger.error("One or more R packages were not installed. Please see the install documentation about installing the associated R packages")
            sys.exit(1)
        #normalized.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)

        # logger.error("False ending of Normalized")
        # logger.debug("final_df should be set as normalized")
        # sys.exit(1)
        final_df = normalized
    elif arguments.method == "Unnormalized":
        #df.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)
        final_df = df
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
            score_df = pd.DataFrame(np.transpose(score_matrix), columns=column_names)

            #score_df.to_csv(sys.stdout, sep=arguments.delimiter, index=arguments.with_index)
            final_df = score_df
        else:
            logger.warning("You must look at '{0}' to decide on the choice of '-n' for specify when generating the dimensionally reduced matrix for clustering function".format(pca_variance_fig_filepath))
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
        tsne_df = pd.DataFrame(np.transpose(tsne), columns=column_names)
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
    if logger.level == logging.DEBUG: # Not working
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

        formatter_result = ("{:9s}\t{:.3f}s\t{:.0f}\t{:.3f}\t{:.3f}"
                            "\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n\n\n")
        sys.stderr.write("\t".join(["Name", "Total_time", "Inertia", "Homogeneity_score", "Completeness_score", "V_measure_score", "Adj_rand_score", "Adj_mutual_info_score", "Silhouette"]) + "\n")
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
    if logger.level == logging.DEBUG:
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
    from kmerdb import fileutil, config, util

    if os.path.splitext(arguments.kdb)[-1] != ".kdb":
        raise IOError("Viewable .kdb filepath does not end in '.kdb'")

    with fileutil.open(arguments.kdb, mode='r') as kdb:
        if kdb.metadata["version"] != config.VERSION:
            logger.warning("KDB version is out of date, may be incompatible with current KDBReader class")
        if arguments.json:
            print(dict(kdb.metadata))
        else:
            yaml.add_representer(OrderedDict, util.represent_ordereddict)
            print(yaml.dump(kdb.metadata))
            print(config.header_delimiter)
            
def view(arguments):
    from kmerdb import fileutil, config, util
    import json
    metadata = None

    def get_header(line, header):
        """
        A little helper recurrence function for grabbing the additions to the header.
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
        
    
    if type(arguments.kdb_in) is None or arguments.kdb_in == "STDIN" or arguments.kdb_in == "/dev/stdin": # Read from STDIN
        logger.warning("Interpreting data from STDIN as uncompressed .kdb input")
        metadata = ''
        kdb_in = None
        if arguments.decompress is True:
            ifile = gzip.open(sys.stdin.buffer, 'rt')
            while type(metadata) is str:
                #logger.debug("looping in gzip for header")
                line = ifile.readline()
                metadata = get_header(line, metadata)
            kdb_in = ifile
        else:
            i = 0
            while type(metadata) is str:
                line = sys.stdin.readline()
                #logger.debug("looping in stdin for header")
                #print(line)

                metadata = get_header(line, metadata)
                i += 1

                if i == 50:
                    raise RuntimeError()
            kdb_in = None
    else:
        assert type(arguments.kdb_in) is str, "kdb_in must be a str"
        if os.path.splitext(arguments.kdb_in)[-1] != ".kdb": # A filepath with invalid suffix
            raise IOError("Viewable .kdb filepath does not end in '.kdb'")
        elif not os.path.exists(arguments.kdb_in):
            raise IOError("Viewable .kdb filepath '{0}' does not exist on the filesystem".format(arguments.kdb_in))
        kdb_in = fileutil.open(arguments.kdb_in, mode='r')
        metadata = kdb_in.metadata
    if metadata["version"] != config.VERSION:
        logger.warning("KDB version is out of date, may be incompatible with current KDBReader class")
    if arguments.kdb_out is None or arguments.kdb_out == "/dev/stdout" or arguments.kdb_out == "STDOUT": # Write to stdout, uncompressed
        if arguments.header:
            yaml.add_representer(OrderedDict, util.represent_ordereddict)
            print(yaml.dump(metadata, sort_keys=False))
            print(config.header_delimiter)
        if kdb_in is None: # Read from STDIN, since there was no kdb_in file.
            logger.info("Reading from STDIN...")
            try:
                for line in sys.stdin:
                    print(line.rstrip())
            except BrokenPipeError as e:
                logger.error(e)
                raise e
        else: # Read from the kdb_in file, and not STDIN
            logger.info("Reading from file...")
            try:
                for line in kdb_in:
                    print(line.rstrip())
            except BrokenPipeError as e:
                logger.error(e)
                raise e
    elif arguments.kdb_out is not None and arguments.compress: # Can't yet write compressed to stdout
        logger.error("Can't write kdb to stdout! We need to use a Bio.bgzf filehandle.")
        sys.exit(1)
    elif type(arguments.kdb_out) is not str:
        raise ValueError("Cannot write a file to an argument that isn't a string")
    elif os.path.exists(arguments.kdb_out):
        logger.warning("Overwriting '{0}'...".format(arguments.kdb_out))
    elif not os.path.exists(arguments.kdb_out):
        logger.debug("Creating '{0}'...".format(arguments.kdb_out))
    else:
        with fileutil.open(arguments.kdb_out, metadata=metadata, mode='wb') as kdb_out:
            try:
                line = None
                if kdb_in is None:
                    while line is None:
                        line = sys.stdin.readline().rstrip()
                        kmer_id, count, kmer_metadata = fileutil.parse_line(line)
                        kdb_out.write("{0}\t{1}\t{2}\n".format(kmer_id, count, kmer_metadata))
                        line = None
                else:
                    for line in kdb_in:
                        line = line.rstrip()
                        kmer_id, count, kmer_metadata = fileutil.parse_line(line)
                        kdb_out.write("{0}\t{1}\t{2}\n".format(kmer_id, count, kmer_metadata))
                    
            except StopIteration as e:
                logger.error(e)
                raise e
            finally:
                #kdb_out._write_block(kdb_out._buffer)
                #kdb_out._handle.flush()
                #kdb_out._handle.close()
                sys.stderr.write(config.DONE)

            
def profile(arguments):
    import math
    from itertools import chain, repeat
    from multiprocessing import Pool
    import json
    import time
    from kmerdb import parse, fileutil, kmer, database
    from kmerdb.config import VERSION

    logger.debug(arguments)


    
    if os.path.splitext(arguments.kdb)[-1] != ".kdb":
        raise IOError("Destination .kdb filepath does not end in '.kdb'")
    
    file_metadata = []
    tempdbs = []
    logger.info("Parsing {0} sequence files to generate a composite k-mer profile...".format(len(list(arguments.seqfile))))
    nullomers = set()
    pool = Pool(processes=arguments.parallel)
    infile = parse.Parseable(arguments) # 
    if fileutil.is_all_fasta(list(arguments.seqfile)):
        logger.info("Processing {0} fasta files...".format(len(list(arguments.seqfile))))
        #logger.info("Processing fastas in parallel")
        logger.debug("Parallel (if specified) mapping the kmerdb.parse.parsefile() method to the seqfile iterable")
        logger.debug("In other words, running the kmerdb.parse.parsefile() method many times on each file multiply specified via the CLI")
        list_of_dbs = pool.map(infile.parsefile, arguments.seqfile)
    else:
        logger.info("Processing (some) fastqs with alternate parallelization schema. WARNING: UNTESTED")
        logger.warning("WARNING: UNTESTED")
        logger.debug("Mapping the kmerdb.parse.parsefile() method to the seqfile iterable")
        list_of_dbs = list(map(infile.parsefile, arguments.seqfile))
    logger.info("Generating list of temporary database handles")
    # list_of_dbs is a list of tuples returned from the mapping (either via pool.map or via regular Pythonic map) of the .parsefile() method
    # each tuple is a triple of the database handle, the metadata, and the set/list of nullomer ids
    for db, m, n in list_of_dbs:
        logger.debug("Appending metadata...")
        file_metadata.append(m)
        db = database.PostgresKdb(arguments.k, arguments.postgres_connection, tablename=db, filename=m["filename"])
        tempdbs.append(db)

        # Calculating nullomer additions
        lenBefore = len(nullomers)
        nullomers = nullomers.union(n)
        lenAfter = len(nullomers)
        newNullomers = lenAfter - lenBefore
        logger.info("{0} contributed {1} unique nullomers to the profile".format(m["filename"], len(n)))
        logger.info("Completed transfer of k-mers from {0} into PostgreSQL database.".format(m["filename"]))
    logger.info("Initial counting process complete, creating BGZF format file (.kdb)...")
    logger.info("Formatting master metadata dictionary...")
    metadata=OrderedDict({
        "version": VERSION,
        "metadata_blocks": 1,
        "k": arguments.k,
        "total_kmers": sum(list(map(lambda x: x["total_kmers"], file_metadata))),
        "unique_kmers": 4**arguments.k - len(nullomers),
        "metadata": False,
        "tags": [],
        "files": file_metadata
    })
        


    
    metadata_bytes = bgzf._as_bytes(yaml.dump(metadata))
    # The number of metadata blocks is the number of bytes of the header block(s) / the number of bytes per block in the BGZF specification
    metadata["metadata_blocks"] = math.ceil( sys.getsizeof(metadata_bytes) / ( 2**16 ) ) # First estimate
    metadata_bytes = bgzf._as_bytes(yaml.dump(metadata))
    metadata["metadata_blocks"] = math.ceil( sys.getsizeof(metadata_bytes) / ( 2**16 ) ) # Second estimate

    logger.info("Collapsing the k-mer counts across the various input files into the final kdb file '{0}'".format(arguments.kdb)) 
    kdb_out = fileutil.open(arguments.kdb, 'wb', metadata=metadata)
    try:
        x = 0
        iterating = True
        while iterating:
            try:
                logger.debug("Collating counts across all files for the {0} k-mer".format(x))
                kmer_dbrecs_per_file = list(map(next, tempdbs)) # Need to rename this variable

                # Unstable code
                #print(kmer_dbrecs_per_file)

                # raise RuntimeError("HOW DOES THIS HAVE NONE")
                if len(kmer_dbrecs_per_file):
                    i = kmer_dbrecs_per_file[0][0] - 1 # Remove 1 for the Sqlite zero-based indexing
                    count = sum([x[1] for x in kmer_dbrecs_per_file]) # The 1th element is the k-mer count, the 0th is the id
                    if arguments.verbose == 2:
                        sys.stderr.write("K-mer counts: {0} = {1}\n".format(list(map(lambda x: x[1], kmer_dbrecs_per_file)), count))
                    if count == 0 and arguments.sparse is True:
                        pass
                    else:
                        seq = kmer.id_to_kmer(i, arguments.k)
                        kmer_metadata = kmer.neighbors(seq, arguments.k) # metadata is initialized by the neighbors
                        # metadata now has three additional properties, based on the total number of times this k-mer occurred. Eventually the dimension of these new properties should match the count.
                        if arguments.all_metadata:


                            seqids = [x[4] for x in kmer_dbrecs_per_file]
                            starts = [x[2] for x in kmer_dbrecs_per_file]
                            reverses = [x[3] for x in kmer_dbrecs_per_file]


                            if len(reverses) == 0:
                                logger.error("REVERSES: {0}".format(reverses[0]))
                                raise RuntimeError("reverses: IS THIS INCORRECT?")
                            elif len(starts) == 0:
                                logger.error("STARTS: {0}".format(starts[0]))
                                raise RuntimeError("starts: IS THIS INCORRECT?")
                            elif len(seqids) == 0:
                                logger.error("SEQIDS: {0}".format(seqids[0]))
                                raise RuntimeError("seqids: IS THIS INCORRECT?")
                            elif len(seqids) == 1 and type(seqids) is list and type(seqids[0]) is list:
                                seqids = seqids[0]
                            elif len(starts) == 1 and type(starts) is list and type(starts[0]) is list:
                                starts = starts[0]
                            elif len(reverses) == 1 and type(reverses) is list and type(reverses[0]) is list:
                                reverses = reverses[0]
                            
                            if "seqids" in kmer_metadata.keys():
                                kmer_metadata["seqids"] += seqids
                            else:
                                kmer_metadata["seqids"] = seqids
                            if "starts" in kmer_metadata.keys():
                                kmer_metadata["starts"] += starts
                            else:
                                kmer_metadata["starts"] = starts
                            if "reverses" in kmer_metadata.keys():
                                kmer_metadata["reverses"] += reverses
                            else:
                                kmer_metadata["reverses"] = reverses
                        else:
                            kmer_metadata["seqids"]   = []
                            kmer_metadata["starts"]   = []
                            kmer_metadata["reverses"] = []
                        if type(kmer_metadata["starts"]) is str:
                            raise TypeError("kmerdb profile could not decode start sites from its preQLite3 database.")
                        elif type(kmer_metadata["starts"]) is list and all(type(x) is int for x in kmer_metadata["starts"]):
                            if arguments.verbose == 2:
                                sys.stderr.write("Parsed {0} k-mer start sites from this sequence.".format(len(kmer_metadata["starts"])))
                        elif type(kmer_metadata["starts"]) is dict:
                            logger.debug("Don't know how starts became a dictionary, but this will not parse correctly. RuntimeError")
                            raise RuntimeError("The implicit type of the Text blob in the Postgres database has changed, and will not parse correctly in kmerdb, rerun with verbose")
                        elif type(kmer_metadata["seqids"]) is list and all(type(x) is str for x in kmer_metadata["seqids"]):
                            if arguments.verbose == 2:
                                sys.stderr.write("Parsed {0} sequence ids associated with this k-mer.".format(len(kmer_metadata["seqids"])))
                        elif type(kmer_metadata["seqids"]) is dict:
                            logger.debug("Don't know how seqids became a dictionary, but this will not parse correctly. RuntimeError")
                            raise RuntimeError("The implicit type of the Text blob in the Postgres database has changed, and will not parse correctly in kmerdb, rerun with verbose")
                        elif type(kmer_metadata["reverses"]) is str:
                            raise TypeError("kmerdb profile could not decode strand information from its PostgreSQL database.")
                        elif type(kmer_metadata["reverses"]) is list and all(type(x) is bool for x in kmer_metadata["reverses"]):
                            if arguments.verbose == 2:
                                sys.stderr.write("Parsed {0} reverse? bools associated with this k-mer.".format(len(kmer_metadata["seqids"])))
                        elif type(kmer_metadata["reverses"]) is dict:
                            logger.debug("Don't know how reverses became a dictionary, but this will not parse correctly. RuntimeError")
                            raise RuntimeError("The implicit type of the Text blob in the Postgres database has changed, and will not parse correctly in kmerdb, rerun with verbose")
                        elif not all(type(x) is bool for x in kmer_metadata["reverses"]):
                            logger.error("kmer metadata: {0}".format(kmer_metadata))
                            logger.error("number of k-mer elements: {0}".format(len(kmer_metadata.values())))
                            logger.error(list(set(type(x) for x in kmer_metadata["reverses"])))
                            raise TypeError("Not all reverse bools were boolean")
                        elif count == 0:
                            n += 1
                        kdb_out.write("{0}\t{1}\t{2}\n".format(i, count, kmer_metadata))
                    x += 1
                else:
                    iterating = False
            except StopIteration as e:
                logger.warn("Shouldn't have encountered StopIteration error...")
                logger.warn("Continuing anyways...")
                iterating = False
        logger.info("Completed the transfer of data from the {0} temporary SQLite3 databases to the kdb file '{1}'".format(len(tempdbs), arguments.kdb))
        logger.info("Done")
    finally:
        import shutil
        from psycopg2 import sql
        kdb_out._write_block(kdb_out._buffer)
        kdb_out._handle.flush()
        kdb_out._handle.close()
        for db in tempdbs:
            if not arguments.keep_db:
                #os.unlink(db.filepath)
                with db.conn.begin():
                    db.conn.execute("DROP TABLE {}".format(db._tablename))
                # else:
            #     kdbfile = arguments.kdb
            #     kdbbasename = arguments.kdb.rstrip(".kdb")
            #     kdbsqlite = kdbbasename + ".sqlite3"
            #     shutil.move(db.filepath, kdbsqlite)
            #     sys.stderr.write("    Database file retained as    '{0}'".format(kdbsqlite))

            if db.conn is not None:
                db.conn.close()
            db._engine.dispose()

    
# def gen_hist(arguments):
#     from kdb import fileutil


#     hist = {}

#     with fileutil.open(arguments.kdb, mode='r') as ifile:
#         i = 0
#         for line in ifile:
#             kmer_id, count = (int(x) for x in line.rstrip().split("\t"))
#             if "\t" in line:
#                 try:
#                     hist[count] += 1
#                 except KeyError as e:
#                     hist[count] = 1
#             else:
#                 print(line)
#                 print(i)
#             i +=1

#     for k, v in hist.items():
#         print(k, v)
        
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
    profile_parser.add_argument("-pg", "--postgres-connection", type=str, help="A postgresql connection string, of the format postgres://user:password@host:port/dbname", required=True)
    profile_parser.add_argument("-p", "--parallel", type=int, default=1, choices=list(range(1, cpu_count()+1)), help="Shred k-mers from reads in parallel")
    profile_parser.add_argument("-pq", "--parallel-fastq", type=int, default=1, help="The number of blocks to read in parallel for reading fastqs")
    profile_parser.add_argument("--batch-size", type=int, default=100000, help="Number of updates to issue per batch to PostgreSQL while counting")
    profile_parser.add_argument("-b", "--fastq-block-size", type=int, default=100000, help="Number of reads to load in memory at once for processing")
    profile_parser.add_argument("-n", type=int, default=1000, help="Number of k-mer metadata records to keep in memory at once before transactions are submitted, this is a space limitation parameter after the initial block of reads is parsed. And during on-disk database generation")
    #profile_parser.add_argument("--keep-S3-file", action="store_true", help="Download S3 file to the current working directory")
    profile_parser.add_argument("--keep-db", action="store_true", help=argparse.SUPPRESS)
    profile_parser.add_argument("--strand-specific", action="store_true", default=False, help="Retain k-mers from the forward strand of the fast(a|q) file only")
    profile_parser.add_argument("--all-metadata", action="store_true", default=False, help="Include read-level k-mer metadata in the .kdb")
    profile_parser.add_argument("--sparse", action="store_true", default=False, help="Whether or not to store the profile as sparse")
    profile_parser.add_argument("-k", default=12, type=int, help="Choose k-mer size (Default: 12)")
    #profile_parser.add_argument("--random-seed", default=sys.)
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
    view_parser.add_argument("-d", "--decompress", action="store_true", help="Decompress input? DEFAULT: ")
    view_parser.add_argument("-c", "--compress", action="store_true", help="Print compressed output")
    view_parser.add_argument("kdb_in", type=str, nargs="?", default=None, help="A k-mer database file (.kdb) to be read (Optional)")
    view_parser.add_argument("kdb_out", type=str, nargs="?", default=None, help="A k-mer database file (.kdb) to be written to (Optional)")
    view_parser.set_defaults(func=view)


    matrix_parser = subparsers.add_parser("matrix", help="Generate a reduced-dimensionality matrix of the n * 4^k (sample x k-mer) data matrix.")
    matrix_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    matrix_parser.add_argument("--with-index", default=False, action="store_true", help="Print the row indices as well")
    matrix_parser.add_argument("-d", "--delimiter", default="\t", type=str, help="The choice of delimiter to parse the input .tsv with. DEFAULT: '\t'")
    matrix_parser.add_argument("--output-delimiter", type=str, default="\t", help="The output delimiter of the final csv/tsv to write. DEFAULT: '\t'")
        
    #matrix_parser.add_argument("--normalize-with", type=str, choices=["ecopy", "DESeq2"], default="DESeq2", help="Normalize with which method? DEFAULT: DESeq2")
    matrix_parser.add_argument("--no-normalized-ints", action="store_true", default=False, help="Don't round normalized counts to the nearest integer")
    matrix_parser.add_argument("-k", default=None, type=int, help="The k-dimension that the files have in common")
    matrix_parser.add_argument("-n", default=None, type=int, help="The number of dimensions to reduce with PCA or t-SNE. DEFAULT: an elbow graph will be generated if -n is not provided to help the user choose -n")

    matrix_parser.add_argument("--perplexity", default=5, type=int, help="The choice of the perplexity for t-SNE based dimensionality reduction")
    matrix_parser.add_argument("method", choices=["PCA", "tSNE", "Normalized", "Unnormalized"], default=None, help="Choice of distance metric between two profiles")
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
    dist_parser.add_argument("-d", "--delimiter", type=str, default="\t", help="The delimiter to use when printing the csv.")
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
    
