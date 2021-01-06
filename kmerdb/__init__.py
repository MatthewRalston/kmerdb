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
import time
from multiprocessing import cpu_count
from collections import OrderedDict


from Bio import bgzf

#import concurrent.futures

global logger
logger = None


def index_file(arguments):
    from kdb import fileutil, index
    with fileutil.open(arguments.kdb, mode='r') as kdb:
        header = kdb.header
    line_index = index.build_line_index_from_kdb(arguments.kdb, header['k'])
    index._write_line_index(arguments.kdbi, line_index)

def distances(arguments):
    import pandas as pd
    import numpy as np
    from scipy.spatial.distance import pdist, squareform


    from kdb import fileutil, distance, config
    n = len(arguments.kdb)


    if len(arguments.kdb) > 1:
        files = list(map(lambda f: fileutil.open(f, 'r'), arguments.kdb))
        if arguments.k is None:
            arguments.k = files[0].k
        if not all(kdbrdr.k == arguments.k for kdbrdr in files):
            logger.error("Files: {0}".format(files))
            logger.error("Choices of k: {0}".format([kdbrdr.k for kdbrdr in files]))
            logger.error("By default the inferred value of k is used when k is not specified at the command line, which was {0}".format(arguments.k))
            raise TypeError("One or more files did not have k set to be equal to {0}".format(arguments.k))
        profiles = np.array(list(map(lambda kdbrdr: kdbrdr.slurp(), files)), dtype="int32")
        # The following does *not* transpose a matrix defined as n x N=4**k
        # n is the number of independent files/samples, (fasta/fastq=>kdb) being assessed
        # N=4**k is the dimensionality of the vector, sparse or not, that makes up the perceived profile, the descriptors is the column dimension.

        # The names for the columns are taken as the basenamed filepath, with all extensions (substrings beginning with '.') stripped.
        logger.info("Converting arrays of k-mer counts into a pandas DataFrame...")
        column_names = list(map(lambda kdbrdr: os.path.basename(kdbrdr._filepath).split(".")[0], files))

        df = pd.DataFrame(profiles, columns=column_names)
    elif len(arguments.kdb) == 1 and (arguments.kdb[0] == "/dev/stdin" or arguments.kdb[0] == "STDIN"):
        logger.info("Hidden: 1 argument. Reading input from stdin")
        try:
            df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        profiles = np.array(df)
        column_names = list(df.columns)
        n = len(column_names)
    elif len(arguments.kdb) == 1 and os.path.splitext(arguments.kdb[0])[-1] == ".tsv":
        logger.info("Hidden: 1 argument. Reading input as tsv")
        try:
            df = pd.read_csv(arguments.kdb[0], sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        profiles = np.array(df)
        column_names = list(df.columns)
        n = len(column_names)
    elif len(arguments.kdb) == 1 and os.path.splitext(arguments.kdb)[-1] == ".kdb":
        logger.error("kdb distance requires more than one .kdb file as positional inputs")
        sys.exit(1)
    else:
        logger.error("bin/kdb.distances() received {0} arguments as input, which were not supported.".format(len(arguments.kdb)))
        sys.exit(1)

    # Masthead stuff
    sys.stderr.write(config.DEFAULT_MASTHEAD)
    if logger.level == logging.DEBUG:
        sys.stderr.write(config.DEBUG_MASTHEAD)
    sys.stderr.write(config.DISTANCE_MASTHEAD)

    logger.info("Calculating a {0}x{0} '{1}' distance matrix...".format(n, arguments.metric))
    
    if arguments.metric in ["spearman"]:
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
                        data[i][j] = distance.correlation(arguments.kdb[i], arguments.kdb[j])
                    elif arguments.metric == "euclidean":
                        data[i][j] = distance.euclidean(arguments.kdb[i], arguments.kdb[j])
                    elif arguments.metric == "spearman":
                        cor, pval = distance.spearman(profiles[:, i], profiles[:, j])
                        # FIXME! also print pval matrices
                        data[i][j] = cor
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

    df = pd.DataFrame(dist, columns=column_names)
    df.to_csv(sys.stdout, sep=arguments.output_delimiter, index=False)

def get_matrix(arguments):
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    
    from kdb import fileutil, config

    
    if len(arguments.kdb) > 1:
        files = list(map(lambda f: fileutil.open(f, 'r'), arguments.kdb))
        if arguments.k is None:
            arguments.k = files[0].k
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
    elif len(arguments.kdb) == 1 and (arguments.kdb == "STDIN" or arguments.kdb == "/dev/stdin"):
        logger.info("Hidden: 1 argument. Reading input as tsv from STDIN")
        try:
            df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        column_names = list(df.columns)
    elif len(arguments.kdb) == 1 and (os.path.splitext(arguments.kdb[0])[-1] == ".tsv"):
        logger.info("Hidden: 1 argument. Reading input as tsv from {0}".format(arguments.kdb[0]))
        try:
            df = pd.read_csv(arguments.kdb[0], sep=arguments.delimiter)
        except pd.errors.EmptyDataError as e:
            logger.error(e)
            logger.error("Pandas error on DataFrame reading. Perhaps a null dataset being read?")
            sys.exit(1)
        column_names = list(df.columns)
    else:
        logger.error("bin/kdb.get_matrix() received {0} arguments as input, and this is not supported.".format(len(arguments.kdb)))
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
            if not arguments.normalized_ints:
                normalized = np.array(np.rint(normalized), dtype="int32")
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
    final_df.to_csv(sys.stdout, sep=arguments.output_delimiter, index=arguments.with_index)
    logger.info("Done printing {0} matrix to STDOUT".format(arguments.method))

    logger.info("Beginning distribution analysis in R...")
    logger.warn("Not implemented in Python...")
    
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

    from kdb import config

    
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

        fig, ax = plt.subplots()
        scatter = ax.scatter(df.iloc[:, 1], df.iloc[:, 2], c=df.iloc[:, 0])
        legend = ax.legend(*scatter.legend_elements(), loc="upper right", title="Cluster")
        ax.add_artist(legend)
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
        ax.add_artist(legend)
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
    
    from kdb import config


    
    if isinstance(arguments.input, argparse.FileType) or isinstance(arguments.input, TextIOWrapper):
        df = pd.read_csv(arguments.input, sep=arguments.delimiter)
        column_names = df.columns
        df = df.transpose()
    elif arguments.input is None:
        df = pd.read_csv(sys.stdin, sep=arguments.delimiter)
        column_names = df.columns
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

    sys.stderr.write("Saving the dendrogram to '{0}'...".format(config.hierarchical_clustering_dendrogram_fig_filepath))
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
    from kdb import fileutil, config
    from kdb.config import VERSION
    from kdb import util

    with fileutil.open(arguments.kdb, mode='r') as kdb:
        if kdb.header["version"] != VERSION:
            logger.warning("KDB version is out of date, may be incompatible with current KDBReader class")
        if arguments.json:
            print(dict(kdb.header))
        else:
            yaml.add_representer(OrderedDict, util.represent_ordereddict)
            print(yaml.dump(kdb.header))
    
            
def view(arguments):
    from kdb import fileutil

    from kdb.config import VERSION

    with fileutil.open(arguments.kdb, mode='r') as kdb:
        if kdb.header["version"] != VERSION:
            logger.warning("KDB version is out of date, may be incompatible with current KDBReader class")
        if arguments.header:
            print(yaml.dump(kdb.header))
        for line in kdb:
            print(line.rstrip())

def profile(arguments):
    import math
    from kdb import parse, fileutil

    from kdb.config import VERSION

    metadata = []
    tempdbs = []
    for f in arguments.seqfile:
        db, m = parse.parsefile(f, arguments.k, p=arguments.parallel, b=arguments.fastq_block_size, stranded=arguments.not_strand_specific)
        metadata.append(m)
        tempdbs.append(db)
    header=OrderedDict({
        "version": VERSION,
        "metadata_blocks": 1,
        "k": arguments.k,
        "metadata": False,
        "tags": [],
        "files": metadata
    })
    
    header_bytes = bgzf._as_bytes(yaml.dump(header))
    # The number of metadata blocks is the number of bytes of the header block(s) / the number of bytes per block in the BGZF specification
    header["metadata_blocks"] = math.ceil( sys.getsizeof(header_bytes) / ( 2**16 ) ) # First estimate
    header_bytes = bgzf._as_bytes(yaml.dump(header))
    header["metadata_blocks"] = math.ceil( sys.getsizeof(header_bytes) / ( 2**16 ) ) # Second estimate
    #header["metadata_blocks"] = 2
    logger.info("Collapsing the k-mer counts across the various input files into the final kdb file '{0}'".format(arguments.kdb))    
    try:
        kdb_out = fileutil.open(arguments.kdb, 'wb', header)
        iterating = True
        while iterating:
            # The 0th element is the count
            try:

                kmer_counts_per_file = list(map(next, tempdbs)) # T
                if len(kmer_counts_per_file):
                    i = kmer_counts_per_file[0][0] - 1 # Remove 1 for the Sqlite zero-based indexing
                    count = sum([x[1] for x in kmer_counts_per_file]) # The 1th element is the k-mer count
                    #sys.stderr.write("\r")
                    if arguments.verbose == 2:
                        sys.stderr.write("K-mer counts: {0} = {1}\n".format(list(map(lambda x: x[1], kmer_counts_per_file)), count))
                    kdb_out.write("{0}\t{1}\n".format(i, count))
                else:
                    iterating = False
            except StopIteration as e:
                logger.warn("Shouldn't have encountered StopIteration error...")
                logger.warn("Continuing anyways...")
                iterating = False
        logger.info("Completed the transfer of data from the {0} temporary SQLite3 databases to the kdb file '{1}'".format(len(tempdbs), arguments.kdb))
        logger.info("Done")
    finally:
        for db in tempdbs:
            db.conn.close()
            db._engine.dispose()
            if not arguments.keep_sqlite:
                os.unlink(db.filepath)
            else:
                logger.debug("    Database file retained:    {0}".format(db.filepath))
        kdb_out._write_block(kdb_out._buffer)
        kdb_out._handle.flush()
        kdb_out._handle.close()


    
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
    logging.basicConfig(level=levels[level], format="%(levelname)s: %(asctime)s %(funcName)s L%(lineno)s| %(message)s", datefmt="%Y/%m%d %I:%M:%S")
    root_logger = logging.getLogger()
    for name in logging.Logger.manager.loggerDict.keys():
        if ('boto' in name) or ('urllib3' in name) or ('s3' in name) or ('findfont' in name):
            logging.getLogger(name).setLevel(logging.WARNING)

    return root_logger

def cli():
    sys.stderr.write("Running kdb script from '{0}'\n".format(__file__))
    sys.stderr.write("Checking installed environment...\n")
    primary_path = sys.path[0]
    ver_info = sys.version_info
    py_version = "python{0}.{1}".format(ver_info.major, ver_info.minor)
    sys.path.append(os.path.join(primary_path, "../lib/{0}/site-packages".format(py_version)))
    #sys.path.append(os.path.join(os.path.dirname(__file__), "../lib"))
    #site_packages = [p for p in sys.path if "kdb" in p]
    
    #sys.stderr.write("PYTHONPATH={0}".format(sys.path))
    #sys.path.remove(os.path.dirname(os.path.abspath(__file__)))


    
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(help="Use --help with sub-commands")


    profile_parser = subparsers.add_parser("profile", help="Parse data into the database from one or more sequence files")
    profile_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    profile_parser.add_argument("-p", "--parallel", type=int, default=1, choices=list(range(1, cpu_count()+1)), help="Shred k-mers from reads in parallel")
    profile_parser.add_argument("-b", "--fastq-block-size", type=int, default=100000, help="Number of reads to load in memory at once for processing")
    #profile_parser.add_argument("--keep-S3-file", action="store_true", help="Download S3 file to the current working directory")
    profile_parser.add_argument("--keep-sqlite", action="store_true", help=argparse.SUPPRESS)
    profile_parser.add_argument("--not-strand-specific", action="store_false", default=True, help="Retain k-mers from the forward strand of the fast(a|q) file only")
    #profile_parser.add_argument("--no-metadata", dest="metadata", action="store_false", default=True, help="Include k-mer metadata in the .kdb")
    profile_parser.add_argument("-k", default=12, type=int, help="Choose k-mer size (Default: 12)")
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
    view_parser.add_argument("kdb", type=str, help="A k-mer database file (.kdb)")
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
    matrix_parser.add_argument("kdb", nargs="+", type=str, metavar="<.kdb>", help="Two or more .kdb files")
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
    hierarchical_parser.add_argument("-m", "--method", type=str, choices=["singe", "complete", "average", "weighted", "centroid", "median", "ward"], default="ward", help="The method for linkage fitting in R to use")
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
    dist_parser.add_argument("kdb", nargs="+", type=str, metavar="<kdbfile1 kdbfile2 ...>", help="Two or more .kdb files")
    dist_parser.set_defaults(func=distances)


    index_parser = subparsers.add_parser("index", help="Create a index file that can be held in memory")
    index_parser.add_argument("-v", "--verbose", help="Prints warnings to the console by default", default=0, action="count")
    index_parser.add_argument("kdb", type=str, help="A k-mer database file (.kdb)")
    index_parser.add_argument("kdbi", type=str, help="Output index file (.kdbi)")
    index_parser.set_defaults(func=index_file)

    
    args=parser.parse_args()
    global logger
    logger = get_root_logger(args.verbose)

    sys.stderr.write("Constructed a logger for the program...\n")
    #logger.debug(sys.path)
    args.func(args)
    
