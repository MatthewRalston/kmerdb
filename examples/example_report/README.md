# k-mer Report

In this example, we will demonstrate how to create an example report for the analysis of k-mer profiles per sample. In this example we will look at several members of the Clostridia, some of the same strain, and how they are related to B. subtilis and E. coli. For those who aren't familiar, Clostridia and B. subtilis are both gram positive, walled bacteria that are typical of the soil. E. coli in contrast is a gram negative enteric bacterium. I expect that the Clostridia all be closely related, while B. subtilis and E. coli would be more distantly related.


The next natural question is... what transformations are needed before analysis can begin? Well, before we even answer that, we should first answer the question: is the dataset worth analyzing? Some questions that follow from this are:

1. Are B. subtilis and E. coli separable in some way directly or indirectly using just k-mer counts?
2. Are the Clostridia kept separate from the other two in a meaningful way as well?
3. How is this result similar or different from a different dimensionality reduction and visualization technique?
4. Can normalization be meaningfully applied?


To answer these questions to some extent, let us look at applying  principal components analysis (PCA) and k-means clustering to see if we can teach k-means to recognize the expected number of labels we're interested in finding in the dataset.

So first we need to apply PCA via the `kdb matrix` command. Let's look at what to do first


## Populate a new working directory

```bash
mkdir examples/example_report2/
cd examples/
cp example_report/index.Rmd example_report2/
cd example_report2/
# Create a new set of profiles for k=$K in the same directory as the input files.
# Please see the manual for GNU parallel if you're not familiar with this usage.
parallel 'kmerdb profile -k $K {} {.}.$K.kdb' ::: $(/bin/ls ../../test/data/*.fasta.gz)
# Generate unnormalized and DESeq2 normalized count matrices
kmerdb matrix Unnormalized test/data/*.$K.kdb > unnormalized_count_matrix.tsv
kmerdb matrix Normalized test/data/*.$K.kdb > normalized_count_matrix.tsv


# Run the PCA and tSNE with k-means clustering
kmerdb matrix PCA # Generate the elbow graph and choose the appropriate version of '-n'
kmerdb matrix PCA -n $N | kdb kmeans -k 3 Biopython
mv kmeans_clustering_of_kmer_profiles.png kmeans_k3_clustering_on_pca3.png
mv kmeans_elbow_graph.png kmeans_k3_elbow_graph_on_pca3.png
kmerdb matrix tSNE -n 2 | kdm kmeans -k 3 Biopython
mv kmeans_clustering_of_kmer_profiles.png kmeans_k3_clustering_on_tsne2.png
mv kmeans_elbow_graph.png kmeans_k3_elbow_graph_on_tsne2.png
# Generate correlation distance matrices for clustering
kmerdb matrix Normalized test/data/*.$K.kdb | kmerdb distance spearman > 8-mer_spearman_dist.tsv
kmerdb kmeans -k 3 Biopython < 8-mer_spearman_dist.tsv
mv kmeans_clustering_of_kmer_profiles.png kmeans_k3_clustering_on_spearman_dist.png
mv kmeans_elbow_graph.png kmeans_k3_elbow_graph_on_spearman_dist.png
kdb matrix Normalized test/data/*.$K.kdb | kmerdb distance correlation > 8-mer_pearson_dist.tsv
kdb kmeans -k 3 Biopython < 8-mer_pearson_dist.tsv
mv kmeans_clustering_of_kmer_profiles.png kmeans_k3_clustering_on_pearson_dist.png
mv kmeans_elbow_graph.png kmeans_k3_elbow_graph_on_spearman_dist.png

# Generate the html report with R. NOTE: not implemented yet.
kmerdb report --unnormalized unnormalized_count_matrix.tsv --normalized normalized_count_matrix.tsv --pca-elbow-graph PCA_variance_accumulation.png --kmeans-pca-clustering kmeans_k3_clustering_on_pca3.png --kmeans-tsne-clustering kmeans_k3_clustering_on_tsne2.png --kmeans-spearman-clustering kmeans_k3_clustering_on_spearman_dist.png --kmeans-on-pca-elbow-graph kmeans_elbow_graph_on_pca3.png --kmeans-on-tsne-elbow-graph kmeans_elbow_graph_on_tsne2.png --kmeans-on-spearman-dist-elbow-graph kmeans_elbow_graph_on_spearman_dist.png --kmeans-on-pearson-dist-elbow-graph kmeans_elbow_graph_on_pearson_dist.png
```


