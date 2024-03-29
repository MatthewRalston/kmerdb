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
# Also be sure to check out the 'kmerdb profile --help' page if you're not familiar with the command
parallel 'kmerdb profile --keep-sqlite -k $K {} {.}.$K.kdb' ::: $(/bin/ls ../../test/data/*.fasta.gz)
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
kmerdb matrix Normalized test/data/*.$K.kdb | kmerdb distance spearman > normalized_spearman_dist.tsv
kmerdb kmeans -k 3 Biopython < 8-mer_spearman_dist.tsv
mv kmeans_clustering_of_kmer_profiles.png kmeans_k3_clustering_on_spearman_dist.png
mv kmeans_elbow_graph.png kmeans_k3_elbow_graph_on_spearman_dist.png
kdb matrix Normalized test/data/*.$K.kdb | kmerdb distance correlation > normalized_pearson_dist.tsv
kdb kmeans -k 3 Biopython < 8-mer_pearson_dist.tsv
mv kmeans_clustering_of_kmer_profiles.png kmeans_k3_clustering_on_pearson_dist.png
mv kmeans_elbow_graph.png kmeans_k3_elbow_graph_on_spearman_dist.png


```


# Monitoring performance

```bash
sudo iotop
iostat -x 1
```


# RStudio

Then customize the report to the needs of your samples and clustering results, and be sure to notify me in an issue if you happen across a case where DESeq2 normalization doesn't make sense with the dataset.

The report should largely speak for itself at this point, as far as what analyses are possible with this framework for assessing fasta or fastq datasets, since clustering was done and a silhouette value is calculated with each clustering.





# References

I grabbed the Rmd template and other things from [dcnorris/preprint-template](https://github.com/dcnorris/preprint-template). Thanks to David Norris for this.



