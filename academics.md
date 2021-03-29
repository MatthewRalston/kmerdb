---
title: For Academics
layout: post
description: 'Quick introduction to Alignment-free methods'
image: assets/images/dna2.gif
nav-menu: true
toc: true
---

* [Introduction](#introduction)
    * [What are k-mer spectra](#what-are-k-mer-spectra)
	* [Alignment-free methods are futureproof](#alignment-free-methods-are-futureproof)
	* [Linear vs quadratic runtime](#linear-vs-quadratic-runtime)
* [Key Research Questions](#key-research-questions)
* [Numbers](#datasets-and-numbers)
* [Technical Specifications](#technical-specifications)
* [Objectives](#objectives)
    * [Graph database](#graph-database)
	* [GPU regression/ML](#tensorflow-support-for-ml-on-gpu)
	* [Algorithm performance](#benchmarking-performance)
	* [Phylogenomics and Metagenomics](#phylogenomics-and-metagenomics)
	* [Profile fidelity](#profile-fidelity) # subsampling/simulated fidelity to the true dataset
	* [Sensitivity analysis](#sensitivity-analysis) # Anvar et al., k = 0-12
	* [Homopolymer treatment](#homopolymer-content) # How to separate or identify k-mer ids that are low-quality
    * [Markov probabilities](#markov-probability)
	* [Dashboards](#dashboards)
* [Future Work](#future-work)

About a 5-10min read.

# Abstract

A brief introduction to the biomathematical topics of k-mer counting, motif abundance, alignment-free quantification and variant analysis, and alignment-free linear runtime are presented to the reader to frame a discussion about common formats for k-mer counts and deBruijn graphs. The k-mer counting technique bypasses the rate-limiting quadratic step in the topic of gene expression to provide faster results for a more rapid, on-line expression measurement experience. Eventually, an on-line real-time gene-expression measurement method may treat genomic sequences as abstract spectra through k-mer counting, much like certain high-throughput chemical and biochemical platforms are capable of analyzing complex, multiplexed spectra, and/or 'unknown' spectra to determine closest known chemical structures.

To this end, we discuss key research topics such as the juxtaposition of alignment and assembly throughout the analytical process, constructive interference of spectra in multiplexed samples, linear runtimes, graph databases, alignment-free sequence similarity, alignment-free sequence distance metrics usable for machine learning, and normalization practices.

To achieve all of these goals, the reader is asked to consider their academic and/or professional networks to locate an advisor who is capable of guideing a student through their first grant-writing and/or fellowship application process. Potentially, I would like to stay for both a masters and a PhD, so this would grant considerable time to manage the scope and timelines of the project.




# Introduction


Here we discuss the underlying bioinformatic techniques behind both alignment and assembly and use this discussion to present a novel alignment-free framework for metagenomic species decomposition, bacterial species inference, metagenome simulation, and Markov probabilities.

A principle challenge in the acceleration and automation of bioinformatic analyses based on Illumina Next Generation Sequencing (NGS) and other platforms is the abstraction of conventional alignment-based methods into an alignment-free analytical space. Bioinformatics is expected to be rate-limiting in the modern biological discovery process, with the recent dramatic reductions in the cost of sequencing leading to mountains of raw data for analysis. Quick, efficient, and accurate methods are needed to identify sequences, quantify abundances, estimate mutational likelihoods, and explore the sequence similarity space.

The first advantage of the so-called "alignment-free" methods is the flexibility to leverage reference sequences when the information is available, but to not restrict the effective utilization rate of the dataset to that which can be aligned with absolute certainty. Now we advocate for data structures, algorithms, and concepts of identity that extends beyond the concepts of reference sequences and alignment.

The second advantage of foregoing the expensive quadratic mapping step (which is required in conventional alignment-based sequencing pipelines) is the comparable accuracy with which expression measurements (Kallisto) and mutations can be called (e.g. Kevlar, Cell 2019). 

Here I provide a brief primer on the two branches of nucleic acid sequencing biology: assembly and alignment. The goal is to familiarize the reader with the essentials of sub-sequence biology and discuss a framework for the investigation of genomic spectra.

## What are K-mer spectra

The first step of many bioinformatic algorithms is to generate the counts of all length k subsequences in the dataset. This vector may be thought of as the k-mer spectra at that choice of k. This k-mer count histogram can be thought of as a kind of genomic or transcriptomic spectrum, a charactristic curve describing the dataet in terms of both the species of interest and the abundance of certain subsequences, which could refelect amplicon targeting in the case of exome sequence, expressed regions in the case of (meta)transcriptomics, or the ratios of composition of microbiomes in the case of metagenomics.

## Alignment-free methods are futureproof

Alignment-free methods are more flexible in their sequencing inputs, in a figurative sense. *Specifically, k-mer spectra are likely to retain their sequence identity in the face of indel-producing sequencing technologies.* Because alignment is not considered, the "errors" from indels are encorporated into inferential errors through a less harmful scoring mechanism than in sequence alignment, and thus a larger quantity of data can be utilized for indel-rich sequencing platforms like Oxford Nanopore, PacBio, and other similar technologies. Additionally, k-mer databases in particular are the backbone of the two major sequencing data applications, sequence alignment (e.g. seed regions, minimizers, etc.) and assembly (e.g. deBruijn graphs). For this reason, the topic of alignment-free methods and in particular, the k-mer spectra abstraction layer, has received considerable attention in the literature for its efficiency.

## Linear vs quadratic runtime

Alignment-free methods can produce quantifications faster and with comparable accuracy to genomic/transcriptomic alignment-based quantification. Several challenges remain including the choice of correct distribution and statistical framewok for modeling and hypothesis testing, Of course, it would be preferable to correctly characterize and then perhaps transform, with appropriate assumptions, the discrete count value into an alternative variable space with conditions suitable for traditional linear modeling, along with the imposed normality assumptions.

How is this more efficient than alignment? Generally, alignment is quadratic in the length of the sequences being aligned. More specifically, the alignment of Illumina based WGS, RNA-Seq, or metagenomic datasets is L*O(nm), where L is the fairly constant read length of the dataset, n is the number of alignable reads, and m is the length of the reference sequence(s). In contrast, counting k-mers is linear, L*O(n), in the number of sequence reads in the case of Illumina data or L*O(m) in the case of fasta input where m is the combined length of chromosomes in the sample's fasta file.

# What key research questions are you trying to answer?

There are multiple candidate research questions that haven't been answered by the existing literature on the topic of k-mer spectra for genomics applications. The first is a general purpose questions, namely "are k-mer spectra useful for calculating sequence similarities?" If we can demonstrate use cases for sequence similarity, then we can provide estimates about the likelihood that an unknown sequence belongs to a given genome, metagenome, transcriptome, etc. using a simple Markov model.

Also, the spectra and their additive nature permits simple deconstruction into known components with linear regression. This could be used to determine the fractions of each species in artificial and/or characterized metagenomes. This would be useful in cases where a researcher wants to quickly identify the fractions of known species that make up the sample. The regression coefficiencts would provide users with estimates of how complete the deconvolution is into the component species. 

Another key finding that I hope to produce with this software, is a graphic from Anvar et al. in the kPAL study that suggests that a choice of k as small as 10 may be sufficient to distinguish between microbial genomes. They reached this conclusion by graphing the median distance between the metagenomes and their corresponding shuffled/randomized sets. The graph displayed a sigmoid shape with the center around k=7. A choice of k=10 was obviously in the plateau and it could be said that this is a type of sensitivity analysis. I'd like to reproduce this figure to confirm what distance is at least sufficient for distinguishing between microbial metagenomes.

The last and perhaps most important question is what to do to best model the multivariate distribution that generates genes or other subsequences from the genomic sequence. Is this really only a bivariate relationship between the k-mer frequencies and the sequence length? Are those the only necessary generating factors?

# Datasets and Numbers

In this next section, I'd like to provide some back-of-the-envelope estimates of memory and space requirements for the hardware used to process the data. Let's start with the assumption that supplementary figure 1 from Anvar et al provides the correct estimate for a minimum choice of k of between k=8:10. 

A 32-bit floating point frequency vector for example, or a 32-bit integer count vector would yield about 0.537Gb of data, per-profile. If you wanted to performa a reasonable regression for some simple microbiological contrasts, you would need hundreds of data points, for potentially as much as > 50 GB of data for a regression model and easily into the terabytes for a machine learning level dataset.

For example Zhang et al. (Microbiome 2020) uses a dataset of 155,810 prokaryotic genomes from refseq. If we were to mmake that many profiles, by my estimate, it would yield about 83.65 TB of data. If we scaled back k to 11, it would take 2.09 TB of memory, more reasonable for our cluser size.


But if we asked a biologist what choice of k they would feel comfortable with, they would choose something similar in specificity to a PCR primer, with k between 22-25 and n=155,810, it would yield between 87 EB and 6 ZB of data.

But do these numbers hold up in the real world? Yes, and no. My old system 'argo' was using only 16 GB of memory and became quickly overwhelmed when analyzing larger choices of k (k=15 was intractable on the old system without significant freezing and system lock-up). My new system 'argo2', the specifications of which are listed below, can tolerate the analysis of 30 x 12-mer profiles with about 60-70% memory usage of 64 GB. It is worth noting that the entire k-mer profile is not held in memory at once, but gradually tallied on-disk. That efficiency aside, each Python process takes about 2-5% of system memory (> 2% on average), and most of this is likely Fasta datasets being read and processed, and active database connections handles used for on-disk k-mer counting.

# Technical specifications

The technical specifications for Argo2 are given below

| Component | Name                                 | Specification                                             |
| --------- | ------------------------------------ | --------------------------------------------------------- |
| CPU       | AMD Ryzen Threadripper 3960X 3.9 GHz | 24 core, 48 thread, 3.9GHz, 128Mb L3 cache                |
| RAM       | G.Skill Ripjaw V Series              | 128 GB (4x32 GB) DDR4 3600 MHz                            |
| GPU       | NVidia RTX A6000                     | 38.7-309.7 TFLOPS, 10752 tensor cores, 48 GB GDDR6 w/ ECC |
| AIC       | Gigabyte AORUS Gen4 PCIe NVMe AIC    | 8TB (4x2 TB) Corsair MP600 M.2 NVMe SSDs in RAID0(21 GB/s)|
| MOBO      | ASUS ROG Strix TRX-40e Gaming        | 3 x PCIex16 gen 4 slots, up to 256 GB RAM, 2 NVMe slots   |
| OS        | Arch Linux 5.11.1                    | Custom Arch linux distribution                            |


# Current state

There have been 3 pre-releases in the codebase thus far, and we are on version number 0.0.7. The codebase has changed into a sophisticated on-disk k-mer counting strategy, with multiple parallelization options. The first of which is native OS parallelism using something like GNU parallels to run the program on many genomes or metagenomes, simultaneously. The second parallelization option use the Python3 multiprocessing library, particularly for processing fastq datasets.

When I say on-disk, I mean the file format I've created is essentially an indexed database. This allows rapid access to sequential transition probabilities, on-disk, during a random walk representing a Markov process. This is a key feature for anyone who wants to seriously traverse a 4<sup>k</sup> dimensional space.

The codebase also currently contains a randomization feature, distance matrices, arbitrary interoperability with Pandas dataframes, clustering features, normalization, standardization, and dimensionality reduction. The suite is currently ready for a regression feature that I've promised, and I'd like to implement this early this Spring. Next, I'd be interested in working on the following features that would make the suite ready for another beta release. 



# Objectives

Let me first define two classes of deliverables. THe first is what I'll call tangibles. Trackable metrics including github commit frequencies, documentation quality, quality of unit testing, reports and progress updates (typically Rmd/latex for characterizations in R, or Jupyter notebooks for Python-based inferences). The second type of deliverable would be referred to as 'other', and includes things like academic updates, collaborations, progress on grant writing and/or manuscripts.

I'd like to communicate this project to an academician who would support my grant writing and potentially support me looking at a PhD program at the University of Delaware. A project of this size would likely see considerable adoption and several users have expressed interest on Reddit. Securing funding is the first step, and a preprint or early manuscript would go a long way to eventually securing a fellowship or similar funding for a PhD project.

I bring the following to the table: advanced study of microbial systems, metabolism, biomathematics, statistics, computer science, software engineering, biochemistry, molecular biology, and genomics. I am adept at building data processing pipelines, storage systems, reproducibility infrastructure, and creating custom reports in standardized formats such as Rmarkdown, LaTeX, and Jupyter. I believe I would be an asset to any laboratory and that this project has obvious merits for microbiologists, bioinformaticians, and mathematicians.

The following are specific goals and milestones for the project.


## Compositional Analysis

One goal of the suite is to produce a type of regression such that exactly n compositional coefficients are obtained. These are the best estimates we can have with the simplest model possible for how the k-mer counts reflect the various genomes that compose the metagenome.

## Graph database

The de Bruijn graph should have a RDF format, which is widely interconvertible, importable, and analyzable by enterprise-grade graph analysis suites (Neptune/SageMaker, Neo4J, etc.). Not only should assembly be done after RDF ingestion, but it should also be an interesting data structure to begin with.

## Tensorflow support for ML on GPU

The suite currently produces a variety of numpy/Pandas ingestible matrices in a pipeline, where the matrices can be normalized, run through PCA etc. on the fly. The end goal is to perform supervised machine learning on the dataset for the purposes of classification. If GPU support was optional, we could see considerable speed-ups on larger datasets.

## Benchmarking performance

It is essential that we compare the performance of this software with its predecessors, using both quantitative and qualitative metrics. For example, multi-profile distance matrices are not supported by any of the predecessor suites. Similarly, it is not posssible to ingest the data into machine learning models or graph engines on even the more sophisticated Jellyfish approach.

Comparing this software with its predecessors both quantitatively and qualitatively is a matter of choosing the most aggressive features and providing the most complete user experience. A deeper and more complete review/summary is underway and is what I would consider my first report. The first report concerns the accuracy of new findings and is considered more important. The second report concerns efficiency of memory and CPU performance.

## Phylogenomics and Metagenomics

In addition to novel characterizations of individual species and metagenomes capable with this approach, by virtue of inter-profile distance metrics, it is possible to deconstruct metagenomes in terms of their component species. If a metagenome contains known species, a regression model may be fit to the available data to determine the fraction of the metagenome that corresponds to each composing species.


## Profile fidelity

As with any study involving sequencing data in the most optimal scenarios, it goes without saying that sequencing errors must be considered. To address the effects of sequencing error on k-mer profiles we suggest two approaches. First, that only relatively high-quality reads should be used. Reads should be trimmed and would ideally align to a reference sequence on their own. We understand departures from this are the norm that we are trying to address, but it is fairly common for most methods, including the predecessors, to ask for high-quality sequencing reads.

The second strategy for characterizing the effect of sequencing error, is to investigate k-mer profile fidelity by comparing the distance between a simulated fastq dataset with basic error rates (substitution and indel errors) and the originating fasta sequence. By selecting an appropriate value of k, sampling fidelity should be high when enough reads are obtained.

## Sensitivity analysis

A key figure in the literature surrounding this topic is from Anvar et al. and is supplemental in the kPAL paper. They perform a sensitivity analysis on the choice of k, comparing simulated metagenomes with randomly shuffled/permuted metagenomes to determine what value of k is sufficient to distinguish the true metagenomes from the shuffled profiles, in terms of Euclidean distance. I'd like to reproduce this figure for a moderately extended value of k.

## Homopolymer content

There are several special types of k-mers that can be considered as part of the profile with particular interest. The first of which, the nullomers (the number of zero-count k-mers), are directly proportional to the theoretical number of observable k-mers (4<sup>k</sup>). An additional type of special k-mer are those which contain high proportions of homopolymer content. These repetitive k-mers are unlikely to be used in many genomes, as they can cause biochemical issues during replication and transcription, such as leaky polymerase activity. Once the choice of k is addressed by sensitivity analysis and profile fidelity above, then homopolymer content would be a naturally interesting subject area to pursue, perhaps as some simple metadata, or as just an area for graphical investigation and characterization of the distributions of homopolymers found in sequenced genomes.

## Markov probability

With each count and accompanying frequency, comes a possible question: what is the likelihood of this sequence given the data? The simplest way to model the nucleic acid sequence in this way is by using a first-order Markov model. A more detailed treatment is available upon request.

Check the following link for details on the Markov formulation. [I will cite thee](#https://matthewralston.github.io/kmerdb/quickstart#kmerdb-probability)


## Dashboards


The final component is a visualization suite, perhaps initially a Rmarkdown knitted as html, with interactive graphics and detailed descriptions of the topics, assumptions, formulae, variables, graphs etc. i.e. a literate program. Literate programming was a popular topic back in BMS, and I seek to practice my own abilities in a considerable way here. This could be a master report, just a simple pdf ready for BioRxiv. Or a more elaborate set of dashboards could be its own web service to analyze kdb output.


# Future Work



