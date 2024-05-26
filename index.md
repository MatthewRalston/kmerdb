---
layout: home
title: Home
landing-title: 'A simple CLI and module for k-mer analysis'
description: A bgzf file format for the analysis and storage of k-mer abundances.
image: null
author: null
show_tile: false
---

The K-mer database profile `kmerdb` is a file format (.kdb) and command-line utility (CLI) for creating k-mer count vector profiles from sequencing data. kmerdb views k-mers and their De Bruijn graphs as fundamental tools for biological sequence abundances and sequence similarities. K-mer methods have been noted as fast and are used in areas from NGS quality control to phylogenomic investigations.

kmerdb is based on the block GNU-zip file (bgzf) standard. Each .kdb file has a metadata header section, much like .bam files. It is essentially a tab-delimited format with a YAML header. Input file checksums, unique and total k-mer counts, and nullomer counts are stored in the metadata block at the top of the file.

Please visit the [Install](#/install) page for details on installation.

See the commands in the [Usage](#/usage) section for an idea of what functionality is built in to kdb.

See the original blog post on the concept [here](https://matthewralston.github.io/blog/kmer-database-format-part-1).

