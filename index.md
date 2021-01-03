---
layout: home
title: Home
landing-title: 'A simple CLI and module for k-mer analysis'
description: A bgzf file format for the analysis and storage of k-mer abundances.
image: null
author: null
show_tile: false
---

The K-mer database (.kdb) is a file format and command-line utility (CLI) for accessing precomputed k-mers from sequencing data. .kdb views k-mer graph databases as a fundamental technology to look at biological sequence identities and abundances. K-mer methods have been noted as fast and are used in areas from NGS quality control to phylogenomic investigations.

.kdb is based on the block GNU-zip file (bgzf) standard. Each .kdb file has a header or metadata section, much like .bam files. It is essentially a tab-delimited format with the last column unstructured for k-mer specific metadata. Input files and total k-mer counts are stored in the metadata block at the top of the file

Please visit the [Install](#/install) page for details on installation.

See the commands in the [Usage](#/usage) section for an idea of what functionality is built in to kdb.

See the original blog post on the concept [here](https://matthewralston.github.io/blog/kmer-database-format-part-1).

