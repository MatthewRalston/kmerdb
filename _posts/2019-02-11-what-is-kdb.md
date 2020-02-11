---
title: 'What is .kdb?'

layout: nil
---

The K-mer database is a file format and command-line utility (CLI) for accessing precomputed k-mers from sequencing data.

The CLI tallies all k-mers from input fasta and fastq files.

The resulting file is a blocked GNU-zip file that can be indexed with the index command.

The Python module and CLI can then retrieve information from the indexed k-mer database quickly to perform calculations.

See the commands below for an idea of what functionality is built in to kdb.
