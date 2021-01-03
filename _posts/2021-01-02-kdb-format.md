---
category: 'What is .kdb?'
title: '.kdb format'
path: '/kdb_format'

layout: nil
---

The k-mer database format is rather simple. It contains a metadata section, followed by a tab delimited format with an unstructured JSON as the final column. Both sections are combined in a compressed layer facilitated by Biopython's bio.bgzf module. 

Each file can be inspected with the view and header commands detailed in the [section below](#/view).
