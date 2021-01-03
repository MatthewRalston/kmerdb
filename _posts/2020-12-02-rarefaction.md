---
category: Usage
title: 'kdb rarefy'
path: '/kdb_rarefy'

layout: nil
---


Suppose you want supply a normalized matrix directly to [ecopy](https://github.com/Auerilas/ecopy) for rarefaction analysis. The rarefaction command requires a tsv and may be piped directly from 'kdb matrix'.


```bash
# You may supply any arbitrary count tsv to rarefy, and it will run that directly.
>./bin/kdb rarefy example.tsv

# Alternatively, you may pipe the result from kdb matrix directly to kdb rarefy.
>./bin/kdb matrix Unnormalized test/data/*.$K.kdb | ./bin/kdb rarefy
```
