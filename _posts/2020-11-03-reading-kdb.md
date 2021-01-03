---
category: API
title: 'Reading .kdb files'
path: '/api_reading_kdb_files'

layout: nil
---


The next obvious topic is reading kdb files. It's very easy to read them line by line, so you should have no problems inventing your own incremental and custom python distance metrics for a reader that is capable of compression but uses a generator style. Note that the lines are parsed raw for speed, and they don't form any memory intensive objects. You are free to make the processing of each line as complex or simple as you want.

```python

from kdb import fileutil

files = ["foo.kdb", "bar.kdb", "baz.kdb"]

objects = [fileutil.open(f, 'r') for f in files]

#full_matrix = [o.slurp() for o in objects]

>for o in objects:
.  for line in o:
.    print(line) # Calculate something with each k-mer count
0    56
1    24

```
