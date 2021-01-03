---
category: API
title: 'Writing .kdb files'
path: '/api_writing_kdb_files'

layout: nil
---

Then we have the topic of writing kdb files. It's very easy again with the built in open method creating and returning instances of class KDBWriter (in this case). It's just a matter of iterating over your data and supplying it to the write method, which forwards to the Bio.bgzf write method through my wrapper class KDBWriter.


```python

from kdb import fileutil

with fileutil.open(f, 'w', header) as kdbfile:
  for x in list:
    kdbfile.write("\t".join(x) + "\n")

```
