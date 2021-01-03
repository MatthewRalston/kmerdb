---
category: Usage
title: 'kdb distance'
path: '/kdb_distance'

layout: nil
---


Suppose you want a distance matrix between profiles; this is made easy with the distance command

```bash
>./bin/kdb distance correlation test/data/*.$K.kdb
```

The result is a symmetric matrix in tsv format with column headers formed from the filenames minus their extensions. It is presumed that to properly analyze the distance matrix, you would name the files after their sample name or their species, or some other identifying features.

