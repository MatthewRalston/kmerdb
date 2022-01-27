#!/bin/bash



rm -rf /home/matt/.pyenv/versions/kdb/lib/python3.10/site-packages/kmerdb /home/matt/.pyenv/versions/kdb/lib/python3.10/site-packages/kmerdb-*.egg-info /ffast2/kdb/kmerdb.egg-info /ffast2/kdb/build /ffast2/kdb/dist
rm -rf /home/matt/.pyenv/versions/3.10.1/envs/kdb/lib/python3.10/site-packages/kmerdb-*
cd /ffast2/kdb/
rm -rf dist build kmerdb.egg-info wheelhouse
cd

