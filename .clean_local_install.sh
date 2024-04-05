#!/bin/bash
cat <<EOF
   Copyright 2020 Matthew Ralston

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
EOF




rm -rf /home/matt/.pyenv/versions/kdb/lib/python3.11/site-packages/kmerdb /home/matt/.pyenv/versions/kdb/lib/python3.11/site-packages/kmerdb-*.egg-info /home/matt/Projects/kdb/kmerdb.egg-info /home/matt/Projects/kdb/build /home/matt/Projects/kdb/dist

rm -rf /home/matt/.pyenv/versions/kdb/lib/python3.12/site-packages/kmerdb /home/matt/.pyenv/versions/kdb/lib/python3.12/site-packages/kmerdb-*.egg-info
rm -rf /home/matt/.pyenv/versions/3.11.7/envs/kdb/lib/python3.11/site-packages/kmerdb-*
rm -rf /home/matt/.pyenv/versions/3.12.2/envs/kdb/lib/python3.12/site-packages/kmerdb-0.7.8.dist-info/
cd /home/matt/Projects/kdb/
rm -rf dist build kmerdb.egg-info wheelhouse
cd

