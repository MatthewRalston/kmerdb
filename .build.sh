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

# cd /ffast2/kdb/
# cd kmerdb
# #cythonize -a -i kmerdb/distance.pyx || echo "Could not Cythonize the Pearson correlation function. Refusing to build." # Requires cython to properly build a manylinux wheel.
# cd ..
# python setup.py sdist bdist_wheel
# auditwheel repair --plat linux_x86_64 dist/kmerdb-*linux_x86_64.whl
# mv wheelhouse/* dist
# #rm dist/*linux_x86_64.whl
# pip install dist/kmerdb-*-linux_x86_64.whl


# #python setup.py install || echo "'python setup.py install' died on your platform. Do you know why?"

# rm -rf wheelhouse





# cd

python -m build
auditwheel repair --plat manylinux2014_x86_64 dist/kmerdb-*linux_x86_64.whl
mv wheelhouse/* dist
rm dist/kmerdb-*linux_x86_64.whl
