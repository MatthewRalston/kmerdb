#!/bin/bash

cd /ffast2/kdb/
cd kmerdb
#cythonize -a -i kmerdb/distance.pyx || echo "Could not Cythonize the Pearson correlation function. Refusing to build." # Requires cython to properly build a manylinux wheel.
cd ..
python setup.py sdist bdist_wheel
/bin/auditwheel repair --plat manylinux2014_x86_64 dist/kmerdb-*linux_x86_64.whl
mv wheelhouse/* dist
rm dist/*linux_x86_64.whl
pip install dist/kmerdb-*-manylinux2014_x86_64.whl


#python setup.py install || echo "'python setup.py install' died on your platform. Do you know why?"

rm -rf wheelhouse





cd
