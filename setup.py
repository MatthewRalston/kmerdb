#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
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

'''

import io
import os
import sys
from shutil import rmtree
sys.path.append(".")

try:
    from setuptools import find_packages
    from setuptools import setup
    from setuptools import Command
    from setuptools.extension import Extension

except ImportError:
    sys.exit(
        "We need the Python library 'setuptools' and 'distutils' to be installed."
        "Try running: python -m ensurepip"
    )

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    sys.exit(
        "We need the Python library 'cython' to be installed."
        "Try running: pip install cython or conda install cython"
        )

try:
    import numpy as np
except ImportError:
    sys.exit("Sadly, the extra feature of the Cython correlation function only works if NumPy is preinstalled, like in a conda environment.")
    # If you're here, welcome to my program. Work with me here. What are you looking for? Pollute my issues.
    # I can get rid of this NumPy requirement if I want, but then the bots will bias my dl numbers. Kthx.
    
if sys.version_info[:2] < (3, 7):
    sys.stderr.write(
        "KMERDB is tested on Python 3.7 or later. "
        "Python %d.%d detected. \n" % sys.version_info[:2]
    )
    sys.exit(1)


# class test_biopython(Command):
#     """Run all of the tests for the package.
#     This is a automatic test run class to make distutils kind of act like
#     perl. With this you can do:
#     python setup.py build
#     python setup.py install
#     python setup.py test
#     """

#     description = "Automatically run the test suite for Biopython."
#     user_options = [("offline", None, "Don't run online tests")]

#     def initialize_options(self):
#         """No-op, initialise options."""
#         self.offline = None

#     def finalize_options(self):
#         """No-op, finalise options."""
#         pass

#     def run(self):
#         """Run the tests."""
#         this_dir = os.getcwd()

#         # change to the test dir and run the tests
#         os.chdir("test")
#         sys.path.insert(0, "")
#         import run_tests
        
#         if self.offline:
#             run_tests.main(["--offline"])
#         else:
#             run_tests.main([])
            
#         # change back to the current directory
#         os.chdir(this_dir)


def can_import(module_name):
    """Check we can import the requested module."""
    try:
        return __import__(module_name)
    except ImportError:
        return None
                
# Package meta-data.
NAME = 'kmerdb'
DESCRIPTION = 'Yet another kmer library for Python'
long_description = 'See README.md for details'
URL = 'https://github.com/MatthewRalston/kmerdb'
CURRENT_RELEASE = "https://github.com/MatthewRalston/kmerdb/archive/v0.7.1.tar.gz"
EMAIL = 'mrals89@gmail.com'
AUTHOR = 'Matthew Ralston'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = "0.7.1"
KEYWORDS = ["bioinformatics", "fastq", "fasta", "k-mer", "kmer", "k-merdb", "kmerdb", "kdb"],
CLASSIFIERS = [
    "Development Status :: 1 - Planning",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Apache Software License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules",
]
# What packages are required for this module to be executed?

#REQUIRED = [l.rstrip() for l in open('./requirements.txt', 'r')]
#REQUIRED.remove("-e git://github.com/MatthewRalston/ecopy.git#egg=ecopy")
#REQUIRED.append("ecopy @ git+https://github.com/MatthewRalston/ecopy@master")


# What packages are optional?
# EXTRAS = {
#     'development': [l.rstrip() for l in open('./requirements-dev.txt', 'r')]

#     # 'fancy feature': ['django'],
# }
if can_import('numpy') is not None:
    import numpy as np
    extensions = [
        Extension("kmerdb.distance", ["kmerdb/distance.pyx"], include_dirs=[np.get_include()], define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],),
    ]
    # Where the magic happens:
    setup(
        name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type='text/markdown',
        author=AUTHOR,
        author_email=EMAIL,
        python_requires=REQUIRES_PYTHON,
        url=URL,
        download_url=CURRENT_RELEASE,
        keywords = KEYWORDS,
        classifiers=CLASSIFIERS,
        packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
        package_dir={'kmerdb': 'kmerdb'},
        package_data={'kmerdb': ['CITATION']},
        # If your package is a single module, use this instead of 'packages':
        #py_modules=['kmerdb'],
        #scripts=['bin/kmerdb', 'bin/kmerdb_report.R'],
        entry_points={
            'console_scripts': ['kmerdb=kmerdb:cli'],
        },
        #install_requires=REQUIRED,#['Cython==0.29.21', 'numpy==1.18.1'],
        #extras_require=EXTRAS,
        include_package_data=True,
        license='Apache-2.0',
        test_suite='test',
        #    tests_require=['mamba', 'expect'],
        #cmdclass={'build_ext': build_ext},
        #ext_modules=cythonize(extensions),
        library_dirs=["."],
        zip_safe=False,
    )


else:
    print("WARNING: Disabling the pyx distances")
    # Where the magic happens:
    setup(
        name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type='text/markdown',
        author=AUTHOR,
        author_email=EMAIL,
        python_requires=REQUIRES_PYTHON,
        url=URL,
        download_url=CURRENT_RELEASE,
        keywords = KEYWORDS,
        classifiers=CLASSIFIERS,
        packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
        package_dir={'kmerdb': 'kmerdb'},
        package_data={'kmerdb': ['CITATION']},
        # If your package is a single module, use this instead of 'packages':
        #py_modules=['kmerdb'],
        #scripts=['bin/kmerdb', 'bin/kmerdb_report.R'],
        entry_points={
            'console_scripts': ['kmerdb=kmerdb:cli'],
        },
        #install_requires=REQUIRED,#['Cython==0.29.21', 'numpy==1.18.1'],
        #extras_require=EXTRAS,
        include_package_data=True,
        license='Apache-2.0',
        test_suite='test',
        #    tests_require=['mamba', 'expect'],
    )
