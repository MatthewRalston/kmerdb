#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from shutil import rmtree
sys.path.append(".")

try:
    from setuptools import find_packages
    from setuptools import setup
    from setuptools import Command
    from setuptools import Extension
except ImportError:
    sys.exit(
        "We need the Python library 'setuptools' to be installed."
        "Try running: python -m ensurepip"
    )

if sys.version_info[:2] < (3, 7):
    sys.stderr.write(
        "KDB is tested on Python 3.7 or later. "
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
NAME = 'kdb'
DESCRIPTION = 'Yet another kmer library for Python'
long_description = 'See README.md for details'
URL = 'https://github.com/MatthewRalston/kdb'
EMAIL = 'mrals89@gmail.com'
AUTHOR = 'Matthew Ralston'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '0.0.1'

# What packages are required for this module to be executed?
#REQUIRED = [l.rstrip() for l in open('./requirements.txt', 'r')]
REQUIRED = [
    'biopython==1.74',
    'boto3==1.10.8',
    'botocore==1.13.8',
    'jsonschema==3.1.1',
    'matplotlib==3.1.3',
    'PyYAML==5.1.2',
    'SQLAlchemy==1.3.13',
    'Sphinx==2.2.0']

# What packages are optional?
EXTRAS = {
    #'development': [l.rstrip() for l in open('./requirements-dev.txt', 'r')]

    # 'fancy feature': ['django'],
}

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
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)"
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    # If your package is a single module, use this instead of 'packages':
    #py_modules=['kdb'],
    scripts=['bin/kdb'],
    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='GPLv3+',
    test_suite='test',
#    tests_require=['mamba', 'expect'],
)
