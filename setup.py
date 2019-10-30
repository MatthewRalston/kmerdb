#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from shutil import rmtree
from config import VERSION

from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = 'kdb'
DESCRIPTION = 'Yet another kmer library for Python'
long_description = 'See README.md for details'
URL = 'https://github.com/MatthewRalston/kdb'
EMAIL = 'mrals89@gmail.com'
AUTHOR = 'Matthew Ralston'
REQUIRES_PYTHON = '>=3.6.0'


# What packages are required for this module to be executed?
REQUIRED = [l.rstrip() for l in open('./requirements.txt', 'r')]
# What packages are optional?
EXTRAS = {
    'development': [l.rstrip() for l in open('./requirements-dev.txt', 'r')]
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
    #packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    # If your package is a single module, use this instead of 'packages':
    py_modules=['mypackage'],
    scripts=['bin/myscript'],
    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='GPLv3+',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy'
    ],

)
