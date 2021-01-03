---
category: 'What is .kdb?'
title: 'Development'
path: '/development'

layout: nil
---


[![PyPI version](https://img.shields.io/pypi/v/kdb.svg)][pip]
[![Python versions](https://img.shields.io/pypi/pyversions/kdb.svg)][Pythons]
[![Travis Build Status](https://travis-ci.org/MatthewRalston/kdb.svg?branch=master)](https://travis-ci.org/MatthewRalston/kdb)
[![Coveralls code coverage](https://img.shields.io/coveralls/MatthewRalston/kdb/master.svg)][Coveralls]
[![ReadTheDocs status](https://readthedocs.org/projects/kdb/badge/?version=stable&style=flat)][RTD]


[pip]: https://pypi.org/project/kdb/
[Pythons]: https://pypi.org/project/kdb/
[Coveralls]: https://coveralls.io/r/MatthewRalston/kdb?branch=master
[RTD]: https://kdb.readthedocs.io/en/latest/

## Travis-CI pytest unit tests
The method for installation and unit tests in a new development environment can be seen in '.travis.yml'. Currently, only unit tests are available for the suite. Acceptance testing has not been implemented. Until then, all acceptance testing is done manually prior to commits, squashes, rebases, and merges. Unit tests may be run with the following:

```bash
python setup.py test
```

