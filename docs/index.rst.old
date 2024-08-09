
kmerdb
===================

- [ x ] [Homepage:](https://matthewralston.github.io/kmerdb)
- [ x ] [Quick Start guide](https://matthewralston.github.io/kmerdb/quickstart)
- [ x ] `kmerdb usage subcommand_name`
  - `profile` - Make k-mer count vectors/profiles, calculate unique k-mer counts, total k-mer counts, nullomer counts. Import to read/write NumPy arrays from profile object attributes.
  - `graph` - Make a weighted edge list of kmer-to-kmer relationships, akin to a De Bruijn graph.
  - `usage` - Display verbose input file/parameter and algorithm details of subcommands.
  - `help` - Display verbose input file/parameter and algorithm details of subcommands.
  - `view` - View .tsv count/frequency vectors with/without preamble.
  - `header` - View YAML formatted header and aggregate counts
  - `matrix` - Collate multiple profiles into a count matrix for dimensionality reduction, etc.
  - `kmeans` - k-means clustering on a distance matrix via Scikit-learn or BioPython with kcluster distances
  - `hierarchical` - hierarchical clustering on a distance matrix via BioPython with linkage choices
  - `distance` - Distance matrices (from kmer count matrices) including SciPy distances, a Pearson correlation coefficient implemented in Cython, and Spearman rank correlation included as additional distances.
  - `index` - Create an index file for the kmer profile (Delayed:)
  - `shuf` - Shuffle a k-mer count vector/profile (Delayed:)
  - `version` - Display kmerdb version number
  - `citation` - Silence citation suggestion
- [ x ] `kmerdb subcommand -h|--help`


k-mer counts from .fa(.gz)/.fq(.gz) sequence data can be computed and stored for access to metadata and count aggregation faculties. For those familiar with `.bam`, a `view` and `header` functions are provided. This file is compatible with `zlib`.

Install with `pip install kmerdb`


Please see the [Quickstart guide](https://matthewralston.github.io/kmerdb/quickstart) for more information about the format, the library, and the project.



Documentation
===================

:README.md: https://github.com/MatthewRalston/kmerdb
	    
:Mainpage: https://matthewralston.github.io/
	   
:QuickStart Guide: https://matthewralston.github.io/kmerdb/quickstart/
		   
:GitHub: https://github.com/MatthewRalston/kmerdb/
	 
:Issue Tracker: https://github.com/MatthewRalston/kmerdb/issues/
		
:TODO: https://github.com/MatthewRalston/kmerdb/blob/main/TODO.org
       
:License URI: https://github.com/MatthewRalston/kmerdb/blob/main/LICENSE.txt
	      
:License: Apache v2.0



IUPAC support
===================



::
  kmerdb profile -k $k -o output input.fa # This will simply discard k-mers containing non-IUPAC characters.



IUPAC residues (ATCG+RYSWKM+BDHV) are kept throughout the k-mer counting. But non-IUPAC residues (N) and characters are trimmed from the sequences prior to k-mer counting. Non-standard IUPAC residues are counted as doublets or triplets.


	  
## Development


::

  python setup.py test


## License

Created by Matthew Ralston - [Scientist, Programmer, Musician](http://matthewralston.github.io) - [Email](mailto:mralston.development@gmail.com)

Distributed under the Apache license. See `LICENSE.txt` for the copy distributed with this project. Open source software is not for everyone, and im the author and maintainer. cheers, on me. You may use and distribute this software, gratis, so long as the original LICENSE.txt is distributed along with the software. This software is distributed AS IS and provides no warranties of any kind.

::
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


## Contributing

1. Fork it (<https://github.com/MatthewRalston/kmerdb/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

## Acknowledgements

Thanks mom and dad and my hometown, awesome hs, and the University of Delaware faculty for support and encouragement. Thanks to my former mentors, bosses, and coworkers. It's been really enjoyable anticipating what the metagenomics community might want from a tool that can handle microbial k-mers well.

Thank you to the authors of kPAL and Jellyfish for the inspiration and bit shifting trick. And thank you to others for the encouragement along the way, who shall remain nameless. 

The intention is that more developers would want to add functionality to the codebase or even just utilize things downstream, but to build out directly with numpy and scipy/scikit as needed to suggest the basic infrastructure for the ML problems and modeling approaches that could be applied to such datasets. This project began under GPL v3.0 and was relicensed with Apache v2. Hopefully this project could gain some interest. I have so much fun working on this project. There's more to it than meets the eye. I'm working on a preprint, and the draft is included in some of the latest versions of the codebase, specifically .Rmd files.

More on the flip-side. It's so complex with technology these days...

.. kmerdb documentation master file, created by
   sphinx-quickstart on Thu Aug  8 00:32:04 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

kmerdb documentation
====================

Add your content using 'reStructuredText' syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: kmerdb.fileutil
   :members:

.. automodule:: kmerdb.graph
   :members:
      
.. automodule:: kmerdb.parse
   :members:

.. automodule:: kmerdb.kmer
   :members:

.. automodule:: kmerdb.util
   :members:

.. automodule:: kmerdb.index
   :members:

.. automodule:: kmerdb.logger
   :members:
  
.. automodule:: kmerdb.distance
   :members:

.. automodule:: kmerdb.python_distances
   :members:

.. automodule:: kmerdb.probability
   :members:

.. autoclass:: kmerdb.fileutil
   :members:

.. autoclass:: kmerdb.graph
   :members:
      
.. autoclass:: kmerdb.parse
   :members:

.. autoclass:: kmerdb.kmer
   :members:

.. autoclass:: kmerdb.util
   :members:

.. autoclass:: kmerdb.index
   :members:

.. autoclass:: kmerdb.logger
   :members:
  
.. autoclass:: kmerdb.distance
   :members:

.. autoclass:: kmerdb.python_distances
   :members:

.. autoclass:: kmerdb.probability
   :members:

.. autoexception:: kmerdb.fileutil
   :members:

.. autoexception:: kmerdb.graph
   :members:
      
.. autoexception:: kmerdb.parse
   :members:

.. autoexception:: kmerdb.kmer
   :members:

.. autoexception:: kmerdb.util
   :members:

.. autoexception:: kmerdb.index
   :members:

.. autoexception:: kmerdb.logger
   :members:
  
.. autoexception:: kmerdb.distance
   :members:

.. autoexception:: kmerdb.python_distances
   :members:

.. autoexception:: kmerdb.probability
   :members:




