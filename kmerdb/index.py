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


import sys
import os
import copy
import array
import gzip

#import jsonschema
sys.path.append('..')

from kmerdb import fileutil


import logging
logger = logging.getLogger(__file__)







# class Index:
#     def __init__(self, k, kdbfile:str=None, indexfile:str=None):
#         # if kdbrdr.__class__.__name__ != fileutil.KDBReader.__name__:
#         #     raise TypeError("kdb.index.IndexBuilder expects a KDBReader object as its first positional argument")
#         if type(k) is not int:
#             raise TypeError("kdb.index.Index expects an int as its first positional argument")
#         self.k = k

#         # 
#         if kdbfile is None and indexfile is None:
#             raise ValueError("kdb.index.Index initialized with no .kdb or index file specified...")            
#         elif kdbfile is not None and type(kdbfile) is not str:
#             raise TypeError("kdb.index.Index expects the keyword argument 'kdbfile', if specified, to be a str")
#         elif indexfile is not None and type(indexfile) is not str:
#             raise TypeError("kdb.index.Index expects the keyword argument 'indexfile', if specified to be a str")
#         else:
#             self.__indexfile = indexfile
#             self.__kdbfile = kdbfile
#             self.index = array.array('Q') # Can be set from the outside


def _write_line_index(indexfile, index):
    if type(indexfile) is not str:
        raise TypeError("kmerdb.index._write_index expects a str as its first positional argument")
    elif not isinstance(index, array.array):
        raise TypeError("kmerdb.index._write_index expects an array.array as its second positional argument")
    with gzip.open(indexfile, 'wt') as ofile:
        for kmer in index:
            ofile.write(str(kmer) + "\n")

                
def _read_line_index(indexfile):
    if type(indexfile) is not str:
        raise TypeError("kmerdb.index._read_index expects a str as its first positional argument")
    line_index = array.array('Q') #, range(4**self.k))
    #i = 0
    with gzip.open(indexfile, 'rt') as ifile:
        for line in ifile:
            #idx[i] = int(line.rstrip()) # FIXME
            idx.append(int(line.rstrip()))
    return idx

        
def build_line_index_from_kdb(kdbfile, k):
    if type(kdbfile) is not str:
        raise TypeError("kmerdb.index.build_line_index_from_kdb expects a str as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.index.build_line_index_from_kdb expects an int as its second positional argument")
    line_index = array.array('Q', range(4**k))
    with fileutil.open(kdbfile, mode='r') as kdbrdr:

        line_index[0] = kdbrdr.tell()
        i = 1
        for line in kdbrdr:
            # ptr = kdbrdr.tell()
            # if i < 10:
            #     print(i, ptr, "| \t", line.rstrip())
            # elif i > 4194300:
            #     print(i, ptr, "| \t", line.rstrip())
                
            try:
                line_index[i] = kdbrdr.tell()
                i+=1
            except IndexError as e: # The final 
                logger.debug("Index generation complete...")
                #logger.warning("{0}:\t{1}".format(i, line.rstrip()))
                #logger.error(e)
    return line_index


