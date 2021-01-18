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
import gzip
import yaml
import numpy as np
#import jsonschema
sys.path.append('..')

from kmerdb import fileutil


import logging
logger = logging.getLogger(__file__)



def open(filepath, mode="r", k=None, idx=None):
    if type(filepath) is not str:
        raise TypeError("kmerdb.index.open expects a str as its first positional argument")
    elif type(mode) is not str:
        raise TypeError("kmerdb.index.open expects a str as its second positional argument")
    elif "w" in mode and type(idx) is not np.ndarray:
        if "x" in mode and type(k) is not int: # Specifically, if its not a
            raise TypeError("kmerdb.index.open in create mode expects an int as its first *args argument")
            
    elif "w" in mode and type(k) is not int:
        raise TypeError("kmerdb.index.open in write mode expects an int as its second *args argument")
    modes = set(list(mode.lower()))
    if modes - set("xrwbt") or len(mode) > len(modes):
        raise ValueError("invalid mode: {}".format(mode))

    creating = "x" in modes
    reading  = "r" in modes
    writing  = "w" in modes
    binary   = "b" in modes
    text     = "t" in modes

    if text and binary:
        raise ValueError("can't have text and binary mode simultaneously")
    elif not (creating or reading or writing):
        raise ValueError("must have exactly one or read/write")

    if "r" in mode.lower():
        return IndexReader(filepath) # Will attempt to read this as an index file
    elif "x" in mode.lower() or set(list("xw")).intersection(modes):
        logger.info("Building index from file and constructing (mode='x') new .kdbi file")
        return IndexBuilder(filepath, k) # Creates a new .kdbi exactly named for the kdb file
    elif "w" in mode.lower():
        logger.info("Writing index to file instead of building index")
        write_index(idx, indexfile, k) # index, indexfile, k


class IndexReader:
    def __init__(self, indexfile:str):
        if type(indexfile) is not str:
            raise TypeError("kmerdb.index.IndexReader expects a str as its first positional argument")
        elif not os.path.exists(indexfile):
            raise TypeError("kmerdb.index.IndexReader expects an existing .kdbi index file as its first positional argument")
        elif not is_gz_file(indexfile):
            raise IOError("kmerdb.index.IndexReader expects a gzip compressed index file as its first positional argument")
        

        idx = np.array([], dtype="int64")
        #i = 0
        with gzip.open(indexfile, 'rt') as ifile:
            first_line = ifile.readline().rstrip()
            try:
                logger.debug("First line of index:\n{0}".format(first_line))
                k = int(first_line)
                idx = np.zeros(4**k, dtype="int64")
            except ValueError as e:
                logger.debug(e)
                raise IOError("Could not determine k by casting the first line of the index as an integer. Improperly formatted index")
            i = 0
            for line in ifile:
                # if i == 0:
                #     logger.warning("i: {0}, line:\n{1}".format(i, line))

                #     sys.exit(1)
                #     pass
                # else:
                kmer_id, line_offset = [int(i) for i in line.rstrip().split("\t")]
                idx[kmer_id] = line_offset
                i += 1
        self.index = idx

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        del self.index
        return

class IndexBuilder:
    def __init__(self, kdbfile:str, k:int):
        if type(kdbfile) is not str:
            raise TypeError("kmerdb.index.IndexBuilder expects a str as its first positional argument")
        elif not os.path.exists(kdbfile):
            logger.error("Could not find {0}".format(kdbfile))
            raise IOError("kmderdb.index.IndexBuilder expects an existing .kdb file as its first positional argument")
        elif os.path.splitext(kdbfile)[-1] != ".kdb":
            raise IOError("kmerdb.index.IndexBuilder was passed a non-kdb filetype as its first positional argument, should be .kdb")
        elif type(k) is not int:
            raise TypeError("kmerdb.index.IndexBuilder expects an int as its second positional argument")
        N = 4**k
        line_index = np.zeros(N, dtype="int64")
        #line_index = np.array([], dtype="int64")
        #line_index = set() # Can only use set transiently as it might sort
        indexfile = kdbfile + "i"
        
        with fileutil.open(kdbfile, mode='r') as kdbrdr:
            pos = kdbrdr.tell() # This could be the byte encoding, since we're reading it as plain text

            logger.info("Reading the .kdb file '{0}' in mode 'r' to generate the index...".format(kdbrdr._filepath))
            logger.info("Top position in the file after metadata parsing is said to be 'with this as kdbrdr:' {0}".format(pos))
            logger.info("KDBReader logs the position in mode 'r' as 'header_offset', after appending blocks: {0}".format(kdbrdr.metadata["header_offset"])) # This could be a compressed offset
            logger.info("Asked for the position one more time in mode 'r' in 'with this' and received: {0}".format(kdbrdr.tell()))
            line_index[0] = kdbrdr.tell()

            # SOmething should be checked here
            #kdbrdr.seek(kdbrdr.metadata["header_offset"])
            #kdbrdr._load_block()
            i = 1
            this_line = line_index[0]
            logger.debug("0th index, points to the first row? : {0}".format(line_index[0]))
            for line in kdbrdr:
                try:
                    next_line = kdbrdr.tell()
                    kmer_id, count = [int(i) for i in line.rstrip().split("\t")[:-1]]
                    line_index[kmer_id] = this_line
                    kdbrdr.seek(this_line)
                    new_line = kdbrdr.readline()
                    if new_line == line:
                        logger.debug("i: {0}, K-mer id: {1}, count: {2}, index line: {3}, offset: {4}".format(i, kmer_id, count, i, this_line))
                    elif next_line == 0 or this_line == 0:
                        logger.warning("i: 0, k-mer id: {1}, this_line: {2}, next_line{3}".format(i, kmer_id, this_line, next_line))
                        raise ValueError("Encountered a line where kdbrdr.tell() is producing 0")
                    else:
                        logger.warning("This line was not reproducibly generated...")
                        logger.warning("i: {0}, K-mer id: {1}, offset: {2}".format(i, kmer_id, this_line))
                        sys.exit(1)
                        
                    this_line = next_line
                    if kdbrdr.tell() == 0:
                        raise IOError("kmerdb.index.IndexBuilder expects the kdbrdr to increment with each line")
                    # elif i > 10:
                    #     logger.error("Position in file: {0}".format(pos))
                    #     sys.exit(1)
                except IndexError as e: # The final 
                    logger.warning("{0}:\t{1}".format(i, line.rstrip()))
                    raise e # STILL DONT KNOW WHATS HAPPENING HERE
                i += 1
        write_index(line_index, indexfile, k)

        
def has_index(kdbfile):
    if type(kdbfile) is not str:
        raise TypeError("kmerdb.index.has_index expects a str as its first positional argument")
    elif not os.path.exists(kdbfile):
        raise IOError("kmerdb.index.has_index received a filename that does not exist")
    elif os.path.splitext(kdbfile)[-1] == ".kdb" and not os.path.exists(kdbfile + "i"):
        return False
    elif os.path.splitext(kdbfile)[-1] == ".kdb" and os.path.exists(kdbfile + "i"):
        return True
    elif os.path.splitext(kdbfile)[-1] != ".kdb":
        raise ValueError("unexpected file extension, not a .kdb file.")
    else:
        raise RuntimeError("kmerdb.index.has_index encountered unexpected user inputs:\nkdbfile: '{0}', exists: {1}".format(kdbfile, os.path.exists(kdbfile)))


def read_line(kdbrdr:fileutil.KDBReader, kdbidx:IndexReader, kmer_id):
    if not isinstance(kdbrdr, fileutil.KDBReader):
        raise TypeError("kmerdb.index.read_line expects a kmerdb.fileutil.KDBReader as its first positional argument")
    elif not isinstance(kdbidx, IndexReader):
        raise TypeError("kmerdb.index.read_line expects a kmerdb.index.IndexReader as its second positional argument")
    elif type(kmer_id) is not int:
        raise TypeError("kmerdb.index.read_line expects an int as its third positional argument")
    if kdbidx.index[kmer_id] == 0:
        logger.warning("Found an index value without an offset, in fact, it's in the header")
        return None, None, None
        #raise ValueError("Invalid index value 0")
    else:
        kdbrdr.seek(kdbidx.index[kmer_id])
        kmer_id, count, metadata = kdbrdr.read_line()
        return kmer_id, count, metadata

        
def is_gz_file(filepath):
    try:
        gzip
        
        with gzip.open(filepath, 'rb') as f:
            f.readline()
        return True
    except OSError as e:
        return False

def write_index(index:np.array, indexfile:str, k:int):
    if not isinstance(index, np.ndarray):
        raise TypeError("kmerdb.index.write_index expects a str as its first positional argument")
    elif type(indexfile) is not str:
        raise TypeError("kmerdb.index.write_index expects a str as its second positional argument")
    elif os.path.exists(indexfile):
        logger.warning("kmerdb.index.write_index is overwriting an existing index '{0}'...".format(indexfile))
        mode = "wt"
    elif type(k) is not int:
        raise TypeError("kmerdb.index.write_index expects an int as its third positional argument")
    else:
        mode = "xt"
    with gzip.open(indexfile, mode) as ofile:
        ofile.write("{0}\n".format(k))
        for i, line_offset in enumerate(index):
            if int(line_offset) == 0:
                logger.debug("K-mer {0} was not observed, not writing to file".format(i))
                pass
            else:
                ofile.write("{0}\t{1}\n".format(i, line_offset))

