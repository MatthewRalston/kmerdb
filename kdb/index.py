import sys
import os
import copy
import json
import time

import jsonschema
from Bio import bgzf
sys.path.append('..')

from kdb import fileutil

import logging
logger = logging.getLogger(__file__)


class IndexBuilder:
    def __init__(self, kdb: fileutil.KDB, kdbrdr: fileutil.KDBReader):
        if kdb.__class__.__name__ != fileutil.KDB.__name__:
            raise TypeError("kdb.index.IndexBuilder expects a KDB object as its first positional argument")
        elif kdbrdr.__class__.__name__ != fileutil.KDBReader.__name__:
            raise TypeError("kdb.index.IndexBuilder expects a KDBReader object as its second positional argument")
        if kdbrdr.loaded is False: # Note: this can be overridden without the data being loaded
            raise TypeError("kdb.index.IndexBuilder expects a fully loaded KDBReader object as its first positional argument")
        self.kdb = kdb
        self.__kdbfile = kdbrdr.filepath

        self.line_index = self._index_lines()
        #self.tag_index  = self._index_tags()

    def _index_lines(self):
        result = []
        with open(self.__kdbfile, 'rb') as ifile:
            self.kdbrdr = fileutil.KDBReader(ifile)

            # Read by block
            offsets = copy.deepcopy(self.kdbrdr._offsets)
            # Omit the header block
            offsets.pop(0)
            # FIXME
            offsets.pop(0)
            i = 0 # Block index
            j = 0 # k-mer  index
            body_start, block1_length, block1_data_start, block1_data_length = offsets.pop(0)
            offset = None
            if block1_data_length > 0: # If the block is not empty, make a virtual offset
                ######      Build+go offset
                offset = bgzf.make_virtual_offset(body_start, 0) # Find the offset of that block
                self.kdbrdr._bgzfrdr.seek(offset)       # Seek the start of the block
                ######      
                line = self.kdbrdr._bgzfrdr.readline().rstrip()  # Read the first line from the block
                while len(line) > 0: #and json.loads(line.split("\t")[2]):
                    #result.append((j, i, raw_start, offset))       # Push (k-mer, block, block_offset, offset) onto the stack
                    #print(i, j, line.rstrip().split("\t")[0:2])
                    #time.sleep(0.5)
                    result.append((j, offset))
                    ######  Build+go offset
                    offset = self.kdbrdr._bgzfrdr.tell()# Find the offset of the next line
                    self.kdbrdr._bgzfrdr.seek(offset)   # Seek the next line
                    ######  
                    line = self.kdbrdr._bgzfrdr.readline().rstrip() # Read a new line
                    j += 1                              # Increment the k-mer/line index
                # If the block and/or the emitted line is empty, skip the block
                if len(line) == 0:
                    logger.debug("\n\nBlock:   {}\nK-mers:   {}".format(i, j))
            else:
                logger.warning("Block {0} starting at {1} has 0 data length".format(i, offset))
                logger.warning("kdb has read {0} k-mers until this point at block {1}".format(j, i))
            i += 1 # Increment the block counter
            # Verify that the number of indexed lines is equal to the number of theoretical k-mers
            if len(result) != self.kdb.num_kmers: 
                logger.error("kdb.index.IndexBuilder encountered a fatal bug: a mismatch of the number of indexed lines and the number of theoretical {}-mers".format(self.kdb.k))
                raise Exception("kdb.index.IndexBuilder.index_lines generated {0} lines offsets for {1} k-mers".format(len(result), self.kdb.num_kmers))

        return tuple(result)
            
    # def _index_tags(self):

    #     tags = self.kdb.header["tags"]
    #     result = {t: [] for t in tags}
    #     for idx, m in enumerate(self.kdb.metadata):
    #         if len(m["tags"]) > 0:
    #             for t in m["tags"]:
    #                 if t not in tags:
    #                     raise Exception("kdb.index.IndexBuilder._index_tags expects the tag '{}' to be in the header dictionary of the KDB file".format(t))
    #                 result[t].append(idx)
    #     return result
