import sys
import os
import copy
import json

import jsonschema
from Bio import bgzf
sys.path.append('..')

from kdb import fileutil

import logging
logger = logging.getLogger(__file__)


class IndexBuilder:
    def __init__(self, kdb: fileutil.KDB, kdbrdr: fileutil.KDBReader):
        if not isinstance(kdb, fileutil.KDB):
            raise TypeError("kdb.index.IndexBuilder expects a KDB object as its first positional argument")
        elif not isinstance(kdbrdr, fileutil.KDBReader):
            raise TypeError("kdb.index.IndexBuilder expects a KDBReader object as its second positional argument")
        if kdbrdr.loaded is False:
            raise TypeError("kdb.index.IndexBuilder expects a fully loaded KDBReader object as its first positional argument")
        self.kdb = kdb
        self.kdbrdr = kdbrdr
        self.line_index = self._index_lines()
        self.tag_index  = self._index_tags()

    def _index_lines(self):
        # Read by block
        offsets = copy.deepcopy(self.kdbrdr._offsets)
        # Omit the header block
        offsets.pop(0)
        result = []
        i = 0
        for raw_start, raw_length, data_start, data_length in offsets:
            if data_length > 0:
                offset = bgzf.make_virtual_offset(raw_start, 0) 
                self.kdbrdr._bgzfrdr.seek(offset) 
                line = self.kdbrdr._bgzfrdr.readline() 
                while line != '': #and json.loads(line.split("\t")[2]):
                    print(line)
                    result.append(offset)
                    offset = self.kdbrdr._bgzfrdr.tell()
                    self.kdbrdr._bgzfrdr.seek(offset)
                    line = self.kdbrdr._bgzfrdr.readline()
                if line == '':
                    logger.debug("Done with block {}...".format(i))
            else:
                logger.warning("Block {0} starting at {1} has 0 data length".format(i, raw_start))
            i += 1
        if len(result) != self.kdb.num_kmers:
            raise Exception("kdb.index.IndexBuilder.index_lines generated {0} line offsets for {1} k-mers")
        return tuple(result)
            
    def _index_tags(self):

        tags = self.kdb.header["tags"]
        result = {t: [] for t in tags}
        for idx, m in enumerate(self.kdb.metadata):
            if len(m["tags"]) > 0:
                for t in m["tags"]:
                    if t not in tags:
                        raise Exception("kdb.index.IndexBuilder._index_tags expects the tag '{}' to be in the header dictionary of the KDB file".format(t))
                    result[t].append(idx)
        return result
