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
import sys
import os
import gzip
import tempfile
import hashlib
import yaml, json
import re
from collections import deque, OrderedDict
import numpy as np

from kmerdb import config

is_integer = re.compile("^[-+]?[0-9]+$")
findall_float = re.compile(r"[-+]?(?:\d*\.\d+|\d+)")

def checksum(filepath:str):
    """Generates md5 and sha256 checksums of a file
    :returns: (md5, sha256)
    :rtype: tuple
    """
    if type(filepath) is not str:
        raise TypeError("kmerdb.util.checksum() expects a str as its argument")
    elif not os.path.exists(filepath):
        raise IOError("kmerdb.parse.SeqParser.__checksum could not find '{}' on the filesystem".format(filepath))
    hash_md5 = hashlib.md5()
    hash_sha256 = hashlib.sha256()
    with open(filepath, 'rb') as ifile:
        for chunk in iter(lambda: ifile.read(4096), b""):
            hash_md5.update(chunk)
            hash_sha256.update(chunk)
    return (hash_md5.hexdigest(), hash_sha256.hexdigest())


def represent_yaml_from_collections_dot_OrderedDict(dumper, data):
    """
    Thanks to Blender and EquipDev on StackOverflow for this handy method to pass to yaml.add_representer
    https://stackoverflow.com/a/16782282

    I use this throughout the metadata representation in the bgzf specification to print out the representer, it's just awesome.

    I really like this answer. Finally got around to liking it this year, in Jan 2022 after using it for like a few years. 

    :param dumper: The OrderedDict_Representer, this faciliatates non-key sorting for optimal metadata block structure.
    :type dumper: 
    :param data:
    :type data: dict
    """

    
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)


def is_gz_file(filepath):
    import gzip
    try:
        
        with gzip.open(filepath, 'rb') as f:
            f.readline()
        return True
    except OSError as e:
        return False



def get_histo(counts):
    """
    Get a histogram from an array.
    """

    if type(counts) is not list: # 10/6/24
        raise TypeError("kmerdb.util.get_histo takes a list as its positional argument")
    
    hist = []
    count_max = int(np.max(np.array(counts, dtype="uint64")))

    for i in range(count_max + 1):
        hist.append(0)

    #print("Count max was : {0} from {1} histogram items".format(count_max, len(hist)))    
    #sys.stderr.write("{0} histogram items\n".format(len(hist)))
    for i in range(len(counts)):
        # 10/7/24
        # print("i: {0}".format(i))
        # print("count array length: {0}".format(len(counts)))
        # print("hist array length: {0}".format(len(hist)))
        # print("histogram item: {0}, {1} in {2} histogram items".format(i, counts[i], len(hist)))
        if counts[i] > 2:
            hist[counts[i]] += 1
    return hist



def is_fasta(fname:str):
    if type(fname) is not str:
        raise TypeError("kmerdb.util.is_fasta() expects a str argument")
    if fname.endswith(".fna") or fname.endswith(".fna.gz") or fname.endswith(".fa.gz") or fname.endswith(".fa") or fname.endswith(".fasta") or fname.endswith(".fasta.gz"):
        return True
    else:
        return False

def is_fastq(fname:str):
    if type(fname) is not str:
        raise TypeError("kmerdb.util.is_fasta() expects a str argument")
    elif fname.endswith(".fastq") or fname.endswith(".fastq.gz") or fname.endswith(".fq.gz") or fname.endswith(".fq"):
        return True
    else:
        return False

    
def is_all_fasta(filenames):
    """Tests if all the strings in a list are fasta format.
    
    :param filenames:
    :type list:
    :returns bool:
    :rtype bool:
    """
    if type(filenames) is not list:
        raise TypeError("kmerdb.util.is_all_fasta() expects a list as its first positional argument")
    elif not all(type(s) is str for s in filenames):
        raise TypeError("kmerdb.util.is_all_fasta() expects a list of str as its first positional argument")
    return all((f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fa.gz') or f.endswith('.fasta.gz') or f.endswith('.fna') or f.endswith('.fna.gz')) for f in filenames)


def isfloat(num):
    """
    Thanks to the author of:
    https://www.programiz.com/python-programming/examples/check-string-number


    :param num: Check if a string is a float, or could be converted to a float
    :type num: str
    :returns: Whether or not the string can be parsed as a float
    :rtype: bool
    """
    if type(num) is str:
        #logger.debug("Type of number being interpreted through pure Python : {0}".format(type(num)))
        findall_float.match(num)
        #logger.debug(re.match(findall_float, num))
        return re.match(findall_float, num) is not None
    elif type(num) is float:
        return True
    elif type(num) is int:
        return False
    else:
        #logger.error(type(num))
        #logger.error(num.dtype)
        raise TypeError("kmerdb.util.isfloat expects a single str as a positional argument")

