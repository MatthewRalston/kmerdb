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
import yaml, json
from collections import deque, OrderedDict



def represent_ordereddict(dumper, data):
    """
    Thanks to Blender and EquipDev on StackOverflow for this handy method to pass to yaml.add_representer
    https://stackoverflow.com/a/16782282
    """
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)


def merge_metadata_lists(k, metadata_across_all_files, new_kmer_metadata):
    """
    Merge two 4**k metadata lists
    """
    if 4**k != len(new_kmer_metadata):
        raise TypeError("kmerdb.util.merge_metadata_lists() expects a new_kmer_metadata list to have 4^{0} elements. Got {1}".format(k, len(new_kmer_metadata)))

    for i, metadata in enumerate(new_kmer_metadata):
        metadata_across_all_files[i] = metadata_across_all_files[i] + metadata

    return metadata_across_all_files
