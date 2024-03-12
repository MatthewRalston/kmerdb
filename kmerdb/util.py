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


from kmerdb import config



def validate_graph_metadata(metadata):
    """
    This function can be abstracted for only using jsonschema non-builtin...

    needs jsonschema

    """
    jsonschema.validate(instance=dict(metadata), schema=config.graph_schema)






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


