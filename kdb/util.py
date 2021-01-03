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



