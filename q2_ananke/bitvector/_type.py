# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType

from ..plugin_setup import plugin
from . import BloomFilterDirectoryFormat

BitVector = SemanticType('BitVector', field_names='type')

BloomFilter = SemanticType('BloomFilter', variant_of=BitVector.field['type'])

plugin.register_semantic_types(BitVector, BloomFilter)

plugin.register_semantic_type_to_format(
    BitVector[BloomFilter],
    artifact_format=BloomFilterDirectoryFormat)

