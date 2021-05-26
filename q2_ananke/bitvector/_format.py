# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import skbio.io

import qiime2.plugin.model as model
from qiime2.plugin import ValidationError
import qiime2

from ..plugin_setup import plugin


class BitVectorFormat(model.BinaryFileFormat):
    metadata_bits = 0

class BloomFilterFormat(BitVectorFormat):
    metadata_bits = 64

class BloomFilterDirectoryFormat(model.DirectoryFormat):
    bloomfilters = model.FileCollection(r'\d+\.bin',
                                        format=BloomFilterFormat)

plugin.register_formats(
        BitVectorFormat, BloomFilterFormat, BloomFilterDirectoryFormat
)

