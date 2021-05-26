# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._format import (
    BloomFilterDirectoryFormat, BitVectorFormat, BloomFilterFormat)
from ._type import (
    BitVector, BloomFilter)

# TODO remove these imports when tests are rewritten. Remove from __all__ too

__all__ = [
    'BloomFilterDirectoryFormat', 'BitVectorFormat', 'BloomFilterFormat',
    'BitVector', 'BloomFilter']
