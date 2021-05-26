# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
from qiime2.plugin import (Plugin, Str, Properties, Choices, Int, Bool, Range,
                           Float, Set, Visualization, Metadata, MetadataColumn,
                           Categorical, Numeric, Citations)

import q2_ananke
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import JoinedSequencesWithQuality
from q2_types.feature_data import FeatureData, Sequence
#from bitvector import BitVector, BloomFilter
citations = Citations.load('citations.bib', package='q2_ananke')

plugin = Plugin(
    name='ananke',
    version=q2_ananke.__version__,
    website='https://github.com/beiko-lab/q2-ananke/',
    package='q2_ananke',
    description=(''),
    short_description='A QIIME2 plugin for the Ananke application',
)

_HASH_OPTS = {'sha1', 'md5'}

plugin.methods.register_function(
   name='dereplicate',
   function=q2_ananke.dereplicate,
   inputs={'sequences': SampleData[JoinedSequencesWithQuality]},
   parameters={'hash_fn': qiime2.plugin.Str %
                       qiime2.plugin.Choices(_HASH_OPTS)},
   outputs=[('dereplicated_table', FeatureTable[Frequency]),
            ('dereplicated_sequences', FeatureData[Sequence])],
   input_descriptions={
       'sequences': 'Input (non-paired) sequence artifact'},
   parameter_descriptions={
       'hash_fn': 'Hash function for feature/sequence names'},
   output_descriptions={
       'dereplicated_table': 'Output dereplicated table',
       'dereplicated_sequences': 'Output dereplicated sequences'},
   description="Dereplicate joined sequences into a table and unique sequences artifact.")

#plugin.methods.register_function(
#    name='compute',
#    function=q2_ananke.compute,
#    inputs={'table': FeatureTable[Frequency]},
#    parameters={'metadata': Metadata,
#                'distance': Str % Choices(['sts', 'euclidean']),
#                'distance_min': Float,
#                'distance_max': Float,
#                'n_steps': Int,
#                'time_col': Str,
#                'series_col': Str,
#                'date_format': Str},
#    description="Compute the distances between features in a table and store the relationships in a Bloom Filter structure",
#    outputs=[('fingerprints', BitVector[BloomFilter])]
#    )
