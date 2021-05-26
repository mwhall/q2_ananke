import os
from collections import Counter, defaultdict
import re
from hashlib import sha1, md5
import gzip
import zlib
from scipy.sparse import lil_matrix
from biom.table import Table
from skbio import Sequence
from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt

HASH_FN_MAP = {"sha1": sha1, "md5": md5}

def dereplicate(sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                hash_fn: str = "sha1") -> (Table, DNAIterator):
    hash_fn = HASH_FN_MAP[hash_fn]
    base_dir = str(sequences)
    sequences = {}
    sequence_counts = defaultdict(Counter)
    count = 0
    sample_ids = set()
    for f in os.listdir(base_dir):
        if f.endswith("_L001_R1_001.fastq.gz"):
            sample_id = re.sub(r'_L001_R[12]_001.fastq.gz', '', f)
            sample_ids.add(sample_id)
            fh = base_dir+"/"+f
            with gzip.open(fh, 'r') as fp:
                line = True
                while line:
                    line = fp.readline() # label
                    if not line:
                        continue
                    line = fp.readline() # sequence
                    seq = line.strip()
                    line = fp.readline() # separator
                    line = fp.readline() # quality
                    seq_hash = hash_fn(seq).hexdigest()
                    if seq_hash not in sequences:
                        sequences[seq_hash] = zlib.compress(seq)
                    sequence_counts[seq_hash][sample_id] += 1
                    count += 1
    seq_hashes = list(sequences.keys())
    seqs = list(sequences.values())
    data = lil_matrix((len(seq_hashes), len(sample_ids)))
    for cidx, sample_id in enumerate(sample_ids):
        for ridx, seq_hash in enumerate(seq_hashes):
            if sample_id in sequence_counts[seq_hash]:
                data[ridx, cidx] = sequence_counts[seq_hash][sample_id]
    del sequence_counts
    biom_table = Table(data, seq_hashes, list(sample_ids))
    del data
    del seq_hashes
    def yield_seqs(sequences):
        for seq_hash in sequences:
            sequence = zlib.decompress(sequences[seq_hash])
            if len(sequence) > 0:
                yield Sequence(sequence, {'id': seq_hash})
    dnait = DNAIterator(yield_seqs(sequences))
    return biom_table, dnait
