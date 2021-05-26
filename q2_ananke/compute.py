import plotly
import plotly.graph_objects as go
import pybloom
from pybloom.utils import range_fn
import numpy as np
from bitarray import bitarray
from struct import pack, unpack
import warnings
import scipy
from qiime2 import Artifact
import biom
import sklearn
import random
import xxhash
from sklearn.metrics import pairwise_distances_chunked
from zipfile import ZipFile
import pandas as pd
from io import BytesIO
import os
sklearn.set_config(working_memory=128)

class DBloomSCAN(object):
    def __init__(self, asv_table, dist_measure='euclidean',
                 dist_range=[0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01], 
                 error_rate = 0.01, max_dist=None, min_abundance=1, start_at=0, bf_edges_per_node=10000):
        self.asv_table = Artifact.load(asv_table).view(biom.Table)
        filtered_asv_table = self.asv_table.filter(self.asv_table.ids(axis="observation")[self.asv_table.sum("observation") > min_abundance],
                                              "observation", inplace=False)
        obs_sorted_abundance = filtered_asv_table.ids('observation')[np.argsort(filtered_asv_table.sum("observation"))[::-1]]
        self.sorted_asv_table = filtered_asv_table.sort_order(obs_sorted_abundance, axis='observation')
        self.matrix = self.sorted_asv_table.matrix_data
        self.min_abundance = min_abundance
        self.bf_edges_per_node=bf_edges_per_node
        self.dist_measure = dist_measure
        self.n_objects = self.matrix.shape[0]
        print("Initializing DBloomSCAN object with %d objects" % (self.n_objects,))
        self.error_rate = error_rate
        self.dist_range = np.array([round(x, 5) for x in dist_range])
        print("Creating Bloom Filters to track distances across %d steps: %s" % (len(self.dist_range), str(self.dist_range),))
        self.bloom_garden = BloomGarden(self.dist_range, bf_edges_per_node*self.n_objects, self.error_rate)
        if max_dist is None:
            self.max_dist = 2*self._sample_distances()
        else:
            self.max_dist = max_dist
        self.started_at = start_at
        self.n_computed = start_at

    def add_distance(self, data):
        #Input distance is expected to be normalized by max_dist already
        i, j, distance = data
        #If the distance is bigger than our largest range, it isn't important, discard ASAP
        if distance > self.dist_range[-1]:
            return
        def bloom_condition(bloom):
            return bloom > distance
        #Given a real distance d, if d <= name, that means that the objects are closer than the threshold,
        #meaning they are neighbours at that distance or less
        i = self.sorted_asv_table.ids('observation')[i]
        j = self.sorted_asv_table.ids('observation')[j]
        pruned_blooms = self.bloom_garden.add((min(i,j),max(i,j)), bloom_condition)
        # If we had to prune a full filter, we kick it out of our dist range
        if pruned_blooms == 0:
            pruned_blooms = None
        else:
            pruned_blooms = -pruned_blooms
        self.dist_range = self.dist_range[:pruned_blooms]
        
    def _sample_distances(self, scale_factor=10):
        print("Sampling distances")
        # Sample some distances to guess the max epsilon, used for scaling
        # In the RAM-hungry version, we know the max up-front and can scale
        # perfectly, but in this case we have to guess and then check
        nrows = self.n_objects + 1
        max_dist = 0.0
        # Sample twice as many distances as we have unique genes
        # TODO: Validate that this is a good enough amount of sampling
        n_iters = scale_factor*int(np.sqrt(nrows))
        print("Sampling %d distances" % (n_iters**2,))
        total_distance = 0
        random_range = random.sample(range(0, nrows-1), n_iters)
        sub_matrix = self.matrix[random_range]
        for chunk in pairwise_distances_chunked(sub_matrix, metric=self.dist_measure, n_jobs=-1):
            distance = chunk.max()
            if distance > max_dist:
                max_dist = distance
        print("Max distance found is %f" % (max_dist,))
        return max_dist

    def are_neighbours(self, i, j, distance, validate = True):
        if distance not in self.bloom_garden.blooms:
            raise ValueError("Distance not found in bloom filters")
        if type(i) == int:
            int_i = i
            i = self.sorted_asv_table.ids('observation')[i]
        else:
            int_i = np.argwhere(self.sorted_asv_table.ids('observation')==i)[0][0]
        if type(j) == int:
            int_j = j
            j = self.sorted_asv_table.ids('observation')[j]
        else:
            int_j = np.argwhere(self.sorted_asv_table.ids('observation')==j)[0][0]
        bf_result = (min(i,j), max(i,j)) in self.bloom_garden.blooms[distance]
        #If the bloom filter says no, they are definitely not neighbours within this distance
        if not bf_result:
            return False
        else:
            if validate:
                #Compute the actual distance
                verified_distance = scipy.spatial.distance.pdist(self.matrix[[int_i,int_j],:].todense(),
                                                                 self.dist_measure)[0]
                verified_distance = verified_distance / self.max_dist
                if verified_distance > distance:
                    return False
                else:
                    return True
            else:
                # Return it as a positive, and put it on the caller to double check
                return True
            
    def _process_distance_chunk(self, D_chunk, start_index):
        D_chunk = D_chunk / self.max_dist
        for index,v in np.ndenumerate(D_chunk):
            i=index[0]+start_index
            j=index[1]
            #Only add the lower triangle since symmetry
            if i<j:
                self.add_distance((i,j,v))
        if i>self.n_computed:
            self.n_computed = i+1 #n is index plus one
        print("%d / %d" % (min(self.n_objects,start_index+D_chunk.shape[0]),self.n_objects))
        return D_chunk
            
    def compute_distances(self):
        print("Beginning memory-chunked computation of distances of full matrix using measure %s" % (self.dist_measure,))
        dist_gen = pairwise_distances_chunked(self.matrix, metric=self.dist_measure, 
                                              n_jobs=-1, reduce_func=self._process_distance_chunk)
        print("0 / %d" % (self.n_objects,))
        for D_chunk in dist_gen:
            pass
        
    def DBSCAN(self, epsilon, min_pts = 2, expand_around=None, max_members=None, warn=True):
        if epsilon not in list(self.dist_range):
            dist_range = self.dist_range
            delta = dist_range - epsilon
            old_epsilon = epsilon
            epsilon = self.dist_range[np.argmin(abs(delta))]
            if warn:
                print("Bloom filter does not exist for this epsilon value, %f. " \
                  "Using the closest precomputed value, %f." % (old_epsilon,epsilon))
        cluster_number = 0
        clustered = set()
        cluster_assignments = {}
        if expand_around is not None:
            if type(expand_around) != int:
                expand_around = int(np.argwhere(self.sorted_asv_table.ids('observation')==expand_around)[0][0])
            index_queue = [ expand_around ]
        else:
            index_queue = range(0, self.n_objects)
            if max_members is not None:
                warnings.warn("max_members ignored, only used with expand_around")
        for i in index_queue:
            if i in clustered:
                continue
            cluster_number += 1
            cluster_assignments[cluster_number] = [i]
            cluster_queue = [i]
            clustered.add(i)
            while cluster_queue:
                k = cluster_queue.pop()
                neighbourhood = []
                for j in range(0, self.n_objects):
                    if (j != k) & (j not in clustered) & \
                       (self.are_neighbours(k, j, epsilon, validate=True)):
                        neighbourhood.append(j)
                        if (expand_around is not None) & (max_members is not None):
                            if len(neighbourhood) + len(clustered) > max_members:
                                raise RuntimeError("Cluster too large, aborting")

                # min_pts neighbourhood size includes the point itself, so we account for that here
                # This means k is a core point
                if len(neighbourhood) >= min_pts - 1:
                    cluster_queue.extend(neighbourhood)
                    #if it is in range of a core point, it's in the cluster
                    cluster_assignments[cluster_number].extend(neighbourhood)
                    clustered.update(neighbourhood)
        return cluster_assignments

    def save_results(self, output_filename):
        print("Writing Ananke results to Zip file...")
        with ZipFile(output_filename,'w') as zf:
            for bloom in self.bloom_garden.blooms:
                bloom_bin = 'bloom_%s.bin'%(str(bloom),)
                with open(bloom_bin,'wb') as ob:
                    self.bloom_garden.blooms[bloom].bitarray.tofile(ob)
                zf.write(bloom_bin)
                os.remove(bloom_bin)
            with open("parameters.txt",'w') as pf:
                for bloom in self.bloom_garden.blooms:
                    pf.write("capacity\t%f\t%d\n" %(bloom,self.bloom_garden.blooms[bloom].capacity))
                for bloom in self.bloom_garden.blooms:
                    pf.write("count\t%f\t%d\n" %(bloom,self.bloom_garden.blooms[bloom].count))
                pf.write("error_rate\t%f\n"%(self.error_rate,))
                pf.write("min_abundance\t%d\n"%(self.min_abundance,))
                pf.write("n_computed\t%d\n"%(self.n_computed,))
                pf.write("bf_edges_per_node\t%d\n"%(self.bf_edges_per_node,))
                pf.write("max_dist\t%f\n"%(self.max_dist,))
            zf.write("parameters.txt")
            os.remove('parameters.txt')
        print("Write complete!")
       
    @classmethod
    def load_results(cls, results_filename, table_artifact_path):
        #Open up zip file
        zf = ZipFile(results_filename)
        #Grab settings
        parameter_table = pd.read_table(BytesIO(zf.read("parameters.txt")), header=None)
        min_abundance = int(parameter_table[parameter_table[0]=="min_abundance"][1])
        error_rate = float(parameter_table[parameter_table[0]=="error_rate"][1])
        n_computed = int(parameter_table[parameter_table[0]=="n_computed"][1])
        bf_edges_per_node = int(parameter_table[parameter_table[0]=="bf_edges_per_node"][1])
        max_dist = float(parameter_table[parameter_table[0]=="max_dist"][1])
        dist_range = list(parameter_table[parameter_table[0]=="count"][1])
        dist_range = [round(x, 5) for x in dist_range]
        #Initialize DBloomSCAN empty object
        dbs = cls(table_artifact_path, bf_edges_per_node=bf_edges_per_node,
                  min_abundance=min_abundance, dist_range=dist_range,
                  error_rate=error_rate, max_dist=max_dist)
        #Replace bitvectors
        for dist in dist_range:
            bloom_array = bitarray(endian="little")
            with BytesIO(zf.read("bloom_%g.bin"%(dist,))) as bff:
                bloom_array.fromfile(bff)
            # Just a little trimming because the saved file sometimes is larger than the in_memory bitarray
            for i in range(0,len(bloom_array)-len(dbs.bloom_garden.blooms[dist].bitarray)):
                bloom_array.pop()
            dbs.bloom_garden.blooms[dist].bitarray = bloom_array 
            counts = parameter_table[parameter_table[0]=="count"]
            count = int(counts[counts[1]==dist][2])
            dbs.bloom_garden.blooms[dist].count = count
        return dbs
        
    def diagnostic_plots(self):
        fig = go.Figure(data=go.Scatter(x=[eps for eps in self.bloom_garden.blooms.keys()],
                                y=[bloom.count for bloom in self.bloom_garden.blooms.values()]))
        fig.update_layout(title="Distance Sampling Plot", xaxis_title="Epsilon", yaxis_title="Number of Neighbours")
        fig.show()
        #Reimport full table
        full_asv_table = self.asv_table
        taxon_abundances = full_asv_table.sum(axis="observation")
        fig = go.Figure()
        fig.add_trace(go.Histogram(x=taxon_abundances[taxon_abundances>=self.min_abundance],
                name="Included in clustering",
                xbins=dict( # bins used for histogram
                start=self.min_abundance,
                end=max(taxon_abundances),
                size=10)))
        fig.add_trace(go.Histogram(x=taxon_abundances[taxon_abundances<self.min_abundance], 
                name="Excluded in clustering",
                xbins=dict( # bins used for histogram
                start=0,
                end=self.min_abundance,
                size=10)))
        fig.update_layout(barmode='overlay')
        fig.update_layout(yaxis_type = "log")
        fig.show()
        unique_seqs_in = taxon_abundances[taxon_abundances>=self.min_abundance]
        unique_seqs_out = taxon_abundances[taxon_abundances<self.min_abundance]
        abundance_in = sum(unique_seqs_in)
        abundance_out = sum(unique_seqs_out)
        unique_seqs_in = len(unique_seqs_in)
        unique_seqs_out = len(unique_seqs_out)

        fig = plotly.subplots.make_subplots(rows=1, cols=2, specs=[[{"type": "pie"}, {"type": "pie"}]])
        fig.update_layout(title={"text":"Unique (left) and Total (right) Sequences Clustered by Ananke","x": 0.5})
        fig.add_trace(go.Pie(labels=["Included","Excluded"], name="Unique", values=[unique_seqs_in,unique_seqs_out]), row=1, col=1)
        fig.add_trace(go.Pie(labels=["Included", "Excluded"], name="Abundance", values=[abundance_in,abundance_out]), row=1, col=2)
        fig.show()
        
class ExternalHashBloom(pybloom.BloomFilter):
    def __init__(self, capacity, error_rate=0.001):
        super().__init__(capacity, error_rate)
        self.make_hashes = make_hashfuncs(self.num_slices, self.bits_per_slice)

    #Overwrite the existing add function, but remove the hash check
    def add_hashes(self, hashes, skip_check = False):
        bitarray = self.bitarray
        bits_per_slice = self.bits_per_slice
        found_all_bits = True
        if self.count > self.capacity:
            raise IndexError("BloomFilter is at capacity")
        offset = 0
        for k in hashes:
            if not skip_check and found_all_bits and not bitarray[offset + k]:
                found_all_bits = False
            self.bitarray[offset + k] = True
            offset += bits_per_slice

        if skip_check:
            self.count += 1
            return False
        elif not found_all_bits:
            self.count += 1
            return False
        else:
            return True
            
class BloomGarden(object):
    def __init__(self, filter_names, capacity, error_rate):
        self.blooms = {}
        for name in filter_names:
            self.blooms[name] = ExternalHashBloom(capacity, error_rate)

    def add(self, key, name_condition):
        prune_list = []
        pruned_blooms = 0
        hashes = None
        for name in self.blooms:
            if not hashes:
                hashes = self.blooms[name].make_hashes(key)
                #This is a generator, so we need to coerce it
                #to something static or the first insert depletes it
                hashes = list(hashes)
            if name_condition(name):
                try:
                    #We can skip the check because each pair we check is unique\
                    self.blooms[name].add(key, skip_check=True)
                except IndexError:
                    print(self.blooms[name].count, self.blooms[name].capacity)
                    print("Bloom filter '%s' hit capacity, closing" % (str(name),))
                    prune_list.append(name)
                    pruned_blooms += 1
        for bloom in prune_list:
            del self.blooms[bloom]
        if not self.blooms:
            raise IndexError("All bloom filters closed. Try using a smaller minimum epsilon value.")
        return pruned_blooms

# This is taken from pybloom, but modified to use xxhash.xxh64()
def make_hashfuncs(num_slices, num_bits):
    if num_bits >= (1 << 31):
        fmt_code, chunk_size = 'Q', 8
    elif num_bits >= (1 << 15):
        fmt_code, chunk_size = 'I', 4
    else:
        fmt_code, chunk_size = 'H', 2
    total_hash_bits = 8 * num_slices * chunk_size
    hashfn = xxhash.xxh64
    fmt = fmt_code * (hashfn().digest_size // chunk_size)
    num_salts, extra = divmod(num_slices, len(fmt))
    if extra:
        num_salts += 1
    salts = tuple(hashfn(hashfn(pack('I', i)).digest()) for i in range_fn(num_salts))
    def _make_hashfuncs(key):
        if isinstance(key, str):
            key = key.encode('utf-8')
        else:
            key = str(key).encode('utf-8')
        i = 0
        for salt in salts:
            h = salt.copy()
            h.update(key)
            for uint in unpack(fmt, h.digest()):
                yield uint % num_bits
                i += 1
                if i >= num_slices:
                    return
    return _make_hashfuncs

