# The following code is adapted from:
# Repository: https://github.com/michael-the1/python-cuckoo/tree/master
# License: MIT

import mmh3
import random
import KmerDecon.bucket as bucket
import gzip
import pickle
from bitarray import bitarray

class CuckooFilter:
    '''
    A Cuckoo filter is a data structure for probablistic set-membership queries.
    We can insert items into the filter and, with a very low false positive probability (FPP), ask whether it contains an item or not.
    Cuckoo filters serve as a drop-in replacement for Bloom filters, but are more space-efficient and support deletion of items.

    Cuckoo filters were originally described in:
        Fan, B., Andersen, D. G., Kaminsky, M., & Mitzenmacher, M. D. (2014, December).
        Cuckoo filter: Practically better than bloom.
        In Proceedings of the 10th ACM International on Conference on emerging Networking Experiments and Technologies (pp. 75-88). ACM.

    Their reference implementation in C++ can be found here: https://github.com/efficient/cuckoofilter
    '''

    def __init__(self, capacity, fingerprint_size, kmer_length,bucket_size=4, max_kicks=500):
        '''
        Initialize Cuckoo filter parameters.

        capacity : size of the filter
            Defines how many buckets the filter contains.
        fingerprint_size: size of the fingerprint in bytes
            A larger fingerprint size results in a lower FPP.
        bucket_size : nr. of entries a bucket can hold
            A bucket can hold multiple entries.
            Default size is 4, which closely approaches the best size for FPP between 0.00001 and 0.002 (see Fan et al.).
            If your targeted FPP is greater than 0.002, a bucket size of 2 is more space efficient.
        max_kicks : nr. of times entries are kicked around before deciding the filter is full
            Defaults to 500. This is an arbitrary number also used by Fan et al. and seems reasonable enough.
        '''
        self.capacity = capacity
        self.fingerprint_size = fingerprint_size
        self.max_kicks = max_kicks
        self.buckets = [bucket.Bucket(size=bucket_size) for _ in range(self.capacity)]
        self.size = 0
        self.bucket_size = bucket_size
        self.kmer_length=kmer_length

    def insert(self, item):
        '''
        Inserts a string into the filter.

        Throws an exception if the insertion fails.
        '''
        self.size = self.size + 1
        fingerprint = self.fingerprint(item)
        if self.contains(item):
            return
        i1, i2 = self.calculate_index_pair(item, fingerprint)

        if self.buckets[i1].insert(fingerprint):
            return i1
        elif self.buckets[i2].insert(fingerprint):
            return i2

        i = random.choice((i1, i2))
        for kick_count in range(self.max_kicks):
            fingerprint = self.buckets[i].swap(fingerprint)
            # i = (i ^ self.index_hash(fingerprint)) % self.capacity
            i = (i ^ self.index_hash(fingerprint.tobytes())) % self.capacity

            if self.buckets[i].insert(fingerprint):
                return i

        self.size = self.size - 1
        raise Exception('Filter is full')

    def contains(self, item):
        '''Checks if a string was inserted into the filter.'''
        fingerprint = self.fingerprint(item)
        i1, i2 = self.calculate_index_pair(item, fingerprint)
        return (fingerprint in self.buckets[i1]) or (fingerprint in self.buckets[i2])

    def delete(self, item):
        '''Removes a string from the filter.'''
        fingerprint = self.fingerprint(item)
        i1, i2 = self.calculate_index_pair(item, fingerprint)
        if self.buckets[i1].delete(fingerprint) or self.buckets[i2].delete(fingerprint):
            self.size = self.size - 1
            return True
        return False

    def index_hash(self, item):
        '''Calculate the (first) index of an item in the filter.'''
        item_hash = mmh3.hash_bytes(item)
        index = int.from_bytes(item_hash, byteorder='big') % self.capacity
        return index

    def calculate_index_pair(self, item, fingerprint):
        '''Calculate both possible indices for the item'''
        i1 = self.index_hash(item)
        # i2 = (i1 ^ self.index_hash(fingerprint)) % self.capacity
        i2 = (i1 ^ self.index_hash(fingerprint.tobytes())) % self.capacity
        return i1, i2

    def fingerprint(self, item):
        '''
        Takes a string and returns its fingerprint in bits.

        The length of the fingerprint is given by fingerprint_size.
        To calculate this fingerprint, we hash the string with MurmurHash3 and truncate the hash.
        '''
        # item_hash = mmh3.hash_bytes(item)
        # return item_hash[:self.fingerprint_size]

        item_hash = bitarray()
        item_hash.frombytes(mmh3.hash_bytes(item))
        return item_hash[:self.fingerprint_size]

    def load_factor(self):
        return self.size / (self.capacity * self.bucket_size)

    def __contains__(self, item):
        return self.contains(item)

    def __repr__(self):
        return '<CuckooFilter: capacity=' + str(self.capacity) + ', fingerprint size=' + str(self.fingerprint_size) + ' byte(s)>'

    def __sizeof__(self):
        return super().__sizeof__() + sum(b.__sizeof__() for b in self.buckets)
    
    def save(self, filename: str) -> None:
        """
        Save the Cuckoo filter to a compressed file using gzip.

        Args:
            filename (str): The filename to save the filter to.
        """
        with gzip.open(filename, 'wb') as f:
             pickle.dump(self, f)
    @classmethod
    def load(cls, filename: str) -> 'CuckooFilter':
        with gzip.open(filename, 'rb') as f:
            return pickle.load(f)
        print(f"Cuckoo Filter load from {filename}")
        return cuckoo_filter
# cf = CuckooFilter(capacity=100000, fingerprint_size=1)
# cf.insert("a")
# a=cf.contains("a")
# print(a)
