# src/KmerDecon/bloom_filter.py

import math
import hashlib
from bitarray import bitarray
from typing import Any

class BloomFilter:
    """
    A Bloom filter for efficient k-mer storage and lookup.
    """

    def __init__(self, expected_elements: int, false_positive_rate: float):
        """
        Initialize the Bloom filter.

        Args:
            expected_elements (int): Estimated number of elements to store.
            false_positive_rate (float): Desired false positive rate.
        """
        self.size = self._get_size(expected_elements, false_positive_rate)
        self.hash_count = self._get_hash_count(self.size, expected_elements)
        self.bit_array = bitarray(self.size)
        self.bit_array.setall(0)

    def _get_size(self, n: int, p: float) -> int:
        """
        Calculate the size of the bit array (m) given expected elements (n) and false positive rate (p).

        Returns:
            int: Size of the bit array.
        """
        m = -(n * math.log(p)) / (math.log(2) ** 2)
        return int(m)

    def _get_hash_count(self, m: int, n: int) -> int:
        """
        Calculate the optimal number of hash functions (k) given bit array size (m) and expected elements (n).

        Returns:
            int: Number of hash functions.
        """
        k = (m / n) * math.log(2)
        return int(k)

    def add(self, item: Any) -> None:
        """
        Add an item to the Bloom filter.

        Args:
            item (Any): The item to add.
        """
        for i in range(self.hash_count):
            digest = hashlib.sha256(f"{item}{i}".encode('utf-8')).hexdigest()
            index = int(digest, 16) % self.size
            self.bit_array[index] = True

    def __contains__(self, item: Any) -> bool:
        """
        Check if an item is in the Bloom filter.

        Args:
            item (Any): The item to check.

        Returns:
            bool: True if the item is probably in the filter, False if definitely not.
        """
        for i in range(self.hash_count):
            digest = hashlib.sha256(f"{item}{i}".encode('utf-8')).hexdigest()
            index = int(digest, 16) % self.size
            if not self.bit_array[index]:
                return False
        return True

    def save(self, filename: str) -> None:
        """
        Save the Bloom filter to a file.

        Args:
            filename (str): The filename to save the filter to.
        """
        with open(filename, 'wb') as f:
            self.bit_array.tofile(f)
        with open(f"{filename}.params", 'w') as f:
            f.write(f"{self.size}\n")
            f.write(f"{self.hash_count}\n")

    @classmethod
    def load(cls, filename: str) -> 'BloomFilter':
        """
        Load a Bloom filter from a file.

        Args:
            filename (str): The filename to load the filter from.

        Returns:
            BloomFilter: The loaded Bloom filter.
        """
        with open(f"{filename}.params", 'r') as f:
            size = int(f.readline())
            hash_count = int(f.readline())
        bloom_filter = cls.__new__(cls)
        bloom_filter.size = size
        bloom_filter.hash_count = hash_count
        bloom_filter.bit_array = bitarray(size)
        with open(filename, 'rb') as f:
            bloom_filter.bit_array.fromfile(f)
        return bloom_filter