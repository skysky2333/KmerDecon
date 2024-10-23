import math
import mmh3
from bitarray import bitarray
from typing import Any

class BloomFilter:
    """
    A Bloom filter for efficient k-mer storage and lookup.
    """

    def __init__(self, expected_elements: int, false_positive_rate: float, kmer_length: int):
        """
        Initialize the Bloom filter.

        Args:
            expected_elements (int): Estimated number of elements to store.
            false_positive_rate (float): Desired false positive rate.
            kmer_length (int): Length of k-mers used.
        """
        self.expected_elements = expected_elements
        self.false_positive_rate = false_positive_rate
        self.kmer_length = kmer_length
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
            digest = mmh3.hash(item, i) % self.size
            self.bit_array[digest] = True

    def __contains__(self, item: Any) -> bool:
        """
        Check if an item is in the Bloom filter.

        Args:
            item (Any): The item to check.

        Returns:
            bool: True if the item is probably in the filter, False if definitely not.
        """
        for i in range(self.hash_count):
            digest = mmh3.hash(item, i) % self.size
            if not self.bit_array[digest]:
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
            f.write(f"{self.expected_elements}\n")
            f.write(f"{self.false_positive_rate}\n")
            f.write(f"{self.kmer_length}\n")
            f.write(f"{len(self.bit_array)}\n")  # Save the actual bit length

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
            expected_elements = int(f.readline())
            false_positive_rate = float(f.readline())
            kmer_length = int(f.readline())
            bitarray_length = int(f.readline())

        bloom_filter = cls.__new__(cls)
        bloom_filter.size = size
        bloom_filter.hash_count = hash_count
        bloom_filter.expected_elements = expected_elements
        bloom_filter.false_positive_rate = false_positive_rate
        bloom_filter.kmer_length = kmer_length
        bloom_filter.bit_array = bitarray()

        with open(filename, 'rb') as f:
            bloom_filter.bit_array.fromfile(f)

        bloom_filter.bit_array = bloom_filter.bit_array[:bitarray_length]

        return bloom_filter