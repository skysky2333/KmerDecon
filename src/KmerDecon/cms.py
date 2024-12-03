import math
import mmh3
import gzip
import numpy as np
from typing import Any, Tuple
import xxhash
class CountMinSketch:
    """
    A Count-Min Sketch for efficient k-mer frequency counting.
    """

    def __init__(self, width: int, depth: int, kmer_length: int):
        """
        Initialize the Count-Min Sketch.

        Args:
            width (int): Number of columns in the sketch (controls accuracy).
            depth (int): Number of rows in the sketch (controls probability of error).
            kmer_length (int): Length of k-mers used.
        """
        self.width = width
        self.depth = depth
        self.kmer_length = kmer_length
        self.table = np.zeros((depth, width), dtype=int)
        self.size=self.width * self.depth

    @staticmethod
    def _optimal_params(expected_elements: int, error_rate: float=0.01) -> Tuple[int, int]:
        """
        Calculate optimal parameters for the Count-Min Sketch.

        Args:
            expected_elements (int): Estimated number of elements to store.
            error_rate (float): Desired error rate.

        Returns:
            Tuple[int, int]: Width and depth of the sketch.
        """
        width = int(math.ceil(math.e / error_rate))*10000
        depth = int(math.ceil(math.log(1 / error_rate)))
        return width, depth

    def add(self, item: Any) -> None:
        """
        Add an item to the Count-Min Sketch.

        Args:
            item (Any): The item to add.
        """
        for i in range(self.depth):
            index = xxhash.xxh32(item.encode('utf-8')).intdigest() % self.width
            self.table[i, index] += 1

    def count(self, item: Any) -> int:
        """
        Get the count of an item in the Count-Min Sketch.
        Args:
            item (Any): The item to query.
        Returns:
            int: Estimated count of the item.
        """
        return min(self.table[i, mmh3.hash(item, i) % self.width] for i in range(self.depth))

    def save(self, filename: str) -> None:
        """
        Save the Count-Min Sketch to a compressed file using gzip.

        Args:
            filename (str): The filename to save the sketch to.
        """
        # Save the table data
        with gzip.open(filename, 'wb') as f:
            np.save(f, self.table)
        
        # Save the parameters (width, depth, kmer_length) to a separate .params file
        with open(f"{filename}.params", 'w') as f:
            f.write(f"{self.width}\n")
            f.write(f"{self.depth}\n")
            f.write(f"{self.kmer_length}\n")

    @classmethod
    def load(cls, filename: str) -> 'CountMinSketch':
        """
        Load a Count-Min Sketch from a compressed file using gzip.

        Args:
            filename (str): The filename to load the sketch from.

        Returns:
            CountMinSketch: The loaded Count-Min Sketch.
        """
        # Load the parameters from the .params file
        with open(f"{filename}.params", 'r') as f:
            width = int(f.readline())
            depth = int(f.readline())
            kmer_length = int(f.readline())

        # Create a new instance and set its parameters
        cms = cls(width, depth, kmer_length)

        # Load the table data from the gzip file
        with gzip.open(filename, 'rb') as f:
            cms.table = np.load(f,allow_pickle=True)

        return cms
# w,d=CountMinSketch._optimal_params(4)
# cms=CountMinSketch(w,d,1)
# cms.add("ATGGTGAAACCTCGTCTCTACTAAAAATACAAAAAAAAATTAGCCGGATGTGGTGGCGGGCACCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCCTGAACCCGGGAGGCGGAGCTT")
# cms.add("ACATCCTGGCTAGCATGGTGAAACCTCGTCTCTACTAAAAATACAAAAAAAAATTAGCCGGATGTGGTGGCGGGCACCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCCTGAACCC")
# a=cms.count("ACATCCTGGCTAGCATGGTGAAACCTCGTCTCTACTAAAAATACAAAAAAAAATTAGCCGGATGTGGTGGCGGGCACCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCCTGAACCC")
# print(a)