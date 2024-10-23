# src/KmerDecon/utils.py

from typing import Iterator

def generate_kmers(sequence: str, k: int) -> Iterator[str]:
    """
    Generate k-mers from a given sequence.

    Args:
        sequence (str): The nucleotide sequence.
        k (int): Length of k-mers to generate.

    Yields:
        Iterator[str]: A generator of k-mers.
    """
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]