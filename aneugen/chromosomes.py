import random
from typing import Tuple
from collections.abc import Sequence

import numpy
import pybedtools
from dinopy import reverse_complement


def simulate_translocate(chromosome1: bytes, chromosome2: bytes, length1: int,
                         length2: int, mode: int=0) -> Tuple[bytes, bytes]:
    """
    Simulate chromosomal translation. This is done by cutting some prefix or
    suffix from one chromosome and and this piece to an other chromosome,
    and vice versa.

    .. seealso::
       https://en.wikipedia.org/wiki/Chromosomal_translocation

    :param chromosome1: DNA sequence of the first chromosome
    :param chromosome2: DNA sequence of the second chromosome
    :param length1: Number of bases to cut from the first chromosome and add
                    this to the second chromosome.
    :param length2: Number of bases to cut from the second chromosome and add
                    this to the first chromosome.
    :param mode: Mode=0 means modifying the suffixes of chromosomes, mode=1
                 means modifying the prefixes.

    """

    if mode == 0:
        new_chromosome1 = chromosome1[:-length1] + chromosome2[-length2:]
        new_chromosome2 = chromosome2[:-length2] + chromosome1[-length1:]
    elif mode == 1:
        new_chromosome1 = chromosome2[:length2] + chromosome1[length1:]
        new_chromosome2 = chromosome1[:length1] + chromosome2[length2:]
    elif mode == 2:
        # Get a suffix from chromosome 1, reverse it, and concatenate it
        # with a suffix of chromosome 2
        new_chromosome1 = (reverse_complement(chromosome1[-length1:]) +
                           chromosome2[length2:])

        # Get a prefix of chromosome 1, and concatenate it with a reversed
        # prefix of chromosome 2
        new_chromosome2 = (chromosome1[:-length1] +
                           reverse_complement(chromosome2[:length2]))
    else:
        raise ValueError("Invalid mode '{}' specified. Possible values are 0, "
                         "1 or 2".format(mode))

    return new_chromosome1, new_chromosome2


def add_mutations(sequence: bytearray, num: int) -> bytearray:
    """
    Randomly mutate a given DNA sequence.

    :param sequence: DNA sequence of the chromosome
    :param num: Number of mutations to apply.
    """

    alphabet_lower = b'actg'
    alphabet = b"ACTG"

    # Determine `num` random positions
    for pos in random.sample(range(len(sequence)), num):
        base = sequence[pos]

        base_index = alphabet_lower.find(base)
        if base_index != -1:
            base = alphabet[base_index]

        new_base = random.choice(list(set(alphabet) - {base}))
        sequence[pos] = new_base

    return sequence


def generate_deletions(sequence: bytearray, num: int, mu: int=20,
                       std: float=6.0):
    for i in range(num):
        size = int(round(numpy.random.normal(mu, std)))
        pos = random.randint(0, len(sequence) - size)

        sequence = sequence[:pos] + sequence[pos+size:]

    return sequence


def find_ty_element_location(annotations: pybedtools.BedTool, start_end=None):
    start, end = -1, -1

    if isinstance(start_end, Sequence):
        start, end = start_end
    elif isinstance(start_end, int):
        start = -1
        end = start_end

    if start > 0:
        annotations = annotations.filter(lambda f: f.start >= start)

    if end > 0:
        annotations = annotations.filter(lambda f: f.stop <= end)

    return annotations.filter(lambda f: f[2] == 'mobile_genetic_element')
