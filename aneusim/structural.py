from typing import Tuple

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
