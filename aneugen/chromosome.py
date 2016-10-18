from typing import Tuple
from skbio import DNA


def simulate_translocate(chromosome1: DNA, chromosome2: DNA, length1: int,
                         length2: int, mode: int=0) -> Tuple[DNA, DNA]:
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
        new_chromosome1 = DNA.concat([chromosome1[:-length1],
                                     chromosome2[-length2:]])
        new_chromosome2 = DNA.concat([chromosome2[:-length2],
                                     chromosome1[-length1:]])
    elif mode == 1:
        new_chromosome1 = DNA.concat([chromosome2[:length2],
                                     chromosome1[length1:]])
        new_chromosome2 = DNA.concat([chromosome1[:length1],
                                     chromosome2[length2:]])
    else:
        raise ValueError("Invalid mode '{}' specified. Possible values are 0 "
                         "or 1".format(mode))

    return new_chromosome1, new_chromosome2
