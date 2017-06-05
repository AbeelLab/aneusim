"""
Utilities to simulate error-free "PacBio" reads from a genome

A lot of code inspired by SimLord:
https://bitbucket.org/genomeinformatics/simlord/
"""

import random
from typing import NamedTuple, List, Tuple

from dinopy import FastaReader, reverse_complement
import numpy
from scipy.stats import lognorm

_ReadLengthDistribution = NamedTuple(
    'ReadLengthDistribution', [
        ('s', float), ('loc', int), ('scale', float)
    ]
)


class ReadLengthDistribution(_ReadLengthDistribution):
    """Represents a log-normal read distribution specified by its
    parameters."""

    def num_reads_for_coverage(self, coverage: int, genome_length: int,
                               min_length: int=0) -> int:
        expected_length = lognorm.expect(
            lambda x: max(x, min_length),
            args=[self.s], loc=self.loc, scale=self.scale,
            lb=-numpy.inf, ub=numpy.inf
        )

        return int(coverage * genome_length / expected_length)

    def get_read_lengths(self, num: int, min_length: int=0):
        sizes = lognorm.rvs(s=self.s, loc=self.loc, scale=self.scale, size=num)

        return [max(min_length, s) for s in sizes]


def read_reference(f: str) -> List[Tuple]:
    fr = FastaReader(f)
    # Use bytearrays because we want to be able to change bases
    # in reads later
    return list(fr.entries(dtype=bytearray))


def generate_reads(chromosomes, dist: ReadLengthDistribution, num_reads: int,
                   min_length: int):
    total_length = sum(c.length for c in chromosomes)
    max_chromosome_length = max(c.length for c in chromosomes)
    weights = [c.length / total_length for c in chromosomes]

    read_lengths = dist.get_read_lengths(num_reads, min_length)
    for curr_read_length in read_lengths:
        curr_read_length = int(min(max_chromosome_length, curr_read_length))

        # Weigh the selection of chromosome by its length, and make sure the
        # chromosome is long enough to be able a produce a read of this size.
        # By weighing the chromosome selection by its length, we ensure equal
        # coverage across the whole genome.
        index = -1
        chromosome = None
        while not chromosome or chromosome.length < curr_read_length:
            index = numpy.random.choice(len(chromosomes), p=weights)
            chromosome = chromosomes[index]

        # Pick a random start location
        start_pos = random.randint(0, chromosome.length-curr_read_length)
        read = chromosome.sequence[start_pos:start_pos+curr_read_length]

        # Replace N's with random bases
        for i, base in enumerate(read):
            if base in {b'n', b'N'}:
                read[i] = random.choice([b'A', b'C', b'T', b'G'])

        # Convert half of the reads to its reverse complement
        if random.choice([0, 1]) == 1:
            read = reverse_complement(read)

        yield read, chromosome.name, start_pos
