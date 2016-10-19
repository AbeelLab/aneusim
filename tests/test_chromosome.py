import pytest

from skbio import DNA
from aneugen import chromosome


def test_translocation():
    chromosome1 = DNA('AAAAAAAACCCC')
    chromosome2 = DNA('TTTTGGGGGGGG')

    new1, new2 = chromosome.simulate_translocate(chromosome1, chromosome2,
                                                 4, 8)
    assert new1 == DNA('AAAAAAAAGGGGGGGG')
    assert new2 == DNA('TTTTCCCC')

    new1, new2 = chromosome.simulate_translocate(chromosome1, chromosome2,
                                                 8, 4, 1)
    assert new1 == DNA('TTTTCCCC')
    assert new2 == DNA('AAAAAAAAGGGGGGGG')


def test_mutations():
    sequence = DNA('AAAAAAAA')
    mutated = chromosome.mutate(sequence, 2)

    # Calculate hamming distance
    assert sequence.distance(mutated) == 0.25
