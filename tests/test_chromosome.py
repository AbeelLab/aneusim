import pytest

from skbio import DNA
from aneugen.chromosomes import (simulate_translocate, add_mutations,
                                 generate_deletions)


def test_translocation():
    chromosome1 = DNA('GAAAAATGCCCT')
    chromosome2 = DNA('CTTTACGGGGGA')

    new1, new2 = simulate_translocate(chromosome1, chromosome2, 4, 8)
    assert new1 == DNA('GAAAAATGACGGGGGA')
    assert new2 == DNA('CTTTCCCT')

    new1, new2 = simulate_translocate(chromosome1, chromosome2, 8, 4, 1)
    assert new1 == DNA('CTTTCCCT')
    assert new2 == DNA('GAAAAATGACGGGGGA')

    new1, new2 = simulate_translocate(chromosome1, chromosome2, 4, 4, 2)
    assert new1 == DNA('AGGGACGGGGGA')
    assert new2 == DNA('GAAAAATGAAAG')


def test_mutations():
    sequence = DNA('AAAAAAAA')
    mutated = add_mutations(sequence, 2)

    # Calculate hamming distance
    assert sequence.distance(mutated) == 0.25


def test_deletions():
    sequence = DNA('AAAAAAAA')
    mutated = generate_deletions(sequence, 1, 3, 0.1)

    assert len(mutated) == len(sequence) - 3

    sequence = DNA('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    mutated = generate_deletions(sequence, 3, 3, 0.1)

    assert len(mutated) == len(sequence) - 9


