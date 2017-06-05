import pytest  # noqa

from aneugen.structural import simulate_translocate, generate_deletions


def test_translocation():
    chromosome1 = b'GAAAAATGCCCT'
    chromosome2 = b'CTTTACGGGGGA'

    new1, new2 = simulate_translocate(chromosome1, chromosome2, 4, 8)
    assert new1 == b'GAAAAATGACGGGGGA'
    assert new2 == b'CTTTCCCT'

    new1, new2 = simulate_translocate(chromosome1, chromosome2, 8, 4, 1)
    assert new1 == b'CTTTCCCT'
    assert new2 == b'GAAAAATGACGGGGGA'

    new1, new2 = simulate_translocate(chromosome1, chromosome2, 4, 4, 2)
    assert new1 == b'AGGGACGGGGGA'
    assert new2 == b'GAAAAATGAAAG'


def test_deletions():
    sequence = b'AAAAAAAA'
    mutated = generate_deletions(sequence, 1, 3, 0.1)

    assert len(mutated) == len(sequence) - 3

    sequence = b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    mutated = generate_deletions(sequence, 3, 3, 0.1)

    assert len(mutated) == len(sequence) - 9
