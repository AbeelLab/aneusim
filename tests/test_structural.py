import pytest  # noqa

from aneusim.structural import simulate_translocate


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
