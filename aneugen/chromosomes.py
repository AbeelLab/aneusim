import random
from typing import List, NamedTuple, Tuple
from collections.abc import Sequence

import numpy
import pybedtools
from skbio import DNA

Annotation = NamedTuple('Annotation', [
    ('name', str), ('start', int), ('length', int), ('trackIndex', 0)
])

CentromereMeta = NamedTuple('CentromereMeta', [
    ('start', int), ('len', int)
])

ChromosomeMeta = NamedTuple('ChromosomeMeta', [
    ('name', str), ('type', str), ('accession', str),
    ('centromere', CentromereMeta), ('annots', List[Annotation])
])


def get_chromosome_metadata(bt: pybedtools.BedTool,
                            chromosome_id: str) -> ChromosomeMeta:
    """
    Retrieves chromosome name, length and centromere metadata from a
    BedTool object.
    """
    output = {
        'accession': chromosome_id,
    }

    metadata = bt.filter(
        lambda e: e[0] == chromosome_id
    ).filter(
        lambda e: e[2] == 'region' or e[2] == 'centromere'
    )

    for feat in metadata:
        if feat[2] == 'region':
            output['length'] = len(feat)
            output['name'] = feat.attrs['Name']
            output['type'] = ('nuclear' if feat.attrs['genome'] ==
                              'chromosome' else 'mitochondrial')
        elif feat[2] == 'centromere':
            if 'Dbxref' in feat.attrs:
                output['centromere'] = CentromereMeta(feat.start, len(feat))

    return ChromosomeMeta(output['name'], output['type'],
                          output['accession'], output['centromere'], [])


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
    elif mode == 2:
        # Get a suffix from chromosome 1, reverse it, and concatenate it
        # with a suffix of chromosome 2
        new_chromosome1 = DNA.concat([
            chromosome1[-length1:].reverse_complement(),
            chromosome2[length2:]
        ])

        # Get a prefix of chromosome 1, and concatenate it with a reversed
        # prefix of chromosome 2
        new_chromosome2 = DNA.concat([
            chromosome1[:-length1],
            chromosome2[:length2].reverse_complement()
        ])
    else:
        raise ValueError("Invalid mode '{}' specified. Possible values are 0, "
                         "1 or 2".format(mode))

    return new_chromosome1, new_chromosome2


def add_mutations(sequence: DNA, num: int) -> DNA:
    """
    Randomly mutate a given DNA sequence.

    :param sequence: DNA sequence of the chromosome
    :param num: Number of mutations to apply.
    """

    alphabet = {'A', 'C', 'T', 'G'}
    new_sequence = sequence.copy(deep=True)

    # Determine `num` random positions
    for pos in random.sample(range(len(sequence)), num):
        base = str(sequence[pos])

        new_base = random.choice(list(alphabet - {base}))
        new_sequence = new_sequence.replace([pos], new_base)

    return new_sequence


def generate_deletions(sequence: DNA, num: int, mu: int=20, std: float=6.0):
    i = 0
    new_sequence = sequence.copy(deep=True)

    while i < num:
        size = int(round(numpy.random.normal(mu, std)))
        pos = random.randint(0, len(new_sequence) - size)

        new_sequence = DNA.concat([new_sequence[:pos],
                                   new_sequence[pos+size:]])
        i += 1

    return new_sequence


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
