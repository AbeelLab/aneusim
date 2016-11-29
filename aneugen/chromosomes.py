import random
from typing import List, NamedTuple, Dict, Tuple

import pybedtools
from skbio import DNA

Annotation = NamedTuple('Annotation', [
    ('name', str), ('start', int), ('length', int), ('trackIndex', 0)
])

CentromereMeta = NamedTuple('CentromereMeta', [
    ('start', int), ('length', int)
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

    The output dict can be exported to JSON to be used with Ideogram.js.
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
    else:
        raise ValueError("Invalid mode '{}' specified. Possible values are 0 "
                         "or 1".format(mode))

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
