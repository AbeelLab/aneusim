#!/usr/bin/env python3
import os
import sys
import random
import logging
import argparse
from configparser import ConfigParser
from typing import Tuple

from skbio import io, DNA


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


#
# Command line interface entry functions
# --------------------------------------
#
def list_records(args):
    for seq in io.read(args.input, format='fasta'):
        print('ID:', seq.metadata['id'], 'Description:',
              seq.metadata['description'])


def contigize(args):
    if not os.path.isdir(args.output_dir):
        raise FileNotFoundError("{} is not an existing directory.".format(
            args.output_dir
        ))

    copy_num_spec = ConfigParser()
    copy_num_spec.read_file(args.spec_file)

    copy_numbers = copy_num_spec['chromosomes']

    for seq in io.read(args.input, format='fasta'):
        chromosome_id = seq.metadata['id']
        if chromosome_id in copy_numbers:
            try:
                copy_number = int(copy_numbers[chromosome_id])
            except ValueError:
                logging.error("Invalid copy number value '{}' for chromosome "
                              "'{}'".format(copy_numbers[chromosome_id],
                                            chromosome_id))
                copy_number = 1
        else:
            logging.warning('No copy number specification found for '
                            'chromosome {}. Assuming a single copy.'.format(
                                chromosome_id))
            copy_number = 1

        for i in range(copy_number):
            output_path = os.path.join(
                args.output_dir,
                "{}.copy{}.fasta".format(chromosome_id, i)
            )

            seq.metadata['id'] = "{}.copy{}".format(chromosome_id, i)

            with open(output_path, 'w') as f:
                io.write(seq, format='fasta', into=f)


def translocate(args):
    chromosome1 = None
    chromosome2 = None

    # Chromosomes each in their own file
    with open(args.chromosome[0]) as f:
        chromosome1 = DNA.read(f, format='fasta', lowercase=True)

    with open(args.chromosome[1]) as f:
        chromosome2 = DNA.read(f, format='fasta', lowercase=True)

    if not chromosome1 or not chromosome2:
        raise ValueError("At least one of the chromosome data could not be "
                         "read.")

    new1, new2 = simulate_translocate(chromosome1, chromosome2,
                                               args.lengths[0],
                                               args.lengths[1])

    if args.in_place:
        file1 = args.chromosome[0]
        file2 = args.chromosome[1]
    else:
        ext_pos = args.chromosome[0].rfind('.')
        file1 = "{}.translocated{}".format(
            args.chromosome[0][:ext_pos], args.chromosome[0][ext_pos:])
        ext_pos = args.chromosome[1].rfind('.')
        file2 = "{}.translocated{}".format(
            args.chromosome[1][:ext_pos], args.chromosome[1][ext_pos:])

    with open(file1, 'w') as f:
        io.write(new1, format='fasta', into=f)

    with open(file2, 'w') as f:
        io.write(new2, format='fasta', into=f)


def mutate(args):
    filename = args.file.name
    sequence = DNA.read(args.file, format='fasta', lowercase=True)
    args.file.close()
    sequence = add_mutations(sequence, args.num)

    if args.in_place:
        with open(filename, 'w') as f:
            io.write(sequence, format='fasta', into=f)
    else:
        io.write(sequence, format='fasta', into=args.output)


def main():
    parser = argparse.ArgumentParser(
        description="A script to help generate synthetic aneuploid genomes."
    )

    subparsers = parser.add_subparsers()

    list_parser = subparsers.add_parser(
        'list', help="List records in a FASTA file."
    )

    list_parser.set_defaults(func=list_records)
    list_parser.add_argument(
        '-i', '--input', type=argparse.FileType('r'), required=True,
        help="The FASTA file to read"
    )

    contigize_parser = subparsers.add_parser(
        'contigize',
        help="Generate individual FASTA files for each chromosome haplotype."
    )

    contigize_parser.set_defaults(func=contigize)
    contigize_parser.add_argument(
        '-i', '--input', type=argparse.FileType('r'), required=True,
        help="Base genome as FASTA file. "
             "Each record should represent a chromosome."
    )
    contigize_parser.add_argument(
        '-s', '--spec-file', type=argparse.FileType('r'), required=True,
        help="Copy number specification for each chromosome, formatted as an "
             "INI file."
    )
    contigize_parser.add_argument(
        'output_dir',
        help="Specify the output directory where all FASTA files will be "
             "stored."
    )

    translocate_parser = subparsers.add_parser(
        'translocate', help="Simulate chromosome translocation."
    )

    translocate_parser.set_defaults(func=translocate)
    translocate_parser.add_argument(
        'chromosome', nargs=2, metavar='chromosome',
        help="Specify the two chromosomes as separate FASTA files."
    )
    translocate_parser.add_argument(
        '-l', '--lengths', nargs=2, metavar='LEN', required=True, type=int,
        help="Number of basepairs translocated for chromosome 1 and 2 "
             "respectively."
    )
    translocate_parser.add_argument(
        '-i', '--in-place', action="store_true", default=False,
        help="Modify the chromosome fasta files in place. Otherwise output "
             "the modified chromosomes as new files in the same directory as "
             "the original files."
    )

    mutate_parser = subparsers.add_parser(
        'mutate', help="Randomly add synthetic mutations to a chromosome."
    )

    mutate_parser.set_defaults(func=mutate)
    mutate_parser.add_argument(
        '-n', '--num', type=int, required=True,
        help="Specify the number of mutations to generate."
    )
    mutate_parser.add_argument(
        'file', type=argparse.FileType('r'), default=sys.stdin,
        help="The FASTA file with the chromosome sequence to read. Defaults "
             "to stdin."
    )
    mutate_parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        required=False,
        help="Output file of the new mutated chromosome, defaults to stdout."
    )
    mutate_parser.add_argument(
        '-i', '--in-place', action="store_true", default=False,
        help="Modify the file in place. Does not work if reading from stdin, "
             "and supersedes the --output option."
    )

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
