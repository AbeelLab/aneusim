#!/usr/bin/env python3
import os
import logging
import argparse
from configparser import ConfigParser

from skbio import io, DNA

from aneugen.chromosome import simulate_translocate


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
        if seq.metadata['id'] in copy_numbers:
            try:
                copy_number = int(copy_numbers[seq.metadata['id']])
            except ValueError:
                logging.error("Invalid copy number value '{}' for chromosome "
                              "'{}'".format(copy_numbers[seq.metadata['id']],
                                            seq.metadata['id']))
                copy_number = 1
        else:
            logging.warning('No copy number specification found for '
                            'chromosome {}. Assuming a single copy.'.format(
                                seq.metadata['id']))
            copy_number = 1

        for i in range(copy_number):
            output_path = os.path.join(
                args.output_dir,
                "{}.copy{}.fasta".format(seq.metadata['id'], i)
            )

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
                                      args.lengths[0], args.lengths[1])

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
        '-i', '--in-place', type=bool, action="store_true", default=False,
        help="Modify the chromosome fasta files in place. Otherwise output "
             "the modified chromosomes as new files in the same directory as "
             "the original files."
    )

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
