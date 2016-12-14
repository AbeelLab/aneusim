#!/usr/bin/env python3
import argparse
import logging
import os
import sys
from configparser import ConfigParser

from skbio import io, DNA
from pybedtools import BedTool, Interval

from aneugen.chromosomes import (simulate_translocate, add_mutations,
                                 generate_deletions, find_ty_element_location)


#
# Command line interface entry functions
# --------------------------------------
#
def list_records(args):
    for seq in io.read(args.input, format='fasta'):
        print('ID:', seq.metadata['id'], 'Description:',
              seq.metadata['description'])


def haplotise(args):
    """Generate separate FASTA files for each chromosome copy."""

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

    if not args.lengths and not args.pos:
        raise argparse.ArgumentError(
            'lengths', "Neither translocation lengths or break positions "
                       "given."
        )

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

    if args.lengths:
        lengths = args.lengths
    else:
        lengths = [-1, -1]

    if args.pos:
        if args.mode == 0:
            lengths[0] = len(chromosome1) - args.pos[0]
            lengths[1] = len(chromosome2) - args.pos[1]
        elif args.mode == 1:
            lengths[0] = args.pos[0]
            lengths[1] = args.pos[1]
        elif args.mode == 2:
            lengths[0] = len(chromosome1) - args.pos[0]
            lengths[1] = args.pos[1]

    new1, new2 = simulate_translocate(chromosome1, chromosome2,
                                      lengths[0], lengths[1])

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


def deletions(args):
    filename = args.file.name
    sequence = DNA.read(args.file, format='fasta', lowercase=True)
    args.file.close()
    sequence = generate_deletions(sequence, args.num, args.mu, args.sigma)

    if args.in_place:
        with open(filename, 'w') as f:
            io.write(sequence, format='fasta', into=f)
    else:
        io.write(sequence, format='fasta', into=args.output)


def find_ty_elems(args):
    annotations = BedTool(args.annotations_file).filter(
        lambda f: f[0] == args.chromosome
    )

    ty_elements = find_ty_element_location(annotations, (args.start, args.end))
    for interval in ty_elements:
        print(interval.start, interval.end)


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

    haplotise_parser = subparsers.add_parser(
        'haplotise',
        help="Generate individual FASTA files for each chromosome haplotype."
    )

    haplotise_parser.set_defaults(func=haplotise)
    haplotise_parser.add_argument(
        '-i', '--input', type=argparse.FileType('r'), required=True,
        help="Base genome as FASTA file. "
             "Each record should represent a chromosome."
    )
    haplotise_parser.add_argument(
        '-s', '--spec-file', type=argparse.FileType('r'), required=True,
        help="Copy number specification for each chromosome, formatted as an "
             "INI file."
    )
    haplotise_parser.add_argument(
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
        '-m', '--mode', choices=(0, 1, 2), default=0, type=int,
        help="Choose translocation mode."
    )
    translocate_parser.add_argument(
        '-l', '--lengths', nargs=2, metavar='LEN', type=int,
        help="Number of basepairs translocated for chromosome 1 and 2 "
             "respectively."
    )
    translocate_parser.add_argument(
        '-p', '--pos', nargs=2, type=int,
        help="Breaking positions for chromosome 1 and 2. Overrides the "
             "--lengths option."
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

    deletions_parser = subparsers.add_parser(
        'deletions', help="Randomly generate deletions throughout the genome."
    )

    deletions_parser.set_defaults(func=deletions)
    deletions_parser.add_argument(
        '-u', '--mu', type=int, default=20,
        help="Mean size of a deletion. Defaults to 20."
    )
    deletions_parser.add_argument(
        '-s', '--sigma', type=float, default=6.0,
        help="Standard deviation of the size of a deletion. Defaults to 6"
    )
    deletions_parser.add_argument(
        '-n', '--num', type=int, required=True,
        help="Total number of deletions to generate."
    )
    deletions_parser.add_argument(
        'file', type=argparse.FileType('r'), default=sys.stdin,
        help="The FASTA file with the chromosome sequence to read. Defaults "
             "to stdin."
    )
    deletions_parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        required=False,
        help="Output file of the new mutated chromosome, defaults to stdout."
    )
    deletions_parser.add_argument(
        '-i', '--in-place', action="store_true", default=False,
        help="Modify the file in place. Does not work if reading from stdin, "
             "and supersedes the --output option."
    )

    find_tyelems_parser = subparsers.add_parser(
        'find_tyelems', help="Find the positions of Ty-elements in an "
                            "chromosome."
    )

    find_tyelems_parser.set_defaults(func=find_ty_elems)
    find_tyelems_parser.add_argument(
        'annotations_file', help="BED or GFF file with genome annotations."
    )
    find_tyelems_parser.add_argument(
        'chromosome', help="Specify chromosome ID"
    )
    find_tyelems_parser.add_argument(
        '-s', '--start', type=int, default=-1,
        help="Start position of the range to search in",

    )
    find_tyelems_parser.add_argument(
        '-e', '--end', type=int, default=-1,
        help="End position of the range to search in",
    )

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
