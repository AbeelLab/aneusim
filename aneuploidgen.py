#!/usr/bin/env python3
import argparse
import logging
import os
import sys
from configparser import ConfigParser

from dinopy import FastaReader, FastaWriter
from pybedtools import BedTool

from aneugen.chromosomes import (simulate_translocate, add_mutations,
                                 generate_deletions, find_ty_element_location)


#
# Command line interface entry functions
# --------------------------------------
#
def list_records(args):
    fr = FastaReader(args.input)
    for c in fr.entries():
        print('ID:', c.name, 'Length:', c.length)


def haplotise(args):
    """Generate separate FASTA files for each chromosome copy."""

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exists_ok=True)

    copy_num_spec = ConfigParser()
    copy_num_spec.read_file(args.spec_file)

    copy_numbers = copy_num_spec['chromosomes']

    fr = FastaReader(args.input)
    for c in fr.entries():
        parts = c.name.decode('utf-8').split()
        chromosome_id = parts[0]
        if chromosome_id in copy_numbers:
            try:
                copy_number = int(copy_numbers[chromosome_id])
            except ValueError:
                logging.error(
                    "Invalid copy number value '{}' for chromosome "
                    "'{}'".format(copy_numbers[chromosome_id], chromosome_id)
                )
                copy_number = 1
        else:
            logging.warning(
                'No copy number specification found for '
                'chromosome {}. Assuming a single copy.'.format(chromosome_id)
            )
            copy_number = 1

        for i in range(copy_number):
            output_path = os.path.join(
                args.output_dir,
                "{}.copy{}.fasta".format(chromosome_id, i)
            )

            copy_name = "{}.copy{} {}".format(chromosome_id, i,
                                              " ".join(parts[1:]))

            with FastaWriter(output_path, force_overwrite=True) as fw:
                fw.write_chromosome((c.sequence, copy_name.encode('utf-8')))


def translocate(args):
    if not args.lengths and not args.pos:
        raise argparse.ArgumentError(
            'lengths', "Neither translocation lengths or break positions "
                       "given."
        )

    chromosome1 = None
    chromosome2 = None

    # Chromosomes each in their own file
    f = FastaReader(args.chromosome[0])
    chromosome1 = next(f.entries())

    f = FastaReader(args.chromosome[1])
    chromosome2 = next(f.entries())

    if not chromosome1 or not chromosome2:
        raise ValueError("At least one of the chromosome data could not be "
                         "read.")

    if args.lengths:
        lengths = args.lengths
    else:
        lengths = [-1, -1]

    if args.pos:
        if args.mode == 0:
            lengths[0] = chromosome1.length - args.pos[0]
            lengths[1] = chromosome2.length - args.pos[1]
        elif args.mode == 1:
            lengths[0] = args.pos[0]
            lengths[1] = args.pos[1]
        elif args.mode == 2:
            lengths[0] = chromosome1.length - args.pos[0]
            lengths[1] = args.pos[1]

    new1, new2 = simulate_translocate(
        chromosome1.sequence, chromosome2.sequence, lengths[0], lengths[1])

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

    with FastaWriter(file1, force_overwrite=True) as fw:
        fw.write_chromosome((new1, chromosome1.name))

    with FastaWriter(file2, force_overwrite=True) as fw:
        fw.write_chromosome((new2, chromosome2.name))


def mutate(args):
    filename = args.file.name
    fr = FastaReader(args.file)
    entry = next(fr.entries(dtype=bytearray))
    new_sequence = add_mutations(entry.sequence, args.num)
    args.file.close()

    into_file = filename if args.in_place else args.output
    with FastaWriter(into_file, force_overwrite=True) as fw:
        fw.write_chromosome((new_sequence, entry.name), dtype=bytearray)


def deletions(args):
    filename = args.file.name

    fr = FastaReader(args.file)
    entry = next(fr.entries())
    args.file.close()

    sequence = generate_deletions(entry.sequence, args.num, args.mu,
                                  args.sigma)

    out_file = filename if args.in_place else args.output
    with FastaWriter(out_file, force_overwrite=True) as fw:
        fw.write_chromosome((sequence, entry.name))


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
        '-i', '--input', type=argparse.FileType('rb'), required=True,
        help="The FASTA file to read"
    )

    haplotise_parser = subparsers.add_parser(
        'haplotise',
        help="Generate individual FASTA files for each chromosome haplotype."
    )

    haplotise_parser.set_defaults(func=haplotise)
    haplotise_parser.add_argument(
        '-i', '--input', type=argparse.FileType('rb'), required=True,
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
        'file', type=argparse.FileType('rb'), default=sys.stdin,
        help="The FASTA file with the chromosome sequence to read. Defaults "
             "to stdin."
    )
    mutate_parser.add_argument(
        '-o', '--output', type=argparse.FileType('wb'), default=sys.stdout,
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
        'file', type=argparse.FileType('rb'), default=sys.stdin,
        help="The FASTA file with the chromosome sequence to read. Defaults "
             "to stdin."
    )
    deletions_parser.add_argument(
        '-o', '--output', type=argparse.FileType('wb'), default=sys.stdout,
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
