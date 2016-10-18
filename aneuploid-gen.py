#!/usr/bin/env python3
import sys
import argparse

from skbio import io, DNA

from aneugen.chromosome import simulate_translocate


def list_records(args):
    for seq in io.read(args.input, format='fasta'):
        print('ID:', seq.metadata['id'], 'Description:',
              seq.metadata['description'])


def copy(args):
    pass


def translocate(args):
    chromosome1 = None
    chromosome2 = None

    if args.input:
        # Chromosomes are individual records in a single FASTA file
        for seq in io.read(args.input, 'fasta', lowercase=True):
            seq = DNA(seq)
            if seq.metadata['id'] == args.chromosome[0]:
                chromosome1 = seq

            if seq.metadata['id'] == args.chromosome[1]:
                chromosome2 = seq

            if chromosome1 and chromosome2:
                break
    else:
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

    io.write(new1, format='fasta', into=args.output)
    io.write(new2, format='fasta', into=args.output)


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

    copy_parser = subparsers.add_parser(
        'copy', help="Copy a specified chromosome and possibly add synthetic "
                     "SNP's"
    )

    copy_parser.set_defaults(func=copy)
    copy_parser.add_argument(
        '-i', '--input', type=argparse.FileType('r'), required=True,
        help="Base genome as FASTA file. "
             "Each record should represent a chromosome."
    )
    copy_parser.add_argument(
        '-c', '--chromosome', required=True,
        help="Specify the FASTA record ID of the chromosome to copy"
    )
    copy_parser.add_argument(
        '-s', '--snp', type=float, required=False, metavar='NUM',
        help="Generate NUM random SNP's in the new chromosome copy."
    )

    translocate_parser = subparsers.add_parser(
        'translocate', help="Simulate chromosome translocation."
    )

    translocate_parser.set_defaults(func=translocate)
    translocate_parser.add_argument(
        'chromosome', nargs=2, metavar='chromosome',
        help="Specify the two chromosomes either as files or FASTA ID's, "
             "depening on wether the -i flag is used."
    )
    translocate_parser.add_argument(
        '-i', '--input', type=argparse.FileType('r'), required=False,
        metavar='FASTA_FILE',
        help="Base genome as FASTA file. "
             "Each record should represent a chromosome."
    )
    translocate_parser.add_argument(
        '-o', '--output', type=argparse.FileType('a'), default=sys.stdout,
        metavar='FILE',
        help="Specify the output file. Defaults to stdout."
    )
    translocate_parser.add_argument(
        '-l', '--lengths', nargs=2, metavar='LEN', required=True, type=int,
        help="Number of basepairs translocated for chromosome 1 and 2 "
             "respectively."
    )

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
