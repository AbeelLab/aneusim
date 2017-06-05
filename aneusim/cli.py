import os
import sys
import json
import logging
import argparse
from configparser import ConfigParser

from dinopy import FastaReader, FastaWriter
from pybedtools import BedTool

from aneusim.structural import (simulate_translocate, generate_deletions,
                                find_ty_element_location)
from aneusim.haplogen import MutationGenerator, MutationType
from aneusim.haplo_spec import get_distance_model, get_dosage
from aneusim.reads import (ReadLengthDistribution, generate_reads,
                           read_reference)

logger = logging.getLogger(__name__)


#
# Command line interface entry functions
# --------------------------------------
#
def list_records(args):
    fr = FastaReader(args.input)
    total_length = 0
    for c in fr.entries():
        print('ID:', c.name, 'Length:', c.length)
        total_length += c.length

    print("Total length:")
    print(total_length, "bp")
    print(total_length/1000, "kbp")
    print(total_length/1000/1000, "Mbp")


def haplogen(args):
    """Generate separate FASTA files for each chromosome copy."""

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exists_ok=True)

    haplotypes_spec = ConfigParser()
    haplotypes_spec.read_file(args.spec_file)

    fr = FastaReader(args.input)
    for c in fr.entries(dtype=bytearray):
        parts = c.name.decode('utf-8').split()
        chromosome_id = parts[0]

        if chromosome_id not in haplotypes_spec:
            logger.warning("No specification on how to generate haplotypes "
                           "for chromosome '{}', ignoring".format(
                               chromosome_id))
            continue

        chromosome_spec = haplotypes_spec[chromosome_id]

        if 'ploidy' not in chromosome_spec:
            logger.error("No ploidy given for chromosome '{}', "
                         "ignoring".format(chromosome_id))
            continue

        ploidy = chromosome_spec.getint('ploidy')
        logger.info("Processing reference chromosome %s with ploidy %d",
                    chromosome_id, ploidy)

        # Obtain distance models
        substitution_dist_model = get_distance_model(chromosome_spec,
                                                     MutationType.SUBSTITUTION)
        insertion_dist_model = get_distance_model(chromosome_spec,
                                                  MutationType.INSERTION)
        deletion_dist_model = get_distance_model(chromosome_spec,
                                                 MutationType.DELETION)

        # Get dosage distribution
        dosage_dist = get_dosage(chromosome_spec)
        biallelic = chromosome_spec.getboolean('biallelic', True)
        ld_reassign_prob = chromosome_spec.getfloat('ld_reassign_prob', 0.4)

        # TODO: allow for configuration of substitution and insertion rates
        mut_gen = MutationGenerator(c.sequence, ploidy, dosage_dist,
                                    biallelic=biallelic,
                                    ld_reassign_prob=ld_reassign_prob)

        haplotypes, mutation_record = mut_gen.generate(
            substitution_dist_model, insertion_dist_model, deletion_dist_model)

        for i, haplotype in enumerate(haplotypes):
            output_path = os.path.join(
                args.output_dir,
                "{}.copy{}.fasta".format(chromosome_id, i)
            )

            copy_name = "{}.copy{} {}".format(chromosome_id, i,
                                              " ".join(parts[1:]))

            with FastaWriter(output_path, force_overwrite=True) as fw:
                fw.write_chromosome((haplotype, copy_name.encode('utf-8')),
                                    dtype=bytearray)

        if args.record_mutations:
            output_file = os.path.join(
                args.output_dir, "{}_mutations.json".format(chromosome_id))

            with open(output_file, "w") as f:
                json.dump(mutation_record.mutations, f, indent=2)


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


def reads(args):
    if not args.coverage and not args.num:
        print("No number of reads given.")
        return

    chromosomes = read_reference(args.input_genome)
    total_length = sum(c.length for c in chromosomes)

    dist = ReadLengthDistribution(*args.lognorm_params)
    num = args.num
    if args.coverage:
        num = dist.num_reads_for_coverage(args.coverage, total_length)
        print(num, "reads required for coverage of", args.coverage,
              file=sys.stderr)

    metadata = {}
    with FastaWriter(args.output, force_overwrite=True) as fw:
        for i, read_data in enumerate(generate_reads(chromosomes, dist, num,
                                                     args.min_length)):
            read, chromosome_name, start_pos = read_data
            chromosome_id = chromosome_name.decode('utf-8').split()[0]

            read_name = "read{}".format(i)
            metadata[read_name] = {
                'chromosome': chromosome_id,
                'start_pos': start_pos,
                'len': len(read)
            }

            fw.write_entry((read, read_name.encode('utf-8')))

    if args.metadata:
        json.dump(metadata, args.metadata)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        description="A script to help generate synthetic aneuploid genomes."
    )
    parser.set_defaults(func=None)

    subparsers = parser.add_subparsers()

    list_parser = subparsers.add_parser(
        'list', help="List records in a FASTA file."
    )

    list_parser.set_defaults(func=list_records)
    list_parser.add_argument(
        '-i', '--input', type=argparse.FileType('rb'), required=True,
        help="The FASTA file to read"
    )

    haplogen_parser = subparsers.add_parser(
        'haplogen',
        help="Generate individual FASTA files for each chromosome haplotype."
    )

    haplogen_parser.set_defaults(func=haplogen)
    haplogen_parser.add_argument(
        'input', type=argparse.FileType('rb'),
        help="Base reference genome as FASTA file. "
             "Each record should represent a chromosome."
    )
    haplogen_parser.add_argument(
        '-s', '--spec-file', type=argparse.FileType('r'), required=True,
        help="Haplotype specification for each chromosome, formatted as an "
             "INI file."
    )
    haplogen_parser.add_argument(
        'output_dir',
        help="Specify the output directory where all FASTA files will be "
             "stored."
    )
    haplogen_parser.add_argument(
        '-r', '--record-mutations', action="store_true", default=False,
        help="Output an additional JSON file for each chromosome specifying "
             "which mutations where made on each haplotype."
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

    reads_parser = subparsers.add_parser(
        'reads', help="Simulate a long-read sequencing experiment, and "
                      "generate a set of reads sampled from the genome."
                      " This read simulator does not contain an error model."
    )

    reads_parser.set_defaults(func=reads)
    reads_parser.add_argument(
        '-ln', '--lognorm-params', type=float,
        default=(0.2001, -10075.4364, 17922.611), nargs=3,
        help="log-normal distribution parameters to sample read lengths from "
             " (sigma, loc, scale)."
    )
    reads_parser.add_argument(
        '-n', '--num', type=int, required=False,
        help="Give the number of reads to generate."
    )
    reads_parser.add_argument(
        '-c', '--coverage', type=int, required=False,
        help="Automatically calculate the number of reads to generate to "
             "ensure the given coverage across the genome. Overrides the "
             "parameter -n/--num"
    )
    reads_parser.add_argument(
        '-l', '--min-length', type=int, default=1000,
        help="Minimum read length. Defaults to 1000."
    )
    reads_parser.add_argument(
        '-m', '--metadata', type=argparse.FileType('w'), default=None,
        required=False,
        help="Specify the filename to output read metadata in JSON format."
    )
    reads_parser.add_argument(
        '-o', '--output', type=argparse.FileType('wb'), default=sys.stdout,
        help="The filename to store the reads in. Defaults to stdout."
    )
    reads_parser.add_argument(
        'input_genome', type=argparse.FileType('rb'), default=sys.stdin,
        help="The reference genome to sample reads from."
    )

    args = parser.parse_args()
    if not args.func:
        parser.print_help()
    else:
        args.func(args)
