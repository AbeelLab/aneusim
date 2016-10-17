#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser(
        description="A script to help generate synthetic aneuploid genomes."
    )

    parser.add_argument(
        '-i', '--input', type=argparse.FileType('r'),
        help="FASTA input file containing the base genome. "
             "Each record should represent a chromosome."
    )


