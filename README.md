`aneusim` - A tool to generate synthetic aneuploid/polyploid genomes
==================================================================

[![Build 
Status](https://travis-ci.org/lrvdijk/aneusim.svg?branch=master)](https://travis-ci.org/lrvdijk/aneusim)

Using an existing (haploid) reference genome as base, this tools helps you to 
build a synthetic aneuploid or polyploid genome derived from your reference. 
This tool has features to generate random insertions, deletions and 
substitutions, using different generation models (log-normal, exponential, 
...), and it is possible to simulate a translocation between two chromosomes. 
Furthermore, it is also able to simulate a sequencing experiment by generating 
reads with a given log-normal read length distribution. However, this tool does 
not contain an error model for the sequencing experiment, and generates 
error-free reads.

Requirements
------------

* Python >= 3.5
* dinopy >= 1.2
* numpy >= 1.10
* scipy >= 1.17
* (tests) pytest

Installation
------------

The tool is not on PyPI yet, so to install this tool clone the repository and 
run:

    pip install -r requirements.txt
    python setup.py install

Documentation
-------------

This package provides the command `aneusim`, which has four subcommands.

### `aneusim list`

Simply list records in a FASTA file and provide some other data.

    usage: aneusim list [-h] -i INPUT

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            The FASTA file to read

### `aneusim haplogen`

Based on a given reference genome, generate different copies for each 
chromosome with possible random mutations.

TODO: document specification file.

    usage: aneusim haplogen [-h] -s SPEC_FILE [-r] input output_dir

    positional arguments:
      input                 Base reference genome as FASTA file. Each record
                            should represent a chromosome.
      output_dir            Specify the output directory where all FASTA files
                            will be stored.

    optional arguments:
      -h, --help            show this help message and exit
      -s SPEC_FILE, --spec-file SPEC_FILE
                            Haplotype specification for each chromosome, formatted
                            as an INI file.
      -r, --record-mutations
                            Output an additional JSON file for each chromosome
                            specifying which mutations where made on each
                            haplotype.

### `aneusim translocate` 

Simulate translocation between two chromosomes.

    usage: aneusim translocate [-h] [-m {0,1,2}] [-l LEN LEN] [-p POS POS] [-i]
                               chromosome chromosome

    positional arguments:
      chromosome            Specify the two chromosomes as separate FASTA files.

    optional arguments:
      -h, --help            show this help message and exit
      -m {0,1,2}, --mode {0,1,2}
                            Choose translocation mode.
      -l LEN LEN, --lengths LEN LEN
                            Number of basepairs translocated for chromosome 1 and
                            2 respectively.
      -p POS POS, --pos POS POS
                            Breaking positions for chromosome 1 and 2. Overrides
                            the --lengths option.
      -i, --in-place        Modify the chromosome fasta files in place. Otherwise
                            output the modified chromosomes as new files in the
                            same directory as the original files.

### `aneusim reads`

Generate error-free reads from a given genome with a given read length 
distribution.

    usage: aneusim reads [-h] [-ln LOGNORM_PARAMS LOGNORM_PARAMS 
    LOGNORM_PARAMS]
                         [-n NUM] [-c COVERAGE] [-l MIN_LENGTH] [-m METADATA]
                         [-o OUTPUT]
                         input_genome

    positional arguments:
      input_genome          The reference genome to sample reads from.

    optional arguments:
      -h, --help            show this help message and exit
      -ln LOGNORM_PARAMS LOGNORM_PARAMS LOGNORM_PARAMS, --lognorm-params LOGNORM_PARAMS LOGNORM_PARAMS LOGNORM_PARAMS
                            log-normal distribution parameters to sample read
                            lengths from (sigma, loc, scale).
      -n NUM, --num NUM     Give the number of reads to generate.
      -c COVERAGE, --coverage COVERAGE
                            Automatically calculate the number of reads to
                            generate to ensure the given coverage across the
                            genome. Overrides the parameter -n/--num
      -l MIN_LENGTH, --min-length MIN_LENGTH
                            Minimum read length. Defaults to 1000.
      -m METADATA, --metadata METADATA
                            Specify the filename to output read metadata in JSON
                            format.
      -o OUTPUT, --output OUTPUT
                            The filename to store the reads in. Defaults to
                            stdout.

