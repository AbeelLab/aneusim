`aneusim` - A tool to generate synthetic aneuploid/polyploid genomes
==================================================================

[![Build 
Status](https://travis-ci.org/AbeelLab/aneusim.svg?branch=master)](https://travis-ci.org/AbeelLab/aneusim)

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

#### CLI Usage

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

#### Haplotype specification configuration

The specification file given with the `-s` option needs to be an INI-style 
configuration file describing how to generate haplotypes from the reference. 
Add a configuration section in this file for each chromosome (or contig) that 
will be used to generate haplotypes. Example:

```ini
[DEFAULT]
substitutions.model = lognormal
substitutions.mean = 3.578
substitutions.sigma = 1.34
substitutions.loc = 2500

insertions.model = 0

deletions.model = lognormal
deletions.mean = 7.34
deletions.sigma = 1.123

deletion_size.model = exponential
deletion_size.lambda = 0.2

dosage_dist.1 = 1.0
dosage_dist.2 = 3/4, 1/4
dosage_dist.3 = 1/2, 2/6, 1/6
dosage_dist.4 = 1/2, 2/8, 1/8, 1/8
dosage_dist.5 = 2/5, 3/10, 1/10, 1/10, 1/10
dosage_dist.6 = 4/12, 3/12, 2/12, 1/12, 1/12, 1/12
dosage_dist.7 = 6/16, 4/16, 2/16, 1/16, 1/16, 1/16, 1/16
dosage_dist.8 = 6/16, 3/16, 2/16, 1/16, 1/16, 1/16, 1/16, 1/16

[BK006935.2]
ploidy = 2

[BK006947.3]
ploidy = 3
```

In this file each configuration section should match the ID of one of your 
FASTA entries. For each matching ID, it will generate `ploidy` amount of 
haplotypes, and randomly mutate each haplotype following the models configured 
by `substitutions`, `insertions`, `deletions`, `deletion_size` and 
`dosage_dist` (described in more detail below). The `DEFAULT` configuration 
section is a special one: this will provide default values for any option not 
set in other configuration sections. It is possible to override any default 
option in a specific, FASTA record matching, configuration section. If the 
FASTA reference file has entries (chromosomes/contigs) with no matching 
configuration section, then these entries will be ignored.

With the options `substitutions`, `insertions`, and `deletions` you can 
configure the distribution of the *distance between two consecutive mutations 
of the same type*. Take a look at the following lines in the above example:

```ini
substitutions.model = lognormal
substitutions.mean = 3.578
substitutions.sigma = 1.34
substitutions.loc = 2500
```

This configures our haplotype generator to generate base substitutions where 
the distance between two consecutive substitutions follow a log-normal 
distribution, with mean=3.578, sigma=1.34 and shifted to the right by 2500 
bases. The same holds for `insertions` and `deletions`. To disable a certain 
kind of mutation, set `type.model` to 0 (see `insertions` in the above 
example).

The haplotype generator has support the following models:

* `lognormal`, parameters `mean`, `sigma`, `loc` (optional). See 
  [SciPy][scipy-lognormal] docs and [Wikipedia][wiki-lognormal] for more 
  information on these parameters.
* `exponential`, parameter `lambda`. See [SciPy][scipy-expon] and 
  [Wikipedia][wiki-expon] for more information on these parameters.
* `fixed`, parameter `value`. Always the same value instead of a value drawn 
  from a probability distribution.

Furthermore, with `deletion_size` you can configure the size distribution of 
the generated deletions. It can follow any of the above described models.

[scipy-lognormal]:https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
[scipy-expon]:https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.expon.html 
[wiki-lognormal]:https://en.wikipedia.org/wiki/Log-normal_distribution
[wiki-expon]:https://en.wikipedia.org/wiki/Exponential_distribution

Now that we know *where* to put mutations, we can now decide *which* haplotypes 
should get the alternative allele. This is configured using `dosage_dist`. The 
*dosage* is defined as the number of haplotypes that get the alternative 
allele. The number after `dosage_dist` reflects the ploidy of the chromosome, 
so you can specify the discrete probability distribution for the *dosage* per 
ploidy. Example:

```ini
dosage_dist.3 = 1/2, 2/6, 1/6
```

For any mutation generated for a chromosome with ploidy 3, the probability that 
only one of the haplotypes gets the alternative allele is 1/2, the probability 
that two out of the three haplotypes get the alternative allele is 2/6 and so 
on.

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

