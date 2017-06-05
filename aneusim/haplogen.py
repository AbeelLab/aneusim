"""
:mod:`aneugen.mutate` - Models for generating random mutations in a genome
==========================================================================

This module contains functions to generate random mutations across a genome,
according to one of our stochastic models.

"""

import enum
import math
import random
from typing import List, Callable, Set, Mapping, Tuple
from collections import defaultdict

import numpy
from numpy import random as nprnd

DistanceModelT = Callable[[], int]
DosageDistributionT = List[float]
InsertionRatesT = List[float]
SubstitutionRatesT = Mapping[bytes, Mapping[bytes, float]]
HaplotypeReturnT = Tuple[List[bytearray], 'MutationRecord']

INSERTION_ALPHABET = (b'A', b'C', b'T', b'G')


class MutationType(enum.Enum):
    INSERTION = 1
    DELETION = 2
    SUBSTITUTION = 3


SPEC_MUTTYPE_KEY = {
    MutationType.SUBSTITUTION: "substitutions",
    MutationType.INSERTION: "insertions",
    MutationType.DELETION: "deletions"
}


class MutationRecord:
    def __init__(self, ploidy: int):
        self.ploidy = ploidy
        self.mutations = []
        for i in range(ploidy):
            self.mutations.append({})

    def record(self, haplotype: int, ref_pos: int, haplo_pos: int,
               alt_assignment: List[int], new_value: bytearray,
               mut_type: MutationType):

        self.mutations[haplotype][haplo_pos] = {
            'ref_pos': ref_pos,
            'alt_assignment': list(alt_assignment),
            'new_value': new_value.decode('ascii'),
            'mut_type': str(mut_type)
        }


class LogNormalDistanceModel:
    def __init__(self, mean: float, sigma: float):
        self.mean = mean
        self.sigma = sigma

    @classmethod
    def from_spec(cls, chromosome_spec, mut_type: MutationType):
        key = SPEC_MUTTYPE_KEY[mut_type]
        mean_key = "{}.mean".format(key)
        sigma_key = "{}.sigma".format(key)

        if mean_key not in chromosome_spec or sigma_key not in chromosome_spec:
            raise KeyError("Chromosome specification does not contain the "
                           "parameters for the log-normal distribution.")

        return cls(chromosome_spec.getfloat(mean_key),
                   chromosome_spec.getfloat(sigma_key))

    def __call__(self) -> int:
        return int(numpy.round(nprnd.lognormal(self.mean, self.sigma, 1)))


class MutationGenerator:
    def __init__(self, reference: bytearray, ploidy: int,
                 dosage_dist: DosageDistributionT,
                 substitution_rates: SubstitutionRatesT=None,
                 insertion_rates: InsertionRatesT=None,
                 biallelic: bool=False,
                 ld_reassign_prob: float=0.4):
        self.reference = reference.upper()
        self.ploidy = ploidy
        self.biallelic = biallelic

        if not math.isclose(sum(dosage_dist), 1.0):
            raise ValueError("Dosage distribution should sum to one. Sum is "
                             "now {:.2f}, dosage dist: {}".format(
                                sum(dosage_dist), dosage_dist))

        self.dosage_dist = dosage_dist
        self.ld_reassign_prob = ld_reassign_prob

        if substitution_rates:
            self.substitution_rates = substitution_rates
        else:
            self.substitution_rates = {}
            alphabet = set(INSERTION_ALPHABET)
            for base in alphabet:
                self.substitution_rates[base[0]] = {
                    b[0]: 1/3 for b in (alphabet - {base})
                }

        if insertion_rates:
            self.insertion_rates = insertion_rates
        else:
            self.insertion_rates = [0.25] * 4  # All bases equally likely

    def _gen_positions(self, distance_model: DistanceModelT):
        pos = 0
        while pos < len(self.reference):
            distance = distance_model()
            pos += distance

            if pos >= len(self.reference):
                break

            dosage_max = min(self.ploidy, len(self.dosage_dist))
            dosage = nprnd.choice(dosage_max, p=self.dosage_dist[:dosage_max])
            dosage += 1

            yield pos, dosage

    def generate(self, substitution_dist_model: DistanceModelT=None,
                 insert_dist_model: DistanceModelT=None,
                 deletion_dist_model: DistanceModelT=None) -> HaplotypeReturnT:

        # Merge all possible mutations by position
        all_mutations = defaultdict(list)

        if substitution_dist_model:
            for pos, dosage in self._gen_positions(substitution_dist_model):
                all_mutations[pos].append((MutationType.SUBSTITUTION, dosage))

        if insert_dist_model:
            for pos, dosage in self._gen_positions(insert_dist_model):
                all_mutations[pos].append((MutationType.INSERTION, dosage))

        if deletion_dist_model:
            for pos, dosage in self._gen_positions(deletion_dist_model):
                all_mutations[pos].append((MutationType.DELETION, dosage))

        # Nothing to mutate
        if len(all_mutations) == 0:
            return [self.reference] * self.ploidy, MutationRecord(self.ploidy)

        # Generate the haplotypes, if there are multiple mutations possible,
        # randomly pick one
        haplotypes = []
        mutation_record = MutationRecord(self.ploidy)
        for i in range(self.ploidy):
            haplotypes.append(bytearray())

        previous_pos = 0
        previous_alt_assignment = None
        for pos in sorted(all_mutations):
            mut_type, dosage = random.choice(all_mutations[pos])

            # Determine which haplotypes get the alternative mutation
            # ---------
            haplo_indices = set(range(self.ploidy))

            # With a given probability, we just take the assignment of the
            # previous mutation position, to account for linkage disequilibrium
            weights = [1-self.ld_reassign_prob, self.ld_reassign_prob]
            assign_to_prev = nprnd.choice(2, p=weights) == 1
            if previous_alt_assignment and assign_to_prev:
                alt_assignment = previous_alt_assignment

                # But we do need to take into account that the dosage may be
                # different. If the previous dosage is lower than the current
                # dosage, then we keep the previous one. Effectively, this
                # means that the dosage = min(previous, current)
                if len(alt_assignment) > dosage:
                    alt_assignment = set(random.sample(alt_assignment, dosage))
            else:
                alt_assignment = set(random.sample(haplo_indices, dosage))

            # Extend haplotypes up until pos
            for i in haplo_indices:
                haplotypes[i].extend(self.reference[previous_pos:pos])

            # Get the possibly additional bases corresponding to the current
            # mutation type
            extra_bases = self._get_bases_for_mutation(pos, alt_assignment,
                                                       mut_type)
            for i in range(self.ploidy):
                mutation_record.record(i, pos, len(haplotypes[i]),
                                       alt_assignment, extra_bases[i],
                                       mut_type)

                haplotypes[i].extend(extra_bases[i])

            previous_alt_assignment = alt_assignment

            # Reference base at `pos` already incorperated through
            # `_get_bases_for_mutation`: therefore, increase pos with 1.
            previous_pos = pos + 1

        return haplotypes, mutation_record

    def _get_bases_for_mutation(self, pos: int, alt_assignment: Set[int],
                                mut_type: MutationType):
        bases = []
        for i in range(self.ploidy):
            bases.append(bytearray())

        ref_base = self.reference[pos]

        if mut_type == MutationType.SUBSTITUTION:
            if self.biallelic:
                alt_base = self._get_substitution(ref_base)

                for i in range(self.ploidy):
                    if i in alt_assignment:
                        bases[i].append(alt_base)
                    else:
                        bases[i].append(ref_base)
            else:
                for i in range(self.ploidy):
                    if i in alt_assignment:
                        alt_base = self._get_substitution(ref_base)
                        bases[i].append(alt_base)
                    else:
                        bases[i].append(ref_base)
        elif mut_type == MutationType.INSERTION:
            if self.biallelic:
                ins_base = self._get_insertion()
                for i in alt_assignment:
                    bases[i].append(ins_base)
            else:
                for i in alt_assignment:
                    ins_base = self._get_insertion()
                    bases[i].append(ins_base)

            # It's an insertion, so add reference base to all haplotypes
            for i in range(self.ploidy):
                bases[i].append(ref_base)
        elif mut_type == MutationType.DELETION:
            for i in alt_assignment:
                bases[i].append(ref_base)

        return bases

    def _get_substitution(self, ref_base):
        possible_sub = list(self.substitution_rates[ref_base].keys())
        weights = list(self.substitution_rates[ref_base].values())

        subst_index = nprnd.choice(len(possible_sub), p=weights)
        return possible_sub[subst_index]

    def _get_insertion(self):
        index = nprnd.choice(len(INSERTION_ALPHABET), p=self.insertion_rates)

        return INSERTION_ALPHABET[index]


DISTANCE_MODELS = {
    'lognormal': LogNormalDistanceModel
}
