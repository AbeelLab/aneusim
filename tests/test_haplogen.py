import pytest  # noqa

from aneusim.haplogen import MutationGenerator


@pytest.fixture
def reference():
    return b"GATCATCGTAGCTACATCGTACGTAGCTACGTAC"


@pytest.fixture
def ploidy():
    return 3


@pytest.fixture
def dosage_dist():
    return [1/2, 2/6, 1/6]


def test_substitutions(reference, ploidy, dosage_dist):
    mut_gen = MutationGenerator(reference, ploidy, dosage_dist, biallelic=True)

    haplotypes, mut_records = mut_gen.generate(
        substitution_dist_model=lambda: 5)

    for haplotype, mut_record in zip(haplotypes, mut_records.mutations):
        assert len(haplotype) == len(reference)

        for record in mut_record:
            assert record['ref_pos'] == record['haplo_pos']
            pos = record['haplo_pos']

            assert haplotype[pos] != reference[pos]


def test_deletions(reference, ploidy, dosage_dist):
    mut_gen = MutationGenerator(reference, ploidy, dosage_dist, biallelic=True)

    haplotypes, mut_records = mut_gen.generate(
        deletion_dist_model=lambda: 5, deletion_size_model=lambda: 1)

    for haplotype, mut_record in zip(haplotypes, mut_records.mutations):
        total_deleted = -sum(r['size'] for r in mut_record)

        assert len(reference) == len(haplotype) + total_deleted

        for record in mut_record:
            ref_pos = record['ref_pos']
            hap_pos = record['haplo_pos']

            assert reference[ref_pos-1] == haplotype[hap_pos-1]
            assert reference[ref_pos+1] == haplotype[hap_pos]


def test_insertions(reference, ploidy, dosage_dist):
    mut_gen = MutationGenerator(reference, ploidy, dosage_dist, biallelic=True)

    haplotypes, mut_records = mut_gen.generate(
        insert_dist_model=lambda: 5)

    for haplotype, mut_record in zip(haplotypes, mut_records.mutations):
        total_inserted = sum(r['size'] for r in mut_record)
        assert len(reference) + total_inserted == len(haplotype)

        for record in mut_record:
            ref_pos = record['ref_pos']
            hap_pos = record['haplo_pos']

            assert reference[ref_pos-1] == haplotype[hap_pos-1]
            assert reference[ref_pos] == haplotype[hap_pos+1]
