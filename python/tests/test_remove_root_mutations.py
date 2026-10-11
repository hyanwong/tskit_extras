"""Tests for removing mutations above local roots."""

import pytest
import tskit

import tskit_extras as tsx


def sample_states(ts):
    return [
        tuple(
            None if genotype == tskit.MISSING_DATA else variant.alleles[genotype]
            for genotype in variant.genotypes
        )
        for variant in ts.variants()
    ]


def single_root_ts(root_mutations=True, descendant_mutation=True):
    tables = tskit.TableCollection(sequence_length=1)
    sample_a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    sample_b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    root = tables.nodes.add_row(time=1)
    tables.edges.add_row(0, 1, root, sample_a)
    tables.edges.add_row(0, 1, root, sample_b)
    site = tables.sites.add_row(
        position=0.5, ancestral_state="A", metadata=b"site"
    )
    parent = tskit.NULL
    if root_mutations:
        parent = tables.mutations.add_row(
            site, root, "G", metadata=b"root-1"
        )
        parent = tables.mutations.add_row(
            site, root, "C", parent=parent, metadata=b"root-2"
        )
    if descendant_mutation:
        tables.mutations.add_row(
            site, sample_a, "T", parent=parent, metadata=b"descendant"
        )
    tables.sort()
    return tables.tree_sequence(), sample_a


def multi_root_ts(mutate_both_roots):
    tables = tskit.TableCollection(sequence_length=1)
    sample_a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    sample_b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    root_a = tables.nodes.add_row(time=1)
    root_b = tables.nodes.add_row(time=1)
    tables.edges.add_row(0, 1, root_a, sample_a)
    tables.edges.add_row(0, 1, root_b, sample_b)
    site = tables.sites.add_row(position=0.5, ancestral_state="A")
    tables.mutations.add_row(site, root_a, "G")
    if mutate_both_roots:
        tables.mutations.add_row(site, root_b, "G")
    tables.sort()
    return tables.tree_sequence()


class TestRemoveRootMutations:
    def test_replaces_root_mutations_with_ancestral_state(self):
        ts, sample_a = single_root_ts()
        result = tsx.remove_root_mutations(ts)

        assert sample_states(result) == sample_states(ts)
        assert ts.site(0).ancestral_state == "A"
        assert ts.num_mutations == 3
        assert result.site(0).ancestral_state == "C"
        assert result.site(0).metadata == b"site"
        assert result.num_mutations == 1
        mutation = result.mutation(0)
        assert mutation.node == sample_a
        assert mutation.parent == tskit.NULL
        assert mutation.derived_state == "T"
        assert mutation.metadata == b"descendant"

    def test_removes_mutation_exposed_as_a_root_by_simplification(self):
        tables = tskit.TableCollection(sequence_length=1)
        sample_a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        sample_b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        outgroup = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        clade = tables.nodes.add_row(time=1)
        root = tables.nodes.add_row(time=2)
        tables.edges.add_row(0, 1, clade, sample_a)
        tables.edges.add_row(0, 1, clade, sample_b)
        tables.edges.add_row(0, 1, root, clade)
        tables.edges.add_row(0, 1, root, outgroup)
        site = tables.sites.add_row(position=0.5, ancestral_state="A")
        tables.mutations.add_row(site, clade, "G")
        tables.sort()
        ts = tables.tree_sequence()
        simplified = ts.simplify(samples=[sample_a, sample_b])
        root_mutation = simplified.mutation(0)
        assert root_mutation.node in simplified.first().roots
        assert root_mutation.edge == tskit.NULL

        result = tsx.remove_root_mutations(simplified)

        assert sample_states(result) == sample_states(simplified)
        assert result.site(0).ancestral_state == "G"
        assert result.num_mutations == 0

    def test_preserves_sites_without_root_mutations(self):
        ts, _ = single_root_ts(root_mutations=False)
        result = tsx.remove_root_mutations(ts)

        assert sample_states(result) == sample_states(ts)
        assert result.site(0).ancestral_state == "A"
        assert result.num_mutations == 1
        assert result.mutation(0).derived_state == "T"

    def test_removes_mutations_on_matching_multiple_roots(self):
        ts = multi_root_ts(mutate_both_roots=True)
        result = tsx.remove_root_mutations(ts)

        assert sample_states(result) == sample_states(ts)
        assert result.site(0).ancestral_state == "G"
        assert result.num_mutations == 0

    def test_raises_when_roots_have_different_inherited_states(self):
        ts = multi_root_ts(mutate_both_roots=False)
        original_states = sample_states(ts)

        with pytest.raises(ValueError, match="Multiple roots"):
            tsx.remove_root_mutations(ts)

        assert sample_states(ts) == original_states
        assert ts.num_mutations == 1

    def test_processes_sites_in_each_local_tree(self):
        tables = tskit.TableCollection(sequence_length=2)
        sample_a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        sample_b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        left_root = tables.nodes.add_row(time=1)
        right_root = tables.nodes.add_row(time=2)
        for sample in (sample_a, sample_b):
            tables.edges.add_row(0, 1, left_root, sample)
            tables.edges.add_row(1, 2, right_root, sample)
        left_site = tables.sites.add_row(position=0.5, ancestral_state="A")
        right_site = tables.sites.add_row(position=1.5, ancestral_state="A")
        tables.mutations.add_row(left_site, left_root, "G")
        tables.mutations.add_row(right_site, right_root, "T")
        tables.sort()
        ts = tables.tree_sequence()

        result = tsx.remove_root_mutations(ts)

        assert sample_states(result) == sample_states(ts)
        assert [site.ancestral_state for site in result.sites()] == ["G", "T"]
        assert result.num_mutations == 0

    def test_keeps_mutations_on_isolated_roots(self):
        tables = tskit.TableCollection(sequence_length=1)
        sample = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        site = tables.sites.add_row(position=0.5, ancestral_state="A")
        tables.mutations.add_row(site, sample, "G")
        ts = tables.tree_sequence()

        result = tsx.remove_root_mutations(ts)

        assert sample_states(result) == sample_states(ts)
        assert result.site(0).ancestral_state == "A"
        assert result.num_mutations == 1

    def test_handles_tree_sequences_without_sites(self):
        ts = tskit.TableCollection(sequence_length=1).tree_sequence()

        result = tsx.remove_root_mutations(ts)

        assert result.num_sites == 0
        assert result.num_mutations == 0

    def test_top_level_alias_matches_tree_sequence_method(self):
        ts, _ = single_root_ts()

        top_level = tsx.remove_root_mutations(ts)
        method = tsx.TreeSequence.remove_root_mutations(ts)

        assert top_level.tables.equals(method.tables, ignore_provenance=True)
