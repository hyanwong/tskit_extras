"""Simulation-based integration tests for OTR parsimony."""

import itertools

import msprime
import numpy as np
import pytest
import tskit

import tskit_extras as tsx


def binary_resolution_score(tree, genotypes, alleles, ancestral_state=None):
    return min(
        len(
            resolution.map_mutations(
                genotypes, alleles, ancestral_state=ancestral_state
            )[1]
        )
        for resolution in tsx.Tree.all_binary_resolutions(
            tree, record_provenance=False
        )
    )


class TestOptimalTreeRefinementParsimonyScore:
    """Tests for parsimony scores optimized over binary refinements."""

    def test_binary_states_match_all_resolutions(self):
        tree = tskit.Tree.generate_star(4)
        alleles = ("A", "G")
        for genotypes in itertools.product(range(2), repeat=4):
            for ancestral_state in (None, "A", "G"):
                assert tsx.Tree.optimal_tree_refinement_parsimony_score(
                    tree, genotypes, alleles, ancestral_state=ancestral_state
                ) == binary_resolution_score(
                    tree, genotypes, alleles, ancestral_state=ancestral_state
                )

    @pytest.mark.parametrize(
        ("genotypes", "ancestral_state"),
        [
            ((0, 1, 2, 0), None),
            ((0, 1, 2, 0), 0),
            ((0, 1, 2, 0), 1),
            ((2, 2, 1, 0), 2),
        ],
    )
    def test_multistate_matches_all_resolutions(self, genotypes, ancestral_state):
        tree = tskit.Tree.generate_star(4)
        alleles = ("A", "C", "G")
        assert tsx.Tree.optimal_tree_refinement_parsimony_score(
            tree, genotypes, alleles, ancestral_state=ancestral_state
        ) == binary_resolution_score(
            tree, genotypes, alleles, ancestral_state=ancestral_state
        )

    def test_multiple_polytomies_match_all_resolutions(self):
        tree = tskit.Tree.generate_balanced(8, arity=3)
        genotypes = [0, 1, 2, 0, 1, 2, 0, 1]
        alleles = ("A", "C", "G")
        for ancestral_state in (None, 0, 1, 2):
            assert tsx.Tree.optimal_tree_refinement_parsimony_score(
                tree, genotypes, alleles, ancestral_state=ancestral_state
            ) == binary_resolution_score(
                tree, genotypes, alleles, ancestral_state=ancestral_state
            )

    def test_internal_samples_and_missing_data(self):
        tables = tskit.TableCollection(sequence_length=1)
        leaves = [
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
            for _ in range(4)
        ]
        internal_sample = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=1)
        root = tables.nodes.add_row(time=2)
        for child in leaves[:2]:
            tables.edges.add_row(0, 1, internal_sample, child)
        for child in (internal_sample, *leaves[2:]):
            tables.edges.add_row(0, 1, root, child)
        tables.sort()
        tree = tables.tree_sequence().first()
        sample_states = {
            leaves[0]: 0,
            leaves[1]: 1,
            leaves[2]: tskit.MISSING_DATA,
            leaves[3]: 0,
            internal_sample: 1,
        }
        genotypes = [sample_states[u] for u in tree.tree_sequence.samples()]
        alleles = ("A", "G")

        assert tsx.Tree.optimal_tree_refinement_parsimony_score(
            tree, genotypes, alleles, ancestral_state="A"
        ) == binary_resolution_score(
            tree, genotypes, alleles, ancestral_state="A"
        )

    def test_multiple_roots_and_fixed_ancestral_state(self):
        tables = tskit.TableCollection(sequence_length=1)
        for _ in range(2):
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        tree = tables.tree_sequence().first()
        genotypes = [1, 1]
        alleles = ("A", "G")

        assert tsx.Tree.optimal_tree_refinement_parsimony_score(
            tree, genotypes, alleles, ancestral_state=0
        ) == binary_resolution_score(
            tree, genotypes, alleles, ancestral_state=0
        )

    def test_matches_map_mutations_on_binary_tree(self):
        tree = tskit.Tree.generate_balanced(8, arity=2)
        genotypes = [0, 1, 1, 0, 1, 0, 0, 1]
        alleles = ("A", "G")
        for ancestral_state in (None, 0, 1):
            assert tsx.Tree.optimal_tree_refinement_parsimony_score(
                tree, genotypes, alleles, ancestral_state=ancestral_state
            ) == len(
                tree.map_mutations(
                    genotypes, alleles, ancestral_state=ancestral_state
                )[1]
            )

    def test_rejects_invalid_genotypes(self):
        tree = tskit.Tree.generate_star(3)
        with pytest.raises(ValueError, match="one entry for each"):
            tsx.Tree.optimal_tree_refinement_parsimony_score(
                tree, [0, 1], ("A", "G")
            )
        with pytest.raises(ValueError, match="non-missing"):
            tsx.Tree.optimal_tree_refinement_parsimony_score(
                tree, [tskit.MISSING_DATA] * 3, ("A", "G")
            )
        with pytest.raises(ValueError, match="valid allele indexes"):
            tsx.Tree.optimal_tree_refinement_parsimony_score(
                tree, [0, 1, 2], ("A", "G")
            )


@pytest.mark.parametrize("removal_fraction", [0.25, 0.5, 0.75])
def test_otr_scores_removed_infinite_sites_mutations(removal_fraction):
    """Check that OTR recovers a single change after unsupported edges are removed.

    Genotypes are saved from the original sequence before masking mutation rows:
    the mutated sequence no longer contains the derived states for those sites.
    Removing unsupported edges then creates polytomies that can make ordinary
    ``map_mutations`` require multiple changes, while OTR should still score one.
    """
    ancestry = msprime.sim_ancestry(
        samples=12,
        sequence_length=100,
        recombination_rate=0.02,
        random_seed=1,
    )
    ts = msprime.sim_mutations(
        ancestry,
        rate=0.03,
        model=msprime.InfiniteSites(),
        random_seed=101,
        discrete_genome=False,
    )

    # Keep only informative, single-mutation sites whose mutation lies on an edge.
    # Copy the observations now, before removing any mutation records.
    focal_mutations = []
    for variant in ts.variants():
        site = variant.site
        if (
            len(site.mutations) == 1
            and site.mutations[0].edge != tskit.NULL
            and len(variant.alleles) == 2
            and len(np.unique(variant.genotypes)) == 2
        ):
            focal_mutations.append(
                (
                    site.position,
                    variant.genotypes.copy(),
                    variant.alleles,
                    site.ancestral_state,
                    site.mutations[0].id,
                )
            )

    keep_mutations = (
        np.random.default_rng(23).random(ts.num_mutations) >= removal_fraction
    )
    tables = ts.dump_tables()
    tables.mutations.keep_rows(keep_mutations)
    tables.sort()
    mutation_pruned_ts = tables.tree_sequence()
    polytomy_ts = tsx.remove_edges(
        mutation_pruned_ts, tsx.unsupported_edges(mutation_pruned_ts)
    )

    num_scored = 0
    num_nonparsimonious = 0
    for position, genotypes, alleles, ancestral_state, mutation_id in focal_mutations:
        if keep_mutations[mutation_id]:
            continue
        tree = polytomy_ts.at(position)
        hartigan_score = len(
            tree.map_mutations(
                genotypes, alleles, ancestral_state=ancestral_state
            )[1]
        )
        otr_score = tsx.Tree.optimal_tree_refinement_parsimony_score(
            tree, genotypes, alleles, ancestral_state=ancestral_state
        )
        assert otr_score == 1
        num_scored += 1
        num_nonparsimonious += hartigan_score > 1

    assert num_scored > 0
    assert num_nonparsimonious > 0
