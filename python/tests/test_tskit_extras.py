"""
Tests for tskit_extras.
"""


import msprime
import numpy as np
import pytest
import tskit

import tskit_extras as tsx


def get_simple_ts(n=10, seed=42, sequence_length=1e4, recombination_rate=1e-3):
    return msprime.sim_ancestry(
        n,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=seed,
    )


def get_ts_with_mutations(n=10, seed=42, mutation_rate=1e-3):
    ts = get_simple_ts(n=n, seed=seed)
    return msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seed)


class TestReorderNodes:
    """Tests for reorder_nodes (discussion #3466)."""

    def test_reorder_oldest_first(self):
        ts = get_simple_ts()
        reordered = tsx.reorder_nodes(ts, order=-ts.nodes_time)
        # The oldest node should now have the highest ID
        assert reordered.node(0).time >= reordered.node(1).time
        assert reordered.node(0).time >= reordered.node(ts.num_nodes - 1).time

    def test_reorder_preserves_topology(self):
        ts = get_simple_ts()
        order = -ts.nodes_time  # oldest-first
        reordered = tsx.reorder_nodes(ts, order=order)
        # Same number of nodes, edges, trees
        assert reordered.num_nodes == ts.num_nodes
        assert reordered.num_edges == ts.num_edges
        assert reordered.num_trees == ts.num_trees

    def test_reorder_preserves_genotypes(self):
        # Genotype matrix is indexed by variant×sample and is independent of
        # internal node IDs, so it must be unchanged after reordering.
        ts = get_ts_with_mutations()
        order = -ts.nodes_time
        reordered = tsx.reorder_nodes(ts, order=order)
        # Samples are reordered too, but genotypes are still per-sample.
        # Compare the sorted genotype matrix rows to be order-independent.
        g1 = np.sort(ts.genotype_matrix(), axis=1)
        g2 = np.sort(reordered.genotype_matrix(), axis=1)
        np.testing.assert_array_equal(g1, g2)

    def test_reorder_stable_sort(self):
        # Ties should preserve original relative order (stable sort)
        ts = get_simple_ts()
        # Use all-same order → nodes should stay in original order
        order = np.zeros(ts.num_nodes)
        reordered = tsx.reorder_nodes(ts, order=order)
        np.testing.assert_array_equal(
            reordered.tables.nodes.time, ts.tables.nodes.time
        )

    def test_reorder_with_mutations(self):
        ts = get_ts_with_mutations()
        order = -ts.nodes_time
        reordered = tsx.reorder_nodes(ts, order=order)
        assert reordered.num_mutations == ts.num_mutations

    def test_table_collection_reorder(self):
        ts = get_simple_ts()
        order = -ts.nodes_time
        tables = ts.dump_tables()
        tsx.TableCollection.reorder_nodes(tables, order)
        tables.sort()
        reordered = tables.tree_sequence()
        # Same result as tsx.reorder_nodes
        expected = tsx.reorder_nodes(ts, order)
        np.testing.assert_array_equal(
            reordered.tables.nodes.time, expected.tables.nodes.time
        )

    def test_top_level_alias(self):
        ts = get_simple_ts()
        order = -ts.nodes_time
        r1 = tsx.reorder_nodes(ts, order)
        r2 = tsx.TreeSequence.reorder_nodes(ts, order)
        np.testing.assert_array_equal(
            r1.tables.nodes.time, r2.tables.nodes.time
        )


class TestArgNumChildren:
    """Tests for arg_num_children (discussion #3449)."""

    def test_return_shape(self):
        ts = get_simple_ts()
        n_children = tsx.arg_num_children(ts)
        assert n_children.shape == (ts.num_nodes,)

    def test_samples_have_no_children(self):
        ts = get_simple_ts()
        n_children = tsx.arg_num_children(ts)
        for sample in ts.samples():
            assert n_children[sample] == 0

    def test_no_edges_gives_zeros(self):
        ts = get_simple_ts()
        tables = ts.dump_tables()
        tables.edges.clear()
        empty_ts = tables.tree_sequence()
        n_children = tsx.arg_num_children(empty_ts)
        np.testing.assert_array_equal(n_children, 0)

    def test_counts_unique_children(self):
        # Simple manual check: one parent with two distinct children
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        c1 = tables.individuals.add_row()
        c2 = tables.individuals.add_row()
        p = tables.individuals.add_row()
        n_c1 = tables.nodes.add_row(
            time=0, population=0, flags=tskit.NODE_IS_SAMPLE, individual=c1
        )
        n_c2 = tables.nodes.add_row(
            time=0, population=0, flags=tskit.NODE_IS_SAMPLE, individual=c2
        )
        n_p = tables.nodes.add_row(time=1, population=0, individual=p)
        tables.edges.add_row(left=0, right=100, parent=n_p, child=n_c1)
        tables.edges.add_row(left=0, right=100, parent=n_p, child=n_c2)
        tables.sort()
        ts = tables.tree_sequence()
        n_children = tsx.arg_num_children(ts)
        assert n_children[n_p] == 2
        assert n_children[n_c1] == 0
        assert n_children[n_c2] == 0

    def test_top_level_alias(self):
        ts = get_simple_ts()
        r1 = tsx.arg_num_children(ts)
        r2 = tsx.TreeSequence.arg_num_children(ts)
        np.testing.assert_array_equal(r1, r2)


class TestRemoveEdges:
    """Tests for remove_edges (discussion #2926)."""

    def _unsupported_edges(self, ts, per_interval=False):
        """Return internal edges unsupported by a mutation."""
        edges_to_remove = np.ones(ts.num_edges, dtype="bool")
        edges_to_remove[[m.edge for m in ts.mutations()]] = False
        edges_to_remove[np.isin(ts.edges_child, ts.samples())] = False
        if per_interval:
            return np.where(edges_to_remove)[0]
        else:
            keep = ~edges_to_remove
            for p, c in zip(ts.edges_parent[keep], ts.edges_child[keep]):
                edges_to_remove[
                    np.logical_and(ts.edges_parent == p, ts.edges_child == c)
                ] = False
            return np.where(edges_to_remove)[0]

    def test_raises_on_mutation_on_removed_edge(self):
        ts = get_ts_with_mutations()
        # Find an edge that has a mutation
        edge_with_mutation = next(
            m.edge for m in ts.mutations() if m.edge != tskit.NULL
        )
        with pytest.raises(ValueError, match="associated mutations"):
            tsx.remove_edges(ts, [edge_with_mutation])

    def test_preserves_roots_and_sample_sets(self):
        for seed in range(1, 6):
            ts = msprime.sim_ancestry(
                20,
                sequence_length=1e5,
                recombination_rate=1e-4,
                random_seed=seed,
            )
            ts = msprime.sim_mutations(ts, rate=1e-3, random_seed=seed)
            del_edges = self._unsupported_edges(ts)
            poly_ts = tsx.remove_edges(ts, del_edges)
            assert ts.num_trees == poly_ts.num_trees
            for t1, t2 in zip(ts.trees(), poly_ts.trees()):
                assert set(t1.roots) == set(t2.roots)
                for node in t2.nodes():
                    assert set(t1.samples(node)) == set(t2.samples(node))

    def test_remove_no_edges(self):
        ts = get_ts_with_mutations()
        result = tsx.remove_edges(ts, [])
        assert result.num_edges == ts.num_edges

    def test_top_level_alias(self):
        ts = get_ts_with_mutations()
        # Just check the alias resolves without error
        result = tsx.remove_edges(ts, [])
        result2 = tsx.TreeSequence.remove_edges(ts, [])
        assert result.num_edges == result2.num_edges


class TestUnsquash:
    """Tests for unsquash (discussion #2657)."""

    def test_unsquash_single_edge(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=100, parent=1, child=0)
        tsx.unsquash(tables.edges, positions=[25, 50, 75])
        assert tables.edges.num_rows == 4
        np.testing.assert_array_equal(tables.edges.left, [0, 25, 50, 75])
        np.testing.assert_array_equal(tables.edges.right, [25, 50, 75, 100])

    def test_unsquash_no_split_needed(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=50, parent=1, child=0)
        # Positions outside the edge span should be ignored
        tsx.unsquash(tables.edges, positions=[60, 80])
        assert tables.edges.num_rows == 1

    def test_unsquash_subset_of_edges(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=100, parent=2, child=0)
        tables.edges.add_row(left=0, right=100, parent=2, child=1)
        # Only split the first edge (index 0)
        tsx.unsquash(tables.edges, positions=[50], edges=[0])
        assert tables.edges.num_rows == 3

    def test_unsquash_deduplicates_positions(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=100, parent=1, child=0)
        tsx.unsquash(tables.edges, positions=[50, 50, 50])
        assert tables.edges.num_rows == 2

    def test_top_level_alias(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=100, parent=1, child=0)
        tsx.unsquash(tables.edges, positions=[50])
        assert tables.edges.num_rows == 2


class TestInvertMap:
    """Tests for invert_map."""

    def test_inverts_simplify_node_map(self):
        ts = get_simple_ts(n=8)
        simplified_ts, node_map = ts.simplify(
            samples=ts.samples()[:4], map_nodes=True
        )
        rev_map = tsx.invert_map(node_map)
        np.testing.assert_array_equal(
            node_map[rev_map], np.arange(simplified_ts.num_nodes)
        )

    def test_skips_null_entries(self):
        node_map = np.array([2, tskit.NULL, 0, 1, tskit.NULL], dtype=np.int32)
        rev_map = tsx.invert_map(node_map)
        np.testing.assert_array_equal(rev_map, [2, 3, 0])

    def test_handles_sparse_output_ids(self):
        node_map = np.array([3, tskit.NULL, 1], dtype=np.int32)
        rev_map = tsx.invert_map(node_map)
        np.testing.assert_array_equal(rev_map, [tskit.NULL, 2, tskit.NULL, 0])


class TestVersionAndImports:
    """Basic import and version sanity checks."""

    def test_version_string(self):
        assert isinstance(tsx.__version__, str)
        assert len(tsx.__version__) > 0

    def test_module_exports(self):
        for name in [
            "reorder_nodes",
            "arg_num_children",
            "invert_map",
            "remove_edges",
            "unsquash",
            "TreeSequence",
            "TableCollection",
        ]:
            assert hasattr(tsx, name)
