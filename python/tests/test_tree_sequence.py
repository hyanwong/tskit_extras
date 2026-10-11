"""Tests for TreeSequence utilities."""

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
        assert reordered.node(0).time >= reordered.node(1).time
        assert reordered.node(0).time >= reordered.node(ts.num_nodes - 1).time

    def test_reorder_preserves_topology(self):
        ts = get_simple_ts()
        order = -ts.nodes_time
        reordered = tsx.reorder_nodes(ts, order=order)
        assert reordered.num_nodes == ts.num_nodes
        assert reordered.num_edges == ts.num_edges
        assert reordered.num_trees == ts.num_trees

    def test_reorder_preserves_genotypes(self):
        ts = get_ts_with_mutations()
        order = -ts.nodes_time
        reordered = tsx.reorder_nodes(ts, order=order)
        g1 = np.sort(ts.genotype_matrix(), axis=1)
        g2 = np.sort(reordered.genotype_matrix(), axis=1)
        np.testing.assert_array_equal(g1, g2)

    def test_reorder_stable_sort(self):
        ts = get_simple_ts()
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

    def test_top_level_alias(self):
        ts = get_simple_ts()
        order = -ts.nodes_time
        r1 = tsx.reorder_nodes(ts, order)
        r2 = tsx.TreeSequence.reorder_nodes(ts, order)
        np.testing.assert_array_equal(r1.tables.nodes.time, r2.tables.nodes.time)


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
        return tsx.unsupported_edges(ts, per_interval=per_interval)

    def test_raises_on_mutation_on_removed_edge(self):
        ts = get_ts_with_mutations()
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
            np.testing.assert_array_equal(
                ts.genotype_matrix(), poly_ts.genotype_matrix()
            )
            for t1, t2 in zip(ts.trees(), poly_ts.trees()):
                assert set(t1.roots) == set(t2.roots)
                for node in t2.nodes():
                    assert set(t1.samples(node)) == set(t2.samples(node))

    def test_removing_edge_above_sample_can_change_genotypes(self):
        tables = tskit.TableCollection(sequence_length=1)
        sample_a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        sample_b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        parent = tables.nodes.add_row(time=1)
        root = tables.nodes.add_row(time=2)
        tables.edges.add_row(0, 1, parent, sample_a)
        tables.edges.add_row(0, 1, parent, sample_b)
        tables.edges.add_row(0, 1, root, parent)
        site = tables.sites.add_row(position=0.5, ancestral_state="A")
        tables.mutations.add_row(site=site, node=parent, derived_state="G")
        tables.sort()
        ts = tables.tree_sequence()
        edge_to_sample = next(
            edge.id for edge in ts.edges() if edge.child == sample_a
        )

        assert edge_to_sample not in tsx.unsupported_edges(ts)
        assert ts.genotype_matrix().tolist() == [[1, 1]]
        collapsed = tsx.remove_edges(ts, [edge_to_sample])

        assert collapsed.genotype_matrix().tolist() == [[tskit.MISSING_DATA, 1]]

    def test_remove_no_edges(self):
        ts = get_ts_with_mutations()
        result = tsx.remove_edges(ts, [])
        assert result.num_edges == ts.num_edges

    def test_top_level_alias(self):
        ts = get_ts_with_mutations()
        result = tsx.remove_edges(ts, [])
        result2 = tsx.TreeSequence.remove_edges(ts, [])
        assert result.num_edges == result2.num_edges

    def test_remove_stacked_edges_over_partial_intervals(self):
        tables = tskit.TableCollection(sequence_length=10)
        a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        child = tables.nodes.add_row(time=1)
        parent = tables.nodes.add_row(time=2)
        root = tables.nodes.add_row(time=3)
        for sample in (a, b):
            tables.edges.add_row(0, 10, child, sample)
        tables.edges.add_row(0, 10, parent, child)
        tables.edges.add_row(0, 5, root, parent)
        tables.edges.add_row(5, 10, root, parent)
        tables.sort()
        ts = tables.tree_sequence()

        reverse_order = np.arange(ts.num_nodes)[::-1]
        inverse_node_map = np.empty(ts.num_nodes, dtype=int)
        inverse_node_map[np.argsort(reverse_order)] = np.arange(ts.num_nodes)
        for input_ts, node_map in (
            (ts, np.arange(ts.num_nodes)),
            (tsx.reorder_nodes(ts, reverse_order), inverse_node_map),
        ):
            remove = [
                edge.id
                for edge in input_ts.edges()
                if (edge.parent, edge.child)
                == (node_map[parent], node_map[child])
                or (
                    (edge.parent, edge.child)
                    == (node_map[root], node_map[parent])
                    and edge.left == 0
                )
            ]
            result = tsx.remove_edges(input_ts, remove)
            expected = {
                (node_map[root], node_map[a], 0, 5),
                (node_map[root], node_map[b], 0, 5),
                (node_map[root], node_map[parent], 5, 10),
                (node_map[parent], node_map[a], 5, 10),
                (node_map[parent], node_map[b], 5, 10),
            }
            assert {
                (e.parent, e.child, e.left, e.right) for e in result.edges()
            } == expected
            for position in (2.5, 7.5):
                tree = result.at(position)
                assert set(tree.samples(tree.root)) == {
                    node_map[a],
                    node_map[b],
                }

    def test_remove_edges_empty_and_validate_edge_ids(self):
        tables = tskit.TableCollection(sequence_length=1)
        ts = tables.tree_sequence()
        assert tsx.remove_edges(ts, []).num_edges == 0
        assert len(tsx.unsupported_edges(ts)) == 0

        ts = get_simple_ts()
        with pytest.raises(ValueError, match="out of bounds"):
            tsx.remove_edges(ts, [-1])
        with pytest.raises(ValueError, match="out of bounds"):
            tsx.remove_edges(ts, [ts.num_edges])
        with pytest.raises(TypeError, match="must be integers"):
            tsx.remove_edges(ts, [1.5])

    def test_unsupported_edges_per_interval_or_relationship(self):
        tables = tskit.TableCollection(sequence_length=10)
        sample_a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        sample_b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        parent = tables.nodes.add_row(time=1)
        root = tables.nodes.add_row(time=2)
        tables.edges.add_row(0, 5, parent, sample_a)
        tables.edges.add_row(5, 10, parent, sample_a)
        tables.edges.add_row(0, 10, parent, sample_b)
        tables.edges.add_row(0, 5, root, parent)
        tables.edges.add_row(5, 10, root, parent)
        site = tables.sites.add_row(position=2, ancestral_state="A")
        tables.mutations.add_row(site=site, node=sample_a, derived_state="G")
        tables.sort()
        ts = tables.tree_sequence()

        by_interval = set(tsx.unsupported_edges(ts, per_interval=True))
        by_relationship = set(tsx.unsupported_edges(ts))
        assert len(by_interval) == 2
        assert len(by_relationship) == 2
        assert {
            (ts.edge(edge_id).parent, ts.edge(edge_id).child)
            for edge_id in by_relationship
        } == {(root, parent)}
        assert {
            (ts.edge(edge_id).parent, ts.edge(edge_id).child,
             ts.edge(edge_id).left, ts.edge(edge_id).right)
            for edge_id in by_interval
        } == {
            (root, parent, 0, 5),
            (root, parent, 5, 10),
        }

    def test_unsupported_edges_alias(self):
        ts = get_ts_with_mutations()
        np.testing.assert_array_equal(
            tsx.unsupported_edges(ts),
            tsx.TreeSequence.unsupported_edges(ts),
        )
