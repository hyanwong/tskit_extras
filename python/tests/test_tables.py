"""Tests for table utilities."""

import numpy as np
import tskit

import tskit_extras as tsx


class TestTableCollection:
    def test_reorder_nodes(self):
        ts = tskit.Tree.generate_balanced(8).tree_sequence
        order = -ts.nodes_time
        tables = ts.dump_tables()
        tsx.TableCollection.reorder_nodes(tables, order)
        tables.sort()
        reordered = tables.tree_sequence()
        expected = tsx.reorder_nodes(ts, order)
        np.testing.assert_array_equal(
            reordered.tables.nodes.time, expected.tables.nodes.time
        )


class TestEdgeTable:
    """Tests for unsquash (discussion #2657)."""

    def test_unsquash_single_edge(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=100, parent=1, child=0)
        tsx.EdgeTable.unsquash(tables.edges, positions=[25, 50, 75])
        assert tables.edges.num_rows == 4
        np.testing.assert_array_equal(tables.edges.left, [0, 25, 50, 75])
        np.testing.assert_array_equal(tables.edges.right, [25, 50, 75, 100])

    def test_unsquash_no_split_needed(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=50, parent=1, child=0)
        tsx.EdgeTable.unsquash(tables.edges, positions=[60, 80])
        assert tables.edges.num_rows == 1

    def test_unsquash_subset_of_edges(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=100, parent=2, child=0)
        tables.edges.add_row(left=0, right=100, parent=2, child=1)
        tsx.EdgeTable.unsquash(tables.edges, positions=[50], edges=[0])
        assert tables.edges.num_rows == 3

    def test_unsquash_deduplicates_positions(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=100, parent=1, child=0)
        tsx.EdgeTable.unsquash(tables.edges, positions=[50, 50, 50])
        assert tables.edges.num_rows == 2

    def test_edge_table_namespace(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()
        tables.nodes.add_row(time=0, population=0, flags=tskit.NODE_IS_SAMPLE)
        tables.nodes.add_row(time=1, population=0)
        tables.edges.add_row(left=0, right=100, parent=1, child=0)
        tsx.EdgeTable.unsquash(tables.edges, positions=[50])
        assert tables.edges.num_rows == 2
