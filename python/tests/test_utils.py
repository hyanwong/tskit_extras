"""Tests for node-map utilities."""

import msprime
import numpy as np
import pytest
import tskit

import tskit_extras as tsx


class TestInvertMap:
    """Tests for inverting node maps returned by simplify."""

    def test_inverts_simplify_node_map(self):
        ts = msprime.sim_ancestry(samples=8, random_seed=42)
        simplified_ts, node_mapping = ts.simplify(
            samples=ts.samples()[:4], map_nodes=True
        )
        reverse_map = tsx.utils.invert_map(node_mapping)

        np.testing.assert_array_equal(
            node_mapping[reverse_map], np.arange(simplified_ts.num_nodes)
        )

    def test_null_entries_and_empty_reverse_map(self):
        node_mapping = np.array([1, tskit.NULL, 0], dtype=np.int32)
        np.testing.assert_array_equal(
            tsx.utils.invert_map(node_mapping), [2, 0]
        )

        all_null = np.array([tskit.NULL, tskit.NULL], dtype=np.int32)
        assert tsx.utils.invert_map(all_null).size == 0

    def test_rejects_malformed_node_maps(self):
        with pytest.raises(ValueError, match="one-dimensional"):
            tsx.utils.invert_map(np.zeros((3, 1), dtype=np.int32))

        sparse_map = np.array([2, tskit.NULL], dtype=np.int32)
        with pytest.raises(ValueError, match="consecutive from zero"):
            tsx.utils.invert_map(sparse_map)

        duplicate_map = np.array([0, 0], dtype=np.int32)
        with pytest.raises(ValueError, match="unique and consecutive"):
            tsx.utils.invert_map(duplicate_map)
