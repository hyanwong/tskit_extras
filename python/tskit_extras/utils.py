"""Standalone utilities that do not operate on a specific tskit object."""

import numpy as np
import tskit


def invert_map(node_mapping):
    """
    Invert a node map returned by ``TreeSequence.simplify(map_nodes=True)``.

    Non-NULL entries in a simplify node map are consecutive IDs starting at
    zero. The returned array is indexed by simplified node ID and contains the
    corresponding original node ID.

    :param array_like node_mapping: The node map returned by simplifying a
        tree sequence with ``map_nodes=True``.
    :return: The reverse node map.
    :rtype: numpy.ndarray

    See https://github.com/hyanwong/tskit_extras/issues/3
    """
    node_mapping = tskit.util.safe_np_int_cast(node_mapping, np.int32)
    if node_mapping.ndim != 1:
        raise ValueError("node_mapping must be one-dimensional")

    kept = node_mapping != tskit.NULL
    new_ids = node_mapping[kept]
    num_kept = len(new_ids)
    if np.any(new_ids < 0) or np.any(new_ids >= num_kept):
        raise ValueError(
            "non-NULL node IDs in node_mapping must be consecutive from zero"
        )

    reverse_map = np.full(num_kept, tskit.NULL, dtype=node_mapping.dtype)
    reverse_map[new_ids] = np.flatnonzero(kept)
    if np.any(reverse_map == tskit.NULL):
        raise ValueError(
            "non-NULL node IDs in node_mapping must be unique and consecutive"
        )
    return reverse_map
