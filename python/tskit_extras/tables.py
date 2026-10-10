"""Utilities that operate on tskit tables."""

import itertools

import numpy as np
import tskit


class TableCollection:
    """
    Methods that operate on a :class:`tskit.TableCollection`.
    """

    @staticmethod
    def reorder_nodes(tables, order):
        """
        Reorder the nodes in a :class:`tskit.TableCollection`, keeping all other
        data identical.

        The passed-in ``order`` is an array of the same length as the number of nodes.
        Each node will be ordered based on its value in ``order``; nodes with the
        lowest value will receive the lowest ID. Order of ties is preserved (stable
        sort).

        :param tskit.TableCollection tables: The table collection to reorder in-place.
        :param array_like order: A 1-D array of numeric values used to determine the
            new ordering of nodes.

        See https://github.com/tskit-dev/tskit/discussions/3466
        """
        node_map = np.argsort(order, kind="stable").astype(
            tables.mutations.node.dtype
        )
        inverse_node_map = np.empty_like(node_map)
        inverse_node_map[node_map] = np.arange(len(node_map), dtype=node_map.dtype)
        tables.nodes.replace_with(tables.nodes[node_map])
        tables.edges.parent = inverse_node_map[tables.edges.parent]
        tables.edges.child = inverse_node_map[tables.edges.child]
        tables.mutations.node = inverse_node_map[tables.mutations.node]
        tables.migrations.node = inverse_node_map[tables.migrations.node]


class EdgeTable:
    """
    Methods that operate on a :class:`tskit.EdgeTable`.
    """

    @staticmethod
    def unsquash(edge_table, positions, edges=None):
        """
        Chop a set of edges into sub-edges at every position in ``positions``.

        For a set of positions and a given set of edges (or all edges if ``edges``
        is ``None``), create a new set of edges which cover the same span but which
        are divided into separate edges at every specified position. This is the
        opposite of :meth:`tskit.EdgeTable.squash`.

        :param tskit.EdgeTable edge_table: The edge table to modify in-place.
        :param array_like positions: Positions at which to split edges.
        :param list edges: Indices of edges to split (``None`` means all edges).

        See https://github.com/tskit-dev/tskit/discussions/2657
        """
        new_edges = tskit.EdgeTable()
        positions = np.unique(positions)  # sort and uniquify
        skip = []
        if edges is not None:
            skip = np.ones(edge_table.num_rows, dtype=bool)
            skip[edges] = False
        for edge, do_skip in itertools.zip_longest(
            edge_table, skip, fillvalue=False
        ):
            if do_skip:
                new_edges.append(edge)
                continue
            for left, right in itertools.pairwise(
                itertools.chain(
                    [edge.left],
                    positions[
                        np.logical_and(
                            positions > edge.left, positions < edge.right
                        )
                    ],
                    [edge.right],
                )
            ):
                new_edges.append(edge.replace(left=left, right=right))
        edge_table.replace_with(new_edges)
