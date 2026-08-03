"""
tskit_extras: Extra utility functions on top of tskit.

Provides tested routines considered too complex or not fully tested enough
for inclusion in core tskit. Functions can be accessed as:

    import tskit_extras as tsx

    tsx.reorder_nodes(ts, order)
    tsx.TreeSequence.reorder_nodes(ts, order)
    tsx.TableCollection.reorder_nodes(tables, order)
"""

import collections
import itertools

import numpy as np
import tskit

from . import _version

__version__ = _version.tskit_extras_version

__all__ = [
    "reorder_nodes",
    "arg_num_children",
    "unsquash",
    "remove_edges",
    "TableCollection",
    "TreeSequence",
]


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


class TreeSequence:
    """
    Methods that operate on a :class:`tskit.TreeSequence`.
    """

    @staticmethod
    def reorder_nodes(ts, order):
        """
        Return a new :class:`tskit.TreeSequence` in which the nodes are reordered
        according to ``order``.

        The passed-in ``order`` is an array of the same length as the number of nodes.
        Each node will be ordered based on its value in ``order``; nodes with the
        lowest value will receive the lowest ID. Order of ties is preserved (stable
        sort).

        :param tskit.TreeSequence ts: The tree sequence to reorder.
        :param array_like order: A 1-D array of numeric values used to determine the
            new ordering of nodes.
        :return: A new tree sequence with reordered nodes.
        :rtype: tskit.TreeSequence

        See https://github.com/tskit-dev/tskit/discussions/3466
        """
        tables = ts.dump_tables()
        TableCollection.reorder_nodes(tables, order)
        tables.sort()
        return tables.tree_sequence()

    @staticmethod
    def arg_num_children(ts):
        """
        Count the number of unique child IDs for each node anywhere in the ARG.

        Relies on edges being sorted by parent ID then by child ID (as guaranteed
        by :meth:`tskit.TreeSequence.sort`).

        :param tskit.TreeSequence ts: The tree sequence to analyse.
        :return: An array of length ``ts.num_nodes`` where element ``i`` is the
            number of unique children of node ``i`` across all edges.
        :rtype: numpy.ndarray

        See https://github.com/tskit-dev/tskit/discussions/3449
        """
        if ts.num_edges == 0:
            return np.zeros(ts.num_nodes, dtype=np.intp)
        same_parent = np.concatenate(
            (ts.edges_parent[1:] == ts.edges_parent[:-1], [False])
        )
        same_child = np.concatenate(
            (ts.edges_child[1:] == ts.edges_child[:-1], [False])
        )
        # last occurrence of each unique child per parent
        is_last_unique = ~same_parent | ~same_child
        return np.bincount(
            ts.edges_parent[is_last_unique], minlength=ts.num_nodes
        )

    @staticmethod
    def remove_edges(ts, edge_id_remove_list):
        """
        Create polytomies by removing edges from a tree sequence.

        When an edge is removed, the children of the removed edge are grafted onto
        the parent of the removed edge over the same genomic span. Mutations associated
        with the edges to be removed must be moved or deleted beforehand (a
        :exc:`ValueError` is raised if any such mutations exist).

        :param tskit.TreeSequence ts: The tree sequence to modify.
        :param list edge_id_remove_list: A list (or other iterable) of integer edge IDs
            to remove.
        :return: A new tree sequence with the specified edges removed.
        :rtype: tskit.TreeSequence

        See https://github.com/tskit-dev/tskit/discussions/2926
        """
        edges_to_remove_by_child = collections.defaultdict(list)
        edge_id_remove_list = set(edge_id_remove_list)
        for m in ts.mutations():
            if m.edge in edge_id_remove_list:
                raise ValueError(
                    "Cannot remove edges that have associated mutations"
                )
        for remove_edge in edge_id_remove_list:
            e = ts.edge(remove_edge)
            edges_to_remove_by_child[e.child].append(e)

        # sort left-to-right for each child
        for k, v in edges_to_remove_by_child.items():
            edges_to_remove_by_child[k] = sorted(v, key=lambda e: e.left)
            # check no overlaps
            for e1, e2 in zip(
                edges_to_remove_by_child[k], edges_to_remove_by_child[k][1:]
            ):
                assert e1.right <= e2.left

        # Sanity check: this means the topmost node will deal with modified edges
        # left at the end
        assert ts.edge(-1).parent not in edges_to_remove_by_child

        new_edges = collections.defaultdict(list)
        tables = ts.dump_tables()
        tables.edges.clear()
        # Edges are sorted by parent time, youngest first, so we can iterate over
        # nodes-as-parents visiting children before parents by using itertools.groupby
        for parent_id, ts_edges in itertools.groupby(
            ts.edges(), lambda e: e.parent
        ):
            # Iterate through the ts edges *plus* the polytomy edges we created in
            # previous steps. This allows re-editing polytomy edges when edges_to_remove
            # are stacked.
            edges = list(ts_edges)
            if parent_id in new_edges:
                edges += new_edges.pop(parent_id)
            if parent_id in edges_to_remove_by_child:
                for e in edges:
                    assert parent_id == e.parent
                    l = -1
                    if e.id in edge_id_remove_list:
                        continue
                    for target_edge in edges_to_remove_by_child[parent_id]:
                        assert target_edge.left > l
                        l = target_edge.left
                        if e.left >= target_edge.right:
                            # target edge is entirely to the LHS of edge e
                            continue
                        elif e.right <= target_edge.left:
                            # target edge is entirely to the RHS of edge e
                            tables.edges.append(e)
                            e = None
                            break
                        else:
                            # edge e overlaps with current target edge
                            if e.left < target_edge.left:
                                tables.edges.add_row(
                                    left=e.left,
                                    right=target_edge.left,
                                    parent=e.parent,
                                    child=e.child,
                                )
                                e = e.replace(left=target_edge.left)
                            if e.right > target_edge.right:
                                assert e.left < target_edge.right
                                new_edges[target_edge.parent].append(
                                    e.replace(
                                        right=target_edge.right,
                                        parent=target_edge.parent,
                                    )
                                )
                                e = e.replace(left=target_edge.right)
                            else:
                                assert e.left < e.right
                                new_edges[target_edge.parent].append(
                                    e.replace(parent=target_edge.parent)
                                )
                                e = None
                                break
                    if e is not None:
                        tables.edges.append(e)
            else:
                for e in edges:
                    if e.id not in edge_id_remove_list:
                        tables.edges.append(e)
        assert len(new_edges) == 0
        tables.sort()
        return tables.tree_sequence()


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

        Requires Python 3.10 or later.

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


# Top-level convenience aliases so that ``tsx.method(...)`` works.

def reorder_nodes(ts, order):
    """Alias for :meth:`TreeSequence.reorder_nodes`."""
    return TreeSequence.reorder_nodes(ts, order)


def arg_num_children(ts):
    """Alias for :meth:`TreeSequence.arg_num_children`."""
    return TreeSequence.arg_num_children(ts)


def remove_edges(ts, edge_id_remove_list):
    """Alias for :meth:`TreeSequence.remove_edges`."""
    return TreeSequence.remove_edges(ts, edge_id_remove_list)


def unsquash(edge_table, positions, edges=None):
    """Alias for :meth:`EdgeTable.unsquash`."""
    return EdgeTable.unsquash(edge_table, positions, edges)
