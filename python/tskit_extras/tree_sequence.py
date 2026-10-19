"""Utilities that operate on tskit tree sequences."""

import collections
import itertools
import operator

import numpy as np
import tskit

from .tables import TableCollection


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
    def _graph_descendants(ts, node):
        """
        Return a boolean array marking the descendants of ``node`` anywhere
        in the graph, regardless of whether those descendants have ``node``
        as a genetic ancestor in any local tree.

        Note that this builds arrays of starts and ends of the child edges for each
        node, which can also be done by using tskit.jit.numba.NumbaTreeSequence
        for potentially faster descendant queries.
        """
        if node < 0:
            raise ValueError(
                "Negative node IDs are disallowed; did you pass tskit.NULL?"
            )

        is_descendant = np.zeros(ts.num_nodes, dtype=bool)
        parent = ts.edges_parent
        # Edges for a given parent are contiguous in the edge table, so the child
        # edges of each node are found from the boundaries between runs of parents.
        child_start = np.zeros(ts.num_nodes, dtype=np.intp)
        child_end = np.zeros(ts.num_nodes, dtype=np.intp)
        if ts.num_edges > 0:
            run_starts = np.flatnonzero(np.r_[True, parent[1:] != parent[:-1]])
            run_ends = np.r_[run_starts[1:], ts.num_edges]
            run_parents = parent[run_starts]
            child_start[run_parents] = run_starts
            child_end[run_parents] = run_ends
        child = ts.edges_child
        # Code from https://github.com/tskit-dev/tskit/discussions/2869

        stack = [node]
        while stack:
            u = stack.pop()
            for e in range(child_start[u], child_end[u]):
                v = child[e]
                if not is_descendant[v]:
                    is_descendant[v] = True
                    stack.append(v)
        return is_descendant

    @staticmethod
    def graph_samples(ts, node):
        """
        Return the sample nodes that are graph descendants of ``node``.

        Graph descendants are reachable by following edges, regardless of
        whether the edges overlap in genomic coordinates. If ``node`` is a
        sample, it is included in the result.

        :param tskit.TreeSequence ts: The tree sequence to analyse.
        :param int node: The node ID.
        :return: A sorted array of sample node IDs.
        :rtype: numpy.ndarray
        """
        is_descendant = TreeSequence._graph_descendants(ts, node)
        if ts.nodes_flags[node] & tskit.NODE_IS_SAMPLE:
            is_descendant[node] = True
        samples = ts.samples()
        return samples[is_descendant[samples]]

    @staticmethod
    def num_graph_samples(ts, node):
        """
        Return the number of sample nodes graph-reachable from ``node``,
        including ``node`` if it is a sample.

        :param tskit.TreeSequence ts: The tree sequence to analyse.
        :param int node: The node ID.
        :return: The number of descendant samples.
        :rtype: int
        """
        return len(TreeSequence.graph_samples(ts, node))

    @staticmethod
    def remove_edges(ts, edge_id_remove_list):
        """
        Create polytomies by removing edges from a tree sequence.

        For each removed edge, its child is disconnected from its parent and the
        child's outgoing edges are grafted onto that parent over the same genomic
        span. If multiple removed edges are stacked, edges are grafted to the nearest
        ancestor that remains connected. Edges carrying mutations cannot be removed;
        move or delete those mutations first. Removing an edge whose child is a sample
        node, including an internal sample, can change that sample's genotype by
        disconnecting it from mutations on its ancestral path. Removing only a subset
        of edges returned by :meth:`unsupported_edges` is intended to preserve sample
        genotypes; arbitrary edge selections are not guaranteed to do so.

        :param tskit.TreeSequence ts: The tree sequence to modify.
        :param iterable edge_id_remove_list: Integer edge IDs to remove.
        :return: A new tree sequence with the specified edges removed.
        :rtype: tskit.TreeSequence

        See https://github.com/tskit-dev/tskit/discussions/2926
        """
        edge_ids_to_remove = set()
        for edge_id in edge_id_remove_list:
            try:
                edge_id = operator.index(edge_id)
            except TypeError as error:
                raise TypeError("edge IDs must be integers") from error
            if edge_id < 0 or edge_id >= ts.num_edges:
                raise ValueError(f"Edge ID {edge_id} is out of bounds")
            edge_ids_to_remove.add(edge_id)

        if len(edge_ids_to_remove) == 0:
            return ts.dump_tables().tree_sequence()

        for m in ts.mutations():
            if m.edge in edge_ids_to_remove:
                raise ValueError(
                    "Cannot remove edges that have associated mutations"
                )

        edges_to_remove_by_child = collections.defaultdict(list)
        for edge_id in edge_ids_to_remove:
            edge = ts.edge(edge_id)
            edges_to_remove_by_child[edge.child].append(edge)
        for edges in edges_to_remove_by_child.values():
            edges.sort(key=lambda edge: edge.left)
            for previous, current in itertools.pairwise(edges):
                assert previous.right <= current.left

        # Edge order is bottom-up in parent time: descendants are processed before
        # ancestors, so grafted edges are queued before their new parent is visited.
        assert ts.edge(-1).parent not in edges_to_remove_by_child

        tables = ts.dump_tables()
        tables.edges.clear()
        new_edges = collections.defaultdict(list)
        for parent_id, parent_edges in itertools.groupby(
            ts.edges(), key=lambda edge: edge.parent
        ):
            edges = list(parent_edges)
            if parent_id in new_edges:
                edges += new_edges.pop(parent_id)
            if parent_id in edges_to_remove_by_child:
                target_edges = edges_to_remove_by_child[parent_id]
                for edge in edges:
                    assert parent_id == edge.parent
                    if edge.id in edge_ids_to_remove:
                        continue
                    remaining_edge = edge
                    previous_target_left = -1
                    for target_edge in target_edges:
                        assert target_edge.left > previous_target_left
                        previous_target_left = target_edge.left
                        if remaining_edge is None:
                            break
                        if remaining_edge.left >= target_edge.right:
                            continue
                        if remaining_edge.right <= target_edge.left:
                            tables.edges.append(remaining_edge)
                            remaining_edge = None
                            break
                        if remaining_edge.left < target_edge.left:
                            tables.edges.append(
                                remaining_edge.replace(right=target_edge.left)
                            )
                            remaining_edge = remaining_edge.replace(
                                left=target_edge.left
                            )
                        if remaining_edge.right > target_edge.right:
                            assert remaining_edge.left < target_edge.right
                            new_edges[target_edge.parent].append(
                                remaining_edge.replace(
                                    right=target_edge.right,
                                    parent=target_edge.parent,
                                )
                            )
                            remaining_edge = remaining_edge.replace(
                                left=target_edge.right
                            )
                        else:
                            assert remaining_edge.left < remaining_edge.right
                            new_edges[target_edge.parent].append(
                                remaining_edge.replace(parent=target_edge.parent)
                            )
                            remaining_edge = None
                            break
                    if remaining_edge is not None:
                        tables.edges.append(remaining_edge)
            else:
                for edge in edges:
                    if edge.id not in edge_ids_to_remove:
                        tables.edges.append(edge)

        assert len(new_edges) == 0

        tables.sort()
        return tables.tree_sequence()

    @staticmethod
    def unsupported_edges(ts, *, per_interval=False):
        """
        Return edge IDs without a mutation, excluding edges above samples.

        Edges carrying mutations are considered supported because removing them could
        change the states inherited by their descendant samples. Edges whose child is a
        sample node are also considered supported, even if they carry no mutation,
        because removing them could disconnect the sample from mutations on its
        ancestral path.

        By default, if any genomic segment of a parent-child relationship has a
        mutation, all segments for that relationship are considered supported. If
        ``per_interval`` is True, each segment of an edge is considered independently.

        :param tskit.TreeSequence ts: The tree sequence to inspect.
        :param bool per_interval: Whether support is assessed per genomic edge interval
            rather than across the whole parent-child relationship.
        :return: The IDs of unsupported edges.
        :rtype: numpy.ndarray

        See https://github.com/tskit-dev/tskit/discussions/2926
        """
        unsupported = np.ones(ts.num_edges, dtype=bool)
        sample_nodes = set(ts.samples())
        if ts.num_edges == 0:
            return np.array([], dtype=np.intp)

        mutation_edges = {m.edge for m in ts.mutations() if m.edge != tskit.NULL}
        if mutation_edges:
            unsupported[list(mutation_edges)] = False
        unsupported[np.isin(ts.edges_child, list(sample_nodes))] = False

        if not per_interval:
            supported_relationships = {
                (ts.edge(edge_id).parent, ts.edge(edge_id).child)
                for edge_id in np.flatnonzero(~unsupported)
            }
            unsupported = np.array(
                [
                    (parent, child) not in supported_relationships
                    for parent, child in zip(ts.edges_parent, ts.edges_child)
                ],
                dtype=bool,
            )

        return np.flatnonzero(unsupported)

    @staticmethod
    def remove_root_mutations(ts):
        """
        Remove mutations on local roots by changing site ancestral states.

        For each site, mutations on non-isolated local roots are removed and their
        final inherited state is used as the site's ancestral state. Mutations on
        other nodes are retained, and mutation parents are recomputed. If local roots
        at a site inherit different states, a single ancestral state cannot represent
        the result and a :class:`ValueError` is raised.

        Sample genotypes are preserved, but removed mutations (including their
        metadata) are not retained.

        :param tskit.TreeSequence ts: The tree sequence to modify.
        :return: A new tree sequence without mutations above local roots.
        :rtype: tskit.TreeSequence

        See https://tskit.dev/tutorials/advanced_simplification.html#removing-mutations-above-the-root
        """
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        for tree in ts.trees():
            for site in tree.sites():
                ancestral_state = site.ancestral_state
                root_states = {
                    root: ancestral_state
                    for root in tree.roots
                    if not tree.is_isolated(root)
                }
                for mutation in site.mutations:
                    if mutation.node in root_states:
                        root_states[mutation.node] = mutation.derived_state
                    else:
                        tables.mutations.append(mutation.replace(parent=tskit.NULL))

                if root_states:
                    states = set(root_states.values())
                    if len(states) != 1:
                        raise ValueError(
                            "Multiple roots with different inherited states "
                            f"exist for the site at position {site.position}"
                        )
                    ancestral_state = states.pop()
                tables.sites.append(site.replace(ancestral_state=ancestral_state))

        tables.compute_mutation_parents()
        return tables.tree_sequence()
