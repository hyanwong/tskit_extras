"""Utilities that operate on tskit trees and tree topologies."""

import collections.abc
import itertools
import json

import numpy as np
import tskit


class TreeNode:
    """
    Utilities for constructing tree topologies as
    :class:`tskit.combinatorics.TreeNode` objects.
    """

    @staticmethod
    def all_binary_trees(leaf_labels):
        """
        Yield every leaf-labelled rooted binary tree for the given leaf labels.

        Each unordered binary topology is yielded once. Trees are returned as
        :class:`tskit.combinatorics.TreeNode` objects, with children in the same
        canonical order as :meth:`tskit.combinatorics.TreeNode.random_binary_tree`.

        This uses the same insertion construction as
        :meth:`tskit.combinatorics.TreeNode.random_binary_tree`: each new leaf is
        attached as the sibling of an existing node under a new internal node.
        Instead of choosing one insertion point at random, this method recursively
        tries each existing node. With distinct leaf labels, removing the last
        added leaf uniquely identifies the preceding tree and insertion point, so
        each topology is yielded once.

        :param iterable leaf_labels: Labels to assign to the leaves. At least one
            label is required.
        :yields: A rooted binary tree with the given leaf labels.
        """
        leaf_labels = tuple(leaf_labels)
        if len(leaf_labels) == 0:
            raise ValueError("At least one leaf label is required")

        tree_node_class = tskit.combinatorics.TreeNode

        def canonical_copy(node, parent=None):
            copied = tree_node_class(parent=parent, label=node.label)
            children_with_keys = [
                canonical_copy(child, copied) for child in node.children
            ]
            children_with_keys.sort(key=lambda item: item[1])
            copied.children = [child for child, _ in children_with_keys]
            if len(children_with_keys) == 0:
                return copied, (1, node.label)
            return copied, (
                sum(key[0] for _, key in children_with_keys),
                min(key[1] for _, key in children_with_keys),
            )

        root = tree_node_class(label=leaf_labels[0])
        nodes = [root]

        def generate(root, next_label):
            if next_label == len(leaf_labels):
                yield canonical_copy(root)[0]
                return

            num_existing_nodes = len(nodes)
            for node_index in range(num_existing_nodes):
                node = nodes[node_index]
                parent = node.parent
                if parent is None:
                    child_index = None
                else:
                    child_index = parent.children.index(node)
                new_leaf = tree_node_class(label=leaf_labels[next_label])
                new_parent = tree_node_class(parent=parent, children=[node, new_leaf])
                node.parent = new_parent
                new_leaf.parent = new_parent
                new_root = root
                if parent is None:
                    new_root = new_parent
                else:
                    parent.children[child_index] = new_parent

                nodes.extend([new_leaf, new_parent])
                yield from generate(new_root, next_label + 1)

                nodes.pop()
                nodes.pop()
                node.parent = parent
                if parent is not None:
                    parent.children[child_index] = node

        yield from generate(root, 1)


def _polytomy_node_costs(child_costs, num_states):
    """
    Return the minimum score below a polytomy root for each root state.

    Each child's cost vector is split into ``mins`` (its best score in any state)
    and ``deltas`` (the extra cost of being in each state). Given the root state
    ``r``, a binary refinement costs one change for each other state ``t`` in a
    set ``T`` of states used by the children, so the cost is
    ``len(T) + sum_i min(deltas[i, T + [r]])``. If ``r`` is not used by any child,
    the root needs two children, so a single-state ``T`` costs one more change.

    Only states that are optimal for at least one child can usefully be in ``T``,
    and identical children are grouped, so the cost is exponential only in the
    number of such states rather than in the number of children.
    """
    costs = np.array(child_costs, dtype=float)
    mins = costs.min(axis=1)
    deltas, counts = np.unique(costs - mins[:, None], axis=0, return_counts=True)
    base = mins.sum()
    candidates = np.flatnonzero((deltas == 0).any(axis=0))

    node_costs = []
    for root_state in range(num_states):
        root_deltas = deltas[:, root_state]
        best = counts @ root_deltas
        others = [s for s in candidates if s != root_state]
        for size in range(1, len(others) + 1):
            for subset in itertools.combinations(others, size):
                subset_deltas = deltas[:, subset].min(axis=1)
                total = size + counts @ np.minimum(subset_deltas, root_deltas)
                if size == 1 and not np.any(root_deltas <= subset_deltas):
                    total += 1
                best = min(best, total)
        node_costs.append(base + best)
    return node_costs


# Above this many states, enumerating every state subset for all sites at once is
# less efficient than the per-site search over only the optimal states.
_MAX_VECTORISED_STATES = 5

# Approximate maximum number of float entries held in cost arrays at once.
_CHUNK_BUDGET = 4_000_000


def _add_edge_cost(costs):
    """
    Return the costs of choosing each state above a subtree, so that a change
    between the chosen state and the subtree's root state costs one.

    ``costs`` has shape ``(num_sites, num_states)``.
    """
    return np.minimum(costs, costs.min(axis=1, keepdims=True) + 1)


def _polytomy_costs_batch(child_costs):
    """
    A version of ``_polytomy_node_costs`` that handles many sites at once.

    ``child_costs`` has shape ``(num_children, num_sites, num_states)``; the result
    has shape ``(num_sites, num_states)``. Every subset of states is considered,
    rather than only those that are optimal for some child, so that the same code
    applies to all sites.
    """
    num_children, num_sites, num_states = child_costs.shape
    mins = child_costs.min(axis=2)
    deltas = child_costs - mins[:, :, None]
    node_costs = np.empty((num_sites, num_states))
    base = mins.sum(axis=0)
    for root_state in range(num_states):
        root_deltas = deltas[:, :, root_state]
        best = root_deltas.sum(axis=0)
        others = [s for s in range(num_states) if s != root_state]
        for size in range(1, len(others) + 1):
            for subset in itertools.combinations(others, size):
                if size == 1:
                    subset_deltas = deltas[:, :, subset[0]]
                else:
                    subset_deltas = deltas[:, :, list(subset)].min(axis=2)
                total = size + np.minimum(subset_deltas, root_deltas).sum(axis=0)
                if size == 1:
                    # No child uses the root state, so the root needs two children
                    total += (root_deltas > subset_deltas).all(axis=0)
                best = np.minimum(best, total)
        node_costs[:, root_state] = base + best
    return node_costs


def _score_sites(tree, nodes, children, sample_columns, genotypes, num_states, anc):
    """
    Return the OTR parsimony score for each row of ``genotypes``.

    ``nodes`` is the postorder node list, ``children`` maps each node to a tuple of
    its children, and ``sample_columns`` maps sample nodes to genotype columns.
    ``anc`` is an integer array with the ancestral state index for each site, or
    ``-1`` for sites where it should be chosen to minimise the score.
    """
    num_sites = genotypes.shape[0]
    states = np.arange(num_states)
    costs = {}
    for u in nodes:
        child_nodes = children[u]
        if len(child_nodes) == 0:
            node_costs = np.zeros((num_sites, num_states))
        elif len(child_nodes) == 1:
            node_costs = _add_edge_cost(costs.pop(child_nodes[0]))
        elif len(child_nodes) == 2:
            node_costs = _add_edge_cost(costs.pop(child_nodes[0]))
            node_costs += _add_edge_cost(costs.pop(child_nodes[1]))
        else:
            child_costs = np.stack([costs.pop(v) for v in child_nodes])
            if num_states <= _MAX_VECTORISED_STATES:
                node_costs = _polytomy_costs_batch(child_costs)
            else:
                node_costs = np.array(
                    [
                        _polytomy_node_costs(child_costs[:, site], num_states)
                        for site in range(num_sites)
                    ]
                )

        column = sample_columns.get(u)
        if column is not None:
            observed = genotypes[:, column]
            excluded = (observed != tskit.MISSING_DATA)[:, None] & (
                observed[:, None] != states[None, :]
            )
            node_costs = np.where(excluded, np.inf, node_costs)
        costs[u] = node_costs

    # All roots share the same state immediately above the forest. If fixed,
    # charge each root edge against it; otherwise minimize over that state.
    total = sum(_add_edge_cost(costs[root]) for root in tree.roots)
    # Negative entries of ``anc`` mean that the state is to be chosen.
    fixed = total[np.arange(num_sites), np.maximum(anc, 0)]
    return np.where(anc >= 0, fixed, total.min(axis=1))


def _check_genotypes(genotypes, num_samples, num_alleles):
    """
    Check a two-dimensional ``(num_sites, num_samples)`` genotype array.
    ``num_alleles`` is ``None``, or an integer or array with the number of alleles
    for each site.
    """
    if genotypes.shape[1] != num_samples:
        raise ValueError("genotypes must have one entry for each tree sequence sample")
    if np.any(genotypes < tskit.MISSING_DATA):
        raise ValueError("genotypes must be valid allele indexes or MISSING_DATA")
    if num_alleles is not None and np.any(
        genotypes >= np.reshape(num_alleles, (-1, 1))
    ):
        raise ValueError("genotypes must be valid allele indexes or MISSING_DATA")
    if not np.all(np.any(genotypes != tskit.MISSING_DATA, axis=1)):
        raise ValueError(
            "At least one non-missing observation must be provided for each site"
        )


def _normalise_alleles(alleles):
    alleles = tuple(alleles)
    if len(alleles) == 0:
        raise ValueError("At least one allele is required")
    if len(alleles) > 64:
        raise ValueError("A maximum of 64 states is supported")
    return alleles


def _split_alleles(alleles):
    """
    Return ``(shared, per_site)`` for the ``alleles`` of a multi-site call. A flat
    sequence of labels is shared by all sites; an iterable of sequences gives the
    alleles of each site in turn.
    """
    if alleles is None:
        return None, None
    iterator = iter(alleles)
    try:
        first = next(iterator)
    except StopIteration:
        return _normalise_alleles(()), None
    iterator = itertools.chain([first], iterator)
    if isinstance(first, str):
        return _normalise_alleles(iterator), None
    return None, iterator


def _resolve_ancestral_state(ancestral_state, alleles):
    """
    Return the ancestral state index, or ``-1`` if it is to be chosen.
    """
    if ancestral_state is None:
        return -1
    if isinstance(ancestral_state, str):
        if alleles is None:
            raise ValueError("alleles must be provided when ancestral_state is a string")
        return alleles.index(ancestral_state)
    if not isinstance(ancestral_state, int | np.integer):
        raise TypeError("ancestral_state must be an integer or allele label")
    ancestral_state = int(ancestral_state)
    if ancestral_state < 0 or (alleles is not None and ancestral_state >= len(alleles)):
        raise ValueError("ancestral_state not between 0 and (num_alleles-1)")
    return ancestral_state

class Tree:
    """
    Methods that operate on a :class:`tskit.Tree`.
    """

    @staticmethod
    def optimal_tree_refinement_parsimony_score(
        tree, genotypes, ancestral_state=None, alleles=None
    ):
        """
        Return the minimum Optimal Tree Refinement (OTR) parsimony score over all
        binary refinements of ``tree``.

        Genotypes are interpreted in the order of ``tree.tree_sequence.samples()``,
        as in :meth:`tskit.Tree.map_mutations`. Genotype entries equal to
        ``tskit.MISSING_DATA`` (``-1``) are treated as missing and do not constrain
        the score. The score counts equal-cost state changes, and observations on
        internal sample nodes are included.

        If ``ancestral_state`` is specified, it is the state immediately ancestral
        to the root(s); a change between it and a root state contributes to the
        score. Otherwise (``ancestral_state=None``), the ancestral state is
        chosen to minimize the score. ``tskit.MISSING_DATA`` is not a valid
        ancestral state. The input tree is not modified.

        The algorithm is linear in the number of children at each polytomy, but
        exponential in the number of distinct states that are optimal for at least
        one child of that polytomy. This is small unless many different states are
        observed below a single polytomy.

        :param tskit.Tree tree: The tree to score.
        :param array_like genotypes: State indexes for the tree sequence samples. This
            can be a one-dimensional array for a single site, or a two-dimensional
            array of shape ``(num_sites, num_samples)`` (including a list of
            per-site sequences), such as that returned by
            :meth:`tskit.TreeSequence.genotype_matrix`, to score many sites on the
            same tree at once. This is much faster than scoring each site
            separately, as the work done for each node is shared by all sites.
            It can also be an iterator that yields one one-dimensional array per
            site, such as ``(v.genotypes for v in ts.variants())``; the sites are
            then read and scored in blocks, so they need not all be held in memory.
        :param ancestral_state: An optional ancestral state label or state index.
            ``alleles`` must be provided when this is a string label. For many
            sites, this can also be an iterable (such as a list, or
            ``itertools.repeat(0)``) with one value for each site, in the same
            order as ``genotypes``, where ``None`` means that the state for that
            site is chosen to minimise the score. Values beyond the number of
            sites are ignored; too few values is an error.
        :param iterable alleles: Optional state labels indexed by ``genotypes``.
            If omitted, the state space is inferred from the genotype indexes
            and any integer ``ancestral_state``. For many sites, a flat sequence of
            labels is used for every site, but an iterable of sequences (such as
            ``(v.alleles for v in ts.variants())``) gives the alleles of each site
            in turn.
        :return: The minimum number of state changes over all binary refinements:
            an integer for one-dimensional ``genotypes``, or an integer array with
            one score per site for two-dimensional or iterator ``genotypes``.
        :rtype: int or numpy.ndarray
        """
        streaming = isinstance(genotypes, collections.abc.Iterator)
        single_site = False
        if not streaming:
            genotypes = tskit.util.safe_np_int_cast(genotypes, np.int8)
            if genotypes.ndim not in (1, 2):
                raise ValueError("genotypes must be a one- or two-dimensional array")
            single_site = genotypes.ndim == 1
            genotypes = np.atleast_2d(genotypes)

        # For many sites, ancestral_state and alleles can each be shared by all
        # sites or be iterables with one value per site.
        ancestral_states = None
        if single_site:
            shared_alleles = None if alleles is None else _normalise_alleles(alleles)
            alleles_by_site = None
        else:
            shared_alleles, alleles_by_site = _split_alleles(alleles)
            if isinstance(ancestral_state, collections.abc.Iterable) and not isinstance(
                ancestral_state, str
            ):
                ancestral_states = iter(ancestral_state)
        shared_ancestral_state = None
        if ancestral_states is None and alleles_by_site is None:
            shared_ancestral_state = _resolve_ancestral_state(
                ancestral_state, shared_alleles
            )

        sample_nodes = tuple(tree.tree_sequence.samples())
        nodes = list(tree.nodes(order="postorder"))
        children = {u: tuple(tree.children(u)) for u in nodes}
        sample_columns = {u: column for column, u in enumerate(sample_nodes)}

        def take(iterator, num_sites, name):
            values = list(itertools.islice(iterator, num_sites))
            if len(values) < num_sites:
                raise ValueError(f"{name} must provide a value for each site")
            return values

        def score_block(block):
            num_sites = len(block)
            if num_sites == 0:
                return np.zeros(0, dtype=np.int64)
            if alleles_by_site is None:
                site_alleles = [shared_alleles] * num_sites
                num_alleles = None if shared_alleles is None else len(shared_alleles)
            else:
                site_alleles = [
                    _normalise_alleles(a)
                    for a in take(alleles_by_site, num_sites, "alleles")
                ]
                num_alleles = np.array([len(a) for a in site_alleles])
            if shared_ancestral_state is not None:
                anc = np.full(num_sites, shared_ancestral_state)
            else:
                if ancestral_states is None:
                    values = [ancestral_state] * num_sites
                else:
                    values = take(ancestral_states, num_sites, "ancestral_state")
                anc = np.array(
                    [
                        _resolve_ancestral_state(value, site_alleles_i)
                        for value, site_alleles_i in zip(values, site_alleles)
                    ]
                )
            _check_genotypes(block, len(sample_nodes), num_alleles)

            # States that are neither observed nor the ancestral state can never
            # lower the score, so they are left out of the calculation.
            num_states = max(int(block.max()), int(anc.max())) + 1
            if num_states > 64:
                raise ValueError("A maximum of 64 states is supported")
            # Bound the memory used by the cost arrays alive at any one time
            chunk_size = max(1, _CHUNK_BUDGET // (num_states * len(nodes)))
            return np.concatenate(
                [
                    _score_sites(
                        tree,
                        nodes,
                        children,
                        sample_columns,
                        block[start : start + chunk_size],
                        num_states,
                        anc[start : start + chunk_size],
                    )
                    for start in range(0, num_sites, chunk_size)
                ]
            ).astype(np.int64)

        if streaming:
            # Stream the sites in blocks, so that they never need to all be in memory
            block_size = max(1, _CHUNK_BUDGET // (4 * len(nodes)))
            results = [np.zeros(0, dtype=np.int64)]
            while rows := list(itertools.islice(genotypes, block_size)):
                block = tskit.util.safe_np_int_cast(rows, np.int8)
                if block.ndim != 2:
                    raise ValueError(
                        "An iterator of genotypes must yield one-dimensional arrays"
                    )
                results.append(score_block(block))
            return np.concatenate(results)

        scores = score_block(genotypes)
        if single_site:
            return int(scores[0])
        return scores

    @staticmethod
    def all_binary_resolutions(tree, *, epsilon=None, record_provenance=True, **kwargs):
        """
        Yield all binary resolutions of a tree with polytomies.

        Each yielded tree is a new :class:`tskit.Tree` with the same interval as
        ``tree``. For a tree with multiple polytomies, all combinations of their
        binary resolutions are yielded.

        :param tskit.Tree tree: The tree to resolve.
        :param float epsilon: The length of each new branch. If ``None``, new
            nodes are placed at the next representable time below their parent.
        :param bool record_provenance: If True, record this operation in each
            yielded tree sequence.
        :param \\**kwargs: Further arguments passed to the returned
            :class:`tskit.Tree`.
        :yields: A tree with all polytomies resolved into binary nodes.
        """
        tables_template = tree.tree_sequence.dump_tables()
        tables_template.keep_intervals([tree.interval], simplify=False)
        polytomies = [
            (u, tuple(tree.children(u)))
            for u in tree.nodes()
            if tree.num_children(u) > 2
        ]
        selected_resolutions = {}

        def build_tree():
            tables = tables_template.copy()
            tables.edges.clear()
            for u in tree.nodes():
                if u in selected_resolutions:
                    root = selected_resolutions[u]
                    stack = [(child, u, tree.time(u)) for child in root.children]
                    while stack:
                        node, parent, parent_time = stack.pop()
                        if node.children:
                            if epsilon is None:
                                child_time = np.nextafter(parent_time, -np.inf)
                            else:
                                child_time = parent_time - epsilon
                            child = tables.nodes.add_row(time=child_time)
                            if parent_time <= child_time:
                                min_child_time = min(
                                    tree.time(v) for v in tree.children(u)
                                )
                                message = (
                                    f"Cannot resolve the degree "
                                    f"{tree.num_children(u)} polytomy rooted at node "
                                    f"{u} with minimum time difference of "
                                    f"{tree.time(u) - min_child_time} to the resolved "
                                    "leaves."
                                )
                                if epsilon is None:
                                    message += (
                                        " The time difference between nodes is so "
                                        "small that more nodes cannot be inserted "
                                        "between within the limits of floating point "
                                        "precision."
                                    )
                                elif parent_time == child_time:
                                    message += (
                                        f" The fixed epsilon value of {epsilon} is "
                                        "too small, resulting in the parent and child "
                                        "times being equal within the limits of "
                                        "numerical precision."
                                    )
                                else:
                                    message += (
                                        f" The fixed epsilon value of {epsilon} is "
                                        "too large, resulting in the parent time "
                                        "being less than the child time."
                                    )
                                raise tskit.LibraryError(message)
                            tables.edges.add_row(*tree.interval, parent, child)
                            for grandchild in node.children:
                                stack.append((grandchild, child, child_time))
                        else:
                            child = node.label
                            child_time = tree.time(child)
                            if parent_time <= child_time:
                                min_child_time = min(
                                    tree.time(v) for v in tree.children(u)
                                )
                                message = (
                                    f"Cannot resolve the degree "
                                    f"{tree.num_children(u)} polytomy rooted at node "
                                    f"{u} with minimum time difference of "
                                    f"{tree.time(u) - min_child_time} to the resolved "
                                    "leaves."
                                )
                                if epsilon is None:
                                    message += (
                                        " The time difference between nodes is so "
                                        "small that more nodes cannot be inserted "
                                        "between within the limits of floating point "
                                        "precision."
                                    )
                                else:
                                    message += (
                                        f" The fixed epsilon value of {epsilon} is "
                                        "too large, resulting in the parent time "
                                        "being less than the child time."
                                    )
                                raise tskit.LibraryError(message)
                            tables.edges.add_row(*tree.interval, parent, child)
                else:
                    for child in tree.children(u):
                        tables.edges.add_row(*tree.interval, u, child)

            if record_provenance:
                parameters = {"command": "all_binary_resolutions"}
                tables.provenances.add_row(
                    record=json.dumps(
                        tskit.provenance.get_provenance_dict(parameters)
                    )
                )
            tables.sort()
            return tables.tree_sequence().at(tree.interval.left, **kwargs)

        def generate(index):
            if index == len(polytomies):
                yield build_tree()
                return
            u, children = polytomies[index]
            for root in TreeNode.all_binary_trees(children):
                selected_resolutions[u] = root
                yield from generate(index + 1)
            del selected_resolutions[u]

        yield from generate(0)
