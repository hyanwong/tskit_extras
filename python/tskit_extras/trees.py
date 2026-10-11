"""Utilities that operate on tskit trees and tree topologies."""

import json

import numpy as np
import tskit


class TreeNode:
    """
    Utilities for constructing tree topologies as :class:`tskit.TreeNode` objects.
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


class Tree:
    """
    Methods that operate on a :class:`tskit.Tree`.
    """

    @staticmethod
    def optimal_tree_refinement_parsimony_score(
        tree, genotypes, alleles, ancestral_state=None
    ):
        """
        Return the minimum Optimal Tree Refinement (OTR) parsimony score over all
        binary refinements of ``tree``.

        Genotypes are interpreted in the order of ``tree.tree_sequence.samples()``,
        as in :meth:`tskit.Tree.map_mutations`. Missing genotypes
        (``tskit.MISSING_DATA``) do not constrain a node. The score counts
        equal-cost state changes, and observations on internal sample nodes are
        included.

        If ``ancestral_state`` is specified, it is the state immediately ancestral
        to the root(s); a change between it and a root state contributes to the
        score. Otherwise, the ancestral state is chosen to minimize the score.
        The input tree is not modified.

        The binary-state algorithm is linear in the number of children at each
        polytomy. For more than two states, an exact subset dynamic program is used;
        its cost is exponential in the degree of each polytomy.

        :param tskit.Tree tree: The tree to score.
        :param array_like genotypes: State indexes for the tree sequence samples.
        :param iterable alleles: The state labels indexed by ``genotypes``.
        :param ancestral_state: An optional ancestral state label or state index.
        :return: The minimum number of state changes over all binary refinements.
        :rtype: int
        """
        alleles = tuple(alleles)
        if len(alleles) == 0:
            raise ValueError("At least one allele is required")
        if len(alleles) > 64:
            raise ValueError("A maximum of 64 states is supported")

        genotypes = tskit.util.safe_np_int_cast(genotypes, np.int8)
        if genotypes.ndim != 1:
            raise ValueError("genotypes must be a one-dimensional array")
        sample_nodes = tuple(tree.tree_sequence.samples())
        if len(genotypes) != len(sample_nodes):
            raise ValueError(
                "genotypes must have one entry for each tree sequence sample"
            )
        if np.any(genotypes < tskit.MISSING_DATA) or np.any(
            genotypes >= len(alleles)
        ):
            raise ValueError("genotypes must be valid allele indexes or MISSING_DATA")
        if not np.any(genotypes != tskit.MISSING_DATA):
            raise ValueError("At least one non-missing observation must be provided")

        if isinstance(ancestral_state, str):
            ancestral_state = alleles.index(ancestral_state)
        elif ancestral_state is not None:
            if not isinstance(ancestral_state, int | np.integer):
                raise TypeError("ancestral_state must be an integer or allele label")
            ancestral_state = int(ancestral_state)
            if ancestral_state < 0 or ancestral_state >= len(alleles):
                raise ValueError("ancestral_state not between 0 and (num_alleles-1)")

        num_states = len(alleles)
        sample_state = dict(zip(sample_nodes, genotypes))
        costs = [None] * tree.tree_sequence.num_nodes

        def cost_to_state(state_costs, state):
            # Add the edge cost for choosing `state` above this child subtree.
            return min(
                cost + (child_state != state)
                for child_state, cost in enumerate(state_costs)
            )

        # costs[u][s] is the minimum score below u conditional on u having state s.
        for u in tree.nodes(order="postorder"):
            children = tuple(tree.children(u))
            child_costs = [costs[v] for v in children]

            if len(children) == 0:
                node_costs = [0] * num_states
            elif len(children) == 1:
                node_costs = [
                    cost_to_state(child_costs[0], state)
                    for state in range(num_states)
                ]
            elif len(children) == 2:
                node_costs = [
                    sum(cost_to_state(vector, state) for vector in child_costs)
                    for state in range(num_states)
                ]
            elif num_states == 2:
                node_costs = []
                for parent_state in range(num_states):
                    # In a binary refinement, child roots can all use the
                    # parent's state, all use the opposite state, or use both.
                    # Those cases cost 0, 2, or 1 changes at the new root,
                    # respectively; the assignment DP finds the cheapest mixed
                    # assignment without enumerating child subsets.
                    same_state_cost = sum(
                        vector[parent_state] for vector in child_costs
                    )
                    opposite_state_cost = (
                        sum(vector[1 - parent_state] for vector in child_costs) + 2
                    )

                    assignment_costs = [0, np.inf, np.inf, np.inf]
                    for vector in child_costs:
                        next_assignment_costs = [np.inf] * 4
                        for category, partial_cost in enumerate(assignment_costs):
                            for state in range(num_states):
                                next_category = category | (1 << state)
                                next_assignment_costs[next_category] = min(
                                    next_assignment_costs[next_category],
                                    partial_cost + vector[state],
                                )
                        assignment_costs = next_assignment_costs

                    mixed_state_cost = assignment_costs[3] + 1
                    node_costs.append(
                        min(
                            same_state_cost,
                            opposite_state_cost,
                            mixed_state_cost,
                        )
                    )
            else:
                num_children = len(children)
                num_subsets = 1 << num_children
                # subset_costs[S][s] is the minimum score for a binary subtree
                # rooted in state s and containing exactly the child subtrees in S.
                # Singleton subsets are the original child nodes; larger subsets
                # are joined by a newly introduced binary node.
                subset_costs = [None] * num_subsets
                costs_to_parent_state = [None] * num_subsets
                for i, vector in enumerate(child_costs):
                    subset = 1 << i
                    subset_costs[subset] = vector
                    costs_to_parent_state[subset] = [
                        cost_to_state(vector, state) for state in range(num_states)
                    ]

                for subset in range(1, num_subsets):
                    if subset.bit_count() == 1:
                        continue
                    best_costs = [np.inf] * num_states
                    first_bit = subset & -subset
                    left_subset = (subset - 1) & subset
                    while left_subset:
                        right_subset = subset ^ left_subset
                        # Each bipartition is symmetric, so require the lowest
                        # indexed child to be in the left subset.
                        if right_subset and left_subset & first_bit:
                            for state in range(num_states):
                                best_costs[state] = min(
                                    best_costs[state],
                                    costs_to_parent_state[left_subset][state]
                                    + costs_to_parent_state[right_subset][state],
                                )
                        left_subset = (left_subset - 1) & subset

                    subset_costs[subset] = best_costs
                    costs_to_parent_state[subset] = [
                        cost_to_state(best_costs, state) for state in range(num_states)
                    ]

                node_costs = subset_costs[-1]

            if tree.is_sample(u):
                observed_state = sample_state[u]
                if observed_state != tskit.MISSING_DATA:
                    node_costs = [
                        cost if state == observed_state else np.inf
                        for state, cost in enumerate(node_costs)
                    ]
            costs[u] = node_costs

        # All roots share the same state immediately above the forest. If fixed,
        # charge each root edge against it; otherwise minimize over that state.
        root_costs = [costs[root] for root in tree.roots]
        if ancestral_state is None:
            return min(
                sum(cost_to_state(vector, state) for vector in root_costs)
                for state in range(num_states)
            )
        return sum(
            cost_to_state(vector, ancestral_state) for vector in root_costs
        )

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
