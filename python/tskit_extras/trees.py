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
