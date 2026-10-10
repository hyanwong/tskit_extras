"""Tests for enumerating binary tree topologies and resolutions."""

import inspect

import pytest
import tskit

import tskit_extras as tsx


class TestAllBinaryTrees:
    """Tests for enumerating binary resolutions of a polytomy."""

    @pytest.mark.parametrize(
        ("num_leaves", "expected_count"),
        [(2, 1), (3, 3), (4, 15), (5, 105)],
    )
    def test_number_of_resolutions(self, num_leaves, expected_count):
        trees = list(tsx.TreeNode.all_binary_trees(range(num_leaves)))
        assert len(trees) == expected_count

    def test_yields_all_unique_binary_trees(self):
        trees = list(tsx.TreeNode.all_binary_trees([3, 1, 4, 2]))
        assert len(trees) == 15

        def signature(node):
            if not node.children:
                return ("leaf", node.label)
            return tuple(signature(child) for child in node.children)

        assert len({signature(tree) for tree in trees}) == 15
        for tree in trees:
            assert tree.parent is None
            leaves = []
            stack = [tree]
            while stack:
                node = stack.pop()
                if node.children:
                    assert len(node.children) == 2
                    for child in node.children:
                        assert child.parent is node
                        stack.append(child)
                else:
                    leaves.append(node.label)
            assert sorted(leaves) == [1, 2, 3, 4]

    def test_newick_and_rank_are_unique(self):
        num_leaves = 5
        trees = [
            root.as_tables(
                num_leaves=num_leaves, span=1, branch_length=1
            ).tree_sequence().first()
            for root in tsx.TreeNode.all_binary_trees(range(num_leaves))
        ]
        assert len(trees) == 105
        assert len({tree.as_newick() for tree in trees}) == 105
        assert len({tree.rank() for tree in trees}) == 105

    def test_all_binary_resolutions_of_star_tree(self):
        tree = tskit.Tree.generate_star(4)
        resolutions = tsx.Tree.all_binary_resolutions(
            tree, record_provenance=False
        )
        assert inspect.isgenerator(resolutions)
        resolutions = list(resolutions)
        assert len(resolutions) == 15
        assert all(
            resolution.num_children(u) <= 2
            for resolution in resolutions
            for u in resolution.nodes()
        )
        assert len({resolution.as_newick() for resolution in resolutions}) == 15
        assert len({resolution.rank() for resolution in resolutions}) == 15
        assert all(resolution.interval == tree.interval for resolution in resolutions)

    def test_all_binary_resolutions_combine_multiple_polytomies(self):
        tables = tskit.TableCollection(sequence_length=1)
        leaves = [
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
            for _ in range(6)
        ]
        left = tables.nodes.add_row(time=1)
        right = tables.nodes.add_row(time=1)
        root = tables.nodes.add_row(time=2)
        for parent, children in [(left, leaves[:3]), (right, leaves[3:])]:
            for child in children:
                tables.edges.add_row(0, 1, parent, child)
        tables.edges.add_row(0, 1, root, left)
        tables.edges.add_row(0, 1, root, right)
        tables.sort()
        tree = tables.tree_sequence().first()

        resolutions = list(
            tsx.Tree.all_binary_resolutions(tree, record_provenance=False)
        )
        assert len(resolutions) == 9
        assert len({resolution.rank() for resolution in resolutions}) == 9

    def test_all_binary_resolutions_mixed_arity(self):
        tree = tskit.Tree.generate_balanced(8, arity=3)
        num_tritomies = sum(tree.num_children(u) == 3 for u in tree.nodes())
        assert num_tritomies == 2
        assert any(tree.num_children(u) == 2 for u in tree.nodes())

        resolutions = list(
            tsx.Tree.all_binary_resolutions(tree, record_provenance=False)
        )
        assert len(resolutions) == 9
        assert len({resolution.rank() for resolution in resolutions}) == len(
            resolutions
        )
        assert all(
            resolution.num_children(u) <= 2
            for resolution in resolutions
            for u in resolution.nodes()
        )

    def test_is_lazy_and_handles_one_leaf(self):
        trees = tsx.TreeNode.all_binary_trees(iter([42]))
        assert inspect.isgenerator(trees)
        tree = next(trees)
        assert tree.label == 42
        assert tree.children == []
        with pytest.raises(StopIteration):
            next(trees)

    def test_requires_at_least_one_leaf(self):
        with pytest.raises(ValueError, match="At least one leaf label"):
            list(tsx.TreeNode.all_binary_trees([]))
