tskit_extras
============

Extra utility functions on top of `tskit <https://tskit.dev/tskit>`_.

These routines are considered too complex or not fully tested enough for
inclusion in core tskit, but are tested and ready to use. They can also serve
as a repository of useful "snippets" that could be reused elsewhere.

Only functions that operate on a :class:`tskit.TreeSequence` have direct
top-level aliases such as ``tsx.reorder_nodes``. Utilities for other tskit
objects are accessed through their corresponding namespace, such as
``tsx.TableCollection`` or ``tsx.EdgeTable``. The
``tsx.TreeNode.all_binary_trees(leaf_labels)`` method lazily yields every
rooted binary resolution of the labelled polytomy. For a full
:class:`tskit.Tree`, ``tsx.Tree.all_binary_resolutions(tree)`` lazily yields
trees for all combinations of its polytomies. The
``tsx.Tree.optimal_tree_refinement_parsimony_score(tree, genotypes, alleles)``
method returns the minimum parsimony score over those binary refinements.

Usage
-----

.. code-block:: python

    import tskit_extras as tsx

    # Reorder nodes in a tree sequence (oldest-first)
    reordered_ts = tsx.reorder_nodes(ts, order=-ts.nodes_time)

    # Or equivalently:
    reordered_ts = tsx.TreeSequence.reorder_nodes(ts, order=-ts.nodes_time)

    # Reorder nodes in a TableCollection in-place
    tsx.TableCollection.reorder_nodes(tables, order=-tables.nodes.time)

    # Count children per parent node across the ARG
    n_children = tsx.arg_num_children(ts)

    # Invert a simplify node map (index by new ID to obtain the original ID)
    simplified_ts, node_map = ts.simplify(samples=ts.samples()[:4], map_nodes=True)
    reverse_map = tsx.utils.invert_map(node_map)

    # Remove edges (create polytomies)
    polytomy_ts = tsx.remove_edges(ts, edge_ids_to_remove)

    # Find edges without mutations, excluding edges above samples
    unsupported = tsx.unsupported_edges(ts, per_interval=True)
    polytomy_ts = tsx.remove_edges(ts, unsupported)

    # Remove mutations above local roots, preserving sample genotypes
    ts_no_root_mutations = tsx.remove_root_mutations(ts)

    # Unsquash edges in an EdgeTable
    tsx.EdgeTable.unsquash(ts.dump_tables().edges, positions=[0.3, 0.7])

    # Lazily enumerate binary resolutions of a polytomy
    for tree in tsx.TreeNode.all_binary_trees([0, 1, 2, 3]):
        process(tree)

    # Enumerate all binary resolutions of a tskit Tree
    for resolved_tree in tsx.Tree.all_binary_resolutions(tree):
        process(resolved_tree)

    # Score observed sample states over all binary refinements
    score = tsx.Tree.optimal_tree_refinement_parsimony_score(
        tree, genotypes, alleles
    )

``remove_edges`` rejects edges carrying mutations. Removing an edge above a sample
node (including an internal sample) can also change that sample's genotype by
disconnecting it from mutations on its ancestral path. Removing only a subset of
edges returned by ``unsupported_edges`` is intended to preserve sample genotypes;
arbitrary edge selections are not guaranteed to do so. ``unsupported_edges`` excludes
both mutation-bearing edges and edges whose child is a sample. By default, support is
assessed across all genomic segments connecting a parent-child pair; set
``per_interval=True`` to assess each segment independently.
``remove_root_mutations`` moves the inherited state at each local root into the
site's ancestral state and removes those root mutations. It raises an error if the
local roots at a site inherit different states.

Functions
---------

Standalone helpers that do not operate on a particular tskit object are
available through the ``tsx.utils`` module. For example,
``tsx.utils.invert_map(node_mapping)`` inverts a simplify node map. This
namespace is intended for such utilities, rather than as a general catch-all.

- ``reorder_nodes(ts, order)`` — reorder nodes in a tree sequence
  (`discussion #3466 <https://github.com/tskit-dev/tskit/discussions/3466>`_)
- ``arg_num_children(ts)`` — count unique children per parent in the ARG
  (`discussion #3449 <https://github.com/tskit-dev/tskit/discussions/3449>`_)
- ``utils.invert_map(node_mapping)`` — invert a simplify node map
  (`issue #3 <https://github.com/hyanwong/tskit_extras/issues/3>`_)
- ``remove_edges(ts, edge_id_remove_list)`` — create polytomies by removing edges
  (`discussion #2926 <https://github.com/tskit-dev/tskit/discussions/2926>`_)
- ``remove_root_mutations(ts)`` — remove mutations above local roots
  (`advanced simplification tutorial <https://tskit.dev/tutorials/advanced_simplification.html#removing-mutations-above-the-root>`_)
- ``unsupported_edges(ts, per_interval=False)`` — find edges without mutations,
  excluding edges above samples
  (`discussion #2926 <https://github.com/tskit-dev/tskit/discussions/2926>`_)
- ``EdgeTable.unsquash(edge_table, positions)`` — opposite of ``EdgeTable.squash()``
  (`discussion #2657 <https://github.com/tskit-dev/tskit/discussions/2657>`_)
- ``TreeNode.all_binary_trees(leaf_labels)`` — yield all rooted binary resolutions of a polytomy
- ``Tree.all_binary_resolutions(tree)`` — yield all combinations resolving a tskit Tree
- ``Tree.optimal_tree_refinement_parsimony_score(
  tree, genotypes, alleles, ancestral_state=None)`` —
  minimum parsimony score over all binary refinements
