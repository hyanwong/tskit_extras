tskit_extras
============

Extra utility functions on top of `tskit <https://tskit.dev/tskit>`_.

These routines are considered too complex or not fully tested enough for
inclusion in core tskit, but are tested and ready to use. They can also serve
as a repository of useful "snippets" that could be reused elsewhere.

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

    # Remove edges (create polytomies)
    polytomy_ts = tsx.remove_edges(ts, edge_ids_to_remove)

    # Unsquash edges in an EdgeTable
    tsx.unsquash(ts.dump_tables().edges, positions=[0.3, 0.7])

Functions
---------

- ``reorder_nodes(ts, order)`` — reorder nodes in a tree sequence
  (`discussion #3466 <https://github.com/tskit-dev/tskit/discussions/3466>`_)
- ``arg_num_children(ts)`` — count unique children per parent in the ARG
  (`discussion #3449 <https://github.com/tskit-dev/tskit/discussions/3449>`_)
- ``remove_edges(ts, edge_id_remove_list)`` — create polytomies by removing edges
  (`discussion #2926 <https://github.com/tskit-dev/tskit/discussions/2926>`_)
- ``unsquash(edge_table, positions)`` — opposite of ``EdgeTable.squash()``
  (`discussion #2657 <https://github.com/tskit-dev/tskit/discussions/2657>`_)
