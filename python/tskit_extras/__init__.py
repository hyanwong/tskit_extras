"""
tskit_extras: Extra utility functions on top of tskit.

Provides tested routines considered too complex or not fully tested enough
for inclusion in core tskit. Functions can be accessed as:

    import tskit_extras as tsx

    tsx.reorder_nodes(ts, order)
    tsx.TreeSequence.reorder_nodes(ts, order)
    tsx.TableCollection.reorder_nodes(tables, order)
"""

from . import _version
from .tables import EdgeTable, TableCollection
from .tree_sequence import TreeSequence
from .trees import Tree, TreeNode

__version__ = _version.tskit_extras_version

__all__ = [
    "reorder_nodes",
    "arg_num_children",
    "remove_edges",
    "TableCollection",
    "TreeSequence",
    "EdgeTable",
    "TreeNode",
    "Tree",
]


def reorder_nodes(ts, order):
    """Alias for :meth:`TreeSequence.reorder_nodes`."""
    return TreeSequence.reorder_nodes(ts, order)


def arg_num_children(ts):
    """Alias for :meth:`TreeSequence.arg_num_children`."""
    return TreeSequence.arg_num_children(ts)


def remove_edges(ts, edge_id_remove_list):
    """Alias for :meth:`TreeSequence.remove_edges`."""
    return TreeSequence.remove_edges(ts, edge_id_remove_list)
