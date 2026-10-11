# tskit_extras

[![Tests](https://github.com/hyanwong/tskit_extras/actions/workflows/tests.yml/badge.svg)](https://github.com/hyanwong/tskit_extras/actions/workflows/tests.yml)

Extra utility functions on top of [tskit](https://tskit.dev/tskit).

These routines are considered too complex or not fully tested enough for
inclusion in core tskit, but are tested and ready to use. They can also serve
as a repository of useful "snippets" that could be reused elsewhere.

Only functions that operate on a `TreeSequence` have direct top-level aliases
such as `tsx.reorder_nodes`. Utilities for other tskit objects are accessed
through their corresponding namespace, such as `tsx.TableCollection` or
`tsx.EdgeTable`. The `tsx.TreeNode.all_binary_trees(leaf_labels)` method lazily
yields every rooted binary resolution of the labeled polytomy. For a full
`tskit.Tree`, `tsx.Tree.all_binary_resolutions(tree)` lazily yields trees for all
combinations of its polytomies.

## Installation

```bash
pip install tskit_extras
```

## Usage

```python
import tskit_extras as tsx

# Reorder nodes in a tree sequence (oldest-first)
reordered_ts = tsx.reorder_nodes(ts, order=-ts.nodes_time)

# Or equivalently, using the class-based API:
reordered_ts = tsx.TreeSequence.reorder_nodes(ts, order=-ts.nodes_time)

# Reorder nodes in a TableCollection in-place
tsx.TableCollection.reorder_nodes(tables, order=-tables.nodes.time)

# Count children per parent node across the ARG
n_children = tsx.arg_num_children(ts)

# Remove edges (create polytomies)
polytomy_ts = tsx.remove_edges(ts, edge_ids_to_remove)

# Find edges without mutations (excluding edges above samples)
unsupported = tsx.unsupported_edges(ts, per_interval=True)
polytomy_ts = tsx.remove_edges(ts, unsupported)

# Unsquash edges in an EdgeTable at given positions
tsx.EdgeTable.unsquash(tables.edges, positions=[0.3, 0.7])

# Lazily enumerate binary resolutions of a polytomy
for tree in tsx.TreeNode.all_binary_trees([0, 1, 2, 3]):
    process(tree)

# Enumerate all binary resolutions of a tskit Tree
for resolved_tree in tsx.Tree.all_binary_resolutions(tree):
    process(resolved_tree)
```

`remove_edges` rejects edges carrying mutations. Removing an edge above a sample node
(including an internal sample) can also change that sample's genotype by disconnecting
it from mutations on its ancestral path. Removing only a subset of edges returned by
`unsupported_edges` is intended to preserve sample genotypes; arbitrary edge selections
are not guaranteed to do so. `unsupported_edges` excludes both mutation-bearing edges
and edges whose child is a sample. By default, support is assessed across all genomic
segments connecting a parent-child pair; set `per_interval=True` to assess each segment
independently.

## Functions

| Function | Description | Source discussion |
|---|---|---|
| `reorder_nodes(ts, order)` | Reorder nodes in a tree sequence | [#3466](https://github.com/tskit-dev/tskit/discussions/3466) |
| `arg_num_children(ts)` | Count unique children per parent in the ARG | [#3449](https://github.com/tskit-dev/tskit/discussions/3449) |
| `remove_edges(ts, edge_id_remove_list)` | Create polytomies by removing edges | [#2926](https://github.com/tskit-dev/tskit/discussions/2926) |
| `unsupported_edges(ts, per_interval=False)` | Find edges without mutations, excluding edges above samples | [#2926](https://github.com/tskit-dev/tskit/discussions/2926) |
| `EdgeTable.unsquash(edge_table, positions)` | Opposite of `EdgeTable.squash()` | [#2657](https://github.com/tskit-dev/tskit/discussions/2657) |

The TreeSequence functions are also available through the TreeSequence namespace:
- `tsx.TreeSequence.reorder_nodes(ts, order)`
- `tsx.TreeSequence.arg_num_children(ts)`
- `tsx.TreeSequence.remove_edges(ts, edge_id_remove_list)`
- `tsx.TreeSequence.unsupported_edges(ts, per_interval=False)`

Utilities for other object types are only available through their namespaces:
- `tsx.TableCollection.reorder_nodes(tables, order)`
- `tsx.EdgeTable.unsquash(edge_table, positions)`
- `tsx.TreeNode.all_binary_trees(leaf_labels)`
- `tsx.Tree.all_binary_resolutions(tree)`

## Repository structure

```
python/              # Python package
  tskit_extras/      # Main package source
  tests/             # Test suite (pytest)
  pyproject.toml     # Package metadata and tool configuration
```

## Development

```bash
cd python
pip install -e ".[dev]"
pytest tests/
ruff check tskit_extras/ tests/
```

## Contributing

This repo follows the same conventions as
[tskit](https://github.com/tskit-dev/tskit). Contributions welcome via
pull request.
