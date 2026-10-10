# tskit_extras

[![Tests](https://github.com/hyanwong/tskit_extras/actions/workflows/tests.yml/badge.svg)](https://github.com/hyanwong/tskit_extras/actions/workflows/tests.yml)

Extra utility functions on top of [tskit](https://tskit.dev/tskit).

These routines are considered too complex or not fully tested enough for
inclusion in core tskit, but are tested and ready to use. They can also serve
as a repository of useful "snippets" that could be reused elsewhere.

Only functions that operate on a `TreeSequence` have direct top-level aliases
such as `tsx.reorder_nodes`. Utilities for other tskit objects are accessed
through their corresponding namespace, such as `tsx.TableCollection` or
`tsx.EdgeTable`.

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

# Unsquash edges in an EdgeTable at given positions
tsx.EdgeTable.unsquash(tables.edges, positions=[0.3, 0.7])
```

## Functions

| Function | Description | Source discussion |
|---|---|---|
| `reorder_nodes(ts, order)` | Reorder nodes in a tree sequence | [#3466](https://github.com/tskit-dev/tskit/discussions/3466) |
| `arg_num_children(ts)` | Count unique children per parent in the ARG | [#3449](https://github.com/tskit-dev/tskit/discussions/3449) |
| `remove_edges(ts, edge_id_remove_list)` | Create polytomies by removing edges | [#2926](https://github.com/tskit-dev/tskit/discussions/2926) |
| `EdgeTable.unsquash(edge_table, positions)` | Opposite of `EdgeTable.squash()` | [#2657](https://github.com/tskit-dev/tskit/discussions/2657) |

The TreeSequence functions are also available through the TreeSequence namespace:
- `tsx.TreeSequence.reorder_nodes(ts, order)`
- `tsx.TreeSequence.arg_num_children(ts)`
- `tsx.TreeSequence.remove_edges(ts, edge_id_remove_list)`

Utilities for other object types are only available through their namespaces:
- `tsx.TableCollection.reorder_nodes(tables, order)`
- `tsx.EdgeTable.unsquash(edge_table, positions)`

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
