# tskit_extras

[![Tests](https://github.com/hyanwong/tskit_extras/actions/workflows/tests.yml/badge.svg)](https://github.com/hyanwong/tskit_extras/actions/workflows/tests.yml)

Extra utility functions on top of [tskit](https://tskit.dev/tskit).

These routines are considered too complex or not fully tested enough for
inclusion in core tskit, but are tested and ready to use. They can also serve
as a repository of useful "snippets" that could be reused elsewhere.

### Calling conventions

Functions that operate on a _tskit_ object are invoked as `.TskitObject.function()`,
for example `.Tree.all_binary_resolutions()`, or `.EdgeTable.unsquash()`.
Stand-alone functions e.g. that operate on an array are stores in the `.utils`
namespace, for example `.utils.invert_map`.

We recommend that you `import tskit_extras as tsx`, and as a convenient shortcut,
functions that operate on a `tskit.TreeSequence` object
have direct top-level aliases. E.g. `tsx.TreeSequence.reorder_nodes(ts, order)`
can equally be called as `tsx.reorder_nodes(ts, order)`. 

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

# Invert a simplify node map (index by new ID to obtain the original ID)
simplified_ts, node_map = ts.simplify(samples=ts.samples()[:4], map_nodes=True)
reverse_map = tsx.utils.invert_map(node_map)

# Remove edges (create polytomies)
polytomy_ts = tsx.remove_edges(ts, edge_ids_to_remove)

# Find edges without mutations (excluding edges above samples)
unsupported = tsx.unsupported_edges(ts, per_interval=True)
polytomy_ts = tsx.remove_edges(ts, unsupported)

# Remove mutations above local roots, preserving sample genotypes
ts_no_root_mutations = tsx.remove_root_mutations(ts)

# Unsquash edges in an EdgeTable at given positions
tsx.EdgeTable.unsquash(tables.edges, positions=[0.3, 0.7])

# Lazily enumerate binary resolutions of a polytomy
for tree in tsx.TreeNode.all_binary_trees([0, 1, 2, 3]):
    process(tree)

# Enumerate all binary resolutions of a tskit Tree
for resolved_tree in tsx.Tree.all_binary_resolutions(tree):
    process(resolved_tree)

# Score observed sample states over all binary refinements
score = tsx.Tree.optimal_tree_refinement_parsimony_score(tree, genotypes, alleles)
```

`remove_edges` rejects edges carrying mutations. Removing an edge above a sample node
(including an internal sample) can also change that sample's genotype by disconnecting
it from mutations on its ancestral path. Removing only a subset of edges returned by
`unsupported_edges` is intended to preserve sample genotypes; arbitrary edge selections
are not guaranteed to do so. `unsupported_edges` excludes both mutation-bearing edges
and edges whose child is a sample. By default, support is assessed across all genomic
segments connecting a parent-child pair; set `per_interval=True` to assess each segment
independently.
`remove_root_mutations` moves the inherited state at each local root into the site's
ancestral state and removes those root mutations. It raises an error if the local
roots at a site inherit different states.

## Functions

| Function | Description | Source discussion |
|---|---|---|
| `reorder_nodes(ts, order)` | Reorder nodes in a tree sequence | [#3466](https://github.com/tskit-dev/tskit/discussions/3466) |
| `arg_num_children(ts)` | Count unique children per parent in the ARG | [#3449](https://github.com/tskit-dev/tskit/discussions/3449) |
| `remove_edges(ts, edge_id_remove_list)` | Create polytomies by removing edges | [#2926](https://github.com/tskit-dev/tskit/discussions/2926) |
| `remove_root_mutations(ts)` | Remove mutations above local roots | [advanced simplification tutorial](https://tskit.dev/tutorials/advanced_simplification.html#removing-mutations-above-the-root) |
| `unsupported_edges(ts, per_interval=False)` | Find edges without mutations, excluding edges above samples | [#2926](https://github.com/tskit-dev/tskit/discussions/2926) |
| `EdgeTable.unsquash(edge_table, positions)` | Opposite of `EdgeTable.squash()` | [#2657](https://github.com/tskit-dev/tskit/discussions/2657) |

The TreeSequence functions are also available through the TreeSequence namespace:
- `tsx.TreeSequence.reorder_nodes(ts, order)`
- `tsx.TreeSequence.arg_num_children(ts)`
- `tsx.TreeSequence.remove_edges(ts, edge_id_remove_list)`
- `tsx.TreeSequence.remove_root_mutations(ts)`
- `tsx.TreeSequence.unsupported_edges(ts, per_interval=False)`

Helpers that do not operate on a particular tskit object are grouped in
`tsx.utils`; for example, `tsx.utils.invert_map(node_mapping)` inverts a
simplify node map. This module is for standalone utilities, rather than a
general catch-all. Helpers for specific non-TreeSequence objects are available
through their namespaces:
- `tsx.TableCollection.reorder_nodes(tables, order)`
- `tsx.EdgeTable.unsquash(edge_table, positions)`
- `tsx.TreeNode.all_binary_trees(leaf_labels)`
- `tsx.Tree.all_binary_resolutions(tree)`
- `tsx.Tree.optimal_tree_refinement_parsimony_score(
  tree, genotypes, alleles, ancestral_state=None)`

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
