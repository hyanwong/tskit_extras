# tskit_extras

[![Tests](https://github.com/hyanwong/tskit_extras/actions/workflows/tests.yml/badge.svg)](https://github.com/hyanwong/tskit_extras/actions/workflows/tests.yml)

Extra utility functions on top of [tskit](https://tskit.dev/tskit).

[Online documentation](https://hyanwong.github.io/tskit_extras/) | [Build the book locally](docs/README.md)

These routines are considered too complex or not fully tested enough for
inclusion in core tskit, but are tested and ready to use. They can also serve
as a repository of useful "snippets" that could be reused elsewhere.

### Calling conventions

Functions that operate on a _tskit_ object are invoked as `.TskitObject.function()`,
for example [`.Tree.all_binary_resolutions()`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.Tree.all_binary_resolutions),
or [`.EdgeTable.unsquash()`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.EdgeTable.unsquash).
Standalone functions, such as those that operate on arrays, are stored in the
`.utils` namespace, for example
[`.utils.invert_map`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.utils.invert_map).

We recommend that you `import tskit_extras as tsx`, and as a convenient shortcut,
functions that operate on a
[`tskit.TreeSequence`](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence) object
have direct top-level aliases. E.g.
[`tsx.TreeSequence.reorder_nodes(ts, order)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.reorder_nodes)
can equally be called as
[`tsx.reorder_nodes(ts, order)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.reorder_nodes).

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

## Functions and methods

<table>
<thead>
<tr><th>Function</th><th>Description</th><th>Source discussion</th></tr>
</thead>
<tbody>
<tr><th colspan="3" scope="rowgroup">Top-level functions</th></tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.arg_num_children"><code>arg_num_children(ts)</code></a></td>
<td>Count unique children per parent in the ARG</td>
<td><a href="https://github.com/tskit-dev/tskit/discussions/3449">#3449</a></td>
</tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.remove_edges"><code>remove_edges(ts, edge_id_remove_list)</code></a></td>
<td>Create polytomies by removing edges</td>
<td><a href="https://github.com/tskit-dev/tskit/discussions/2926">#2926</a></td>
</tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.remove_root_mutations"><code>remove_root_mutations(ts)</code></a></td>
<td>Remove mutations above local roots</td>
<td><a href="https://tskit.dev/tutorials/advanced_simplification.html#removing-mutations-above-the-root">Advanced simplification tutorial</a></td>
</tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.reorder_nodes"><code>reorder_nodes(ts, order)</code></a></td>
<td>Reorder nodes in a tree sequence</td>
<td><a href="https://github.com/tskit-dev/tskit/discussions/3466">#3466</a></td>
</tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.unsupported_edges"><code>unsupported_edges(ts, per_interval=False)</code></a></td>
<td>Find edges without mutations, excluding edges above samples</td>
<td><a href="https://github.com/tskit-dev/tskit/discussions/2926">#2926</a></td>
</tr>
<tr><td colspan="3"><hr></td></tr>
<tr><th colspan="3" scope="rowgroup">Utility functions</th></tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.utils.invert_map"><code>utils.invert_map(node_mapping)</code></a></td>
<td>Invert a simplify node map</td>
<td><a href="https://github.com/hyanwong/tskit_extras/issues/3">#3</a></td>
</tr>
<tr><td colspan="3"><hr></td></tr>
<tr><th colspan="3" scope="rowgroup">Functions on other tskit objects</th></tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.EdgeTable.unsquash"><code>EdgeTable.unsquash(edge_table, positions)</code></a></td>
<td>Opposite of <a href="https://tskit.dev/tskit/docs/stable/python-api.html#tskit.EdgeTable.squash"><code>EdgeTable.squash()</code></a></td>
<td><a href="https://github.com/tskit-dev/tskit/discussions/2657">#2657</a></td>
</tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TableCollection.reorder_nodes"><code>TableCollection.reorder_nodes(tables, order)</code></a></td>
<td>Reorder nodes in a TableCollection</td>
<td>—</td>
</tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.Tree.all_binary_resolutions"><code>Tree.all_binary_resolutions(tree)</code></a></td>
<td>Enumerate all binary resolutions of a tskit Tree</td>
<td>—</td>
</tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.Tree.optimal_tree_refinement_parsimony_score"><code>Tree.optimal_tree_refinement_parsimony_score(tree, genotypes, alleles, ancestral_state=None)</code></a></td>
<td>Minimum parsimony score over all binary refinements</td>
<td>—</td>
</tr>
<tr>
<td><a href="https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeNode.all_binary_trees"><code>TreeNode.all_binary_trees(leaf_labels)</code></a></td>
<td>Enumerate all rooted binary resolutions of a polytomy</td>
<td>—</td>
</tr>
</tbody>
</table>

### TreeSequence functions

The TreeSequence functions are also available through the TreeSequence namespace:
- [`tsx.TreeSequence.reorder_nodes(ts, order)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.reorder_nodes)
- [`tsx.TreeSequence.arg_num_children(ts)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.arg_num_children)
- [`tsx.TreeSequence.remove_edges(ts, edge_id_remove_list)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.remove_edges)
- [`tsx.TreeSequence.remove_root_mutations(ts)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.remove_root_mutations)
- [`tsx.TreeSequence.unsupported_edges(ts, per_interval=False)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeSequence.unsupported_edges)

### General utilities

Helpers that do not operate on a particular tskit object are grouped in
`tsx.utils`; for example,
[`tsx.utils.invert_map(node_mapping)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.utils.invert_map)
inverts a simplify node map. 

### Functions that operate on other _tskit_ objects

Helpers for specific non-TreeSequence objects are available through their namespaces:
- [`tsx.TableCollection.reorder_nodes(tables, order)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TableCollection.reorder_nodes)
- [`tsx.EdgeTable.unsquash(edge_table, positions)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.EdgeTable.unsquash)
- [`tsx.TreeNode.all_binary_trees(leaf_labels)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.TreeNode.all_binary_trees)
- [`tsx.Tree.all_binary_resolutions(tree)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.Tree.all_binary_resolutions)
- [`tsx.Tree.optimal_tree_refinement_parsimony_score(tree, genotypes, alleles, ancestral_state=None)`](https://hyanwong.github.io/tskit_extras/api.html#tskit_extras.Tree.optimal_tree_refinement_parsimony_score)

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
