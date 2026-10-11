---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Examples

These executable examples use [msprime](https://tskit.dev/msprime/) to create
tree sequences and demonstrate the public `tskit_extras` API.

## Simulate and simplify

```{code-cell} ipython3
import msprime
import tskit_extras as tsx

ts = msprime.sim_ancestry(
    samples=4,
    sequence_length=1_000,
    recombination_rate=1e-8,
    population_size=10_000,
    random_seed=13,
)
ts.num_trees
```

Simplification can return a map from original node IDs to simplified node IDs.
Invert it with the standalone helper in `tsx.utils`:

```{code-cell} ipython3
simplified_ts, node_map = ts.simplify(samples=ts.samples()[:4], map_nodes=True)
reverse_map = tsx.utils.invert_map(node_map)
len(reverse_map), simplified_ts.num_nodes
```

The top-level functions are convenient aliases for TreeSequence methods:

```{code-cell} ipython3
reordered_ts = tsx.reorder_nodes(ts, order=-ts.nodes_time)
reordered_ts.num_nodes
```
