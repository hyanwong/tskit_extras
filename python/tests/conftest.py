"""
Shared test fixtures and helpers for tskit_extras tests.
"""

import msprime


def get_simple_ts(n=10, seed=42, sequence_length=1e4, recombination_rate=1e-3):
    """Return a simple tree sequence for testing."""
    return msprime.sim_ancestry(
        n,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=seed,
    )


def get_ts_with_mutations(n=10, seed=42, mutation_rate=1e-3):
    """Return a tree sequence with mutations for testing."""
    ts = get_simple_ts(n=n, seed=seed)
    return msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seed)
