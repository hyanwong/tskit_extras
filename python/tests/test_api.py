"""Tests for the public tskit_extras API."""

import inspect
from importlib.metadata import metadata

import tskit_extras as tsx


def test_minimum_python_version():
    assert metadata("tskit_extras")["Requires-Python"] == ">=3.11"


def test_version_string():
    assert isinstance(tsx.__version__, str)
    assert len(tsx.__version__) > 0


def test_module_exports():
    for name in [
        "reorder_nodes",
        "arg_num_children",
        "remove_edges",
        "unsupported_edges",
        "TreeSequence",
        "TableCollection",
        "EdgeTable",
        "TreeNode",
        "Tree",
    ]:
        assert hasattr(tsx, name)
    assert not hasattr(tsx, "unsquash")


def test_top_level_functions_take_tree_sequence_first():
    functions = [
        (name, function)
        for name, function in inspect.getmembers(tsx, inspect.isfunction)
        if function.__module__ == tsx.__name__
    ]
    assert functions
    for name, function in functions:
        parameters = inspect.signature(function).parameters
        assert parameters and next(iter(parameters)) == "ts", (
            f"Top-level function {name} must take a TreeSequence as its "
            "first argument"
        )
