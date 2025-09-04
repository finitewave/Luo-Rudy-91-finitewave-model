"""
ops.py — mathematical core (“ops”) for a Finitewave model.

Contract (TL;DR):
- This module is the single source of truth for the model equations.
- Provide pure Python functions with scalar inputs/outputs (no NumPy arrays, no classes, no globals).
- Do NOT add numba/jax/torch here — the Finitewave runtime will wrap these functions for you.
- Stimulus and time integration are handled outside of the model. Here you only return time derivatives.
- At minimum implement:
    * get_variables() -> dict[str, float]
    * get_parameters() -> dict[str, float]
    * calc_* functions that return d(state)/dt (one per state variable).
"""

__all__ = (
    "get_variables",
    "get_parameters",
    "calc_rhs",  # add other calc_* functions as needed
    # e.g. "calc_dv",
)


def get_variables() -> dict[str, float]:
    """
    Returns default initial values for state variables.
    Example (ionic model): {"u": -84.0, "m": 0.01, "h": 0.99}
    Example (phenomenological): {"u": 0.0, "v": 0.1}
    """
    raise NotImplementedError("The get_variables method must be implemented in a subclass.")


def get_parameters() -> dict[str, float]:
    """
    Returns default parameter values for the model.
    Example (ionic model): {"g_Na": 120.0, "E_Na": 50.0, "C_m": 1.0}
    Example (phenomenological): {"a": 0.1, "b": 0.5}
    """
    raise NotImplementedError("The get_parameters method must be implemented in a subclass.")


def calc_rhs() -> float:
    """
    Computes the right-hand side of the model.
    """
    raise NotImplementedError("The calc_rhs method must be implemented in a subclass.")

