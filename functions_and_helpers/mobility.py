# functions_and_helpers/mobility.py
"""
Reusable mobility computation utilities for Robot Geometry toolbox.

Provides:
- JOINT_FI_DEFAULTS: mapping from joint code to default relative freedom fi.
- joint_type_to_fi(typ): returns default fi or raises KeyError if unknown.
- compute_mobility(n, fi_list, spherical=False, deltaM=0.0): returns a dict with
  spatial/spherical mobility values (general and closed-loop when applicable).
- is_single_chain(n, j): convenience check for closed single-chain (n==j).
"""

from typing import Iterable, Dict, Optional

# Default relative freedoms (fi) for standard joint types
# R: revolute (1), P: prismatic (1), C: cylindric (2), B: ball-and-socket (3), PL: planar (3)
JOINT_FI_DEFAULTS = {
    "R": 1.0,
    "P": 1.0,
    "C": 2.0,
    "B": 3.0,   # ball-and-socket
    "PL": 3.0,  # planar
}

def joint_type_to_fi(typ: str) -> float:
    """
    Return the default relative freedom (fi) for a joint type code.
    Valid types: 'R', 'P', 'C', 'B', 'PL'.
    """
    typ = (typ or "").strip().upper()
    if typ not in JOINT_FI_DEFAULTS:
        raise KeyError(f"Unknown joint type '{typ}'. Valid: {', '.join(JOINT_FI_DEFAULTS.keys())}")
    return JOINT_FI_DEFAULTS[typ]


def is_single_chain(n: int, j: int) -> bool:
    """Returns True if the mechanism is a single-chain closed loop (n == j)."""
    return int(n) == int(j)


def compute_mobility(n: int, fi_list: Iterable[float], spherical: bool = False, deltaM: float = 0.0) -> Dict[str, Optional[float]]:
    """
    Compute mobility values for spatial/spherical mechanisms.

    Args:
        n: number of bodies (include ground).
        fi_list: iterable of per-joint relative freedoms (fi).
        spherical: if True, interpret 'primary' mechanism as spherical; values for both
                   spatial and spherical are returned regardless.
        deltaM: special-geometry correction term (added after base mobility).

    Returns:
        dict with keys:
            'sum_fi'             : sum of fi
            'spatial_general'    : 6(n-1) - Σ(6 - fi) + ΔM
            'spatial_closed'     : Σfi - 6 + ΔM  (only if n == j, else None)
            'spherical_general'  : 3(n-1) - Σ(3 - fi) + ΔM
            'spherical_closed'   : Σfi - 3 + ΔM  (only if n == j, else None)
            'group'              : Σfi - 3        (only if n == j, else None)
    """
    n = int(n)
    fi = [float(x) for x in fi_list]
    sum_fi = sum(fi)

    spatial_general = 6 * (n - 1) - sum((6 - x) for x in fi) + float(deltaM)
    spherical_general = 3 * (n - 1) - sum((3 - x) for x in fi) + float(deltaM)

    closed = False
    spatial_closed = None
    spherical_closed = None
    group = None
    # We treat 'closed' purely as n == j
    j = len(fi)
    if n == j:
        closed = True
        spatial_closed = sum_fi - 6 + float(deltaM)
        spherical_closed = sum_fi - 3 + float(deltaM)
        group = sum_fi - 3

    return {
        "sum_fi": sum_fi,
        "spatial_general": spatial_general,
        "spatial_closed": spatial_closed,
        "spherical_general": spherical_general,
        "spherical_closed": spherical_closed,
        "group": group,
    }
