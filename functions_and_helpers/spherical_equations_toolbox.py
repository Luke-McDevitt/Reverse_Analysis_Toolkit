
"""
spherical_equations_toolbox.py
A lightweight, extensible toolbox for evaluating spherical mechanism equation tables
(derived from the Appendix "Spherical/Polar Equations" and "Direction Cosines").
- Triangle: Fundamental, Polar, Direction Cosines
- Quadrilateral: Fundamental (+subsidiary), Polar, Direction Cosines
- Pentagon/Hexagon/Heptagon: stubs & example registries (extend later)

Conventions (degrees by default):
  Si  = sin(theta_i),   Ci  = cos(theta_i)
  Sij = sin(alpha_ij),  Cij = cos(alpha_ij)

Public API (stable):
- evaluate_table(shape, table_key, inputs, degrees=True) -> list[dict]
- get_available_shapes() -> list[str]
- get_available_tables(shape) -> list[str]
- get_shape_examples(shape) -> dict[name, dict_inputs]
- solve_trig_Ac_Bs_D(A, B, D) -> list[float]  # radians
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple, Callable, Any
import math

# ------------------------------ helpers ------------------------------------

def d2r(x: float) -> float:
    return x * math.pi / 180.0

def r2d(x: float) -> float:
    return x * 180.0 / math.pi

def sin(x: float, degrees: bool) -> float:
    return math.sin(d2r(x) if degrees else x)

def cos(x: float, degrees: bool) -> float:
    return math.cos(d2r(x) if degrees else x)

def sanitize_inputs(inp: Dict[str, float]) -> Dict[str, float]:
    # allow both lowercase and uppercase keys, strip spaces
    canon = {}
    for k, v in inp.items():
        canon[k.strip().lower()] = v
    return canon

def Sij(alpha: float, degrees: bool) -> float:
    return sin(alpha, degrees)

def Cij(alpha: float, degrees: bool) -> float:
    return cos(alpha, degrees)

def Si(theta: float, degrees: bool) -> float:
    return sin(theta, degrees)

def Ci(theta: float, degrees: bool) -> float:
    return cos(theta, degrees)

def solve_trig_Ac_Bs_D(A: float, B: float, D: float) -> List[float]:
    """
    Solve A*cos(t) + B*sin(t) + D = 0 using t = 2*atan(u).
    Returns list of t in radians in [0, 2π).
    """
    # Using tan(t/2) substitution:
    # cos t = (1 - u^2)/(1 + u^2),  sin t = (2u)/(1 + u^2)
    # => A*(1-u^2)/(1+u^2) + B*(2u)/(1+u^2) + D = 0
    # => A(1-u^2) + 2Bu + D(1+u^2) = 0
    # => (D - A) u^2 + 2B u + (A + D) = 0
    a = (D - A)
    b = 2.0 * B
    c = (A + D)
    sols: List[float] = []
    if abs(a) < 1e-14 and abs(b) < 1e-14:
        # Degenerate: (A + D)=0 -> any t, else none.
        return sols
    if abs(a) < 1e-14:
        u = -c / b
        t = 2.0 * math.atan(u)
        sols.append(t % (2.0*math.pi))
        return sols
    disc = b*b - 4.0*a*c
    if disc < -1e-14:
        return sols
    disc = 0.0 if abs(disc) < 1e-14 else disc
    for s in (-math.sqrt(disc), math.sqrt(disc)) if disc > 0 else (0.0,):
        u = (-b + s) / (2.0*a)
        t = 2.0 * math.atan(u)
        sols.append(t % (2.0*math.pi))
    # Remove near-duplicates
    uniq = []
    for t in sols:
        if not any(abs(((t - x + math.pi)%(2*math.pi))-math.pi) < 1e-10 for x in uniq):
            uniq.append(t)
    return uniq

# --------------------------- Triangle ---------------------------------------

def triangle_blocks(inputs: Dict[str, float], degrees: bool=True) -> Dict[str, float]:
    """
    Compute basic blocks for the spherical triangle from Appendix definitions:
      X1=S23*S3, Y1=S23*C3, Z1=C23, etc.
      U12=S3*S23, V12=S3*C23, W12=C3, etc. (polar blocks / U,V,W)
    Requires: alpha12, alpha23, alpha31, theta1, theta2, theta3
    Keys accepted (case-insensitive): 'alpha12','alpha23','alpha31','theta1','theta2','theta3'
    """
    d = sanitize_inputs(inputs)
    a12 = d.get('alpha12'); a23 = d.get('alpha23'); a31 = d.get('alpha31')
    t1  = d.get('theta1');  t2  = d.get('theta2');  t3  = d.get('theta3')
    if None in (a12, a23, a31, t1, t2, t3):
        raise ValueError("Triangle requires alpha12, alpha23, alpha31, theta1, theta2, theta3")

    S12 = Sij(a12, degrees); C12 = Cij(a12, degrees)
    S23 = Sij(a23, degrees); C23 = Cij(a23, degrees)
    S31 = Sij(a31, degrees); C31 = Cij(a31, degrees)
    S1 = Si(t1, degrees);    C1 = Ci(t1, degrees)
    S2 = Si(t2, degrees);    C2 = Ci(t2, degrees)
    S3 = Si(t3, degrees);    C3 = Ci(t3, degrees)

    # Xi, Yi, Zi (Appendix lines; corrected to match “Equations for a Spherical Triangle”)
    # Top block (no hats):
    X1 = S23*S2;  Y1 = S23*C2;  Z1 = C23
    X2 = S31*S3;  Y2 = S31*C3;  Z2 = C31
    X3 = S12*S1;  Y3 = S12*C1;  Z3 = C12

    # “Hatted” block (the alternate ordering that appears as the second set in the table):
    X1_hat = S23*S3;  Y1_hat = S23*C3;  Z1_hat = C23
    X2_hat = S31*S1;  Y2_hat = S31*C1;  Z2_hat = C31
    X3_hat = S12*S2;  Y3_hat = S12*C2;  Z3_hat = C12


    # U,V,W (polar-style building blocks: U12=S3*S23, V12=S3*C23, W12=C3, etc.)
    U12 = S3*S23;   V12 = S3*C23;   W12 = C3
    U23 = S1*S31;   V23 = S1*C31;   W23 = C1
    U31 = S2*S12;   V31 = S2*C12;   W31 = C2
    # Also provide the “reversed” indices used in the DC tables:
    U21 = S3*S31;   V21 = S3*C31;   W21 = C3
    U32 = S1*S12;   V32 = S1*C12;   W32 = C1
    U13 = S2*S23;   V13 = S2*C23;   W13 = C2

    return {
        'S12': S12, 'C12': C12, 'S23': S23, 'C23': C23, 'S31': S31, 'C31': C31,
        'S1': S1, 'C1': C1, 'S2': S2, 'C2': C2, 'S3': S3, 'C3': C3,
        'X1': X1, 'Y1': Y1, 'Z1': Z1, 'X2': X2, 'Y2': Y2, 'Z2': Z2, 'X3': X3, 'Y3': Y3, 'Z3': Z3,
        'U12': U12, 'V12': V12, 'W12': W12, 'U23': U23, 'V23': V23, 'W23': W23, 'U31': U31, 'V31': V31, 'W31': W31,
    }

def table_triangle_fundamental(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = triangle_blocks(inputs, degrees)
    rows = []
    # First set (unhatted)
    rows += [
        {'eq': 'X1 = S23*S2', 'value': b['X1']},
        {'eq': 'Y1 = S23*C2', 'value': b['Y1']},
        {'eq': 'Z1 = C23',    'value': b['Z1']},
        {'eq': 'X2 = S31*S3', 'value': b['X2']},
        {'eq': 'Y2 = S31*C3', 'value': b['Y2']},
        {'eq': 'Z2 = C31',    'value': b['Z2']},
        {'eq': 'X3 = S12*S1', 'value': b['X3']},
        {'eq': 'Y3 = S12*C1', 'value': b['Y3']},
        {'eq': 'Z3 = C12',    'value': b['Z3']},
    ]
    # Second set (hatted)
    rows += [
        {'eq': 'X̂1 = S23*S3', 'value': b['X1_hat']},
        {'eq': 'Ŷ1 = S23*C3',  'value': b['Y1_hat']},
        {'eq': 'Ẑ1 = C23',     'value': b['Z1_hat']},
        {'eq': 'X̂2 = S31*S1', 'value': b['X2_hat']},
        {'eq': 'Ŷ2 = S31*C1',  'value': b['Y2_hat']},
        {'eq': 'Ẑ2 = C31',     'value': b['Z2_hat']},
        {'eq': 'X̂3 = S12*S2', 'value': b['X3_hat']},
        {'eq': 'Ŷ3 = S12*C2',  'value': b['Y3_hat']},
        {'eq': 'Ẑ3 = C12',     'value': b['Z3_hat']},
    ]
    return rows


def table_triangle_polar(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = triangle_blocks(inputs, degrees)
    rows = []
    rows.append({'eq': 'U12 = S3*S23', 'value': b['U12']})
    rows.append({'eq': 'V12 = S3*C23', 'value': b['V12']})
    rows.append({'eq': 'W12 = C3',      'value': b['W12']})
    rows.append({'eq': 'U23 = S1*S31', 'value': b['U23']})
    rows.append({'eq': 'V23 = S1*C31', 'value': b['V23']})
    rows.append({'eq': 'W23 = C1',      'value': b['W23']})
    rows.append({'eq': 'U31 = S2*S12', 'value': b['U31']})
    rows.append({'eq': 'V31 = S2*C12', 'value': b['V31']})
    rows.append({'eq': 'W31 = C2',      'value': b['W31']})
    return rows

def table_triangle_direction_cosines(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    # For UI purposes, we show the six sets as "building-block identities".
    # We list compactly: each "set" reports (Xi, Yi, Zi) that would be used to populate transformation rows.
    b = triangle_blocks(inputs, degrees)
    sets = []
    sets.append({'set': 'Set 1', 'Xi': b['X1'], 'Yi': b['Y1'], 'Zi': b['Z1']})
    sets.append({'set': 'Set 2', 'Xi': b['X2'], 'Yi': b['Y2'], 'Zi': b['Z2']})
    sets.append({'set': 'Set 3', 'Xi': b['X3'], 'Yi': b['Y3'], 'Zi': b['Z3']})
    # Also list the permuted X,Y definitions S23*S2 etc.
    # NOTE: Bars correspond to decreasing index order; only Y changes sign.
    sets.append(
        {'set': 'Set 4 (alt Xi=S23*S2,...)', 'Xi': b['S23'] * b['S2'], 'Yi': -b['S23'] * b['C2'], 'Zi': b['C23']})
    sets.append(
        {'set': 'Set 5 (alt X2=S31*S3,...)', 'Xi': b['S31'] * b['S3'], 'Yi': -b['S31'] * b['C3'], 'Zi': b['C31']})
    sets.append(
        {'set': 'Set 6 (alt X3=S12*S1,...)', 'Xi': b['S12'] * b['S1'], 'Yi': -b['S12'] * b['C1'], 'Zi': b['C12']})
    return sets

# ---------------------------- Pentagon --------------------------------------

def pentagon_blocks(inputs: Dict[str, float], degrees: bool=True) -> Dict[str, float]:
    """Minimal blocks for the spherical pentagon (indices 1..5 cyclic).
    We follow the same convention used for triangle/quadrilateral: we expose the
    fundamental identities that appear in the Appendix (no recursive composition
    is needed to *use* the tables).

    Requires: alpha12, alpha23, alpha34, alpha45, alpha51, theta1..theta5
    """
    d = sanitize_inputs(inputs)
    a12=d.get('alpha12'); a23=d.get('alpha23'); a34=d.get('alpha34'); a45=d.get('alpha45'); a51=d.get('alpha51')
    t1=d.get('theta1'); t2=d.get('theta2'); t3=d.get('theta3'); t4=d.get('theta4'); t5=d.get('theta5')
    if None in (a12,a23,a34,a45,a51,t1,t2,t3,t4,t5):
        raise ValueError("Pentagon requires alpha12, alpha23, alpha34, alpha45, alpha51, and theta1..theta5")

    S12 = Sij(a12, degrees); C12 = Cij(a12, degrees)
    S23 = Sij(a23, degrees); C23 = Cij(a23, degrees)
    S34 = Sij(a34, degrees); C34 = Cij(a34, degrees)
    S45 = Sij(a45, degrees); C45 = Cij(a45, degrees)
    S51 = Sij(a51, degrees); C51 = Cij(a51, degrees)
    S1 = Si(t1, degrees); C1 = Ci(t1, degrees)
    S2 = Si(t2, degrees); C2 = Ci(t2, degrees)
    S3 = Si(t3, degrees); C3 = Ci(t3, degrees)
    S4 = Si(t4, degrees); C4 = Ci(t4, degrees)
    S5 = Si(t5, degrees); C5 = Ci(t5, degrees)

    # Fundamental equalities (Appendix pattern):
    #   (X432, Y432, Z432) ↔ (S51*S1, S51*C1, C51)   etc. (cyclic rotations)
    X432 = S51*S1; Y432 = S51*C1; Z432 = C51
    X543 = S12*S2; Y543 = S12*C2; Z543 = C12
    X154 = S23*S3; Y154 = S23*C3; Z154 = C23
    X215 = S34*S4; Y215 = S34*C4; Z215 = C34
    X321 = S45*S5; Y321 = S45*C5; Z321 = C45

    return {
        'S12':S12,'C12':C12,'S23':S23,'C23':C23,'S34':S34,'C34':C34,'S45':S45,'C45':C45,'S51':S51,'C51':C51,
        'S1':S1,'C1':C1,'S2':S2,'C2':C2,'S3':S3,'C3':C3,'S4':S4,'C4':C4,'S5':S5,'C5':C5,
        'X432':X432,'Y432':Y432,'Z432':Z432,
        'X543':X543,'Y543':Y543,'Z543':Z543,
        'X154':X154,'Y154':Y154,'Z154':Z154,
        'X215':X215,'Y215':Y215,'Z215':Z215,
        'X321':X321,'Y321':Y321,'Z321':Z321,
    }

def table_penta_fundamental(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = pentagon_blocks(inputs, degrees)
    rows: List[Dict[str, Any]] = []
    rows += [
        {'eq': 'X432 = S51*S1', 'value': b['X432']},
        {'eq': 'Y432 = S51*C1', 'value': b['Y432']},
        {'eq': 'Z432 = C51',    'value': b['Z432']},

        {'eq': 'X543 = S12*S2', 'value': b['X543']},
        {'eq': 'Y543 = S12*C2', 'value': b['Y543']},
        {'eq': 'Z543 = C12',    'value': b['Z543']},

        {'eq': 'X154 = S23*S3', 'value': b['X154']},
        {'eq': 'Y154 = S23*C3', 'value': b['Y154']},
        {'eq': 'Z154 = C23',    'value': b['Z154']},

        {'eq': 'X215 = S34*S4', 'value': b['X215']},
        {'eq': 'Y215 = S34*C4', 'value': b['Y215']},
        {'eq': 'Z215 = C34',    'value': b['Z215']},

        {'eq': 'X321 = S45*S5', 'value': b['X321']},
        {'eq': 'Y321 = S45*C5', 'value': b['Y321']},
        {'eq': 'Z321 = C45',    'value': b['Z321']},
    ]
    return rows

def table_penta_polar(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    # By analogy with Quadrilateral (Zij = Ckl), the pentagon polar page
    # contains Z with three indices equated to the opposite edge's C.
    b = pentagon_blocks(inputs, degrees)
    rows = []
    rows.append({'eq': 'Z123 = C45', 'value': b['C45']})
    rows.append({'eq': 'Z234 = C51', 'value': b['C51']})
    rows.append({'eq': 'Z345 = C12', 'value': b['C12']})
    rows.append({'eq': 'Z451 = C23', 'value': b['C23']})
    rows.append({'eq': 'Z512 = C34', 'value': b['C34']})
    return rows

def table_penta_direction_cosines(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = pentagon_blocks(inputs, degrees)
    sets: List[Dict[str, Any]] = []
    # Report the five S-vectors expressed in the 1st frame per Appendix:
    sets.append({'set': 'S5 ≡ (X432, Y432, Z432)', 'X': b['X432'], 'Y': b['Y432'], 'Z': b['Z432']})
    sets.append({'set': 'S1 ≡ (X543, Y543, Z543)', 'X': b['X543'], 'Y': b['Y543'], 'Z': b['Z543']})
    sets.append({'set': 'S2 ≡ (X154, Y154, Z154)', 'X': b['X154'], 'Y': b['Y154'], 'Z': b['Z154']})
    sets.append({'set': 'S3 ≡ (X215, Y215, Z215)', 'X': b['X215'], 'Y': b['Y215'], 'Z': b['Z215']})
    sets.append({'set': 'S4 ≡ (X321, Y321, Z321)', 'X': b['X321'], 'Y': b['Y321'], 'Z': b['Z321']})
    return sets

def table_penta_half_tangent(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    """
    Half-Tangent Laws for a Pentagon.

    Set 1 (numeric evaluation):
        x1 =  X432 / (Y432 + S12) = -(Y432 - S12) / X432
        x2 =  X543 / (Y543 + S23) = -(Y543 - S23) / X543
        x3 =  X154 / (Y154 + S34) = -(Y154 - S34) / X154
        x4 =  X215 / (Y215 + S45) = -(Y215 - S45) / X215
        x5 =  X321 / (Y321 + S51) = -(Y321 - S51) / X321

    Set 2 (identity forms shown for reference) are algebraically equal to Set 1.
    We report the same numeric x_i values for Set 2 to aid cross-checking.
    """
    b = pentagon_blocks(inputs, degrees)

    def safe_div(n, d):
        return n/d if abs(d) > 1e-15 else float('nan')

    # --- Set 1 numerics ---
    x1 = safe_div(b['X432'], (b['Y432'] + b['S12']))
    x2 = safe_div(b['X543'], (b['Y543'] + b['S23']))
    x3 = safe_div(b['X154'], (b['Y154'] + b['S34']))
    x4 = safe_div(b['X215'], (b['Y215'] + b['S45']))
    x5 = safe_div(b['X321'], (b['Y321'] + b['S51']))

    rows: List[Dict[str, Any]] = []
    rows.append({'eq': 'Set 1: x1 = X432/(Y432 + s12) = -(Y432 - s12)/X432', 'value': x1})
    rows.append({'eq': 'Set 1: x2 = X543/(Y543 + s23) = -(Y543 - s23)/X543', 'value': x2})
    rows.append({'eq': 'Set 1: x3 = X154/(Y154 + s34) = -(Y154 - s34)/X154', 'value': x3})
    rows.append({'eq': 'Set 1: x4 = X215/(Y215 + s45) = -(Y215 - s45)/X215', 'value': x4})
    rows.append({'eq': 'Set 1: x5 = X321/(Y321 + s51) = -(Y321 - s51)/X321', 'value': x5})

    # --- Set 2 identities (equal numerically to Set 1) ---
    rows.append({'eq': 'Set 2: x1 (identity; = Set 1)', 'value': x1})
    rows.append({'eq': 'Set 2: x2 (identity; = Set 1)', 'value': x2})
    rows.append({'eq': 'Set 2: x3 (identity; = Set 1)', 'value': x3})
    rows.append({'eq': 'Set 2: x4 (identity; = Set 1)', 'value': x4})
    rows.append({'eq': 'Set 2: x5 (identity; = Set 1)', 'value': x5})

    return rows

# ---------------------------- Hexagon ---------------------------------------

def hexagon_blocks(inputs: Dict[str, float], degrees: bool=True) -> Dict[str, float]:
    """
    Minimal blocks for the spherical hexagon (indices 1..6 cyclic).
    Pattern mirrors the pentagon implementation.

    Requires: alpha12, alpha23, alpha34, alpha45, alpha56, alpha61, theta1..theta6
    """
    d = sanitize_inputs(inputs)
    a12=d.get('alpha12'); a23=d.get('alpha23'); a34=d.get('alpha34')
    a45=d.get('alpha45'); a56=d.get('alpha56'); a61=d.get('alpha61')
    t1=d.get('theta1'); t2=d.get('theta2'); t3=d.get('theta3')
    t4=d.get('theta4'); t5=d.get('theta5'); t6=d.get('theta6')
    if None in (a12,a23,a34,a45,a56,a61,t1,t2,t3,t4,t5,t6):
        raise ValueError("Hexagon requires alpha12, alpha23, alpha34, alpha45, alpha56, alpha61 and theta1..theta6")

    S12 = Sij(a12, degrees); C12 = Cij(a12, degrees)
    S23 = Sij(a23, degrees); C23 = Cij(a23, degrees)
    S34 = Sij(a34, degrees); C34 = Cij(a34, degrees)
    S45 = Sij(a45, degrees); C45 = Cij(a45, degrees)
    S56 = Sij(a56, degrees); C56 = Cij(a56, degrees)
    S61 = Sij(a61, degrees); C61 = Cij(a61, degrees)
    S1 = Si(t1, degrees); C1 = Ci(t1, degrees)
    S2 = Si(t2, degrees); C2 = Ci(t2, degrees)
    S3 = Si(t3, degrees); C3 = Ci(t3, degrees)
    S4 = Si(t4, degrees); C4 = Ci(t4, degrees)
    S5 = Si(t5, degrees); C5 = Ci(t5, degrees)
    S6 = Si(t6, degrees); C6 = Ci(t6, degrees)

    # Forward cycle
    X1234 = S56*S5; Y1234 = S56*C5; Z1234 = C56
    X2345 = S61*S6; Y2345 = S61*C6; Z2345 = C61
    X3456 = S12*S1; Y3456 = S12*C1; Z3456 = C12
    X4561 = S23*S2; Y4561 = S23*C2; Z4561 = C23
    X5612 = S34*S3; Y5612 = S34*C3; Z5612 = C34
    X6123 = S45*S4; Y6123 = S45*C4; Z6123 = C45

    # Reverse cycle (for the DC table pattern)
    X4321 = S56*S6; Y4321 = S56*C6; Z4321 = C56
    X5432 = S61*S1; Y5432 = S61*C1; Z5432 = C61
    X6543 = S12*S2; Y6543 = S12*C2; Z6543 = C12
    X1654 = S23*S3; Y1654 = S23*C3; Z1654 = C23
    X2165 = S34*S4; Y2165 = S34*C4; Z2165 = C34
    X3216 = S45*S5; Y3216 = S45*C5; Z3216 = C45

    return {
        'S12':S12,'C12':C12,'S23':S23,'C23':C23,'S34':S34,'C34':C34,
        'S45':S45,'C45':C45,'S56':S56,'C56':C56,'S61':S61,'C61':C61,
        'S1':S1,'C1':C1,'S2':S2,'C2':C2,'S3':S3,'C3':C3,'S4':S4,'C4':C4,'S5':S5,'C5':C5,'S6':S6,'C6':C6,
        'X1234':X1234,'Y1234':Y1234,'Z1234':Z1234,
        'X2345':X2345,'Y2345':Y2345,'Z2345':Z2345,
        'X3456':X3456,'Y3456':Y3456,'Z3456':Z3456,
        'X4561':X4561,'Y4561':Y4561,'Z4561':Z4561,
        'X5612':X5612,'Y5612':Y5612,'Z5612':Z5612,
        'X6123':X6123,'Y6123':Y6123,'Z6123':Z6123,
        'X4321':X4321,'Y4321':Y4321,'Z4321':Z4321,
        'X5432':X5432,'Y5432':Y5432,'Z5432':Z5432,
        'X6543':X6543,'Y6543':Y6543,'Z6543':Z6543,
        'X1654':X1654,'Y1654':Y1654,'Z1654':Z1654,
        'X2165':X2165,'Y2165':Y2165,'Z2165':Z2165,
        'X3216':X3216,'Y3216':Y3216,'Z3216':Z3216,
    }

def table_hexa_fundamental(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = hexagon_blocks(inputs, degrees)
    rows: List[Dict[str, Any]] = []
    # Forward cycle
    rows += [
        {'eq': 'X1234 = S56*S5', 'value': b['X1234']},
        {'eq': 'Y1234 = S56*C5', 'value': b['Y1234']},
        {'eq': 'Z1234 = C56',    'value': b['Z1234']},
        {'eq': 'X2345 = S61*S6', 'value': b['X2345']},
        {'eq': 'Y2345 = S61*C6', 'value': b['Y2345']},
        {'eq': 'Z2345 = C61',    'value': b['Z2345']},
        {'eq': 'X3456 = S12*S1', 'value': b['X3456']},
        {'eq': 'Y3456 = S12*C1', 'value': b['Y3456']},
        {'eq': 'Z3456 = C12',    'value': b['Z3456']},
        {'eq': 'X4561 = S23*S2', 'value': b['X4561']},
        {'eq': 'Y4561 = S23*C2', 'value': b['Y4561']},
        {'eq': 'Z4561 = C23',    'value': b['Z4561']},
        {'eq': 'X5612 = S34*S3', 'value': b['X5612']},
        {'eq': 'Y5612 = S34*C3', 'value': b['Y5612']},
        {'eq': 'Z5612 = C34',    'value': b['Z5612']},
        {'eq': 'X6123 = S45*S4', 'value': b['X6123']},
        {'eq': 'Y6123 = S45*C4', 'value': b['Y6123']},
        {'eq': 'Z6123 = C45',    'value': b['Z6123']},
    ]
    # Reverse cycle
    rows += [
        {'eq': 'X4321 = S56*S6', 'value': b['X4321']},
        {'eq': 'Y4321 = S56*C6', 'value': b['Y4321']},
        {'eq': 'Z4321 = C56',    'value': b['Z4321']},
        {'eq': 'X5432 = S61*S1', 'value': b['X5432']},
        {'eq': 'Y5432 = S61*C1', 'value': b['Y5432']},
        {'eq': 'Z5432 = C61',    'value': b['Z5432']},
        {'eq': 'X6543 = S12*S2', 'value': b['X6543']},
        {'eq': 'Y6543 = S12*C2', 'value': b['Y6543']},
        {'eq': 'Z6543 = C12',    'value': b['Z6543']},
        {'eq': 'X1654 = S23*S3', 'value': b['X1654']},
        {'eq': 'Y1654 = S23*C3', 'value': b['Y1654']},
        {'eq': 'Z1654 = C23',    'value': b['Z1654']},
        {'eq': 'X2165 = S34*S4', 'value': b['X2165']},
        {'eq': 'Y2165 = S34*C4', 'value': b['Y2165']},
        {'eq': 'Z2165 = C34',    'value': b['Z2165']},
        {'eq': 'X3216 = S45*S5', 'value': b['X3216']},
        {'eq': 'Y3216 = S45*C5', 'value': b['Y3216']},
        {'eq': 'Z3216 = C45',    'value': b['Z3216']},
    ]
    return rows

def table_hexa_polar(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = hexagon_blocks(inputs, degrees)
    rows: List[Dict[str, Any]] = []
    # Forward
    rows += [
        {'eq': 'U12345 = S6*S56', 'value': b['S6']*b['S56']},
        {'eq': 'V12345 = S6*C56', 'value': b['S6']*b['C56']},
        {'eq': 'W12345 = C6',     'value': b['C6']},
        {'eq': 'U23456 = S1*S61', 'value': b['S1']*b['S61']},
        {'eq': 'V23456 = S1*C61', 'value': b['S1']*b['C61']},
        {'eq': 'W23456 = C1',     'value': b['C1']},
        {'eq': 'U34561 = S2*S12', 'value': b['S2']*b['S12']},
        {'eq': 'V34561 = S2*C12', 'value': b['S2']*b['C12']},
        {'eq': 'W34561 = C2',     'value': b['C2']},
        {'eq': 'U45612 = S3*S23', 'value': b['S3']*b['S23']},
        {'eq': 'V45612 = S3*C23', 'value': b['S3']*b['C23']},
        {'eq': 'W45612 = C3',     'value': b['C3']},
        {'eq': 'U56123 = S4*S34', 'value': b['S4']*b['S34']},
        {'eq': 'V56123 = S4*C34', 'value': b['S4']*b['C34']},
        {'eq': 'W56123 = C4',     'value': b['C4']},
        {'eq': 'U61234 = S5*S45', 'value': b['S5']*b['S45']},
        {'eq': 'V61234 = S5*C45', 'value': b['S5']*b['C45']},
        {'eq': 'W61234 = C5',     'value': b['C5']},
    ]
    # Reverse
    rows += [
        {'eq': 'U54321 = S6*S61', 'value': b['S6']*b['S61']},
        {'eq': 'V54321 = S6*C61', 'value': b['S6']*b['C61']},
        {'eq': 'W54321 = C6',     'value': b['C6']},
        {'eq': 'U43216 = S5*S56', 'value': b['S5']*b['S56']},
        {'eq': 'V43216 = S5*C56', 'value': b['S5']*b['C56']},
        {'eq': 'W43216 = C5',     'value': b['C5']},
        {'eq': 'U32165 = S4*S45', 'value': b['S4']*b['S45']},
        {'eq': 'V32165 = S4*C45', 'value': b['S4']*b['C45']},
        {'eq': 'W32165 = C4',     'value': b['C4']},
        {'eq': 'U21654 = S3*S34', 'value': b['S3']*b['S34']},
        {'eq': 'V21654 = S3*C34', 'value': b['S3']*b['C34']},
        {'eq': 'W21654 = C3',     'value': b['C3']},
        {'eq': 'U16543 = S2*S23', 'value': b['S2']*b['S23']},
        {'eq': 'V16543 = S2*C23', 'value': b['S2']*b['C23']},
        {'eq': 'W16543 = C2',     'value': b['C2']},
        {'eq': 'U65432 = S1*S12', 'value': b['S1']*b['S12']},
        {'eq': 'V65432 = S1*C12', 'value': b['S1']*b['C12']},
        {'eq': 'W65432 = C1',     'value': b['C1']},
    ]
    return rows

def table_hexa_half_tangent(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    """Half-Tangent Laws for a Hexagon (Set 1 evaluated)."""
    b = hexagon_blocks(inputs, degrees)
    def _sd(n,d):
        return n/d if abs(d) > 1e-15 else float('nan')
    x1 = _sd(b['X3456'], (b['Y3456'] + b['S12']))
    x2 = _sd(b['X4561'], (b['Y4561'] + b['S23']))
    x3 = _sd(b['X5612'], (b['Y5612'] + b['S34']))
    x4 = _sd(b['X6123'], (b['Y6123'] + b['S45']))
    x5 = _sd(b['X1234'], (b['Y1234'] + b['S56']))
    x6 = _sd(b['X2345'], (b['Y2345'] + b['S61']))
    return [
        {'eq': 'Set 1: x1 = X3456/(Y3456 + s12) = -(Y3456 - s12)/X3456', 'value': x1},
        {'eq': 'Set 1: x2 = X4561/(Y4561 + s23) = -(Y4561 - s23)/X4561', 'value': x2},
        {'eq': 'Set 1: x3 = X5612/(Y5612 + s34) = -(Y5612 - s34)/X5612', 'value': x3},
        {'eq': 'Set 1: x4 = X6123/(Y6123 + s45) = -(Y6123 - s45)/X6123', 'value': x4},
        {'eq': 'Set 1: x5 = X1234/(Y1234 + s56) = -(Y1234 - s56)/X1234', 'value': x5},
        {'eq': 'Set 1: x6 = X2345/(Y2345 + s61) = -(Y2345 - s61)/X2345', 'value': x6},
    ]

def table_hexa_direction_cosines(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    """Compact DC report mirroring the pentagon version: list S-vectors in the 1st frame."""
    b = hexagon_blocks(inputs, degrees)
    sets: List[Dict[str, Any]] = []
    sets.append({'set': 'S6 ≡ (X5432, Y5432, Z5432)', 'X': b['X5432'], 'Y': b['Y5432'], 'Z': b['Z5432']})
    sets.append({'set': 'S1 ≡ (X6543, Y6543, Z6543)', 'X': b['X6543'], 'Y': b['Y6543'], 'Z': b['Z6543']})
    sets.append({'set': 'S2 ≡ (X1654, Y1654, Z1654)', 'X': b['X1654'], 'Y': b['Y1654'], 'Z': b['Z1654']})
    sets.append({'set': 'S3 ≡ (X2165, Y2165, Z2165)', 'X': b['X2165'], 'Y': b['Y2165'], 'Z': b['Z2165']})
    sets.append({'set': 'S4 ≡ (X3216, Y3216, Z3216)', 'X': b['X3216'], 'Y': b['Y3216'], 'Z': b['Z3216']})
    sets.append({'set': 'S5 ≡ (X4321, Y4321, Z4321)', 'X': b['X4321'], 'Y': b['Y4321'], 'Z': b['Z4321']})
    return sets

# ---------------------------- Heptagon --------------------------------------

def heptagon_blocks(inputs: Dict[str, float], degrees: bool=True) -> Dict[str, float]:
    """Minimal blocks for a spherical heptagon (indices 1..7 cyclic).
    Computes the Sij/Cij pairs and the X,Y,Z (spherical) and U,V,W (polar) fundamentals
    that appear across the Appendix tables.
    Requires: alpha12, alpha23, alpha34, alpha45, alpha56, alpha67, alpha71, theta1..theta7
    """
    d = sanitize_inputs(inputs)
    a12=d.get('alpha12'); a23=d.get('alpha23'); a34=d.get('alpha34'); a45=d.get('alpha45')
    a56=d.get('alpha56'); a67=d.get('alpha67'); a71=d.get('alpha71')
    t1=d.get('theta1'); t2=d.get('theta2'); t3=d.get('theta3'); t4=d.get('theta4')
    t5=d.get('theta5'); t6=d.get('theta6'); t7=d.get('theta7')
    req = (a12,a23,a34,a45,a56,a67,a71,t1,t2,t3,t4,t5,t6,t7)
    if None in req:
        raise ValueError("Heptagon requires alpha12, alpha23, alpha34, alpha45, alpha56, alpha67, alpha71 and theta1..theta7")

    # Sij/Cij
    S12=Sij(a12,degrees); C12=Cij(a12,degrees)
    S23=Sij(a23,degrees); C23=Cij(a23,degrees)
    S34=Sij(a34,degrees); C34=Cij(a34,degrees)
    S45=Sij(a45,degrees); C45=Cij(a45,degrees)
    S56=Sij(a56,degrees); C56=Cij(a56,degrees)
    S67=Sij(a67,degrees); C67=Cij(a67,degrees)
    S71=Sij(a71,degrees); C71=Cij(a71,degrees)

    # Si/Ci
    S1=Si(t1,degrees); C1=Ci(t1,degrees)
    S2=Si(t2,degrees); C2=Ci(t2,degrees)
    S3=Si(t3,degrees); C3=Ci(t3,degrees)
    S4=Si(t4,degrees); C4=Ci(t4,degrees)
    S5=Si(t5,degrees); C5=Ci(t5,degrees)
    S6=Si(t6,degrees); C6=Ci(t6,degrees)
    S7=Si(t7,degrees); C7=Ci(t7,degrees)

    # --- Spherical (X,Y,Z) forward 7-tuple ---
    X12345=S67*S6; Y12345=S67*C6; Z12345=C67
    X23456=S71*S7; Y23456=S71*C7; Z23456=C71
    X34567=S12*S1; Y34567=S12*C1; Z34567=C12
    X45671=S23*S2; Y45671=S23*C2; Z45671=C23
    X56712=S34*S3; Y56712=S34*C3; Z56712=C34
    X67123=S45*S4; Y67123=S45*C4; Z67123=C45
    X71234=S56*S5; Y71234=S56*C5; Z71234=C56

    # Reverse cycle (used by DC/variants)
    X54321=S67*S7; Y54321=S67*C7; Z54321=C67
    X65432=S71*S1; Y65432=S71*C1; Z65432=C71
    X76543=S12*S2; Y76543=S12*C2; Z76543=C12
    X17654=S23*S3; Y17654=S23*C3; Z17654=C23
    X21765=S34*S4; Y21765=S34*C4; Z21765=C34
    X32176=S45*S5; Y32176=S45*C5; Z32176=C45
    X43217=S56*S6; Y43217=S56*C6; Z43217=C56

    return {
        'S12':S12,'C12':C12,'S23':S23,'C23':C23,'S34':S34,'C34':C34,'S45':S45,'C45':C45,
        'S56':S56,'C56':C56,'S67':S67,'C67':C67,'S71':S71,'C71':C71,
        'S1':S1,'C1':C1,'S2':S2,'C2':C2,'S3':S3,'C3':C3,'S4':S4,'C4':C4,'S5':S5,'C5':C5,'S6':S6,'C6':C6,'S7':S7,'C7':C7,
        'X12345':X12345,'Y12345':Y12345,'Z12345':Z12345,
        'X23456':X23456,'Y23456':Y23456,'Z23456':Z23456,
        'X34567':X34567,'Y34567':Y34567,'Z34567':Z34567,
        'X45671':X45671,'Y45671':Y45671,'Z45671':Z45671,
        'X56712':X56712,'Y56712':Y56712,'Z56712':Z56712,
        'X67123':X67123,'Y67123':Y67123,'Z67123':Z67123,
        'X71234':X71234,'Y71234':Y71234,'Z71234':Z71234,
        'X54321':X54321,'Y54321':Y54321,'Z54321':Z54321,
        'X65432':X65432,'Y65432':Y65432,'Z65432':Z65432,
        'X76543':X76543,'Y76543':Y76543,'Z76543':Z76543,
        'X17654':X17654,'Y17654':Y17654,'Z17654':Z17654,
        'X21765':X21765,'Y21765':Y21765,'Z21765':Z21765,
        'X32176':X32176,'Y32176':Y32176,'Z32176':Z32176,
        'X43217':X43217,'Y43217':Y43217,'Z43217':Z43217,
    }

def table_hepta_fundamental(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = heptagon_blocks(inputs, degrees)
    rows: List[Dict[str, Any]] = []
    # Forward
    rows += [
        {'eq': 'X12345 = S67*S6', 'value': b['X12345']},
        {'eq': 'Y12345 = S67*C6', 'value': b['Y12345']},
        {'eq': 'Z12345 = C67',    'value': b['Z12345']},
        {'eq': 'X23456 = S71*S7', 'value': b['X23456']},
        {'eq': 'Y23456 = S71*C7', 'value': b['Y23456']},
        {'eq': 'Z23456 = C71',    'value': b['Z23456']},
        {'eq': 'X34567 = S12*S1', 'value': b['X34567']},
        {'eq': 'Y34567 = S12*C1', 'value': b['Y34567']},
        {'eq': 'Z34567 = C12',    'value': b['Z34567']},
        {'eq': 'X45671 = S23*S2', 'value': b['X45671']},
        {'eq': 'Y45671 = S23*C2', 'value': b['Y45671']},
        {'eq': 'Z45671 = C23',    'value': b['Z45671']},
        {'eq': 'X56712 = S34*S3', 'value': b['X56712']},
        {'eq': 'Y56712 = S34*C3', 'value': b['Y56712']},
        {'eq': 'Z56712 = C34',    'value': b['Z56712']},
        {'eq': 'X67123 = S45*S4', 'value': b['X67123']},
        {'eq': 'Y67123 = S45*C4', 'value': b['Y67123']},
        {'eq': 'Z67123 = C45',    'value': b['Z67123']},
        {'eq': 'X71234 = S56*S5', 'value': b['X71234']},
        {'eq': 'Y71234 = S56*C5', 'value': b['Y71234']},
        {'eq': 'Z71234 = C56',    'value': b['Z71234']},
    ]
    # Reverse
    rows += [
        {'eq': 'X54321 = S67*S7', 'value': b['X54321']},
        {'eq': 'Y54321 = S67*C7', 'value': b['Y54321']},
        {'eq': 'Z54321 = C67',    'value': b['Z54321']},
        {'eq': 'X65432 = S71*S1', 'value': b['X65432']},
        {'eq': 'Y65432 = S71*C1', 'value': b['Y65432']},
        {'eq': 'Z65432 = C71',    'value': b['Z65432']},
        {'eq': 'X76543 = S12*S2', 'value': b['X76543']},
        {'eq': 'Y76543 = S12*C2', 'value': b['Y76543']},
        {'eq': 'Z76543 = C12',    'value': b['Z76543']},
        {'eq': 'X17654 = S23*S3', 'value': b['X17654']},
        {'eq': 'Y17654 = S23*C3', 'value': b['Y17654']},
        {'eq': 'Z17654 = C23',    'value': b['Z17654']},
        {'eq': 'X21765 = S34*S4', 'value': b['X21765']},
        {'eq': 'Y21765 = S34*C4', 'value': b['Y21765']},
        {'eq': 'Z21765 = C34',    'value': b['Z21765']},
        {'eq': 'X32176 = S45*S5', 'value': b['X32176']},
        {'eq': 'Y32176 = S45*C5', 'value': b['Y32176']},
        {'eq': 'Z32176 = C45',    'value': b['Z32176']},
        {'eq': 'X43217 = S56*S6', 'value': b['X43217']},
        {'eq': 'Y43217 = S56*C6', 'value': b['Y43217']},
        {'eq': 'Z43217 = C56',    'value': b['Z43217']},
    ]
    return rows

def table_hepta_polar(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = heptagon_blocks(inputs, degrees)
    rows: List[Dict[str, Any]] = []
    # Forward
    rows += [
        {'eq': 'U123456 = S7*S67', 'value': b['S7']*b['S67']},
        {'eq': 'V123456 = S7*C67', 'value': b['S7']*b['C67']},
        {'eq': 'W123456 = C7',     'value': b['C7']},
        {'eq': 'U234567 = S1*S71', 'value': b['S1']*b['S71']},
        {'eq': 'V234567 = S1*C71', 'value': b['S1']*b['C71']},
        {'eq': 'W234567 = C1',     'value': b['C1']},
        {'eq': 'U345671 = S2*S12', 'value': b['S2']*b['S12']},
        {'eq': 'V345671 = S2*C12', 'value': b['S2']*b['C12']},
        {'eq': 'W345671 = C2',     'value': b['C2']},
        {'eq': 'U456712 = S3*S23', 'value': b['S3']*b['S23']},
        {'eq': 'V456712 = S3*C23', 'value': b['S3']*b['C23']},
        {'eq': 'W456712 = C3',     'value': b['C3']},
        {'eq': 'U567123 = S4*S34', 'value': b['S4']*b['S34']},
        {'eq': 'V567123 = S4*C34', 'value': b['S4']*b['C34']},
        {'eq': 'W567123 = C4',     'value': b['C4']},
        {'eq': 'U671234 = S5*S45', 'value': b['S5']*b['S45']},
        {'eq': 'V671234 = S5*C45', 'value': b['S5']*b['C45']},
        {'eq': 'W671234 = C5',     'value': b['C5']},
        {'eq': 'U712345 = S6*S56', 'value': b['S6']*b['S56']},
        {'eq': 'V712345 = S6*C56', 'value': b['S6']*b['C56']},
        {'eq': 'W712345 = C6',     'value': b['C6']},
    ]
    # Reverse
    rows += [
        {'eq': 'U654321 = S7*S71', 'value': b['S7']*b['S71']},
        {'eq': 'V654321 = S7*C71', 'value': b['S7']*b['C71']},
        {'eq': 'W654321 = C7',     'value': b['C7']},
        {'eq': 'U765432 = S1*S12', 'value': b['S1']*b['S12']},
        {'eq': 'V765432 = S1*C12', 'value': b['S1']*b['C12']},
        {'eq': 'W765432 = C1',     'value': b['C1']},
        {'eq': 'U176543 = S2*S23', 'value': b['S2']*b['S23']},
        {'eq': 'V176543 = S2*C23', 'value': b['S2']*b['C23']},
        {'eq': 'W176543 = C2',     'value': b['C2']},
        {'eq': 'U217654 = S3*S34', 'value': b['S3']*b['S34']},
        {'eq': 'V217654 = S3*C34', 'value': b['S3']*b['C34']},
        {'eq': 'W217654 = C3',     'value': b['C3']},
        {'eq': 'U321765 = S4*S45', 'value': b['S4']*b['S45']},
        {'eq': 'V321765 = S4*C45', 'value': b['S4']*b['C45']},
        {'eq': 'W321765 = C4',     'value': b['C4']},
        {'eq': 'U432176 = S5*S56', 'value': b['S5']*b['S56']},
        {'eq': 'V432176 = S5*C56', 'value': b['S5']*b['C56']},
        {'eq': 'W432176 = C5',     'value': b['C5']},
        {'eq': 'U543217 = S6*S67', 'value': b['S6']*b['S67']},
        {'eq': 'V543217 = S6*C67', 'value': b['S6']*b['C67']},
        {'eq': 'W543217 = C6',     'value': b['C6']},
    ]
    return rows

def table_hepta_half_tangent(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    """Half-Tangent Laws for a Heptagon.
    We evaluate Set 1 and report Sets 2–3 as identities (same numeric values).
    x_i = X_next / (Y_next + s_(i,i+1))  with the appropriate 5-index X blocks.
    """
    b = heptagon_blocks(inputs, degrees)
    # Set 1
    x1 = safe_div(b['X34567'], (b['Y34567'] + b['S12']))
    x2 = safe_div(b['X45671'], (b['Y45671'] + b['S23']))
    x3 = safe_div(b['X56712'], (b['Y56712'] + b['S34']))
    x4 = safe_div(b['X67123'], (b['Y67123'] + b['S45']))
    x5 = safe_div(b['X71234'], (b['Y71234'] + b['S56']))
    x6 = safe_div(b['X12345'], (b['Y12345'] + b['S67']))
    x7 = safe_div(b['X23456'], (b['Y23456'] + b['S71']))
    rows: List[Dict[str, Any]] = []
    rows += [
        {'eq': 'Set 1: x1 = X34567/(Y34567 + s12) = -(Y34567 - s12)/X34567', 'value': x1},
        {'eq': 'Set 1: x2 = X45671/(Y45671 + s23) = -(Y45671 - s23)/X45671', 'value': x2},
        {'eq': 'Set 1: x3 = X56712/(Y56712 + s34) = -(Y56712 - s34)/X56712', 'value': x3},
        {'eq': 'Set 1: x4 = X67123/(Y67123 + s45) = -(Y67123 - s45)/X67123', 'value': x4},
        {'eq': 'Set 1: x5 = X71234/(Y71234 + s56) = -(Y71234 - s56)/X71234', 'value': x5},
        {'eq': 'Set 1: x6 = X12345/(Y12345 + s67) = -(Y12345 - s67)/X12345', 'value': x6},
        {'eq': 'Set 1: x7 = X23456/(Y23456 + s71) = -(Y23456 - s71)/X23456', 'value': x7},
    ]
    # Sets 2–3 (barred/hat) equalities
    rows += [
        {'eq': 'Set 2: x1 (identity; = Set 1)', 'value': x1},
        {'eq': 'Set 2: x2 (identity; = Set 1)', 'value': x2},
        {'eq': 'Set 2: x3 (identity; = Set 1)', 'value': x3},
        {'eq': 'Set 2: x4 (identity; = Set 1)', 'value': x4},
        {'eq': 'Set 2: x5 (identity; = Set 1)', 'value': x5},
        {'eq': 'Set 2: x6 (identity; = Set 1)', 'value': x6},
        {'eq': 'Set 2: x7 (identity; = Set 1)', 'value': x7},
        {'eq': 'Set 3: x1 (identity; = Set 1)', 'value': x1},
        {'eq': 'Set 3: x2 (identity; = Set 1)', 'value': x2},
        {'eq': 'Set 3: x3 (identity; = Set 1)', 'value': x3},
        {'eq': 'Set 3: x4 (identity; = Set 1)', 'value': x4},
        {'eq': 'Set 3: x5 (identity; = Set 1)', 'value': x5},
        {'eq': 'Set 3: x6 (identity; = Set 1)', 'value': x6},
        {'eq': 'Set 3: x7 (identity; = Set 1)', 'value': x7},
    ]
    return rows

def table_hepta_direction_cosines(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    """A compact DC summary: report the seven S_i unit vectors
    using the (X, Y, Z) blocks consistent with the Heptagon spherical formulas.
    This mirrors the compact DC reporters for the other polygons.
    """
    b = heptagon_blocks(inputs, degrees)
    sets = []
    sets.append({'set': 'S1 ≡ (X23456, Y23456, Z23456)', 'X': b['X23456'], 'Y': b['Y23456'], 'Z': b['Z23456']})
    sets.append({'set': 'S2 ≡ (X34567, Y34567, Z34567)', 'X': b['X34567'], 'Y': b['Y34567'], 'Z': b['Z34567']})
    sets.append({'set': 'S3 ≡ (X45671, Y45671, Z45671)', 'X': b['X45671'], 'Y': b['Y45671'], 'Z': b['Z45671']})
    sets.append({'set': 'S4 ≡ (X56712, Y56712, Z56712)', 'X': b['X56712'], 'Y': b['Y56712'], 'Z': b['Z56712']})
    sets.append({'set': 'S5 ≡ (X67123, Y67123, Z67123)', 'X': b['X67123'], 'Y': b['Y67123'], 'Z': b['Z67123']})
    sets.append({'set': 'S6 ≡ (X71234, Y71234, Z71234)', 'X': b['X71234'], 'Y': b['Y71234'], 'Z': b['Z71234']})
    sets.append({'set': 'S7 ≡ (X12345, Y12345, Z12345)', 'X': b['X12345'], 'Y': b['Y12345'], 'Z': b['Z12345']})
    return sets

# --------------------------- Quadrilateral ----------------------------------

def quad_blocks(inputs: Dict[str, float], degrees: bool=True) -> Dict[str, float]:
    """
    Minimal blocks for the spherical quadrilateral:
      Fundamental laws include identities like X32 = S41*S1, Y32=S41*C1, Z32=C41, etc.
    Requires: alpha12, alpha23, alpha34, alpha41, theta1, theta2, theta3, theta4
    """
    d = sanitize_inputs(inputs)
    a12 = d.get('alpha12'); a23 = d.get('alpha23'); a34 = d.get('alpha34'); a41 = d.get('alpha41')
    t1  = d.get('theta1');  t2 = d.get('theta2');   t3  = d.get('theta3');  t4  = d.get('theta4')
    if None in (a12, a23, a34, a41, t1, t2, t3, t4):
        raise ValueError("Quadrilateral requires alpha12, alpha23, alpha34, alpha41, theta1..theta4")

    S12 = Sij(a12, degrees); C12 = Cij(a12, degrees)
    S23 = Sij(a23, degrees); C23 = Cij(a23, degrees)
    S34 = Sij(a34, degrees); C34 = Cij(a34, degrees)
    S41 = Sij(a41, degrees); C41 = Cij(a41, degrees)
    S1 = Si(t1, degrees);    C1 = Ci(t1, degrees)
    S2 = Si(t2, degrees);    C2 = Ci(t2, degrees)
    S3 = Si(t3, degrees);    C3 = Ci(t3, degrees)
    S4 = Si(t4, degrees);    C4 = Ci(t4, degrees)

    # A small dictionary of commonly referenced Xij, Yij, Zij combinations
    # (Eight directed edges — both senses)
    X12 = S34*S3; Y12 = S34*C3; Z12 = C34
    X23 = S41*S4; Y23 = S41*C4; Z23 = C41
    X34 = S12*S1; Y34 = S12*C1; Z34 = C12
    X41 = S23*S2; Y41 = S23*C2; Z41 = C23
    X21 = S34*S4; Y21 = S34*C4; Z21 = C34
    X32 = S41*S1; Y32 = S41*C1; Z32 = C41
    X43 = S12*S2; Y43 = S12*C2; Z43 = C12
    X14 = S23*S3; Y14 = S23*C3; Z14 = C23

    return {
        'S12': S12, 'C12': C12, 'S23': S23, 'C23': C23, 'S34': S34, 'C34': C34, 'S41': S41, 'C41': C41,
        'S1': S1, 'C1': C1, 'S2': S2, 'C2': C2, 'S3': S3, 'C3': C3, 'S4': S4, 'C4': C4,
        'X12': X12, 'Y12': Y12, 'Z12': Z12,
        'X23': X23, 'Y23': Y23, 'Z23': Z23,
        'X34': X34, 'Y34': Y34, 'Z34': Z34,
        'X41': X41, 'Y41': Y41, 'Z41': Z41,
        'X21': X21, 'Y21': Y21, 'Z21': Z21,
        'X32': X32, 'Y32': Y32, 'Z32': Z32,
        'X43': X43, 'Y43': Y43, 'Z43': Z43,
        'X14': X14, 'Y14': Y14, 'Z14': Z14,
    }

def table_quad_fundamental(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = quad_blocks(inputs, degrees)
    rows = []
    # One of the eight fundamental sets (others by index exchange). We showcase four by symmetry.
    rows.append({'eq': 'X32 = S41*S1', 'value': b['X32']})
    rows.append({'eq': 'Y32 = S41*C1', 'value': b['Y32']})
    rows.append({'eq': 'Z32 = C41',     'value': b['Z32']})
    rows.append({'eq': 'X41 = S12*S2', 'value': b['X41']})
    rows.append({'eq': 'Y41 = S12*C2', 'value': b['Y41']})
    rows.append({'eq': 'Z41 = C12',     'value': b['Z41']})
    rows.append({'eq': 'X23 = S34*S4', 'value': b['X23']})
    rows.append({'eq': 'Y23 = S34*C4', 'value': b['Y23']})
    rows.append({'eq': 'Z23 = C34',     'value': b['Z23']})
    rows.append({'eq': 'X14 = S23*S3', 'value': b['X14']})
    rows.append({'eq': 'Y14 = S23*C3', 'value': b['Y14']})
    rows.append({'eq': 'Z14 = C23',     'value': b['Z14']})
    return rows

def table_quad_polar(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    # For compactness we expose the "polar" Z- and W- equalities listed in the Appendix for quadrilateral.
    b = quad_blocks(inputs, degrees)
    rows = []
    rows.append({'eq': 'Z12 = C34', 'value': b['C34']})
    rows.append({'eq': 'Z23 = C41', 'value': b['C41']})
    rows.append({'eq': 'Z34 = C12', 'value': b['C12']})
    rows.append({'eq': 'Z41 = C23', 'value': b['C23']})
    # W123 = C4, etc., require nested composition; not needed for simple numeric check in this baseline.
    return rows

def table_quad_half_tangent(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    """
    Half-Tangent Laws for a Quadrilateral.

    Set 1 (Appendix):
        x1 =  X34 / (Y34 + S12) = -(Y34 - S12)/X34
        x2 =  X41 / (Y41 + S23) = -(Y41 - S23)/X41
        x3 =  X12 / (Y12 + S34) = -(Y12 - S34)/X12
        x4 =  X23 / (Y23 + S41) = -(Y23 - S41)/X23

    Set 2 (Appendix, bar = decreasing order ⇒ Y flips sign):
        x1 = (X4 - X̄2)/(Y4 - Ȳ2) = -(Y4 + Ȳ2)/(X4 + X̄2), etc.
        (We report Set-2 as identities; numerically they equal Set-1.)
    """
    b = quad_blocks(inputs, degrees)

    def safe_div(n, d):
        return n/d if abs(d) > 1e-15 else float('nan')

    rows: List[Dict[str, Any]] = []
    # --- Set 1 numerics ---
    x1 = safe_div(b['X34'], (b['Y34'] + b['S12']))
    x2 = safe_div(b['X41'], (b['Y41'] + b['S23']))
    x3 = safe_div(b['X12'], (b['Y12'] + b['S34']))
    x4 = safe_div(b['X23'], (b['Y23'] + b['S41']))

    rows.append({'eq': 'Set 1: x1 = X34/(Y34 + s12) = -(Y34 - s12)/X34', 'value': x1})
    rows.append({'eq': 'Set 1: x2 = X41/(Y41 + s23) = -(Y41 - s23)/X41', 'value': x2})
    rows.append({'eq': 'Set 1: x3 = X12/(Y12 + s34) = -(Y12 - s34)/X12', 'value': x3})
    rows.append({'eq': 'Set 1: x4 = X23/(Y23 + s41) = -(Y23 - s41)/X23', 'value': x4})

    # --- Set 2 identities (equal to Set 1) ---
    rows.append({'eq': 'Set 2: x1 = (X4 - X̄2)/(Y4 - Ȳ2) = -(Y4 + Ȳ2)/(X4 + X̄2)', 'value': x1})
    rows.append({'eq': 'Set 2: x2 = (X1 - X̄3)/(Y1 - Ȳ3) = -(Y1 + Ȳ3)/(X1 + X̄3)', 'value': x2})
    rows.append({'eq': 'Set 2: x3 = (X2 - X̄4)/(Y2 - Ȳ4) = -(Y2 + Ȳ4)/(X2 + X̄4)', 'value': x3})
    rows.append({'eq': 'Set 2: x4 = (X3 - X̄1)/(Y3 - Ȳ1) = -(Y3 + Ȳ1)/(X3 + X̄1)', 'value': x4})

    return rows

def table_quad_direction_cosines(inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    b = quad_blocks(inputs, degrees)
    sets = []
    # Report the four "edges" worth of X,Y,Z as representative DC sets
    sets.append({'set': 'Edge (3->2)', 'X': b['X32'], 'Y': b['Y32'], 'Z': b['Z32']})
    sets.append({'set': 'Edge (4->1)', 'X': b['X41'], 'Y': b['Y41'], 'Z': b['Z41']})
    sets.append({'set': 'Edge (2->3)', 'X': b['X23'], 'Y': b['Y23'], 'Z': b['Z23']})
    sets.append({'set': 'Edge (1->4)', 'X': b['X14'], 'Y': b['Y14'], 'Z': b['Z14']})
    return sets

# ------------------------ Pentagon / Hexagon / Heptagon ---------------------
# Stubs with example registries (extend later with full tables)

def stub_table(_inputs: Dict[str, float], _degrees: bool=True) -> List[Dict[str, Any]]:
    return [{'info': 'Not yet implemented. Use Triangle/Quadrilateral tables for now.'}]

# ------------------------ Appendix table names & formula strings -------------
# This block lets the test GUI list *all 19* Appendix tables and print a table
# even when no numeric inputs are provided.  For Triangle and Quadrilateral we
# include the full equation lists; for the other polygons we provide placeholders
# so the names appear in the dropdown now and the formulas can be filled later.

# 19 tables in the Appendix (3 + 4 + 4 + 4 + 4)
_ALL_APPENDIX_TABLES: list[tuple[str,str]] = [
    ("Triangle", "Equations for a Spherical Triangle"),
    ("Triangle", "Equations for a Polar Triangle"),
    ("Triangle", "Direction Cosines - Spatial Triangle"),
    ("Quadrilateral", "Equations for a Spherical Quadrilateral"),
    ("Quadrilateral", "Equations for a Polar Quadrilateral"),
    ("Quadrilateral", "Half-Tangent Laws for a Quadrilateral"),
    ("Quadrilateral", "Direction Cosines - Spatial Quadrilateral"),
    ("Pentagon", "Equations for a Spherical Pentagon"),
    ("Pentagon", "Equations for a Polar Pentagon"),
    ("Pentagon", "Half-Tangent Laws for a Pentagon"),
    ("Pentagon", "Direction Cosines - Spatial Pentagon"),
    ("Hexagon", "Equations for a Spherical Hexagon"),
    ("Hexagon", "Equations for a Polar Hexagon"),
    ("Hexagon", "Half-Tangent Laws for a Hexagon"),
    ("Hexagon", "Direction Cosines - Spatial Hexagon"),
    ("Heptagon", "Equations for a Spherical Heptagon"),
    ("Heptagon", "Equations for a Polar Heptagon"),
    ("Heptagon", "Half-Tangent Laws for a Heptagon"),
    ("Heptagon", "Direction Cosines - Spatial Heptagon"),
]

# Canonical formula strings used by the "Output Table" button
_APP_FORMULAS: dict[str, list[str]] = {
    # Triangle
    "Triangle — Equations for a Spherical Triangle": [
        "X1 = S23*S3", "Y1 = S23*C3", "Z1 = C23",
        "X2 = S31*S1", "Y2 = S31*C1", "Z2 = C31",
        "X3 = S12*S2", "Y3 = S12*C2", "Z3 = C12",
    ],
    "Triangle — Equations for a Polar Triangle": [
        "U12 = S3*S23", "V12 = S3*C23", "W12 = C3",
        "U23 = S1*S31", "V23 = S1*C31", "W23 = C1",
        "U31 = S2*S12", "V31 = S2*C12", "W31 = C2",
    ],
    "Triangle — Direction Cosines - Spatial Triangle": [
        "Set 1: Xi = S23*S3,  Yi = S23*C3,  Zi = C23",
        "Set 2: Xi = S31*S1,  Yi = S31*C1,  Zi = C31",
        "Set 3: Xi = S12*S2,  Yi = S12*C2,  Zi = C12",
        "Set 4: Xi = S23*S2,  Yi = -S23*C2,  Zi = C23",
        "Set 5: Xi = S31*S3,  Yi = -S31*C3,  Zi = C31",
        "Set 6: Xi = S12*S1,  Yi = -S12*C1,  Zi = C12",
    ],

    # Quadrilateral
    "Quadrilateral — Equations for a Spherical Quadrilateral": [
        "Fundamental Formulas:",
        "  X12 = s34 s3      ;  Y12 = s34 c3      ;  Z12 = c34",
        "  X23 = s41 s4      ;  Y23 = s41 c4      ;  Z23 = c41",
        "  X34 = s12 s1      ;  Y34 = s12 c1      ;  Z34 = c12",
        "  X41 = s23 s2      ;  Y41 = s23 c2      ;  Z41 = c23",
        "  X21 = s34 s4      ;  Y21 = s34 c4      ;  Z21 = c34",
        "  X32 = s41 s1      ;  Y32 = s41 c1      ;  Z32 = c41",
        "  X43 = s12 s2      ;  Y43 = s12 c2      ;  Z43 = c12",
        "  X14 = s23 s3      ;  Y14 = s23 c3      ;  Z14 = c23",
        "",
        "Subsidiary Formulas:",
        "  X12 = X̄3      ;   −X̄12 = Ȳ3      ;   Z1 = Z3",
        "  X23 = X̄4      ;   −X̄23 = Ȳ4      ;   Z2 = Z4",
        "  X34 = X̄1      ;   −X̄34 = Ȳ1      ;   Z3 = Z1",
        "  X41 = X̄2      ;   −X̄41 = Ȳ2      ;   Z4 = Z2",
        "  X21 = X̄4      ;   −X̄21 = Ȳ4      ;   Z1 = Z3",
        "  X32 = X̄1      ;   −X̄32 = Ȳ1      ;   Z2 = Z4",
        "  X43 = X̄2      ;   −X̄43 = Ȳ2      ;   Z3 = Z1",
        "  X14 = X̄3      ;   −X̄14 = Ȳ3      ;   Z4 = Z2",
    ],
    "Quadrilateral — Equations for a Polar Quadrilateral": [
        "Z12 = C34", "Z23 = C41", "Z34 = C12", "Z41 = C23",
    ],
    "Quadrilateral — Half-Tangent Laws for a Quadrilateral": [
        "Set 1:",
        "  x1 =  X34 / (Y34 + s12)  =  -(Y34 - s12) / X34",
        "  x2 =  X41 / (Y41 + s23)  =  -(Y41 - s23) / X41",
        "  x3 =  X12 / (Y12 + s34)  =  -(Y12 - s34) / X12",
        "  x4 =  X23 / (Y23 + s41)  =  -(Y23 - s41) / X23",
        "",
        "Set 2:",
        "  x1 = (X4 - X̄2) / (Y4 - Ȳ2)  =  -(Y4 + Ȳ2) / (X4 + X̄2)",
        "  x2 = (X1 - X̄3) / (Y1 - Ȳ3)  =  -(Y1 + Ȳ3) / (X1 + X̄3)",
        "  x3 = (X2 - X̄4) / (Y2 - Ȳ4)  =  -(Y2 + Ȳ4) / (X2 + X̄4)",
        "  x4 = (X3 - X̄1) / (Y3 - Ȳ1)  =  -(Y3 + Ȳ1) / (X3 + X̄1)",
    ],
    "Quadrilateral — Direction Cosines - Spatial Quadrilateral": [
        "Edge (3→2): X32 = S41*S1,  Y32 = S41*C1,  Z32 = C41",
        "Edge (4→1): X41 = S12*S2,  Y41 = S12*C2,  Z41 = C12",
        "Edge (2→3): X23 = S34*S4,  Y23 = S34*C4,  Z23 = C34",
        "Edge (1→4): X14 = S23*S3,  Y14 = S23*C3,  Z14 = C23",
    ],
    "Pentagon — Equations for a Spherical Pentagon": [
        "Fundamental Formulas (cyclic):",
        "  X432 = s51 s1     Y432 = s51 c1     Z432 = c51",
        "  X543 = s12 s2     Y543 = s12 c2     Z543 = c12",
        "  X154 = s23 s3     Y154 = s23 c3     Z154 = c23",
        "  X215 = s34 s4     Y215 = s34 c4     Z215 = c34",
        "  X321 = s45 s5     Y321 = s45 c5     Z321 = c45",
    ],
    "Pentagon — Equations for a Polar Pentagon": [
        "Z123 = c45   Z234 = c51   Z345 = c12   Z451 = c23   Z512 = c34",
        "Note: additional U/V/W rows can be generated similarly if needed.",
    ],
    "Pentagon — Half-Tangent Laws for a Pentagon": [
        "Set 1:",
        "  x1 =  X432 / (Y432 + s12)  =  -(Y432 - s12) / X432",
        "  x2 =  X543 / (Y543 + s23)  =  -(Y543 - s23) / X543",
        "  x3 =  X154 / (Y154 + s34)  =  -(Y154 - s34) / X154",
        "  x4 =  X215 / (Y215 + s45)  =  -(Y215 - s45) / X215",
        "  x5 =  X321 / (Y321 + s51)  =  -(Y321 - s51) / X321",
        "",
        "Set 2:",
        "  Identity forms equal to Set 1 (barred variants shown in the figure).",
    ],
    "Pentagon — Direction Cosines - Spatial Pentagon": [
        "Set 1:",
        "  S1 (0, 0, 1)         S2 (0, -s12, c12)      S3 (X2,  Ȳ2,  Z2)      S4 (X32,  Y32,  Z32)      S5 (X432, Y432, Z432)",
        "  a12 (1, 0, 0)        a23 (c2, s2 c12, 0)    a34 (W32,  -U*321,  U321)   a45 (W432, -U*432,  U432)   a51 (c1, -s1, 0)",
        "Set 2:",
        "  S2 (0, 0, 1)         S3 (0, -s23, c23)      S4 (X3,  Ȳ3,  Z3)      S5 (X43,  Y43,  Z43)      S1 (X543, Y543, Z543)",
        "  a23 (1, 0, 0)        a34 (c3, s3 c23, 0)    a45 (W43,  -U*432,  U432)   a51 (W543, -U*543,  U543)   a12 (c2, -s2, 0)",
        "Set 3:",
        "  S3 (0, 0, 1)         S4 (0, -s34, c34)      S5 (X4,  Ȳ4,  Z4)      S1 (X54,  Y54,  Z54)      S2 (X154, Y154, Z154)",
        "  a34 (1, 0, 0)        a45 (c4, s4 c34, 0)    a51 (W54,  -U*543,  U543)   a12 (W154, -U*154,  U154)   a23 (c3, -s3, 0)",
        "Set 4:",
        "  S4 (0, 0, 1)         S5 (0, -s45, c45)      S1 (X5,  Ȳ5,  Z5)      S2 (X15,  Y15,  Z15)      S3 (X215, Y215, Z215)",
        "  a45 (1, 0, 0)        a51 (c5, s5 c45, 0)    a12 (W15,  -U*154,  U154)   a23 (W215, -U*215,  U215)   a34 (c4, -s4, 0)",
        "Set 5:",
        "  S5 (0, 0, 1)         S1 (0, -s51, c51)      S2 (X1,  Ȳ1,  Z1)      S3 (X21,  Y21,  Z21)      S4 (X321, Y321, Z321)",
        "  a51 (1, 0, 0)        a12 (c1, s1 c51, 0)    a23 (W21,  -U*21,   U21)    a34 (W321, -U*321,  U321)   a45 (c5, -s5, 0)",
        "Set 6:",
        "  S1 (0, 0, 1)         S5 (0,  s51, c51)      S4 (X5, -Y5,  Z5)       S3 (X45, -Y45, Z45)       S2 (X345, -Y345, Z345)",
        "  a51 (1, 0, 0)        a45 (c5, -s5 c51, 0)   a34 (W45,   U*451,  U451)   a23 (W345,  U*345,  U345)   a12 (c1,  s1, 0)",
        "Set 7:",
        "  S5 (0, 0, 1)         S4 (0,  s45, c45)      S3 (X4, -Y4,  Z4)       S2 (X34, -Y34, Z34)       S1 (X234, -Y234, Z234)",
        "  a45 (1, 0, 0)        a34 (c4, -s4 c45, 0)   a23 (W34,   U*345,  U345)   a12 (W234,  U*234,  U234)   a51 (c5,  s5, 0)",
        "Set 8:",
        "  S4 (0, 0, 1)         S3 (0,  s34, c34)      S2 (X3, -Y3,  Z3)       S1 (X23, -Y23, Z23)       S5 (X123, -Y123, Z123)",
        "  a34 (1, 0, 0)        a23 (c3, -s3 c34, 0)   a12 (W23,   U*23,   U23)    a51 (W123,  U*123,  U123)   a45 (c4,  s4, 0)",
        "Set 9:",
        "  S3 (0, 0, 1)         S2 (0,  s23, c23)      S1 (X2, -Y2,  Z2)       S5 (X12, -Y12, Z12)       S4 (X512, -Y512, Z512)",
        "  a23 (1, 0, 0)        a12 (c2, -s2 c23, 0)   a51 (W12,   U*12,   U12)    a45 (W512,  U*512,  U512)   a34 (c3,  s3, 0)",
        "Set 10:",
        "  S2 (0, 0, 1)         S1 (0,  s12, c12)      S5 (X1, -Y1,  Z1)       S4 (X51, -Y51, Z51)       S3 (X451, -Y451, Z451)",
        "  a12 (1, 0, 0)        a51 (c1, -s1 c12, 0)   a45 (W51,   U*51,   U51)    a34 (W451,  U*451,  U451)   a23 (c2,  s2, 0)",
    ],

}

# Add placeholder lists for the remaining shapes so their names appear now.
for poly in ("Pentagon","Hexagon","Heptagon"):
    for (_, tname) in filter(lambda p: p[0]==poly, _ALL_APPENDIX_TABLES):
        key = f"{poly} — {tname}"
        _APP_FORMULAS.setdefault(key, [f"{tname} — not yet implemented; formulas will be added from the Appendix."])

def get_all_appendix_table_names() -> list[str]:
    """Return 'Shape — Table' names for all 19 Appendix tables (for the test GUI)."""
    return [f"{shape} — {tname}" for (shape, tname) in _ALL_APPENDIX_TABLES]

def get_appendix_table(name: str) -> list[str]:
    """Return the raw formula strings for a table name produced by get_all_appendix_table_names()."""
    if name not in _APP_FORMULAS:
        raise KeyError(f"Unknown Appendix table: {name}")
    return _APP_FORMULAS[name][:]

# ------------------------------ Registry ------------------------------------
# ------------------------ Appendix table catalog (text output only) ----------
# Keys match exactly what the Test GUI lists in “All Appendix Tables”.
# Notation: Sij=sin(alpha_ij), Cij=cos(alpha_ij), Si=sin(theta_i), Ci=cos(theta_i).
# The “hatted” set in the Triangle page is rendered with diacritics (X̂, Ŷ, Ẑ).

_APP_TEXT: Dict[str, List[str]] = {
    # --- Triangle pages ---
    "Triangle — Equations for a Spherical Triangle": [
        "X1 = S23*S2    Y1 = S23*C2    Z1 = C23",
        "X2 = S31*S3    Y2 = S31*C3    Z2 = C31",
        "X3 = S12*S1    Y3 = S12*C1    Z3 = C12",
        "X̂1 = S23*S3    Ŷ1 = S23*C3    Ẑ1 = C23",
        "X̂2 = S31*S1    Ŷ2 = S31*C1    Ẑ2 = C31",
        "X̂3 = S12*S2    Ŷ3 = S12*C2    Ẑ3 = C12",
    ],
    "Triangle — Equations for a Polar Triangle": [
        "U12 = S3*S23    V12 = S3*C23    W12 = C3",
        "U23 = S1*S31    V23 = S1*C31    W23 = C1",
        "U31 = S2*S12    V31 = S2*C12    W31 = C2",
        "U21 = S3*S31    V21 = S3*C31    W21 = C3",
        "U32 = S1*S12    V32 = S1*C12    W32 = C1",
        "U13 = S2*S23    V13 = S2*C23    W13 = C2",
    ],
    # The Direction Cosines page is very layout-heavy; we print it as structured rows.
    # (If you want the a_ij rows added verbatim too, say: “add the a_ij rows for DC tables”.)
    "Triangle — Direction Cosines - Spatial Triangle": [
        "Set 1:",
        "  S1 (0, 0, 1)        a12 (1, 0, 0)",
        "  S2 (0, -s12, c12)   a23 (c2,  s2 c12,  U21)",
        "  S3 (X2, Y2, Z2)     a31 (c1, -s1,      0  )",
        "—",
        "Set 2:",
        "  S2 (0, 0, 1)        a23 (1, 0, 0)",
        "  S3 (0, -s23, c23)   a31 (c3,  s3 c23,  U32)",
        "  S1 (X3, Y3, Z3)     a12 (c2, -s2,      0  )",
        "—",
        "Set 3:",
        "  S3 (0, 0, 1)        a31 (1, 0, 0)",
        "  S1 (0, -s31, c31)   a12 (c1,  s1 c31,  U13)",
        "  S2 (X1, Y1, Z1)     a23 (c3, -s3,      0  )",
        "—",
        "Set 4:",
        "  S1 (0, 0, 1)        a31 (1, 0, 0)",
        "  S3 (0,  s31, c31)   a23 (c3, -s3 c31,  U31)",
        "  S2 (X3, -Y3, Z3)    a12 (c1,  s1,      0  )",
        "—",
        "Set 5:",
        "  S3 (0, 0, 1)        a23 (1, 0, 0)",
        "  S2 (0,  s23, c23)   a12 (c2, -s2 c23,  U23)",
        "  S1 (X2, -Y2, Z2)    a31 (c3,  s3,      0  )",
        "—",
        "Set 6:",
        "  S2 (0, 0, 1)        a12 (1, 0, 0)",
        "  S1 (0,  s12, c12)   a31 (c1, -s1 c12,  U12)",
        "  S3 (X1, -Y1, Z1)    a23 (c2,  s2,      0  )",
    ],

    # --- Quadrilateral pages ---
    "Quadrilateral — Equations for a Spherical Quadrilateral": [
        "Fundamental Formulas:",
        "  X12 = s34 s3     Y12 = s34 c3     Z12 = c34",
        "  X23 = s41 s4     Y23 = s41 c4     Z23 = c41",
        "  X34 = s12 s1     Y34 = s12 c1     Z34 = c12",
        "  X41 = s23 s2     Y41 = s23 c2     Z41 = c23",
        "—",
        "  X21 = s34 s4     Y21 = s34 c4     Z21 = c34",
        "  X32 = s41 s1     Y32 = s41 c1     Z32 = c41",
        "  X43 = s12 s2     Y43 = s12 c2     Z43 = c12",
        "  X14 = s23 s3     Y14 = s23 c3     Z14 = c23",
        "—",
        "Subsidiary Formulas:",
        "  X12 = X̄3     −X̄12 =  Ȳ3        Z1 =  Z3",
        "  X23 = X̄4     −X̄23 =  Ȳ4        Z2 =  Z4",
        "  X34 = X̄1     −X̄34 =  Ȳ1        Z3 =  Z1",
        "  X41 = X̄2     −X̄41 =  Ȳ2        Z4 =  Z2",
        "—",
        "  X21 = X̄4     −X̄21 =  Ȳ4        Z1 =  Z4",
        "  X32 = X̄1     −X̄32 =  Ȳ1        Z2 =  Z1",
        "  X43 = X̄2     −X̄43 =  Ȳ2        Z3 =  Z2",
        "  X14 = X̄3     −X̄14 =  Ȳ3        Z1 =  Z3",
    ],
    "Quadrilateral — Equations for a Polar Quadrilateral": [
        "Fundamental Formulas:",
        "  U123 = s4 s34    V123 = s4 c34    W123 = c4",
        "  U234 = s1 s41    V234 = s1 c41    W234 = c1",
        "  U341 = s2 s12    V341 = s2 c12    W341 = c2",
        "  U412 = s3 s23    V412 = s3 c23    W412 = c3",
        "—",
        "  U321 = s4 s41    V321 = s4 c41    W321 = c4",
        "  U214 = s3 s34    V214 = s3 c34    W214 = c3",
        "  U143 = s2 s23    V143 = s2 c23    W143 = c2",
        "  U432 = s1 s12    V432 = s1 c12    W432 = c1",
        "—",
        "Subsidiary Formulas:",
        "  U123 = U43     Ū123 = -V43     W12 = W43",
        "  U234 = U14     Ū234 = -V14     W23 = W14",
        "  U341 = U21     Ū341 = -V21     W34 = W21",
        "  U412 = U32     Ū412 = -V32     W41 = W32",
        "—",
        "  U321 = U41     Ū321 = -V41     W32 = W41",
        "  U214 = U34     Ū214 = -V34     W21 = W34",
        "  U143 = U23     Ū143 = -V23     W14 = W23",
        "  U432 = U12     Ū432 = -V12     W43 = W12",
    ],
    "Quadrilateral — Half-Tangent Laws for a Quadrilateral": [
        "Set 1:",
        "x1 =  X34/(Y34 + S12)  = -(Y34 - S12)/X34    |   x1 =  X32/(Y32 + S41)  = -(Y32 - S41)/X32",
        "x2 =  X41/(Y41 + S23)  = -(Y41 - S23)/X41    |   x2 =  X43/(Y43 + S12)  = -(Y43 - S12)/X43",
        "x3 =  X12/(Y12 + S34)  = -(Y12 - S34)/X12    |   x3 =  X14/(Y14 + S23)  = -(Y14 - S23)/X14",
        "x4 =  X23/(Y23 + S41)  = -(Y23 - S41)/X23    |   x4 =  X21/(Y21 + S34)  = -(Y21 - S34)/X21",
        "Set 2:",
        "x1 = (X4 - X̄2)/(Y4 - Ŷ2) = -(Y4 + Ŷ2)/(X4 + X̄2)",
        "x2 = (X1 - X̄3)/(Y1 - Ŷ3) = -(Y1 + Ŷ3)/(X1 + X̄3)",
        "x3 = (X2 - X̄4)/(Y2 - Ŷ4) = -(Y2 + Ŷ4)/(X2 + X̄4)",
        "x4 = (X3 - X̄1)/(Y3 - Ŷ1) = -(Y3 + Ŷ1)/(X3 + X̄1)",
    ],
    # Direction Cosines – Spatial Quadrilateral (sets only for now; a_ij rows on request)
    "Quadrilateral — Direction Cosines - Spatial Quadrilateral": [
        "Set 1:",
        "  S1 (0, 0, 1)           a12 (1, 0, 0)",
        "  S2 (0, -s12, c12)      a23 (c2,  s2 c12,  U21)",
        "  S3 (X2,  Y2,  Z2)      a34 (W32, -Ū321,  U321)",
        "  S4 (X32, Y32, Z32)     a41 (c1,  -s1,     0   )",
        "—",
        "Set 2:",
        "  S2 (0, 0, 1)           a23 (1, 0, 0)",
        "  S3 (0, -s23, c23)      a34 (c3,  s3 c23,  U32)",
        "  S4 (X3,  Y3,  Z3)      a41 (W43, -Ū432,  U432)",
        "  S1 (X43, Y43, Z43)     a12 (c2,  -s2,     0   )",
        "—",
        "Set 3:",
        "  S3 (0, 0, 1)           a34 (1, 0, 0)",
        "  S4 (0, -s34, c34)      a41 (c4,  s4 c34,  U43)",
        "  S1 (X4,  Y4,  Z4)      a12 (W14, -Ū143,  U143)",
        "  S2 (X14, Y14, Z14)     a23 (c3,  -s3,     0   )",
        "—",
        "Set 4:",
        "  S4 (0, 0, 1)           a41 (1, 0, 0)",
        "  S1 (0, -s41, c41)      a12 (c1,  s1 c41,  U14)",
        "  S2 (X1,  Y1,  Z1)      a23 (W21, -Ū214,  U214)",
        "  S3 (X21, Y21, Z21)     a34 (c4,  -s4,     0   )",
        "—",
        "Set 5:",
        "  S1 (0, 0, 1)           a41 (1, 0, 0)",
        "  S4 (0,  s41, c41)      a34 (c4, -s4 c41,  U41)",
        "  S3 (X4, -Y4, Z4)       a23 (W34,  Ū341,   U341)",
        "  S2 (X34,-Y34, Z34)     a12 (c1,   s1,     0   )",
        "—",
        "Set 6:",
        "  S4 (0, 0, 1)           a34 (1, 0, 0)",
        "  S3 (0,  s34, c34)      a23 (c3, -s3 c34,  U34)",
        "  S2 (X3, -Y3, Z3)       a12 (W23,  Ū234,   U234)",
        "  S1 (X23,-Y23, Z23)     a41 (c4,   s4,     0   )",
        "—",
        "Set 7:",
        "  S3 (0, 0, 1)           a23 (1, 0, 0)",
        "  S2 (0,  s23, c23)      a12 (c2, -s2 c23,  U23)",
        "  S1 (X2, -Y2, Z2)       a41 (W12,  Ū123,   U123)",
        "  S4 (X12,-Y12, Z12)     a34 (c3,   s3,     0   )",
        "—",
        "Set 8:",
        "  S2 (0, 0, 1)           a12 (1, 0, 0)",
        "  S1 (0,  s12, c12)      a41 (c1, -s1 c12,  U12)",
        "  S4 (X1, -Y1, Z1)       a34 (W41,  U*412,  U412)",
        "  S3 (X41,-Y41, Z41)     a23 (c2,   s2,     0   )",
    ],
    # --- Pentagon pages ---
    "Pentagon — Equations for a Spherical Pentagon": [
        "Fundamental Formulas:",
        "  X123 = s45 s4     Y123 = s45 c4     Z123 = c45",
        "  X234 = s51 s5     Y234 = s51 c5     Z234 = c51",
        "  X345 = s12 s1     Y345 = s12 c1     Z345 = c12",
        "  X451 = s23 s2     Y451 = s23 c2     Z451 = c23",
        "  X512 = s34 s3     Y512 = s34 c3     Z512 = c34",
        "—",
        "  X321 = s45 s5     Y321 = s45 c5     Z321 = c45",
        "  X432 = s51 s1     Y432 = s51 c1     Z432 = c51",
        "  X543 = s12 s2     Y543 = s12 c2     Z543 = c12",
        "  X154 = s23 s3     Y154 = s23 c3     Z154 = c23",
        "  X215 = s34 s4     Y215 = s34 c4     Z215 = c34",
        "—",
        "Subsidiary Formulas — Set 1:",
        "  X123 = X̄4     ;  −X̄123 = Ȳ4     ;  Z12 = Z̄4",
        "  X234 = X̄5     ;  −X̄234 = Ȳ5     ;  Z23 = Z̄5",
        "  X345 = X̄1     ;  −X̄345 = Ȳ1     ;  Z34 = Z̄1",
        "  X451 = X̄2     ;  −X̄451 = Ȳ2     ;  Z45 = Z̄2",
        "  X512 = X̄3     ;  −X̄512 = Ȳ3     ;  Z51 = Z̄3",
        "—",
        "  X321 = X5      ;  −X̄321 = Y5      ;  Z32 = Z5",
        "  X432 = X1      ;  −X̄432 = Y1      ;  Z43 = Z1",
        "  X543 = X2      ;  −X̄543 = Y2      ;  Z54 = Z2",
        "  X154 = X3      ;  −X̄154 = Y3      ;  Z15 = Z3",
        "  X215 = X4      ;  −X̄215 = Y4      ;  Z21 = Z4",
        "—",
        "Subsidiary Formulas — Set 2:",
        "  X12 = X43      ;  Y12 = −X̄43      ;  Z12 = Z̄4",
        "  X23 = X54      ;  Y23 = −X̄54      ;  Z23 = Z̄5",
        "  X34 = X15      ;  Y34 = −X̄15      ;  Z34 = Z̄1",
        "  X45 = X21      ;  Y45 = −X̄21      ;  Z45 = Z̄2",
        "  X51 = X32      ;  Y51 = −X̄32      ;  Z51 = Z̄3",
        "—",
        "  X32 = X51      ;  Y32 = −X̄51      ;  Z32 = Z5",
        "  X43 = X12      ;  Y43 = −X̄12      ;  Z43 = Z1",
        "  X54 = X23      ;  Y54 = −X̄23      ;  Z54 = Z2",
        "  X15 = X34      ;  Y15 = −X̄34      ;  Z15 = Z3",
        "  X21 = X45      ;  Y21 = −X̄45      ;  Z21 = Z4",
    ],

    "Pentagon — Equations for a Polar Pentagon": [
        "Fundamental Formulas:",
        "  U1234 = s5 s45    V1234 = s5 c45    W1234 = c5",
        "  U2345 = s1 s51    V2345 = s1 c51    W2345 = c1",
        "  U3451 = s2 s12    V3451 = s2 c12    W3451 = c2",
        "  U4512 = s3 s23    V4512 = s3 c23    W4512 = c3",
        "  U5123 = s4 s34    V5123 = s4 c34    W5123 = c4",
        "—",
        "  U4321 = s5 s51    V4321 = s5 c51    W4321 = c5",
        "  U3215 = s4 s45    V3215 = s4 c45    W3215 = c4",
        "  U2154 = s3 s34    V2154 = s3 c34    W2154 = c3",
        "  U1543 = s2 s23    V1543 = s2 c23    W1543 = c2",
        "  U5432 = s1 s12    V5432 = s1 c12    W5432 = c1",
        "—",
        "Subsidiary Formulas — Set 1:",
        "  U1234 = U54      ;  Ū1234 = −V54      ;  W123 = W54",
        "  U2345 = U15      ;  Ū2345 = −V15      ;  W234 = W15",
        "  U3451 = U21      ;  Ū3451 = −V21      ;  W345 = W21",
        "  U4512 = U32      ;  Ū4512 = −V32      ;  W451 = W32",
        "  U5123 = U43      ;  Ū5123 = −V43      ;  W512 = W43",
        "—",
        "  U4321 = U51      ;  Ū4321 = −V51      ;  W432 = W51",
        "  U3215 = U45      ;  Ū3215 = −V45      ;  W321 = W45",
        "  U2154 = U34      ;  Ū2154 = −V34      ;  W215 = W34",
        "  U1543 = U23      ;  Ū1543 = −V23      ;  W154 = W23",
        "  U5432 = U12      ;  Ū5432 = −V12      ;  W543 = W12",
        "—",
        "Subsidiary Formulas — Set 2:",
        "  U123 = U543      ;  V123 = −Ū543      ;  W123 = W54",
        "  U234 = U154      ;  V234 = −Ū154      ;  W234 = W15",
        "  U345 = U215      ;  V345 = −Ū215      ;  W345 = W21",
        "  U451 = U321      ;  V451 = −Ū321      ;  W451 = W32",
        "  U512 = U432      ;  V512 = −Ū432      ;  W512 = W43",
        "—",
        "  U432 = U512      ;  V432 = −Ū512      ;  W432 = W51",
        "  U321 = U451      ;  V321 = −Ū451      ;  W321 = W45",
        "  U215 = U345      ;  V215 = −Ū345      ;  W215 = W34",
        "  U154 = U234      ;  V154 = −Ū234      ;  W154 = W23",
        "  U543 = U123      ;  V543 = −Ū123      ;  W543 = W12",
    ],

    "Pentagon — Half-Tangent Laws for a Pentagon": [
        "Set 1:",
        "x1 =  X345/(Y345 + s12) = -(Y345 - s12)/X345    |   x1 =  X432/(Y432 + s51) = -(Y432 - s51)/X432",
        "x2 =  X451/(Y451 + s23) = -(Y451 - s23)/X451    |   x2 =  X543/(Y543 + s12) = -(Y543 - s12)/X543",
        "x3 =  X512/(Y512 + s34) = -(Y512 - s34)/X512    |   x3 =  X154/(Y154 + s23) = -(Y154 - s23)/X154",
        "x4 =  X123/(Y123 + s45) = -(Y123 - s45)/X123    |   x4 =  X215/(Y215 + s34) = -(Y215 - s34)/X215",
        "x5 =  X234/(Y234 + s51) = -(Y234 - s51)/X234    |   x5 =  X321/(Y321 + s45) = -(Y321 - s45)/X321",
        "—",
        "Set 2:",
        "x1 = (X45 − X̄2)/(Y45 − Ȳ2) = −(Y45 + Ȳ2)/(X45 + X̄2)    |   x1 = (X32 − X5)/(Y32 − Y5) = −(Y32 + Y5)/(X32 + X5)",
        "x2 = (X51 − X̄3)/(Y51 − Ȳ3) = −(Y51 + Ȳ3)/(X51 + X̄3)    |   x2 = (X43 − X1)/(Y43 − Y1) = −(Y43 + Y1)/(X43 + X1)",
        "x3 = (X12 − X̄4)/(Y12 − Ȳ4) = −(Y12 + Ȳ4)/(X12 + X̄4)    |   x3 = (X54 − X2)/(Y54 − Y2) = −(Y54 + Y2)/(X54 + X2)",
        "x4 = (X23 − X̄5)/(Y23 − Ȳ5) = −(Y23 + Ȳ5)/(X23 + X̄5)    |   x4 = (X15 − X3)/(Y15 − Y3) = −(Y15 + Y3)/(X15 + X3)",
        "x5 = (X34 − X̄1)/(Y34 − Ȳ1) = −(Y34 + Ȳ1)/(X34 + X̄1)    |   x5 = (X21 − X4)/(Y21 − Y4) = −(Y21 + Y4)/(X21 + X4)",
    ],

    # Direction Cosines – 10 set headers (now with a_ij rows)
    "Pentagon — Direction Cosines - Spatial Pentagon": [
        "Set 1:  S1(0,0,1) ; S2(0,−s12,c12) ; S3(X2, Ȳ2, Z2) ; S4(X32, Y32, Z32) ; S5(X432, Y432, Z432)",
        "       a12(1,0,0) ; a23(c2, s2 c12, U21) ; a34(W32, −U*321, U321) ; a45(W432, −U*4321, U4321) ; a51(c1, −s1, 0)",
        "Set 2:  S2(0,0,1) ; S3(0,−s23,c23) ; S4(X3, Ȳ3, Z3) ; S5(X43, Y43, Z43) ; S1(X543, Y543, Z543)",
        "       a23(1,0,0) ; a34(c3, s3 c23, U32) ; a45(W43, −U*432, U432) ; a51(W543, −U*5432, U5432) ; a12(c2, −s2, 0)",
        "Set 3:  S3(0,0,1) ; S4(0,−s34,c34) ; S5(X4, Ȳ4, Z4) ; S1(X54, Y54, Z54) ; S2(X154, Y154, Z154)",
        "       a34(1,0,0) ; a45(c4, s4 c34, U43) ; a51(W54, −U*543, U543) ; a12(W154, −U*1543, U1543) ; a23(c3, −s3, 0)",
        "Set 4:  S4(0,0,1) ; S5(0,−s45,c45) ; S1(X5, Ȳ5, Z5) ; S2(X15, Y15, Z15) ; S3(X215, Y215, Z215)",
        "       a45(1,0,0) ; a51(c5, s5 c45, U54) ; a12(W15, −U*154, U154) ; a23(W215, −U*2154, U2154) ; a34(c4, −s4, 0)",
        "Set 5:  S5(0,0,1) ; S1(0,−s51,c51) ; S2(X1, Ȳ1, Z1) ; S3(X21, Y21, Z21) ; S4(X321, Y321, Z321)",
        "       a51(1,0,0) ; a12(c1, s1 c51, U15) ; a23(W21, −U*215, U215) ; a34(W321, −U*3215, U3215) ; a45(c5, −s5, 0)",
        "Set 6:  S1(0,0,1) ; S5(0, s51,c51) ; S4(X5, −Y5, Z5) ; S3(X45, −Y45, Z45) ; S2(X345, −Y345, Z345)",
        "       a51(1,0,0) ; a45(c5, −s5 c51, U51) ; a34(W45, U*451, U451) ; a23(W345, U*3451, U3451) ; a12(c1, s1, 0)",
        "Set 7:  S5(0,0,1) ; S4(0, s45,c45) ; S3(X4, −Y4, Z4) ; S2(X34, −Y34, Z34) ; S1(X234, −Y234, Z234)",
        "       a45(1,0,0) ; a34(c4, −s4 c45, U45) ; a23(W34, U*345, U345) ; a12(W234, U*2345, U2345) ; a51(c5, s5, 0)",
        "Set 8:  S4(0,0,1) ; S3(0, s34,c34) ; S2(X3, −Y3, Z3) ; S1(X23, −Y23, Z23) ; S5(X123, −Y123, Z123)",
        "       a34(1,0,0) ; a23(c3, −s3 c34, U34) ; a12(W23, U*234, U234) ; a51(W123, U*1234, U1234) ; a45(c4, s4, 0)",
        "Set 9:  S3(0,0,1) ; S2(0, s23,c23) ; S1(X2, −Y2, Z2) ; S5(X12, −Y12, Z12) ; S4(X512, −Y512, Z512)",
        "       a23(1,0,0) ; a12(c2, −s2 c23, U23) ; a51(W12, U*123, U123) ; a45(W512, U*5123, U5123) ; a34(c3, s3, 0)",
        "Set10:  S2(0,0,1) ; S1(0, s12,c12) ; S5(X1, −Y1, Z1) ; S4(X51, −Y51, Z51) ; S3(X451, −Y451, Z451)",
        "       a12(1,0,0) ; a51(c1, −s1 c12, U12) ; a45(W51, U*512, U512) ; a34(W451, U*4512, U4512) ; a23(c2, s2, 0)",
    ],
    "Hexagon — Equations for a Spherical Hexagon": [
        "Fundamental Formulas:",
        "  X1234 = s56 s5     Y1234 = s56 c5     Z1234 = c56",
        "  X2345 = s61 s6     Y2345 = s61 c6     Z2345 = c61",
        "  X3456 = s12 s1     Y3456 = s12 c1     Z3456 = c12",
        "  X4561 = s23 s2     Y4561 = s23 c2     Z4561 = c23",
        "  X5612 = s34 s3     Y5612 = s34 c3     Z5612 = c34",
        "  X6123 = s45 s4     Y6123 = s45 c4     Z6123 = c45",
        "—",
        "  X4321 = s56 s6     Y4321 = s56 c6     Z4321 = c56",
        "  X5432 = s61 s1     Y5432 = s61 c1     Z5432 = c61",
        "  X6543 = s12 s2     Y6543 = s12 c2     Z6543 = c12",
        "  X1654 = s23 s3     Y1654 = s23 c3     Z1654 = c23",
        "  X2165 = s34 s4     Y2165 = s34 c4     Z2165 = c34",
        "  X3216 = s45 s5     Y3216 = s45 c5     Z3216 = c45",
        "—",
        "Subsidiary Formulas:",
        "  Set 1:",
        "    X1234 = X̄5      −X*1234 = Ȳ5      Z123 = Z5",
        "    X2345 = X̄6      −X*2345 = Ȳ6      Z234 = Z6",
        "    X3456 = X̄1      −X*3456 = Ȳ1      Z345 = Z1",
        "    X4561 = X̄2      −X*4561 = Ȳ2      Z456 = Z2",
        "    X5612 = X̄3      −X*5612 = Ȳ3      Z561 = Z3",
        "    X6123 = X̄4      −X*6123 = Ȳ4      Z612 = Z4",
        "  Set 2:",
        "    X123 = X54       Y123 = −X*54      Z123 = Z̄5",
        "    X234 = X65       Y234 = −X*65      Z234 = Z̄6",
        "    X345 = X16       Y345 = −X*16      Z345 = Z̄1",
        "    X456 = X21       Y456 = −X*21      Z456 = Z̄2",
        "    X561 = X32       Y561 = −X*32      Z561 = Z̄3",
        "    X612 = X43       Y612 = −X*43      Z612 = Z̄4",
        "    X432 = X61       Y432 = −X*61      Z432 = Z̄6",
        "    X543 = X12       Y543 = −X*12      Z543 = Z̄1",
        "    X654 = X23       Y654 = −X*23      Z654 = Z̄2",
        "    X165 = X34       Y165 = −X*34      Z165 = Z̄3",
        "    X216 = X45       Y216 = −X*45      Z216 = Z̄4",
        "    X321 = X56       Y321 = −X*56      Z321 = Z̄5",
        "  Set 3:",
        "    X123 = X54       −X*123 = Y54      Z12 = Z54",
        "    X234 = X65       −X*234 = Y65      Z23 = Z65",
        "    X345 = X16       −X*345 = Y16      Z34 = Z16",
        "    X456 = X21       −X*456 = Y21      Z45 = Z21",
        "    X561 = X32       −X*561 = Y32      Z56 = Z32",
        "    X612 = X43       −X*612 = Y43      Z61 = Z43",
        "    X432 = X61       −X*432 = Y61      Z43 = Z61",
        "    X543 = X12       −X*543 = Y12      Z54 = Z12",
        "    X654 = X23       −X*654 = Y23      Z65 = Z23",
        "    X165 = X34       −X*165 = Y34      Z16 = Z34",
        "    X216 = X45       −X*216 = Y45      Z21 = Z45",
        "    X321 = X56       −X*321 = Y56      Z32 = Z56",
    ],
    "Hexagon — Equations for a Polar Hexagon": [
        "Fundamental Formulas:",
        "  U12345 = s6 s56    V12345 = s6 c56    W12345 = c6",
        "  U23456 = s1 s61    V23456 = s1 c61    W23456 = c1",
        "  U34561 = s2 s12    V34561 = s2 c12    W34561 = c2",
        "  U45612 = s3 s23    V45612 = s3 c23    W45612 = c3",
        "  U56123 = s4 s34    V56123 = s4 c34    W56123 = c4",
        "  U61234 = s5 s45    V61234 = s5 c45    W61234 = c5",
        "—",
        "  U54321 = s6 s61    V54321 = s6 c61    W54321 = c6",
        "  U43216 = s5 s56    V43216 = s5 c56    W43216 = c5",
        "  U32165 = s4 s45    V32165 = s4 c45    W32165 = c4",
        "  U21654 = s3 s34    V21654 = s3 c34    W21654 = c3",
        "  U16543 = s2 s23    V16543 = s2 c23    W16543 = c2",
        "  U65432 = s1 s12    V65432 = s1 c12    W65432 = c1",
        "—",
        "Subsidiary Formulas:",
        "  Set 1:",
        "    U12345 = U65     U*12345 = −V65     W12345 = W65",
        "    U23456 = U16     U*23456 = −V16     W23456 = W16",
        "    U34561 = U21     U*34561 = −V21     W34561 = W21",
        "    U45612 = U32     U*45612 = −V32     W45612 = W32",
        "    U56123 = U43     U*56123 = −V43     W56123 = W43",
        "    U61234 = U54     U*61234 = −V54     W61234 = W54",
        "    U54321 = U61     U*54321 = −V61     W54321 = W61",
        "    U43216 = U56     U*43216 = −V56     W43216 = W56",
        "    U32165 = U45     U*32165 = −V45     W32165 = W45",
        "    U21654 = U34     U*21654 = −V34     W21654 = W34",
        "    U16543 = U23     U*16543 = −V23     W16543 = W23",
        "    U65432 = U12     U*65432 = −V12     W65432 = W12",
        "  Set 2:",
        "    U1234 = U654     V1234 = −U*654     W1234 = W654",
        "    U2345 = U165     V2345 = −U*165     W2345 = W16",
        "    U3456 = U216     V3456 = −U*216     W3456 = W21",
        "    U4561 = U321     V4561 = −U*321     W4561 = W32",
        "    U5612 = U432     V5612 = −U*432     W5612 = W43",
        "    U6123 = U543     V6123 = −U*543     W6123 = W54",
        "    U5432 = U612     V5432 = −U*612     W5432 = W61",
        "    U4321 = U561     V4321 = −U*561     W4321 = W56",
        "    U3216 = U456     V3216 = −U*456     W3216 = W45",
        "    U2165 = U345     V2165 = −U*345     W2165 = W34",
        "    U1654 = U234     V1654 = −U*234     W1654 = W23",
        "    U6543 = U123     V6543 = −U*123     W6543 = W12",
        "  Set 3:",
        "    U1234 = U654     U*1234 = −V654     W1234 = W654",
        "    U2345 = U165     U*2345 = −V165     W2345 = W16",
        "    U3456 = U216     U*3456 = −V216     W3456 = W21",
        "    U4561 = U321     U*4561 = −V321     W4561 = W32",
        "    U5612 = U432     U*5612 = −V432     W5612 = W43",
        "    U6123 = U543     U*6123 = −V543     W6123 = W54",
        "    U5432 = U612     U*5432 = −V612     W5432 = W61",
        "    U4321 = U561     U*4321 = −V561     W4321 = W56",
        "    U3216 = U456     U*3216 = −V456     W3216 = W45",
        "    U2165 = U345     U*2165 = −V345     W2165 = W34",
        "    U1654 = U234     U*1654 = −V234     W1654 = W23",
        "    U6543 = U123     U*6543 = −V123     W6543 = W12",
    ],
    "Hexagon — Half-Tangent Laws for a Hexagon": [
        "Set 1 (left):",
        "  x1 =  X3456/(Y3456 + s12)  =  −(Y3456 − s12)/X3456",
        "  x2 =  X4561/(Y4561 + s23)  =  −(Y4561 − s23)/X4561",
        "  x3 =  X5612/(Y5612 + s34)  =  −(Y5612 − s34)/X5612",
        "  x4 =  X6123/(Y6123 + s45)  =  −(Y6123 − s45)/X6123",
        "  x5 =  X1234/(Y1234 + s56)  =  −(Y1234 − s56)/X1234",
        "  x6 =  X2345/(Y2345 + s61)  =  −(Y2345 − s61)/X2345",
        "Set 1 (right):",
        "  x1 =  X5432/(Y5432 + s61)  =  −(Y5432 − s61)/X5432",
        "  x2 =  X6543/(Y6543 + s12)  =  −(Y6543 − s12)/X6543",
        "  x3 =  X1654/(Y1654 + s23)  =  −(Y1654 − s23)/X1654",
        "  x4 =  X2165/(Y2165 + s34)  =  −(Y2165 − s34)/X2165",
        "  x5 =  X3216/(Y3216 + s45)  =  −(Y3216 − s45)/X3216",
        "  x6 =  X4321/(Y4321 + s56)  =  −(Y4321 − s56)/X4321",
        "—",
        "Set 2 (barred forms):",
        "  x1 = (X456 − X̄2)/(Y456 − Ȳ2) = −(Y456 + Ȳ2)/(X456 + X̄2)",
        "  x2 = (X561 − X̄3)/(Y561 − Ȳ3) = −(Y561 + Ȳ3)/(X561 + X̄3)",
        "  x3 = (X612 − X̄4)/(Y612 − Ȳ4) = −(Y612 + Ȳ4)/(X612 + X̄4)",
        "  x4 = (X123 − X̄5)/(Y123 − Ȳ5) = −(Y123 + Ȳ5)/(X123 + X̄5)",
        "  x5 = (X234 − X̄6)/(Y234 − Ȳ6) = −(Y234 + Ȳ6)/(X234 + X̄6)",
        "  x6 = (X345 − X̄1)/(Y345 − Ȳ1) = −(Y345 + Ȳ1)/(X345 + X̄1)",
        "  x1 = (X432 − X6)/(Y432 − Y6)  = −(Y432 + Y6)/(X432 + X6)",
        "  x2 = (X543 − X1)/(Y543 − Y1)  = −(Y543 + Y1)/(X543 + X1)",
        "  x3 = (X654 − X2)/(Y654 − Y2)  = −(Y654 + Y2)/(X654 + X2)",
        "  x4 = (X165 − X3)/(Y165 − Y3)  = −(Y165 + Y3)/(X165 + X3)",
        "  x5 = (X216 − X4)/(Y216 − Y4)  = −(Y216 + Y4)/(X216 + X4)",
        "  x6 = (X321 − X5)/(Y321 − Y5)  = −(Y321 + Y5)/(X321 + X5)",
        "—",
        "Set 3:",
        "  x1 =  (X56 − X32)/(Y56 − Y32)  =  −(Y56 + Y32)/(X56 + X32)",
        "  x2 =  (X61 − X43)/(Y61 − Y43)  =  −(Y61 + Y43)/(X61 + X43)",
        "  x3 =  (X12 − X54)/(Y12 − Y54)  =  −(Y12 + Y54)/(X12 + X54)",
        "  x4 =  (X23 − X65)/(Y23 − Y65)  =  −(Y23 + Y65)/(X23 + X65)",
        "  x5 =  (X34 − X16)/(Y34 − Y16)  =  −(Y34 + Y16)/(X34 + X16)",
        "  x6 =  (X45 − X21)/(Y45 − Y21)  =  −(Y45 + Y21)/(X45 + X21)",
    ],
    "Hexagon — Direction Cosines - Spatial Hexagon": [
        "Set 1:  S1(0,0,1) ; S2(0,−s12,c12) ; S3(X2, Ȳ2, Z2) ; S4(X32, Y32, Z32) ; S5(X432, Y432, Z432) ; S6(X5432, Y5432, Z5432)",
        "       a12(1,0,0) ; a23(c2, s2 c12, U21) ; a34(W32, −U*321, U321) ; a45(W432, −U*4321, U4321) ; a56(W5432, −U*5432, U5432) ; a61(c1, −s1, 0)",
        "Set 2:  S2(0,0,1) ; S3(0,−s23,c23) ; S4(X3, Ȳ3, Z3) ; S5(X43, Y43, Z43) ; S6(X543, Y543, Z543) ; S1(X6543, Y6543, Z6543)",
        "       a23(1,0,0) ; a34(c3, s3 c23, U32) ; a45(W43, −U*432, U432) ; a56(W543, −U*543, U543) ; a61(W6543, −U*6543, U6543) ; a12(c2, −s2, 0)",
        "Set 3:  S3(0,0,1) ; S4(0,−s34,c34) ; S5(X4, Ȳ4, Z4) ; S6(X54, Y54, Z54) ; S1(X654, Y654, Z654) ; S2(X1654, Y1654, Z1654)",
        "       a34(1,0,0) ; a45(c4, s4 c34, U43) ; a56(W54, −U*543, U543) ; a61(W654, −U*654, U654) ; a12(W1654, −U*1654, U1654) ; a23(c3, −s3, 0)",
        "Set 4:  S4(0,0,1) ; S5(0,−s45,c45) ; S6(X5, Ȳ5, Z5) ; S1(X65, Y65, Z65) ; S2(X165, Y165, Z165) ; S3(X2165, Y2165, Z2165)",
        "       a45(1,0,0) ; a56(c5, s5 c45, U54) ; a61(W65, −U*654, U654) ; a12(W165, −U*165, U165) ; a23(W2165, −U*2165, U2165) ; a34(c4, −s4, 0)",
        "Set 5:  S5(0,0,1) ; S6(0,−s56,c56) ; S1(X6, Ȳ6, Z6) ; S2(X16, Y16, Z16) ; S3(X216, Y216, Z216) ; S4(X3216, Y3216, Z3216)",
        "       a56(1,0,0) ; a61(c6, s6 c56, U65) ; a12(W16, −U*165, U165) ; a23(W216, −U*216, U216) ; a34(W3216, −U*3216, U3216) ; a45(c5, −s5, 0)",
        "Set 6:  S6(0,0,1) ; S1(0,−s61,c61) ; S2(X1, Ȳ1, Z1) ; S3(X21, Y21, Z21) ; S4(X321, Y321, Z321) ; S5(X4321, Y4321, Z4321)",
        "       a61(1,0,0) ; a12(c1, s1 c61, U16) ; a23(W21, −U*21, U21) ; a34(W321, −U*321, U321) ; a45(W4321, −U*4321, U4321) ; a56(c6, −s6, 0)",
        "Set 7:  S1(0,0,1) ; S6(0, s61,c61) ; S5(X6, −Y6, Z6) ; S4(X56, −Y56, Z56) ; S3(X456, −Y456, Z456) ; S2(X3456, −Y3456, Z3456)",
        "       a12(1,0,0) ; a61(c1, −s1 c61, −U16) ; a56(W6, U*65, −U65) ; a45(W56, U*56, −U56) ; a34(W456, U*456, −U456) ; a23(c2, s2, 0)",
        "Set 8:  S6(0,0,1) ; S5(0, s56,c56) ; S4(X5, −Y5, Z5) ; S3(X45, −Y45, Z45) ; S2(X345, −Y345, Z345) ; S1(X2345, −Y2345, Z2345)",
        "       a61(1,0,0) ; a56(c6, −s6 c56, −U65) ; a45(W5, U*54, −U54) ; a34(W45, U*45, −U45) ; a23(W345, U*345, −U345) ; a12(c1, s1, 0)",
        "Set 9:  S5(0,0,1) ; S4(0, s45,c45) ; S3(X4, −Y4, Z4) ; S2(X34, −Y34, Z34) ; S1(X234, −Y234, Z234) ; S6(X1234, −Y1234, Z1234)",
        "       a56(1,0,0) ; a45(c5, −s5 c45, −U54) ; a34(W4, U*43, −U43) ; a23(W34, U*34, −U34) ; a12(W234, U*234, −U234) ; a61(c6, s6, 0)",
        "Set10:  S4(0,0,1) ; S3(0, s34,c34) ; S2(X3, −Y3, Z3) ; S1(X23, −Y23, Z23) ; S6(X123, −Y123, Z123) ; S5(X6123, −Y6123, Z6123)",
        "       a45(1,0,0) ; a34(c4, −s4 c34, −U43) ; a23(W3, U*32, −U32) ; a12(W23, U*23, −U23) ; a61(W123, U*123, −U123) ; a56(c5, s5, 0)",
        "Set11:  S3(0,0,1) ; S2(0, s23,c23) ; S1(X2, −Y2, Z2) ; S6(X12, −Y12, Z12) ; S5(X612, −Y612, Z612) ; S4(X5612, −Y5612, Z5612)",
        "       a34(1,0,0) ; a23(c3, −s3 c23, −U32) ; a12(W2, U*21, −U21) ; a61(W12, U*12, −U12) ; a56(W612, U*612, −U612) ; a45(c4, s4, 0)",
        "Set12:  S2(0,0,1) ; S1(0, s12,c12) ; S6(X1, −Y1, Z1) ; S5(X61, −Y61, Z61) ; S4(X561, −Y561, Z561) ; S3(X4561, −Y4561, Z4561)",
        "       a23(1,0,0) ; a12(c2, −s2 c12, −U21) ; a61(W1, U*16, −U16) ; a56(W61, U*61, −U61) ; a45(W561, U*561, −U561) ; a34(c3, s3, 0)",
    ],
    "Heptagon — Equations for a Spherical Heptagon": [
        "Fundamental Formulas:",
        "  X12345 = s67 s6    Y12345 = s67 c6    Z12345 = c67",
        "  X23456 = s71 s7    Y23456 = s71 c7    Z23456 = c71",
        "  X34567 = s12 s1    Y34567 = s12 c1    Z34567 = c12",
        "  X45671 = s23 s2    Y45671 = s23 c2    Z45671 = c23",
        "  X56712 = s34 s3    Y56712 = s34 c3    Z56712 = c34",
        "  X67123 = s45 s4    Y67123 = s45 c4    Z67123 = c45",
        "  X71234 = s56 s5    Y71234 = s56 c5    Z71234 = c56",
        "—",
        "  X54321 = s67 s7    Y54321 = s67 c7    Z54321 = c67",
        "  X65432 = s71 s1    Y65432 = s71 c1    Z65432 = c71",
        "  X76543 = s12 s2    Y76543 = s12 c2    Z76543 = c12",
        "  X17654 = s23 s3    Y17654 = s23 c3    Z17654 = c23",
        "  X21765 = s34 s4    Y21765 = s34 c4    Z21765 = c34",
        "  X32176 = s45 s5    Y32176 = s45 c5    Z32176 = c45",
        "  X43217 = s56 s6    Y43217 = s56 c6    Z43217 = c56",
        "—",
        "Subsidiary Formulas — Set 1:",
        "  X12345 = X̄6      −X*12345 = Ȳ6      Z1234 = Z̄6",
        "  X23456 = X̄7      −X*23456 = Ȳ7      Z2345 = Z̄7",
        "  X34567 = X̄1      −X*34567 = Ȳ1      Z3456 = Z̄1",
        "  X45671 = X̄2      −X*45671 = Ȳ2      Z4567 = Z̄2",
        "  X56712 = X̄3      −X*56712 = Ȳ3      Z5671 = Z̄3",
        "  X67123 = X̄4      −X*67123 = Ȳ4      Z6712 = Z̄4",
        "  X71234 = X̄5      −X*71234 = Ȳ5      Z7123 = Z̄5",
        "  X54321 = X7       −X*54321 = Y7       Z5432 = Z7",
        "  X65432 = X1       −X*65432 = Y1       Z6543 = Z1",
        "  X76543 = X2       −X*76543 = Y2       Z7654 = Z2",
        "  X17654 = X3       −X*17654 = Y3       Z1765 = Z3",
        "  X21765 = X4       −X*26715 = Y4       Z2176 = Z4",
        "  X32176 = X5       −X*32176 = Y5       Z3217 = Z5",
        "  X43217 = X6       −X*43217 = Y6       Z4321 = Z6",
        "—",
        "Subsidiary Formulas — Set 2:",
        "  X1234 = X65       Y1234 = −X*65       Z1234 = Z̄6",
        "  X2345 = X76       Y2345 = −X*76       Z2345 = Z̄7",
        "  X3456 = X17       Y3456 = −X*17       Z3456 = Z̄1",
        "  X4567 = X21       Y4567 = −X*21       Z4567 = Z̄2",
        "  X5671 = X32       Y5671 = −X*32       Z5671 = Z̄3",
        "  X6712 = X43       Y6712 = −X*43       Z6712 = Z̄4",
        "  X7123 = X54       Y7123 = −X*54       Z7123 = Z̄5",
        "  X5432 = X71       Y5432 = −X*71       Z5432 = Z7",
        "  X6543 = X12       Y6543 = −X*12       Z6543 = Z1",
        "  X7654 = X23       Y7654 = −X*23       Z7654 = Z2",
        "  X1765 = X34       Y1765 = −X*34       Z1765 = Z3",
        "  X2176 = X45       Y2176 = −X*45       Z2176 = Z4",
        "  X3217 = X56       Y3217 = −X*56       Z3217 = Z5",
        "  X4321 = X67       Y4321 = −X*67       Z4321 = Z6",
        "—",
        "Subsidiary Formulas — Set 3:",
        "  X1234 = X65       −X*1234 = Y65       Z123 = Z65",
        "  X2345 = X76       −X*2345 = Y76       Z234 = Z76",
        "  X3456 = X17       −X*3456 = Y17       Z345 = Z17",
        "  X4567 = X21       −X*4567 = Y21       Z456 = Z21",
        "  X5671 = X32       −X*5671 = Y32       Z567 = Z32",
        "  X6712 = X43       −X*6712 = Y43       Z671 = Z43",
        "  X7123 = X54       −X*7123 = Y54       Z712 = Z54",
        "  X5432 = X71       −X*5432 = Y71       Z543 = Z71",
        "  X6543 = X12       −X*6543 = Y12       Z654 = Z12",
        "  X7654 = X23       −X*7654 = Y23       Z765 = Z23",
        "  X1765 = X34       −X*1765 = Y34       Z176 = Z34",
        "  X2176 = X45       −X*2176 = Y45       Z217 = Z45",
        "  X3217 = X56       −X*3217 = Y56       Z321 = Z56",
        "  X4321 = X67       −X*4321 = Y67       Z432 = Z67",
        "—",
        "Subsidiary Formulas — Set 4:",
        "  X123 = X654       Y123 = −X*654       Z123 = Z65",
        "  X234 = X765       Y234 = −X*765       Z234 = Z76",
        "  X345 = X176       Y345 = −X*176       Z345 = Z17",
        "  X456 = X217       Y456 = −X*217       Z456 = Z21",
        "  X567 = X321       Y567 = −X*321       Z567 = Z32",
        "  X671 = X432       Y671 = −X*432       Z671 = Z43",
        "  X712 = X543       Y712 = −X*543       Z712 = Z54",
        "  X543 = X712       Y543 = −X*712       Z543 = Z71",
        "  X654 = X123       Y654 = −X*123       Z654 = Z12",
        "  X765 = X234       Y765 = −X*234       Z765 = Z23",
        "  X176 = X345       Y176 = −X*345       Z176 = Z34",
        "  X217 = X456       Y217 = −X*456       Z217 = Z45",
        "  X321 = X567       Y321 = −X*567       Z321 = Z56",
        "  X432 = X671       Y432 = −X*671       Z432 = Z67",
    ],
    "Heptagon — Equations for a Polar Heptagon": [
        "Fundamental Formulas:",
        "  U123456 = s7 s67   V123456 = s7 c67   W123456 = c7",
        "  U234567 = s1 s71   V234567 = s1 c71   W234567 = c1",
        "  U345671 = s2 s12   V345671 = s2 c12   W345671 = c2",
        "  U456712 = s3 s23   V456712 = s3 c23   W456712 = c3",
        "  U567123 = s4 s34   V567123 = s4 c34   W567123 = c4",
        "  U671234 = s5 s45   V671234 = s5 c45   W671234 = c5",
        "  U712345 = s6 s56   V712345 = s6 c56   W712345 = c6",
        "—",
        "  U654321 = s7 s71   V654321 = s7 c71   W654321 = c7",
        "  U765432 = s1 s12   V765432 = s1 c12   W765432 = c1",
        "  U176543 = s2 s23   V176543 = s2 c23   W176543 = c2",
        "  U217654 = s3 s34   V217654 = s3 c34   W217654 = c3",
        "  U321765 = s4 s45   V321765 = s4 c45   W321765 = c4",
        "  U432176 = s5 s56   V432176 = s5 c56   W432176 = c5",
        "  U543217 = s6 s67   V543217 = s6 c67   W543217 = c6",
        "—",
        "Subsidiary Formulas — Set 1:",
        "  U123456 = U76      U*123456 = −V76      W12345 = W76",
        "  U234567 = U17      U*234567 = −V17      W23456 = W17",
        "  U345671 = U21      U*345671 = −V21      W34567 = W21",
        "  U456712 = U32      U*456712 = −V32      W45671 = W32",
        "  U567123 = U43      U*567123 = −V43      W56712 = W43",
        "  U671234 = U54      U*671234 = −V54      W67123 = W54",
        "  U712345 = U65      U*712345 = −V65      W71234 = W65",
        "  U654321 = U71      U*654321 = −V71      W65432 = W71",
        "  U543217 = U67      U*543217 = −V67      W54321 = W67",
        "  U432176 = U56      U*432176 = −V56      W43217 = W56",
        "  U321765 = U45      U*321765 = −V45      W32176 = W45",
        "  U217654 = U34      U*217654 = −V34      W21765 = W34",
        "  U176543 = U23      U*176543 = −V23      W17654 = W23",
        "  U765432 = U12      U*765432 = −V12      W76543 = W12",
        "—",
        "Subsidiary Formulas — Set 2:",
        "  U12345 = U765      V12345 = −U*765      W12345 = W76",
        "  U23456 = U176      V23456 = −U*176      W23456 = W17",
        "  U34567 = U217      V34567 = −U*217      W34567 = W21",
        "  U45671 = U321      V45671 = −U*321      W45671 = W32",
        "  U56712 = U432      V56712 = −U*432      W56712 = W43",
        "  U67123 = U543      V67123 = −U*543      W67123 = W54",
        "  U71234 = U654      V71234 = −U*654      W71234 = W65",
        "  U65432 = U712      V65432 = −U*712      W65432 = W71",
        "  U54321 = U671      V54321 = −U*671      W54321 = W67",
        "  U43217 = U567      V43217 = −U*567      W43217 = W56",
        "  U32176 = U456      V32176 = −U*456      W32176 = W45",
        "  U21765 = U345      V21765 = −U*345      W21765 = W34",
        "  U17654 = U234      V17654 = −U*234      W17654 = W23",
        "  U76543 = U123      V76543 = −U*123      W76543 = W12",
        "—",
        "Subsidiary Formulas — Set 3:",
        "  U12345 = U765      U*12345 = −V765      W12345 = W76",
        "  U23456 = U176      U*23456 = −V176      W23456 = W17",
        "  U34567 = U217      U*34567 = −V217      W34567 = W21",
        "  U45671 = U321      U*45671 = −V321      W45671 = W32",
        "  U56712 = U432      U*56712 = −V432      W56712 = W43",
        "  U67123 = U543      U*67123 = −V543      W67123 = W54",
        "  U71234 = U654      U*71234 = −V654      W71234 = W65",
        "  U65432 = U712      U*65432 = −V712      W65432 = W71",
        "  U54321 = U671      U*54321 = −V671      W54321 = W67",
        "  U43217 = U567      U*43217 = −V567      W43217 = W56",
        "  U32176 = U456      U*32176 = −V456      W32176 = W45",
        "  U21765 = U345      U*21765 = −V345      W21765 = W34",
        "  U17654 = U234      U*17654 = −V234      W17654 = W23",
        "  U76543 = U123      U*76543 = −V123      W76543 = W12",
        "—",
        "Subsidiary Formulas — Set 4:",
        "  U1234  = U7654     V1234  = −U*7654     W1234  = W765",
        "  U2345  = U1765     V2345  = −U*1765     W2345  = W176",
        "  U3456  = U2176     V3456  = −U*2176     W3456  = W217",
        "  U4567  = U3217     V4567  = −U*3217     W4567  = W321",
        "  U5671  = U4321     V5671  = −U*4321     W5671  = W432",
        "  U6712  = U5432     V6712  = −U*5432     W6712  = W543",
        "  U7123  = U6543     V7123  = −U*6543     W7123  = W654",
        "  U6543  = U7123     V6543  = −U*7123     W6543  = W712",
        "  U5432  = U6712     V5432  = −U*6712     W5432  = W671",
        "  U4321  = U5671     V4321  = −U*5671     W4321  = W567",
        "  U3217  = U4567     V3217  = −U*4567     W3217  = W456",
        "  U2176  = U3456     V2176  = −U*3456     W2176  = W345",
        "  U1765  = U2345     V1765  = −U*2345     W1765  = W234",
        "  U7654  = U1234     V7654  = −U*1234     W7654  = W123",
    ],
    "Heptagon — Half-Tangent Laws for a Heptagon": [
        "Set 1 (forward cycle):",
        "  x1 =  X34567 / (Y34567 + s12)  =  -(Y34567 - s12) / X34567",
        "  x2 =  X45671 / (Y45671 + s23)  =  -(Y45671 - s23) / X45671",
        "  x3 =  X56712 / (Y56712 + s34)  =  -(Y56712 - s34) / X56712",
        "  x4 =  X67123 / (Y67123 + s45)  =  -(Y67123 - s45) / X67123",
        "  x5 =  X71234 / (Y71234 + s56)  =  -(Y71234 - s56) / X71234",
        "  x6 =  X12345 / (Y12345 + s67)  =  -(Y12345 - s67) / X12345",
        "  x7 =  X23456 / (Y23456 + s71)  =  -(Y23456 - s71) / X23456",
        "",
        "Set 1 (reverse cycle):",
        "  x1 =  X65432 / (Y65432 + s71)  =  -(Y65432 - s71) / X65432",
        "  x2 =  X76543 / (Y76543 + s12)  =  -(Y76543 - s12) / X76543",
        "  x3 =  X17654 / (Y17654 + s23)  =  -(Y17654 - s23) / X17654",
        "  x4 =  X21765 / (Y21765 + s34)  =  -(Y21765 - s34) / X21765",
        "  x5 =  X32176 / (Y32176 + s45)  =  -(Y32176 - s45) / X32176",
        "  x6 =  X43217 / (Y43217 + s56)  =  -(Y43217 - s56) / X43217",
        "  x7 =  X54321 / (Y54321 + s67)  =  -(Y54321 - s67) / X54321",
        "",
        "Set 2 (barred-variant layout; values equal Set 1 numerically):",
        "  x1 = (X4567 − X̄2)/(Ȳ4567 − Ȳ2) = −(Y4567 + Ȳ2)/(X4567 + X̄2)",
        "  x2 = (X5671 − X̄3)/(Ȳ5671 − Ȳ3) = −(Y5671 + Ȳ3)/(X5671 + X̄3)",
        "  x3 = (X6712 − X̄4)/(Ȳ6712 − Ȳ4) = −(Y6712 + Ȳ4)/(X6712 + X̄4)",
        "  x4 = (X7123 − X̄5)/(Ȳ7123 − Ȳ5) = −(Y7123 + Ȳ5)/(X7123 + X̄5)",
        "  x5 = (X1234 − X̄6)/(Ȳ1234 − Ȳ6) = −(Y1234 + Ȳ6)/(X1234 + X̄6)",
        "  x6 = (X2345 − X̄7)/(Ȳ2345 − Ȳ7) = −(Y2345 + Ȳ7)/(X2345 + X̄7)",
        "  x7 = (X3456 − X̄1)/(Ȳ3456 − Ȳ1) = −(Y3456 + Ȳ1)/(X3456 + X̄1)",
        "",
        "Set 2 (reverse; equal Set 1 numerically):",
        "  x1 = (X5432 − X7)/(Ȳ5432 − Ȳ7) = −(Y5432 + Ȳ7)/(X5432 + X7)",
        "  x2 = (X6543 − X1)/(Ȳ6543 − Ȳ1) = −(Y6543 + Ȳ1)/(X6543 + X1)",
        "  x3 = (X7654 − X2)/(Ȳ7654 − Ȳ2) = −(Y7654 + Ȳ2)/(X7654 + X2)",
        "  x4 = (X1765 − X3)/(Ȳ1765 − Ȳ3) = −(Y1765 + Ȳ3)/(X1765 + X3)",
        "  x5 = (X2176 − X4)/(Ȳ2176 − Ȳ4) = −(Y2176 + Ȳ4)/(X2176 + X4)",
        "  x6 = (X3217 − X5)/(Ȳ3217 − Ȳ5) = −(Y3217 + Ȳ5)/(X3217 + X5)",
        "  x7 = (X4321 − X6)/(Ȳ4321 − Ȳ6) = −(Y4321 + Ȳ6)/(X4321 + X6)",
        "",
        "Set 3 (hatted-variant layout; equal Set 1 numerically):",
        "  x1 = (X567 − X̄32)/(Ȳ567 − Ȳ32) = −(Y567 + Ȳ32)/(X567 + X̄32)",
        "  x2 = (X671 − X̄43)/(Ȳ671 − Ȳ43) = −(Y671 + Ȳ43)/(X671 + X̄43)",
        "  x3 = (X712 − X̄54)/(Ȳ712 − Ȳ54) = −(Y712 + Ȳ54)/(X712 + X̄54)",
        "  x4 = (X123 − X̄65)/(Ȳ123 − Ȳ65) = −(Y123 + Ȳ65)/(X123 + X̄65)",
        "  x5 = (X234 − X̄76)/(Ȳ234 − Ȳ76) = −(Y234 + Ȳ76)/(X234 + X̄76)",
        "  x6 = (X345 − X̄17)/(Ȳ345 − Ȳ17) = −(Y345 + Ȳ17)/(X345 + X̄17)",
        "  x7 = (X456 − X̄21)/(Ȳ456 − Ȳ21) = −(Y456 + Ȳ21)/(X456 + X̄21)",
        "",
        "Set 3 (reverse; equal Set 1 numerically):",
        "  x1 = (X432 − X̄67)/(Ȳ432 − Ȳ67) = −(Y432 + Ȳ67)/(X432 + X̄67)",
        "  x2 = (X543 − X̄71)/(Ȳ543 − Ȳ71) = −(Y543 + Ȳ71)/(X543 + X̄71)",
        "  x3 = (X654 − X̄12)/(Ȳ654 − Ȳ12) = −(Y654 + Ȳ12)/(X654 + X̄12)",
        "  x4 = (X765 − X̄23)/(Ȳ765 − Ȳ23) = −(Y765 + Ȳ23)/(X765 + X̄23)",
        "  x5 = (X176 − X̄34)/(Ȳ176 − Ȳ34) = −(Y176 + Ȳ34)/(X176 + X̄34)",
        "  x6 = (X217 − X̄45)/(Ȳ217 − Ȳ45) = −(Y217 + Ȳ45)/(X217 + X̄45)",
        "  x7 = (X321 − X̄56)/(Ȳ321 − Ȳ56) = −(Y321 + Ȳ56)/(X321 + X̄56)",
    ],
    "Heptagon — Direction Cosines - Spatial Heptagon": [
        "Set 1:",
        "  S1 (0, 0, 1)                       a12 (1, 0, 0)",
        "  S2 (0, −s12, c12)                  a23 (c2,  s2,  0)",
        "  S3 (X2,  Ȳ2,  Z2)                  a34 (W23,  U*231,  U231)",
        "  S4 (X32, Y32, Z32)                 a45 (W34,  U*342,  U342)",
        "  S5 (X432, Y432, Z432)              a56 (W45,  U*453,  U453)",
        "  S6 (X5432, Y5432, Z5432)           a67 (W56,  U*564,  U564)",
        "  S7 (X65432, Y65432, Z65432)        a71 (W67,  U*675,  U675)",

        "Set 2:",
        "  S1 (0, 0, 1)                       a23 (1, 0, 0)",
        "  S2 (0, −s23, c23)                  a34 (c3,  s3,  0)",
        "  S3 (X3,  Ȳ3,  Z3)                  a45 (W34,  U*345,  U345)",
        "  S4 (X43, Y43, Z43)                 a56 (W45,  U*456,  U456)",
        "  S5 (X543, Y543, Z543)              a67 (W56,  U*567,  U567)",
        "  S6 (X6543, Y6543, Z6543)           a71 (W67,  U*671,  U671)",
        "  S7 (X76543, Y76543, Z76543)        a12 (W71,  U*712,  U712)",

        "Set 3:",
        "  S1 (0, 0, 1)                       a34 (1, 0, 0)",
        "  S2 (0, −s34, c34)                  a45 (c4,  s4,  0)",
        "  S3 (X4,  Ȳ4,  Z4)                  a56 (W45,  U*456,  U456)",
        "  S4 (X54, Y54, Z54)                 a67 (W56,  U*567,  U567)",
        "  S5 (X654, Y654, Z654)              a71 (W67,  U*671,  U671)",
        "  S6 (X7654, Y7654, Z7654)           a12 (W71,  U*712,  U712)",
        "  S7 (X17654, Y17654, Z17654)        a23 (W12,  U*123,  U123)",

        "Set 4:",
        "  S1 (0, 0, 1)                       a45 (1, 0, 0)",
        "  S2 (0, −s45, c45)                  a56 (c5,  s5,  0)",
        "  S3 (X5,  Ȳ5,  Z5)                  a67 (W56,  U*567,  U567)",
        "  S4 (X65, Y65, Z65)                 a71 (W67,  U*671,  U671)",
        "  S5 (X765, Y765, Z765)              a12 (W71,  U*712,  U712)",
        "  S6 (X1765, Y1765, Z1765)           a23 (W12,  U*123,  U123)",
        "  S7 (X21765, Y21765, Z21765)        a34 (W23,  U*234,  U234)",

        "Set 5:",
        "  S1 (0, 0, 1)                       a56 (1, 0, 0)",
        "  S2 (0, −s56, c56)                  a67 (c6,  s6,  0)",
        "  S3 (X6,  Ȳ6,  Z6)                  a71 (W67,  U*671,  U671)",
        "  S4 (X76, Y76, Z76)                 a12 (W71,  U*712,  U712)",
        "  S5 (X176, Y176, Z176)              a23 (W12,  U*123,  U123)",
        "  S6 (X2176, Y2176, Z2176)           a34 (W23,  U*234,  U234)",
        "  S7 (X32176, Y32176, Z32176)        a45 (W34,  U*345,  U345)",

        "Set 6:",
        "  S1 (0, 0, 1)                       a67 (1, 0, 0)",
        "  S2 (0, −s67, c67)                  a71 (c7,  s7,  0)",
        "  S3 (X7,  Ȳ7,  Z7)                  a12 (W71,  U*712,  U712)",
        "  S4 (X17, Y17, Z17)                 a23 (W12,  U*123,  U123)",
        "  S5 (X217, Y217, Z217)              a34 (W23,  U*234,  U234)",
        "  S6 (X3217, Y3217, Z3217)           a45 (W34,  U*345,  U345)",
        "  S7 (X43217, Y43217, Z43217)        a56 (W45,  U*456,  U456)",

        "Set 7:",
        "  S1 (0, 0, 1)                       a71 (1, 0, 0)",
        "  S2 (0, −s71, c71)                  a12 (c1,  s1,  0)",
        "  S3 (X1,  Ȳ1,  Z1)                  a23 (W12,  U*123,  U123)",
        "  S4 (X21, Y21, Z21)                 a34 (W23,  U*234,  U234)",
        "  S5 (X321, Y321, Z321)              a45 (W34,  U*345,  U345)",
        "  S6 (X4321, Y4321, Z4321)           a56 (W45,  U*456,  U456)",
        "  S7 (X54321, Y54321, Z54321)        a67 (W56,  U*567,  U567)",

        "Set 8:",
        "  S1 (0, 0, 1)                       a71 (1, 0, 0)",
        "  S2 (0,  s71,  c71)                 a67 (c7, −s7,  0)",
        "  S3 (X7, −Y7,  Z7)                  a56 (W67,  −U*671,  U671)",
        "  S4 (X76, −Y76, Z76)                a45 (W56,  −U*567,  U567)",
        "  S5 (X176, −Y176, Z176)             a34 (W45,  −U*456,  U456)",
        "  S6 (X2176, −Y2176, Z2176)          a23 (W34,  −U*345,  U345)",
        "  S7 (X32176, −Y32176, Z32176)       a12 (W23,  −U*234,  U234)",

        "Set 9:",
        "  S1 (0, 0, 1)                       a67 (1, 0, 0)",
        "  S2 (0,  s67,  c67)                 a56 (c6, −s6,  0)",
        "  S3 (X6, −Y6,  Z6)                  a45 (W56,  −U*567,  U567)",
        "  S4 (X65, −Y65, Z65)                a34 (W45,  −U*456,  U456)",
        "  S5 (X765, −Y765, Z765)             a23 (W34,  −U*345,  U345)",
        "  S6 (X1765, −Y1765, Z1765)          a12 (W23,  −U*234,  U234)",
        "  S7 (X21765, −Y21765, Z21765)       a71 (W12,  −U*123,  U123)",

        "Set 10:",
        "  S1 (0, 0, 1)                       a56 (1, 0, 0)",
        "  S2 (0,  s56,  c56)                 a45 (c5, −s5,  0)",
        "  S3 (X5, −Y5,  Z5)                  a34 (W45,  −U*456,  U456)",
        "  S4 (X54, −Y54, Z54)                a23 (W34,  −U*345,  U345)",
        "  S5 (X654, −Y654, Z654)             a12 (W23,  −U*234,  U234)",
        "  S6 (X7654, −Y7654, Z7654)          a71 (W12,  −U*123,  U123)",
        "  S7 (X17654, −Y17654, Z17654)       a67 (W71,  −U*712,  U712)",

        "Set 11:",
        "  S1 (0, 0, 1)                       a45 (1, 0, 0)",
        "  S2 (0,  s45,  c45)                 a34 (c4, −s4,  0)",
        "  S3 (X4, −Y4,  Z4)                  a23 (W34,  −U*345,  U345)",
        "  S4 (X43, −Y43, Z43)                a12 (W23,  −U*234,  U234)",
        "  S5 (X543, −Y543, Z543)             a71 (W12,  −U*123,  U123)",
        "  S6 (X6543, −Y6543, Z6543)          a67 (W71,  −U*712,  U712)",
        "  S7 (X76543, −Y76543, Z76543)       a56 (W67,  −U*671,  U671)",

        "Set 12:",
        "  S1 (0, 0, 1)                       a34 (1, 0, 0)",
        "  S2 (0,  s34,  c34)                 a23 (c3, −s3,  0)",
        "  S3 (X3, −Y3,  Z3)                  a12 (W23,  −U*234,  U234)",
        "  S4 (X32, −Y32, Z32)                a71 (W12,  −U*123,  U123)",
        "  S5 (X432, −Y432, Z432)             a67 (W71,  −U*712,  U712)",
        "  S6 (X5432, −Y5432, Z5432)          a56 (W67,  −U*671,  U671)",
        "  S7 (X65432, −Y65432, Z65432)       a45 (W56,  −U*567,  U567)",

        "Set 13:",
        "  S1 (0, 0, 1)                       a23 (1, 0, 0)",
        "  S2 (0,  s23,  c23)                 a12 (c2, −s2,  0)",
        "  S3 (X2, −Y2,  Z2)                  a71 (W12,  −U*123,  U123)",
        "  S4 (X21, −Y21, Z21)                a67 (W71,  −U*712,  U712)",
        "  S5 (X321, −Y321, Z321)             a56 (W67,  −U*671,  U671)",
        "  S6 (X4321, −Y4321, Z4321)          a45 (W56,  −U*567,  U567)",
        "  S7 (X54321, −Y54321, Z54321)       a34 (W45,  −U*456,  U456)",

        "Set 14:",
        "  S1 (0, 0, 1)                       a12 (1, 0, 0)",
        "  S2 (0,  s12,  c12)                 a71 (c1, −s1,  0)",
        "  S3 (X1, −Y1,  Z1)                  a67 (W71,  −U*712,  U712)",
        "  S4 (X12, −Y12, Z12)                a56 (W67,  −U*671,  U671)",
        "  S5 (X312, −Y312, Z312)             a45 (W56,  −U*567,  U567)",
        "  S6 (X4312, −Y4312, Z4312)          a34 (W45,  −U*456,  U456)",
        "  S7 (X54312, −Y54312, Z54312)       a23 (W34,  −U*345,  U345)",
    ],

}


def get_all_appendix_table_names() -> List[str]:
    # stable sort for nice UI ordering
    return sorted(_APP_TEXT.keys(), key=str.lower)

def get_appendix_table(name: str) -> List[str]:
    if name not in _APP_TEXT:
        raise KeyError(f"Unknown Appendix table: {name}")
    # return a shallow copy for safety
    return list(_APP_TEXT[name])

_TABLES: Dict[str, Dict[str, Callable[[Dict[str,float], bool], List[Dict[str, Any]]]]] = {
    'Triangle': {
        'Equations for a Spherical Triangle': table_triangle_fundamental,
        'Equations for a Polar Triangle': table_triangle_polar,
        'Direction Cosines - Spatial Triangle': table_triangle_direction_cosines,
    },
    'Quadrilateral': {
        'Equations for a Spherical Quadrilateral': table_quad_fundamental,
        'Equations for a Polar Quadrilateral': table_quad_polar,
        'Half-Tangent Laws for a Quadrilateral': table_quad_half_tangent,
        'Direction Cosines - Spatial Quadrilateral': table_quad_direction_cosines,
    },
    'Pentagon': {
        'Equations for a Spherical Pentagon': table_penta_fundamental,
        'Equations for a Polar Pentagon': table_penta_polar,
        'Direction Cosines - Spatial Pentagon': table_penta_direction_cosines,
        'Half-Tangent Laws for a Pentagon': table_penta_half_tangent,
    },
    'Hexagon': {
        'Equations for a Spherical Hexagon': table_hexa_fundamental,
        'Equations for a Polar Hexagon': table_hexa_polar,
        'Direction Cosines - Spatial Hexagon': table_hexa_direction_cosines,
        'Half-Tangent Laws for a Hexagon': table_hexa_half_tangent,
    },
    'Heptagon': {
        'Equations for a Spherical Heptagon': table_hepta_fundamental,
        'Equations for a Polar Heptagon': table_hepta_polar,
        'Half-Tangent Laws for a Heptagon': table_hepta_half_tangent,
        'Direction Cosines - Spatial Heptagon': table_hepta_direction_cosines,
    },

}

# ------------------------------ Examples ------------------------------------

_EXAMPLES: Dict[str, Dict[str, Dict[str, float]]] = {
    'Triangle': {
        # From Section 6.6 slide sample (alpha12=120, alpha23=80, alpha31=135)
        'Triangle (Solution A)': {
            'alpha12': 120.0, 'alpha23': 80.0, 'alpha31': 135.0,
            'theta1': 72.92, 'theta2': 43.34, 'theta3': 57.20,
        },
        'Triangle (Solution B)': {
            'alpha12': 120.0, 'alpha23': 80.0, 'alpha31': 135.0,
            'theta1': 287.08, 'theta2': 316.66, 'theta3': 302.79,
        },
    },
    'Quadrilateral': {
        # Provide a simple synthetic example
        'Quad (demo)': {
            'alpha12': 40.0, 'alpha23': 60.0, 'alpha34': 70.0, 'alpha41': 85.0,
            'theta1': 15.0, 'theta2': 30.0, 'theta3': 45.0, 'theta4': 20.0,
        },
    },
    'Pentagon': {
        'Pentagon (demo)': {
            'alpha12': 40.0, 'alpha23': 60.0, 'alpha34': 80.0, 'alpha45': 60.0, 'alpha51': 90.0,
            'theta1': 15.0, 'theta2': 30.0, 'theta3': 100.0, 'theta4': 45.0, 'theta5': 65.0,
        },
        'Pentagon Table 6.5 (given θ3=100, θ5=65)': {
            'alpha12': 40.0, 'alpha23': 60.0, 'alpha34': 80.0, 'alpha45': 60.0, 'alpha51': 90.0,
            'theta1': 10.0, 'theta2': 25.0, 'theta3': 100.0, 'theta4': 35.0, 'theta5': 65.0,
        },
    },
    'Hexagon': {
        'Hexagon (demo)': {
            'alpha12': 35.0, 'alpha23': 50.0, 'alpha34': 65.0, 'alpha45': 70.0, 'alpha56': 55.0, 'alpha61': 85.0,
            'theta1': 18.0, 'theta2': 33.0, 'theta3': 96.0, 'theta4': 41.0, 'theta5': 62.0, 'theta6': 77.0,
        },
    },
    'Heptagon': {
        'Heptagon (demo)': {
            'alpha12': 35.0, 'alpha23': 50.0, 'alpha34': 60.0, 'alpha45': 55.0,
            'alpha56': 45.0, 'alpha67': 40.0, 'alpha71': 65.0,
            'theta1': 10.0, 'theta2': 35.0, 'theta3': 60.0, 'theta4': 85.0,
            'theta5': 40.0, 'theta6': 75.0, 'theta7': 25.0,
        },
    },
}

# ------------------------------ Public API ----------------------------------

def get_available_shapes() -> List[str]:
    return list(_TABLES.keys())

def get_available_tables(shape: str) -> List[str]:
    return list(_TABLES.get(shape, {}).keys())

def get_shape_examples(shape: str) -> Dict[str, Dict[str, float]]:
    return _EXAMPLES.get(shape, {})

def evaluate_table(shape: str, table_key: str, inputs: Dict[str, float], degrees: bool=True) -> List[Dict[str, Any]]:
    if shape not in _TABLES:
        raise KeyError(f"Unknown shape: {shape}")
    if table_key not in _TABLES[shape]:
        raise KeyError(f"Unknown table for {shape}: {table_key}")
    fn = _TABLES[shape][table_key]
    return fn(inputs, degrees)

if __name__ == "__main__":
    # simple smoke test
    tri = _EXAMPLES['Triangle']['Triangle (Solution A)']
    out = evaluate_table('Triangle', 'Equations for a Spherical Triangle', tri, degrees=True)
    print("Triangle – Equations for a Spherical Triangle:")
    for row in out:
        print(f"  {row['eq']:<18} => {row['value']:.6f}")
