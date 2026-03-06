#!/usr/bin/env python3
"""
35_stress_test_duality.py — Stress test: spectral duality on Dic_n, n=2..101.

Item 19 from tracker_v2: aggressive test of duality on borderline cases.

  - n odd (unobstructed): duality should FAIL
  - n even (obstructed): duality should HOLD
  - Full check: m(χ,j) + m(χ,j*-j) = 1 for ALL j ∈ [0,j*], not just j=0

Note: elements are built ANALYTICALLY (a^k = (cos kπ/n, sin kπ/n, 0, 0)),
not by iterative qmul. Iterative construction causes float drift at large n
(n≥77), corrupting conjugacy class structure via phantom qkeys.

Tests:
  [1] A₀ duality: Dic_n, n=2..101 (100 groups, ~5000 j-pairs)
  [2] Binary multiplicities + zero fraction (obstructed groups)
  [3] ALL scalar characters: Dic_n, n=2..30
  [4] j* is last zero: m(A₀,j) ≥ 1 for j > j*
  [5] Triple equivalence (A)⟺(B)⟺(C): n=2..50

RAW OUTPUT:
===========

======================================================================
[1] A₀ DUALITY: Dic_n, n = 2..101
    Full check: m(A₀,j) + m(A₀,j*-j) = 1 for ALL j ∈ [0,j*]
======================================================================

  Dic_  2 (|G|=   8, j*=  1): duality HOLD (even→HOLD) ✓
  Dic_  3 (|G|=  12, j*=  2): duality FAIL (odd→FAIL) ✓
  Dic_  4 (|G|=  16, j*=  3): duality HOLD (even→HOLD) ✓
  Dic_  5 (|G|=  20, j*=  4): duality FAIL (odd→FAIL) ✓
  Dic_  6 (|G|=  24, j*=  5): duality HOLD (even→HOLD) ✓
  Dic_  7 (|G|=  28, j*=  6): duality FAIL (odd→FAIL) ✓
  Dic_  8 (|G|=  32, j*=  7): duality HOLD (even→HOLD) ✓
  Dic_  9 (|G|=  36, j*=  8): duality FAIL (odd→FAIL) ✓
  Dic_ 10 (|G|=  40, j*=  9): duality HOLD (even→HOLD) ✓
  Dic_ 11 (|G|=  44, j*= 10): duality FAIL (odd→FAIL) ✓
  Dic_ 12 (|G|=  48, j*= 11): duality HOLD (even→HOLD) ✓
  Dic_ 20 (|G|=  80, j*= 19): duality HOLD (even→HOLD) ✓
  Dic_ 40 (|G|= 160, j*= 39): duality HOLD (even→HOLD) ✓
  Dic_ 60 (|G|= 240, j*= 59): duality HOLD (even→HOLD) ✓
  Dic_ 80 (|G|= 320, j*= 79): duality HOLD (even→HOLD) ✓
  Dic_100 (|G|= 400, j*= 99): duality HOLD (even→HOLD) ✓
  Dic_101 (|G|= 404, j*=100): duality FAIL (odd→FAIL) ✓

  ✓ ALL 100 groups pass (100 asserts). Total j-pairs: 5150
  Numerical: max |m - round(m)| = 2.7e-03 (at Dic_91)

======================================================================
[2] BINARY MULTIPLICITIES + ZERO FRACTION (n even)
======================================================================

  Dic_  2: j*=  1, m ∈ {0,1} ✓, zeros = 1/2 = 0.500 ✓
  Dic_  4: j*=  3, m ∈ {0,1} ✓, zeros = 2/4 = 0.500 ✓
  Dic_  6: j*=  5, m ∈ {0,1} ✓, zeros = 3/6 = 0.500 ✓
  Dic_  8: j*=  7, m ∈ {0,1} ✓, zeros = 4/8 = 0.500 ✓
  Dic_ 10: j*=  9, m ∈ {0,1} ✓, zeros = 5/10 = 0.500 ✓
  Dic_ 20: j*= 19, m ∈ {0,1} ✓, zeros = 10/20 = 0.500 ✓
  Dic_ 30: j*= 29, m ∈ {0,1} ✓, zeros = 15/30 = 0.500 ✓
  Dic_ 50: j*= 49, m ∈ {0,1} ✓, zeros = 25/50 = 0.500 ✓
  Dic_100: j*= 99, m ∈ {0,1} ✓, zeros = 50/100 = 0.500 ✓

======================================================================
[3] ALL SCALAR CHARACTERS: Dic_n, n = 2..30
    χ_{α,β}: χ(a^k) = α^k, χ(a^k·b) = α^k·β
    Scalar ⟺ χ(-1) = α^n = 1
======================================================================

  Dic_ 2: 4 scalar chars, duality all HOLD ✓
  Dic_ 3: 2 scalar chars, duality all FAIL ✓
  Dic_ 4: 4 scalar chars, duality all HOLD ✓
  Dic_ 5: 2 scalar chars, duality all FAIL ✓
  Dic_ 6: 4 scalar chars, duality all HOLD ✓
  Dic_ 7: 2 scalar chars, duality all FAIL ✓
  Dic_ 8: 4 scalar chars, duality all HOLD ✓
  Dic_ 9: 2 scalar chars, duality all FAIL ✓
  Dic_10: 4 scalar chars, duality all HOLD ✓
  Dic_11: 2 scalar chars, duality all FAIL ✓
  Dic_12: 4 scalar chars, duality all HOLD ✓
  Dic_13: 2 scalar chars, duality all FAIL ✓
  Dic_14: 4 scalar chars, duality all HOLD ✓
  Dic_15: 2 scalar chars, duality all FAIL ✓
  Dic_16: 4 scalar chars, duality all HOLD ✓
  Dic_17: 2 scalar chars, duality all FAIL ✓
  Dic_18: 4 scalar chars, duality all HOLD ✓
  Dic_19: 2 scalar chars, duality all FAIL ✓
  Dic_20: 4 scalar chars, duality all HOLD ✓
  Dic_21: 2 scalar chars, duality all FAIL ✓
  Dic_22: 4 scalar chars, duality all HOLD ✓
  Dic_23: 2 scalar chars, duality all FAIL ✓
  Dic_24: 4 scalar chars, duality all HOLD ✓
  Dic_25: 2 scalar chars, duality all FAIL ✓
  Dic_26: 4 scalar chars, duality all HOLD ✓
  Dic_27: 2 scalar chars, duality all FAIL ✓
  Dic_28: 4 scalar chars, duality all HOLD ✓
  Dic_29: 2 scalar chars, duality all FAIL ✓
  Dic_30: 4 scalar chars, duality all HOLD ✓

  ✓ 88 character tests. Total character×j pairs: 1408

======================================================================
[4] j* IS LAST ZERO (n even, check j ∈ (j*, j*+50])
======================================================================

  Dic_ 2: m(j*=1)=0, m(j)≥1 for j ∈ [2, 51] ✓
  Dic_ 4: m(j*=3)=0, m(j)≥1 for j ∈ [4, 53] ✓
  Dic_ 6: m(j*=5)=0, m(j)≥1 for j ∈ [6, 55] ✓
  Dic_ 8: m(j*=7)=0, m(j)≥1 for j ∈ [8, 57] ✓
  Dic_10: m(j*=9)=0, m(j)≥1 for j ∈ [10, 59] ✓
  Dic_20: m(j*=19)=0, m(j)≥1 for j ∈ [20, 69] ✓
  Dic_50: m(j*=49)=0, m(j)≥1 for j ∈ [50, 99] ✓

======================================================================
[5] TRIPLE EQUIVALENCE (A)⟺(B)⟺(C): n=2..50
    (A) -1 ∈ [G,G]
    (B) |H|/ord(h) even for all h ∈ H\{e}
    (C) spectral duality
======================================================================

  Dic_ 2: (A)=True   (B)=True   (C)=True  ✓
  Dic_ 3: (A)=False  (B)=False  (C)=False ✓
  Dic_ 4: (A)=True   (B)=True   (C)=True  ✓
  Dic_ 5: (A)=False  (B)=False  (C)=False ✓
  Dic_ 6: (A)=True   (B)=True   (C)=True  ✓
  Dic_ 7: (A)=False  (B)=False  (C)=False ✓
  Dic_ 8: (A)=True   (B)=True   (C)=True  ✓
  Dic_ 9: (A)=False  (B)=False  (C)=False ✓
  Dic_10: (A)=True   (B)=True   (C)=True  ✓
  Dic_11: (A)=False  (B)=False  (C)=False ✓
  Dic_12: (A)=True   (B)=True   (C)=True  ✓
  Dic_20: (A)=True   (B)=True   (C)=True  ✓
  Dic_30: (A)=True   (B)=True   (C)=True  ✓
  Dic_40: (A)=True   (B)=True   (C)=True  ✓
  Dic_50: (A)=True   (B)=True   (C)=True  ✓

  ✓ All 49 groups: (A)⟺(B)⟺(C) (49 asserts)

======================================================================
SUMMARY
======================================================================

  Stress test: 100 binary dihedral groups Dic_n (n=2..101)

  [1] A₀ duality (n=2..101):      100 asserts, 5150 j-pairs
  [2] Binary mult + zero frac:     239 asserts
  [3] All scalar chars (n≤30):      88 asserts, 1408 char×j pairs
  [4] j* last zero:                357 asserts
  [5] (A)⟺(B)⟺(C) (n≤50):         49 asserts
                                TOTAL:  833 asserts

  No borderline failures. Duality ⟺ n even, exhaustively.

All checks passed.
"""
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))

import numpy as np
from src.quaternion import qmul, qkey
from src.group import compute_conjugacy_classes, compute_commutator_subgroup

MINUS1_KEY = qkey(np.array([-1, 0, 0, 0], dtype=float))


def chi_su2(j, trace):
    """SU(2) spin-j character at given trace = 2cos(θ/2)."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def build_dicyclic(n):
    """Build Dic_n (= binary dihedral 2D_n) of order 4n directly.

    Elements: a^k (k=0..2n-1) and a^k·b (k=0..2n-1)
    where a = (cos π/n, sin π/n, 0, 0), b = (0,0,1,0).

    a^k computed ANALYTICALLY (not iteratively) to avoid float drift.
    a^k = (cos(kπ/n), sin(kπ/n), 0, 0)
    a^k·b = (0, 0, cos(kπ/n), sin(kπ/n))

    Returns (elements, k_vals, types):
      k_vals[i] = power of a for element i
      types[i] = 0 (a-type: a^k) or 1 (b-type: a^k·b)
    """
    elements = []
    k_vals = []
    types = []

    for k in range(2 * n):
        angle = k * np.pi / n
        c, s = np.cos(angle), np.sin(angle)

        # a^k = (cos(kπ/n), sin(kπ/n), 0, 0)
        elements.append(np.array([c, s, 0, 0], dtype=float))
        k_vals.append(k)
        types.append(0)

        # a^k · b = (0, 0, cos(kπ/n), sin(kπ/n))
        elements.append(np.array([0, 0, c, s], dtype=float))
        k_vals.append(k)
        types.append(1)

    return np.array(elements), np.array(k_vals), np.array(types)


def m_A0(j, classes, G):
    """Multiplicity of trivial character A₀ in spin-j rep."""
    return sum(c['size'] * chi_su2(j, c['trace']) for c in classes) / G


# ============================================================
# [1] A₀ duality: Dic_n, n = 2..101
# ============================================================

print("=" * 70)
print("[1] A₀ DUALITY: Dic_n, n = 2..101")
print("    Full check: m(A₀,j) + m(A₀,j*-j) = 1 for ALL j ∈ [0,j*]")
print("=" * 70)
print()

n_total_pairs = 0
n_assert_1 = 0
global_max_err = 0.0

for n in range(2, 102):
    elements, _, _ = build_dicyclic(n)
    G = len(elements)
    assert G == 4 * n, f"Dic_{n}: |G|={G} ≠ {4*n}"

    j_star = n - 1
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])

    # Full duality check with numerical error tracking
    duality_ok = True
    group_max_err = 0.0
    for j in range(j_star + 1):
        raw_j = m_A0(j, classes, G)
        raw_d = m_A0(j_star - j, classes, G)
        group_max_err = max(group_max_err,
                           abs(raw_j - round(raw_j)),
                           abs(raw_d - round(raw_d)))
        m_j = int(round(raw_j))
        m_d = int(round(raw_d))
        if m_j + m_d != 1:
            duality_ok = False
            break

    if group_max_err > global_max_err:
        global_max_err = group_max_err
        worst_n = n
    global_max_err = max(global_max_err, group_max_err)

    expected = (n % 2 == 0)
    assert duality_ok == expected, (
        f"Dic_{n}: duality={'HOLD' if duality_ok else 'FAIL'}, "
        f"expected={'HOLD' if expected else 'FAIL'}")
    n_assert_1 += 1

    n_total_pairs += j_star + 1

    if n <= 12 or n % 20 == 0 or n == 101:
        status = "HOLD" if duality_ok else "FAIL"
        exp_str = "even→HOLD" if expected else "odd→FAIL"
        print(f"  Dic_{n:>3d} (|G|={G:>4d}, j*={j_star:>3d}): "
              f"duality {status} ({exp_str}) ✓")

print()
print(f"  ✓ ALL 100 groups pass ({n_assert_1} asserts). "
      f"Total j-pairs: {n_total_pairs}")
print(f"  Numerical: max |m - round(m)| = {global_max_err:.1e} (at Dic_{worst_n})")
assert global_max_err < 0.01, f"Numerical error too large: {global_max_err:.1e}"

# ============================================================
# [2] Binary multiplicities + zero fraction (obstructed)
# ============================================================

print()
print("=" * 70)
print("[2] BINARY MULTIPLICITIES + ZERO FRACTION (n even)")
print("=" * 70)
print()

n_assert_2 = 0

for n in [2, 4, 6, 8, 10, 20, 30, 50, 100]:
    elements, _, _ = build_dicyclic(n)
    G = len(elements)
    j_star = n - 1
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])

    n_zeros = 0
    for j in range(j_star + 1):
        m = int(round(m_A0(j, classes, G)))
        assert m in (0, 1), f"Dic_{n}: m(A₀,{j}) = {m}"
        n_assert_2 += 1
        if m == 0:
            n_zeros += 1

    frac = n_zeros / (j_star + 1)
    assert abs(frac - 0.5) < 1e-10, f"Dic_{n}: zero fraction = {frac}"
    n_assert_2 += 1

    print(f"  Dic_{n:>3d}: j*={j_star:>3d}, m ∈ {{0,1}} ✓, "
          f"zeros = {n_zeros}/{j_star+1} = 0.500 ✓")

# ============================================================
# [3] ALL scalar characters: Dic_n, n = 2..30
# ============================================================

print()
print("=" * 70)
print("[3] ALL SCALAR CHARACTERS: Dic_n, n = 2..30")
print("    χ_{α,β}: χ(a^k) = α^k, χ(a^k·b) = α^k·β")
print("    Scalar ⟺ χ(-1) = α^n = 1")
print("=" * 70)
print()

# Scalar characters of Dic_n:
#   α = χ(a) ∈ {±1}  (since a² ∈ [G,G])
#   β = χ(b) with β² = α^n
#   Scalar: α^n = 1
#
#   n even: α=1 → β=±1; α=-1 → (-1)^n=1 → β=±1. Total: 4 chars.
#   n odd:  α=1 → β=±1; α=-1 → (-1)^n=-1 → β=±i (not real). Total: 2 chars.
#   (For n odd, only α=1 gives scalar characters.)

n_char_pairs = 0
n_assert_3 = 0

for n in range(2, 31):
    elements, k_vals, types = build_dicyclic(n)
    G = len(elements)
    j_star = n - 1
    expected = (n % 2 == 0)

    scalar_chars = [(1, 1), (1, -1)]
    if n % 2 == 0:
        scalar_chars += [(-1, 1), (-1, -1)]

    all_ok = True
    for alpha, beta in scalar_chars:
        # Character values on each element
        chi_vals = np.zeros(G)
        for i in range(G):
            k = k_vals[i]
            t = types[i]
            chi_vals[i] = (alpha ** k) * (beta ** t)

        # Full duality check for this character
        duality_ok = True
        for j in range(j_star + 1):
            m_j = 0.0
            m_d = 0.0
            for i in range(G):
                trace = 2 * elements[i][0]
                cv = chi_vals[i]
                m_j += cv * chi_su2(j, trace)
                m_d += cv * chi_su2(j_star - j, trace)
            m_j = int(round(m_j / G))
            m_d = int(round(m_d / G))
            if m_j + m_d != 1:
                duality_ok = False
                break

        assert duality_ok == expected, (
            f"Dic_{n} χ({alpha},{beta}): "
            f"duality={'HOLD' if duality_ok else 'FAIL'}, "
            f"expected={'HOLD' if expected else 'FAIL'}")
        n_assert_3 += 1
        n_char_pairs += j_star + 1

        if not (duality_ok == expected):
            all_ok = False

    n_chars = len(scalar_chars)
    status = "all HOLD" if expected else "all FAIL"
    print(f"  Dic_{n:>2d}: {n_chars} scalar chars, duality {status} ✓")

print()
print(f"  ✓ {n_assert_3} character tests. "
      f"Total character×j pairs: {n_char_pairs}")

# ============================================================
# [4] j* is last zero (obstructed, selected n)
# ============================================================

print()
print("=" * 70)
print("[4] j* IS LAST ZERO (n even, check j ∈ (j*, j*+50])")
print("=" * 70)
print()

J_EXTRA = 50
n_assert_4 = 0

for n in [2, 4, 6, 8, 10, 20, 50]:
    elements, _, _ = build_dicyclic(n)
    G = len(elements)
    j_star = n - 1
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])

    # j* is a zero
    assert int(round(m_A0(j_star, classes, G))) == 0
    n_assert_4 += 1

    # All j > j* have m ≥ 1
    for j in range(j_star + 1, j_star + J_EXTRA + 1):
        m = int(round(m_A0(j, classes, G)))
        assert m >= 1, f"Dic_{n}: m(A₀,{j})=0 for j>{j_star}!"
        n_assert_4 += 1

    print(f"  Dic_{n:>2d}: m(j*={j_star})=0, "
          f"m(j)≥1 for j ∈ [{j_star+1}, {j_star+J_EXTRA}] ✓")

# ============================================================
# [5] Triple equivalence (A)⟺(B)⟺(C): n=2..50
# ============================================================

print()
print("=" * 70)
print("[5] TRIPLE EQUIVALENCE (A)⟺(B)⟺(C): n=2..50")
print("    (A) -1 ∈ [G,G]")
print("    (B) |H|/ord(h) even for all h ∈ H\\{e}")
print("    (C) spectral duality")
print("=" * 70)
print()

n_assert_5 = 0

for n in range(2, 51):
    elements, _, _ = build_dicyclic(n)
    G = len(elements)
    j_star = n - 1

    # (A) obstruction
    comm = compute_commutator_subgroup(elements)
    A = MINUS1_KEY in set(comm.keys())

    # (B) divisibility: H = D_n, |H| = 2n
    # Element orders in D_n: divisors of n (rotations) + order 2 (reflections)
    # |H|/d = 2n/d. For d | n: 2n/d. For d = 2: n.
    # (B) holds iff ALL ratios even.
    # Rotations: 2n/d is always even (since 2n/d = 2·(n/d)).
    # Reflections: 2n/2 = n. Even iff n even.
    # So (B) ⟺ n even.
    H = 2 * n
    so3_orders = {2}  # reflections
    for d in range(2, n + 1):
        if n % d == 0:
            so3_orders.add(d)
    B = all((H // d) % 2 == 0 for d in so3_orders)

    # (C) duality (already checked in [1] but verify again)
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    C_ok = True
    for j in range(j_star + 1):
        m_j = int(round(m_A0(j, classes, G)))
        m_d = int(round(m_A0(j_star - j, classes, G)))
        if m_j + m_d != 1:
            C_ok = False
            break

    assert A == B == C_ok, (
        f"Dic_{n}: (A)={A}, (B)={B}, (C)={C_ok}")
    n_assert_5 += 1

    if n <= 12 or n % 10 == 0:
        print(f"  Dic_{n:>2d}: (A)={str(A):5s}  (B)={str(B):5s}  "
              f"(C)={str(C_ok):5s} ✓")

print()
print(f"  ✓ All 49 groups: (A)⟺(B)⟺(C) ({n_assert_5} asserts)")

# ============================================================
# Summary
# ============================================================

total_asserts = n_assert_1 + n_assert_2 + n_assert_3 + n_assert_4 + n_assert_5

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("  Stress test: 100 binary dihedral groups Dic_n (n=2..101)")
print()
print(f"  [1] A₀ duality (n=2..101):     {n_assert_1:>4d} asserts, "
      f"{n_total_pairs} j-pairs")
print(f"  [2] Binary mult + zero frac:    {n_assert_2:>4d} asserts")
print(f"  [3] All scalar chars (n≤30):    {n_assert_3:>4d} asserts, "
      f"{n_char_pairs} char×j pairs")
print(f"  [4] j* last zero:               {n_assert_4:>4d} asserts")
print(f"  [5] (A)⟺(B)⟺(C) (n≤50):       {n_assert_5:>4d} asserts")
print(f"  {'TOTAL':>35s}: {total_asserts:>4d} asserts")
print()
print("  No borderline failures. Duality ⟺ n even, exhaustively.")
print()
print("All checks passed.")
