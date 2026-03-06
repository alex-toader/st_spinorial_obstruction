#!/usr/bin/env python3
"""
C4: Discrete anomaly interpretation + C5: Scope beyond SU(2).

PART A — DISCRETE ANOMALY INTERPRETATION
=========================================

The obstruction admits a cohomological interpretation as a class in
H²(G/Z₂; Z₂), measuring whether the central extension splits.

Setup: For G ⊂ SU(2) finite, the exact sequence
  1 → Z₂ → G → G̅ → 1   (where G̅ = G/{±1} ⊂ SO(3))

defines a class [ξ] ∈ H²(G̅; Z₂). The extension splits iff [ξ] = 0.

The obstruction -1 ∈ [G,G] is equivalent to:
  -1 maps to trivial element in G^{ab} = G/[G,G]

This means: the Z₂ central extension, when restricted to [G,G], is trivial.
But the extension itself is NON-TRIVIAL (as G ≠ G̅ × Z₂ in general).

The "anomaly" interpretation:
  - The Z₂ center symmetry of SU(2) acts on the space of 1D chars of G
  - χ(-1) = ±1 labels the Z₂ representation
  - When -1 ∈ [G,G]: ALL 1D chars have χ(-1) = +1
  - This means: the Z₂ center is "screened" — it acts trivially on ALL
    scalar quantization sectors
  - Analogy with 't Hooft anomaly: a center symmetry that cannot be
    "detected" by local probes (here: 1D characters)

This is NOT a genuine 't Hooft anomaly in the QFT sense (no gauge fields,
no 4-form). It is a STRUCTURAL analogy. The obstruction is better described
as a cohomological property of the central extension.

Concrete numerical check: verify that the extension class corresponds to
-1 ∈ [G,G] for all groups.

PART B — SCOPE BEYOND SU(2)
============================

Question: does Q₈ ⊄ G ⟹ -1 ∉ [G,G] hold beyond finite SU(2) subgroups?

Answer: the implication Q₈ ⊂ G ⟹ -1 ∈ [G,G] is ALWAYS true (for any
group with a central element -1, if Q₈ embeds then [i,j] = -1).

The reverse -1 ∈ [G,G] ⟹ Q₈ ⊂ G uses the SU(2) structure crucially:
  - The proof requires a,b to be UNIT QUATERNIONS
  - If G is not in SU(2), [a,b] = -1 does not force a,b to be pure imaginary

COUNTEREXAMPLE: Take G = SL(2,3) (order 24, = 2T as abstract group).
  Embed in GL(2,C) rather than SU(2).
  -1 ∈ [G,G] still holds (group-theoretic property).
  Q₈ ⊂ G still holds (subgroup structure preserved).
  So the equivalence holds HERE, but for a different reason (abstract group).

The real question is: can we find G with a central element z of order 2,
z ∈ [G,G], but Q₈ ⊄ G?

COUNTEREXAMPLE: G = Z₂ × Z₂ × Z₂, z = (1,1,1).
  [G,G] = {e} since G is abelian. z ∉ [G,G]. No counterexample.

COUNTEREXAMPLE: G = D₄ (dihedral of order 8, NOT in SU(2)).
  [D₄, D₄] = Z₂ = center. So -1 ∈ [D₄, D₄]. But Q₈ ⊄ D₄?
  D₄ has elements of order 1,2,4. Q₈ has elements of order 1,2,4.
  |D₄| = 8 = |Q₈|. If Q₈ ⊂ D₄, then Q₈ = D₄. But D₄ ≇ Q₈.
  D₄ has 5 elements of order 2; Q₈ has only 1 (-1).
  So Q₈ ⊄ D₄, yet -1 ∈ [D₄, D₄]. COUNTEREXAMPLE FOUND.

CONCLUSION: The equivalence Q₈ ⊂ G ⟺ -1 ∈ [G,G] is SPECIFIC to
finite subgroups of SU(2). Beyond SU(2), the Q₈ criterion fails.
The SU(2) structure (quaternion algebra, unit norm) is essential.

RAW OUTPUT:
===========

======================================================================
PART A: CENTRAL EXTENSION AND Z₂ SCREENING
======================================================================

For G ⊂ SU(2) finite, the central extension:
  1 → Z₂ → G → G̅ = G/{±1} → 1

The image of -1 in G^{ab} = G/[G,G] determines whether
the Z₂ center is 'visible' to scalar quantization:
  -1 ∉ [G,G]: image nontrivial → Z₂ detectable → spinorial sectors exist
  -1 ∈ [G,G]: image trivial → Z₂ screened → no spinorial sectors

   Group  |G|  |[G,G]|   -1∈[G,G]  Z₂ screened  spinorial
  --------------------------------------------------------
      Z₆    6        1      False           no        yes  ✓
    Dic₃   12        3      False           no        yes  ✓
    Dic₅   20        5      False           no        yes  ✓
      Q₈    8        2       True          yes         no  ✓
    Dic₄   16        4       True          yes         no  ✓
      2T   24        8       True          yes         no  ✓
      2O   48       24       True          yes         no  ✓
      2I  120      120       True          yes         no  ✓

Pattern: Z₂ screening (no spinorial) ⟺ -1 ∈ [G,G] ⟺ Q₈ ⊂ G.
This is a 'discrete screening' of the Z₂ center symmetry.

======================================================================
PART B: SCOPE — COUNTEREXAMPLE BEYOND SU(2)
======================================================================

D₄ (dihedral group of order 8) is NOT a subgroup of SU(2),
but HAS a center Z₂ = {e, r²} with r² ∈ [D₄, D₄].

  |D₄| = 8
  Center = ['r^0', 'r^2']
  r² = [-1.  0.  0. -1.] (= -I)
  [D₄, D₄] = ['r^0', 'r^2']
  r² ∈ [D₄, D₄]: True

  Element orders in D₄:
    r^0: order 1
    r^2: order 2
    sr^0: order 2
    sr^1: order 2
    sr^2: order 2
    sr^3: order 2

  Elements of order 2 in D₄: 5
  Elements of order 2 in Q₈: 1 (only -1)

  D₄ has 5 elements of order 2, Q₈ has only 1.
  Therefore Q₈ ⊄ D₄ (as abstract groups).

  But r² ∈ [D₄, D₄] (center element in commutator subgroup).

  COUNTEREXAMPLE: -1 ∈ [G,G] but Q₈ ⊄ G.
  The equivalence fails outside SU(2). ✓

======================================================================
PART C: WHY THE EQUIVALENCE IS SPECIFIC TO SU(2)
======================================================================

The classification-free proof (script 10) uses:
  [a,b] = -1 ⟹ ab = -ba ⟹ a₀ = b₀ = 0, ā ⊥ b̄

This requires a,b to be UNIT QUATERNIONS (|a| = |b| = 1).
The key step 'a₀² = -|ā|² is impossible' uses:
  a₀² + |ā|² = 1  (unit norm constraint)

In D₄ ⊂ GL(2,R), the matrices are NOT unit quaternions.
The anticommutation ab = -ba still holds for
  a = r = [[0,-1],[1,0]], b = s = [[1,0],[0,-1]]
  rs = [0. 1. 1. 0.], sr = [ 0. -1. -1.  0.]
  rs = -sr: True

But a and b are NOT pure imaginary unit quaternions:
  r has eigenvalues ±i (okay)
  s has eigenvalues ±1 (real, NOT purely imaginary)
  s is a REFLECTION, not a rotation.

The SU(2) structure excludes reflections. This is why the
Q₈ criterion works for SU(2) subgroups but not in general.

======================================================================
SUMMARY
======================================================================

C4 (Discrete anomaly):
  The obstruction is a Z₂ SCREENING phenomenon.
  When -1 ∈ [G,G]: the Z₂ center acts trivially on all
  scalar quantization sectors. Not a 't Hooft anomaly in the
  QFT sense, but a structural property of the central extension.
  Worth a Remark in the paper.

C5 (Scope beyond SU(2)):
  COUNTEREXAMPLE: D₄ has -1 ∈ [G,G] but Q₈ ⊄ D₄.
  The equivalence Q₈ ⊂ G ⟺ -1 ∈ [G,G] is SPECIFIC to
  finite SU(2) subgroups. The quaternion unit norm constraint
  is essential for the classification-free proof.
  Worth a Remark clarifying the scope.
"""

import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))

import numpy as np
from src.quaternion import qmul, qinv, qkey
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                        build_binary_icosahedral, compute_commutator_subgroup,
                        MINUS_ONE_KEY)


# =============================================================
# Part A: Discrete anomaly — central extension analysis
# =============================================================

print("=" * 70)
print("PART A: CENTRAL EXTENSION AND Z₂ SCREENING")
print("=" * 70)
print()
print("For G ⊂ SU(2) finite, the central extension:")
print("  1 → Z₂ → G → G̅ = G/{±1} → 1")
print()
print("The image of -1 in G^{ab} = G/[G,G] determines whether")
print("the Z₂ center is 'visible' to scalar quantization:")
print("  -1 ∉ [G,G]: image nontrivial → Z₂ detectable → spinorial sectors exist")
print("  -1 ∈ [G,G]: image trivial → Z₂ screened → no spinorial sectors")
print()

# Verify: image of -1 in abelianization for all groups
groups_all = [
    ("Z₆ (cyclic)", None, 6),
    ("2D₃ (Dic₃)", None, 12),
    ("2D₅ (Dic₅)", None, 20),
    ("Q₈ (Dic₂)", None, 8),
    ("2T", build_binary_tetrahedral(), 24),
    ("2O", build_binary_octahedral(), 48),
    ("2I", build_binary_icosahedral(), 120),
]

# Build dicyclic groups
def build_dicyclic(N):
    elems = []
    for k in range(2 * N):
        theta = k * np.pi / N
        elems.append(np.array([np.cos(theta), np.sin(theta), 0, 0]))
    for k in range(2 * N):
        theta = k * np.pi / N
        q = qmul(np.array([np.cos(theta), np.sin(theta), 0, 0]),
                 np.array([0, 0, 1, 0]))
        elems.append(q)
    return np.array(elems)

def build_cyclic_double(n):
    elems = []
    for k in range(2 * n):
        theta = k * np.pi / n
        elems.append(np.array([np.cos(theta), np.sin(theta), 0, 0]))
    return np.array(elems)

concrete = [
    ("Z₆", build_cyclic_double(3), False),
    ("Dic₃", build_dicyclic(3), False),
    ("Dic₅", build_dicyclic(5), False),
    ("Q₈", build_dicyclic(2), True),
    ("Dic₄", build_dicyclic(4), True),
    ("2T", build_binary_tetrahedral(), True),
    ("2O", build_binary_octahedral(), True),
    ("2I", build_binary_icosahedral(), True),
]

print(f"  {'Group':>6} {'|G|':>4} {'|[G,G]|':>8} {'-1∈[G,G]':>10} {'Z₂ screened':>12} {'spinorial':>10}")
print("  " + "-" * 56)

for name, elems, expected_obs in concrete:
    order = len(elems)
    comm = compute_commutator_subgroup(elems)
    m1_in = MINUS_ONE_KEY in comm
    screened = "yes" if m1_in else "no"
    spin = "no" if m1_in else "yes"
    ok = (m1_in == expected_obs)
    m1_str = "True" if m1_in else "False"
    print(f"  {name:>6} {order:>4} {len(comm):>8} {m1_str:>10}"
          f" {screened:>12} {spin:>10}  {'✓' if ok else '✗'}")

print()
print("Pattern: Z₂ screening (no spinorial) ⟺ -1 ∈ [G,G] ⟺ Q₈ ⊂ G.")
print("This is a 'discrete screening' of the Z₂ center symmetry.")

# =============================================================
# Part B: Scope beyond SU(2) — counterexample D₄
# =============================================================

print()
print("=" * 70)
print("PART B: SCOPE — COUNTEREXAMPLE BEYOND SU(2)")
print("=" * 70)
print()
print("D₄ (dihedral group of order 8) is NOT a subgroup of SU(2),")
print("but HAS a center Z₂ = {e, r²} with r² ∈ [D₄, D₄].")
print()

# Build D₄ abstractly as permutation group
# D₄ = ⟨r, s | r⁴ = s² = e, srs = r⁻¹⟩
# Elements: e, r, r², r³, s, sr, sr², sr³

# Compute [D₄, D₄] abstractly
# [r, s] = r·s·r⁻¹·s⁻¹ = r·s·r³·s = r·(sr³)·s
# sr = r⁻¹s = r³s, so s·r³ = r·s (from srs=r⁻¹ → sr=r³s → sr³ = r·s)
# Actually: srs⁻¹ = r⁻¹ = r³. So s·r·s = r³.
# [r,s] = r·s·r⁻¹·s⁻¹ = r·s·r³·s.
# s·r³ = r·s (from s·r·s = r³ → s·r = r³·s → s·r³ = (s·r)·r² = r³·s·r² = r³·r²·s⁻¹... hmm
# Let me just enumerate.

# Represent D₄ as 2×2 matrices:
# r = rotation by 90° = [[0,-1],[1,0]]
# s = reflection = [[1,0],[0,-1]]

def mat_mul(A, B):
    return A @ B

def mat_eq(A, B):
    return np.allclose(A, B)

r = np.array([[0, -1], [1, 0]], dtype=float)
s = np.array([[1, 0], [0, -1]], dtype=float)
e = np.eye(2, dtype=float)

D4_elements = []
names_D4 = []
for i in range(4):
    ri = np.linalg.matrix_power(r, i)
    D4_elements.append(ri)
    names_D4.append(f"r^{i}")
for i in range(4):
    ri = np.linalg.matrix_power(r, i)
    D4_elements.append(s @ ri)
    names_D4.append(f"sr^{i}")

print(f"  |D₄| = {len(D4_elements)}")

# Center: elements that commute with everything
center = []
for i, gi in enumerate(D4_elements):
    central = True
    for gj in D4_elements:
        if not mat_eq(gi @ gj, gj @ gi):
            central = False
            break
    if central:
        center.append(names_D4[i])
print(f"  Center = {center}")

# r² = [[-1,0],[0,-1]] = -I, the "center" element
r2 = np.linalg.matrix_power(r, 2)
print(f"  r² = {r2.flatten()} (= -I)")

# Commutator subgroup
comm_set = set()
for gi in D4_elements:
    for gj in D4_elements:
        comm = gi @ gj @ np.linalg.inv(gi) @ np.linalg.inv(gj)
        # Round and find matching element
        for k, gk in enumerate(D4_elements):
            if mat_eq(comm, gk):
                comm_set.add(k)
                break

# Close under multiplication
changed = True
comm_list = list(comm_set)
while changed:
    changed = False
    new_set = set(comm_set)
    for a in comm_list:
        for b in comm_list:
            prod = D4_elements[a] @ D4_elements[b]
            for k, gk in enumerate(D4_elements):
                if mat_eq(prod, gk):
                    if k not in new_set:
                        new_set.add(k)
                        changed = True
                    break
    comm_set = new_set
    comm_list = list(comm_set)

comm_names = [names_D4[i] for i in sorted(comm_set)]
print(f"  [D₄, D₄] = {comm_names}")

# Check r² ∈ [D₄, D₄]
r2_idx = 2  # r^2
r2_in_comm = r2_idx in comm_set
print(f"  r² ∈ [D₄, D₄]: {r2_in_comm}")

# Check Q₈ ⊂ D₄: count elements of each order
print()
print("  Element orders in D₄:")
for i, gi in enumerate(D4_elements):
    for ord in range(1, 9):
        if mat_eq(np.linalg.matrix_power(gi, ord), e):
            if ord <= 2:
                print(f"    {names_D4[i]}: order {ord}")
            break

n_order2 = sum(1 for gi in D4_elements
               if mat_eq(gi @ gi, e) and not mat_eq(gi, e))
print(f"\n  Elements of order 2 in D₄: {n_order2}")
print(f"  Elements of order 2 in Q₈: 1 (only -1)")
print()
print(f"  D₄ has {n_order2} elements of order 2, Q₈ has only 1.")
print(f"  Therefore Q₈ ⊄ D₄ (as abstract groups).")
print()
print(f"  But r² ∈ [D₄, D₄] (center element in commutator subgroup).")
print()
print(f"  COUNTEREXAMPLE: -1 ∈ [G,G] but Q₈ ⊄ G.")
print(f"  The equivalence fails outside SU(2). ✓")

# =============================================================
# Part C: Why it fails — the quaternion constraint
# =============================================================

print()
print("=" * 70)
print("PART C: WHY THE EQUIVALENCE IS SPECIFIC TO SU(2)")
print("=" * 70)
print()
print("The classification-free proof (script 10) uses:")
print("  [a,b] = -1 ⟹ ab = -ba ⟹ a₀ = b₀ = 0, ā ⊥ b̄")
print()
print("This requires a,b to be UNIT QUATERNIONS (|a| = |b| = 1).")
print("The key step 'a₀² = -|ā|² is impossible' uses:")
print("  a₀² + |ā|² = 1  (unit norm constraint)")
print()
print("In D₄ ⊂ GL(2,R), the matrices are NOT unit quaternions.")
print("The anticommutation ab = -ba still holds for")
print("  a = r = [[0,-1],[1,0]], b = s = [[1,0],[0,-1]]")
ab = r @ s
ba = s @ r
print(f"  rs = {ab.flatten()}, sr = {ba.flatten()}")
print(f"  rs = -sr: {mat_eq(ab, -ba)}")
print()
print("But a and b are NOT pure imaginary unit quaternions:")
print("  r has eigenvalues ±i (okay)")
print("  s has eigenvalues ±1 (real, NOT purely imaginary)")
print("  s is a REFLECTION, not a rotation.")
print()
print("The SU(2) structure excludes reflections. This is why the")
print("Q₈ criterion works for SU(2) subgroups but not in general.")

# =============================================================
# Summary
# =============================================================

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("C4 (Discrete anomaly):")
print("  The obstruction is a Z₂ SCREENING phenomenon.")
print("  When -1 ∈ [G,G]: the Z₂ center acts trivially on all")
print("  scalar quantization sectors. Not a 't Hooft anomaly in the")
print("  QFT sense, but a structural property of the central extension.")
print("  Worth a Remark in the paper.")
print()
print("C5 (Scope beyond SU(2)):")
print("  COUNTEREXAMPLE: D₄ has -1 ∈ [G,G] but Q₈ ⊄ D₄.")
print("  The equivalence Q₈ ⊂ G ⟺ -1 ∈ [G,G] is SPECIFIC to")
print("  finite SU(2) subgroups. The quaternion unit norm constraint")
print("  is essential for the classification-free proof.")
print("  Worth a Remark clarifying the scope.")
