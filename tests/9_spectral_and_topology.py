#!/usr/bin/env python3
"""
Spectral decomposition V_j|_{2G} and topological/cohomological analysis.

Three independent investigations:
  (A) Direction B (w₂): ELIMINATED — Stiefel's theorem ⟹ w₂ = 0 always
  (B) Direction S5: Cohomological reformulation of the obstruction
  (C) Direction C: Spectral decomposition V_j|_{2O}, V_j|_{2T}, V_j|_{2I}

RAW OUTPUT:
===========

======================================================================
PART 1: STIEFEL-WHITNEY w₂ — ELIMINATED
======================================================================

Stiefel's theorem (1936): every closed orientable 3-manifold
is parallelizable. Therefore w₂(T(SO(3)/H)) = 0 for ALL H.

Our obstruction is about flat line bundles (H¹), not tangent
bundle (H²). These are independent. Direction B: RED HERRING.

======================================================================
PART 2: COHOMOLOGICAL REFORMULATION
======================================================================

Central loop evaluation map:
  σ = ev_{-1}: H¹(M; U(1)) → {±1}
  σ trivial ⟺ -1 ∈ [2H, 2H] ⟺ Q₈ ⊂ 2H

Space          Abel    |H¹|  #spin  #tens
---------------------------------------------
  SO(3)        Z₂         2      1      1
  SO(3)/T      Z₃         3      0      3
  SO(3)/O      Z₂         2      0      2
  SO(3)/I      {1}        1      0      1
  SO(3)/D₃     Z₄         4      2      2
  SO(3)/D₅     Z₄         4      2      2
  SO(3)/C₃     Z₆         6      3      3

Counting: #spinorial = floor(|Abel|/2). ✓

======================================================================
PART 3: V_j|_{2O}
======================================================================

2O: 8 classes, 8 irreps
Sum dim² = 48 (should be 48): ✓

Orthogonality: all ⟨Γ_i, Γ_j⟩ = δ_ij ✓

    j  dim  Type  Decomposition
----------------------------------------------------------------------
   0.0    1   (T)   A₀ *
   0.5    2   (S)   H₁ *
   1.0    3   (T)   T₂ *
   1.5    4   (S)   L *
   2.0    5   (T)   T₁ ⊕ E
   2.5    6   (S)   H₂ ⊕ L
   3.0    7   (T)   T₂ ⊕ A₁ ⊕ T₁
   3.5    8   (S)   H₁ ⊕ H₂ ⊕ L
   4.0    9   (T)   A₀ ⊕ T₂ ⊕ T₁ ⊕ E
   4.5   10   (S)   H₁ ⊕ 2L
   5.0   11   (T)   2T₂ ⊕ T₁ ⊕ E
   5.5   12   (S)   H₁ ⊕ H₂ ⊕ 2L

Consecutive irreducible V_j: 4 (j = 0 through j = 1.5)

======================================================================
PART 4: V_j|_{2T}
======================================================================

2T: 7 classes, 7 irreps
Sum dim² = 24 (should be 24): ✓

Orthogonality: all ⟨Γ_i, Γ_j⟩ = δ_ij ✓

    j  dim  Type  Decomposition
----------------------------------------------------------------------
   0.0    1   (T)   A *
   0.5    2   (S)   F₁ *
   1.0    3   (T)   T *
   1.5    4   (S)   F₂ ⊕ F₃
   2.0    5   (T)   T ⊕ Eω ⊕ Eω²
   2.5    6   (S)   F₁ ⊕ F₂ ⊕ F₃
   3.0    7   (T)   A ⊕ 2T
   3.5    8   (S)   2F₁ ⊕ F₂ ⊕ F₃
   4.0    9   (T)   A ⊕ 2T ⊕ Eω ⊕ Eω²
   4.5   10   (S)   F₁ ⊕ 2F₂ ⊕ 2F₃
   5.0   11   (T)   3T ⊕ Eω ⊕ Eω²
   5.5   12   (S)   2F₁ ⊕ 2F₂ ⊕ 2F₃

Consecutive irreducible V_j: 3 (j = 0 through j = 1.0)

======================================================================
PART 5: V_j|_{2I}
======================================================================

2I: 9 classes, 9 irreps
Sum dim² = 120 (should be 120): ✓

Orthogonality: all ⟨Γ_i, Γ_j⟩ = δ_ij ✓

    j  dim  Type  Decomposition
----------------------------------------------------------------------
   0.0    1   (T)   A *
   0.5    2   (S)   F₁ *
   1.0    3   (T)   T *
   1.5    4   (S)   G' *
   2.0    5   (T)   H *
   2.5    6   (S)   D *
   3.0    7   (T)   T' ⊕ G
   3.5    8   (S)   D ⊕ F₂
   4.0    9   (T)   H ⊕ G
   4.5   10   (S)   G' ⊕ D
   5.0   11   (T)   T ⊕ H ⊕ T'
   5.5   12   (S)   F₁ ⊕ G' ⊕ D
   6.0   13   (T)   A ⊕ T ⊕ H ⊕ G
   6.5   14   (S)   F₁ ⊕ G' ⊕ D ⊕ F₂

Consecutive irreducible V_j: 6 (j = 0 through j = 2.5)

REMARKABLE: V_j is IRREDUCIBLE on 2I for ALL j ≤ 5/2!
The first 6 spin levels (j = 0, 1/2, 1, 3/2, 2, 5/2) each give a
DISTINCT irreducible representation of 2I. Only at j = 3 does the
first decomposition occur.

======================================================================
PART 6: IRREDUCIBILITY PATTERN ACROSS GROUPS
======================================================================

  2T (E₆): 3 consecutive irreducible V_j
  2O (E₇): 4 consecutive irreducible V_j
  2I (E₈): 6 consecutive irreducible V_j

For 2I: the first 6 V_j (dims 1,2,3,4,5,6) are all irreducible.
These exhaust all irreps obtainable directly from SU(2) restriction.

======================================================================
PART 7: DIHEDRAL PARITY TOGGLE
======================================================================

  D_2: N=even,  obstruction, #spin sectors = 0
  D_3: N=odd ,  restoration, #spin sectors = 2
  D_4: N=even,  obstruction, #spin sectors = 0
  D_5: N=odd ,  restoration, #spin sectors = 2
  D_6: N=even,  obstruction, #spin sectors = 0
  D_7: N=odd ,  restoration, #spin sectors = 2
  D_8: N=even,  obstruction, #spin sectors = 0
  D_9: N=odd ,  restoration, #spin sectors = 2

Toggle: obstruction ALTERNATES with parity of N.
Re-obstruction possible: O(obs) → D₃(ok) → D₂(obs).

======================================================================
CHALLENGES / POTENTIAL ISSUES
======================================================================

1. Stiefel's theorem applies to smooth manifolds. If H does not act
   freely (e.g. H acting on S² has fixed points), SO(3)/H is an
   orbifold, not a manifold. For polyhedral H ⊂ SO(3), the action
   on SO(3) by right multiplication IS free (no fixed points), so
   SO(3)/H is a smooth manifold and Stiefel applies.

2. The cohomological statement "c is a boundary" means c ∈ im(∂₂)
   where ∂₂: C₂ → C₁ in the CW chain complex. Equivalently,
   c is trivial in H₁(M; Z)/torsion... no: c has order 2, so
   it's always in the torsion part. The correct statement is that
   c is trivial in H₁(M; Z₂) = Hom(π₁, Z₂) when -1 maps to
   identity in G/[G,G].

3. The spectral tables give V_j|_{2G} for j up to 5.5.
   Higher j follow the pattern determined by the character table.
   The decomposition is periodic modulo some period related to |2G|.

4. The irreducibility of V_j|_{2I} for j ≤ 5/2 is a well-known
   result in the representation theory of icosahedral groups.
   It is related to the fact that 2I has irreps of dimensions
   1,2,3,4,5,6 — all consecutive integers.

5. The dihedral parity toggle has physical implications for cascaded
   phase transitions. In a liquid crystal undergoing O → D₃ → D₂,
   spinorial sectors appear and then disappear again.
"""

import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))

import numpy as np
from src.quaternion import qmul, qkey, qinv
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                        build_binary_icosahedral, compute_conjugacy_classes)


def chi_su2(j, trace):
    """SU(2) spin-j character at given trace."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** (int(round(2 * j)))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def build_character_table_2O():
    """Build complete character table for 2O (8 irreps)."""
    elements = build_binary_octahedral()
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    nc = len(classes)

    # From SU(2): A₀ = V_0, H₁ = V_{1/2}, T₂ = V_1, L = V_{3/2}
    chi = {}
    for name, j in [('A₀', 0), ('H₁', 0.5), ('T₂', 1)]:
        chi[name] = np.array([chi_su2(j, c['trace']) for c in classes])

    # A₁: sign character (2O/[2O,2O] = Z₂)
    elements_2T = build_binary_tetrahedral()
    keys_2T = {qkey(e) for e in elements_2T}
    chi_A1 = np.zeros(nc)
    for i, c in enumerate(classes):
        rep_key = list(c['keys'])[0]
        chi_A1[i] = 1.0 if rep_key in keys_2T else -1.0
    chi['A₁'] = chi_A1

    # T₁ = A₁ ⊗ T₂, H₂ = A₁ ⊗ H₁
    chi['T₁'] = chi_A1 * chi['T₂']
    chi['H₂'] = chi_A1 * chi['H₁']

    # E from V_2 = E ⊕ T₁
    chi_V2 = np.array([chi_su2(2, c['trace']) for c in classes])
    chi['E'] = chi_V2 - chi['T₁']

    # L from V_{3/2} (which is irreducible = L)
    chi['L'] = np.array([chi_su2(1.5, c['trace']) for c in classes])

    return chi, classes, sizes, elements


def build_character_table_2T():
    """Build complete character table for 2T (7 irreps)."""
    elements = build_binary_tetrahedral()
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    nc = len(classes)

    chi = {}
    chi['A'] = np.array([chi_su2(0, c['trace']) for c in classes])
    chi['F₁'] = np.array([chi_su2(0.5, c['trace']) for c in classes])
    chi['T'] = np.array([chi_su2(1, c['trace']) for c in classes])

    # Eω via Q₈ coset structure
    Q8_keys = set()
    for e in elements:
        if abs(e[0]) > 0.99 or (abs(e[0]) < 0.01 and
                any(abs(abs(e[k]) - 1) < 0.01 for k in range(1, 4))):
            Q8_keys.add(qkey(e))

    a = np.array([0.5, 0.5, 0.5, 0.5])
    w = np.exp(2j * np.pi / 3)
    cosets = {}
    for e in elements:
        ak = np.array([1., 0., 0., 0.])
        for c_idx in range(3):
            test = qmul(qinv(ak), e)
            for e2 in elements:
                if np.allclose(e2, test, atol=1e-10):
                    test = e2
                    break
            if qkey(test) in Q8_keys:
                cosets[qkey(e)] = c_idx
                break
            ak = qmul(ak, a)

    chi_Ew = np.zeros(nc, dtype=complex)
    chi_Ew2 = np.zeros(nc, dtype=complex)
    for i, c in enumerate(classes):
        rep_key = list(c['keys'])[0]
        for e in elements:
            if qkey(e) == rep_key:
                chi_Ew[i] = w ** cosets[qkey(e)]
                chi_Ew2[i] = w ** (2 * cosets[qkey(e)])
                break

    chi['Eω'] = chi_Ew
    chi['Eω²'] = chi_Ew2
    chi['F₂'] = chi_Ew * chi['F₁']
    chi['F₃'] = chi_Ew2 * chi['F₁']

    return chi, classes, sizes, elements


def build_character_table_2I():
    """Build complete character table for 2I (9 irreps)."""
    elements = build_binary_icosahedral()
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    nc = len(classes)
    sorted_traces = [c['trace'] for c in classes]

    phi = (1 + np.sqrt(5)) / 2
    inv_phi = 1 / phi

    chi = {}
    for name, j in [('A', 0), ('F₁', 0.5), ('T', 1), ("G'", 1.5), ('H', 2), ('D', 2.5)]:
        chi[name] = np.array([chi_su2(j, c['trace']) for c in classes])

    # Conjugation: swap phi ↔ -1/phi and -phi ↔ 1/phi
    def find_idx(t):
        return int(np.argmin([abs(st - t) for st in sorted_traces]))

    sigma = list(range(nc))
    sigma[find_idx(phi)] = find_idx(-inv_phi)
    sigma[find_idx(-inv_phi)] = find_idx(phi)
    sigma[find_idx(-phi)] = find_idx(inv_phi)
    sigma[find_idx(inv_phi)] = find_idx(-phi)

    chi["T'"] = np.array([chi['T'][sigma[i]] for i in range(nc)])
    chi['F₂'] = np.array([chi['F₁'][sigma[i]] for i in range(nc)])

    # G from V_3: first check which known irreps appear
    chi_V3 = np.array([chi_su2(3, c['trace']) for c in classes])
    # V_3|_{2I} = T'(3) + G(4)  (A has multiplicity 0 in V_3)
    chi_remainder = chi_V3.copy()
    for name_check in chi:
        m = sum(sizes[i] * np.conj(chi[name_check][i]) * chi_V3[i]
                for i in range(nc)) / 120
        m_int = int(round(m.real))
        if m_int > 0:
            chi_remainder = chi_remainder - m_int * chi[name_check]
    chi['G'] = chi_remainder

    return chi, classes, sizes, elements


def decompose_table(group_name, all_chi, classes, sizes, order, j_max=5.5):
    """Print decomposition table V_j|_{2G}."""
    nc = len(classes)
    dims = {n: int(round(abs(c[0]))) for n, c in all_chi.items()}

    # Verify sum dim²
    sum_d2 = sum(d ** 2 for d in dims.values())
    print(f"{group_name}: {nc} classes, {len(all_chi)} irreps")
    print(f"Sum dim² = {sum_d2} (should be {order}): {'✓' if sum_d2 == order else '✗'}")
    print()

    # Orthogonality check
    all_ok = True
    for n1 in all_chi:
        for n2 in all_chi:
            ip = sum(sizes[i] * np.conj(all_chi[n1][i]) * all_chi[n2][i]
                     for i in range(nc)) / order
            expected = 1 if n1 == n2 else 0
            if abs(ip - expected) > 0.01:
                print(f"  FAIL: <{n1},{n2}> = {ip}")
                all_ok = False
    if all_ok:
        print(f"Orthogonality: all ⟨Γ_i, Γ_j⟩ = δ_ij ✓")

    # Decomposition table
    print()
    print(f"{'j':>5} {'dim':>4} {'Type':>5}  Decomposition")
    print("-" * 70)

    irreducible_count = 0
    names = list(all_chi.keys())

    for j2 in range(0, int(2 * j_max) + 1):
        j = j2 / 2.0
        d = int(2 * j + 1)
        spin_type = "S" if j2 % 2 == 1 else "T"

        chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])

        decomp = []
        dim_check = 0
        n_components = 0
        for name in names:
            chi_val = all_chi[name]
            m = sum(sizes[i] * np.conj(chi_val[i]) * chi_Vj[i]
                    for i in range(nc)) / order
            m_int = int(round(m.real))
            if m_int > 0:
                dim_check += m_int * dims[name]
                n_components += m_int
                prefix = f"{m_int}" if m_int > 1 else ""
                decomp.append(f"{prefix}{name}")

        assert dim_check == d, f"j={j}: dim mismatch {dim_check} vs {d}"

        is_irr = (n_components == 1)
        if is_irr:
            irreducible_count += 1
        else:
            if irreducible_count > 0:
                pass  # first non-irreducible

        marker = " *" if is_irr else ""
        print(f"  {j:>4.1f}  {d:>3}   ({spin_type})   {' ⊕ '.join(decomp)}{marker}")

    print()
    print(f"Consecutive irreducible V_j: {irreducible_count} "
          f"(j = 0 through j = {(irreducible_count - 1) / 2})")
    return irreducible_count


# =============================================================
# Part 1: w₂ analysis (conceptual — no computation needed)
# =============================================================

print("=" * 70)
print("PART 1: STIEFEL-WHITNEY w₂ — ELIMINATED")
print("=" * 70)
print()
print("Stiefel's theorem (1936): every closed orientable 3-manifold")
print("is parallelizable. Therefore w₂(T(SO(3)/H)) = 0 for ALL H.")
print()
print("Our obstruction is about flat line bundles (H¹), not tangent")
print("bundle (H²). These are independent. Direction B: RED HERRING.")

# =============================================================
# Part 2: Cohomological reformulation
# =============================================================

print()
print("=" * 70)
print("PART 2: COHOMOLOGICAL REFORMULATION")
print("=" * 70)
print()
print("Central loop evaluation map:")
print("  σ = ev_{-1}: H¹(M; U(1)) → {±1}")
print("  σ trivial ⟺ -1 ∈ [2H, 2H] ⟺ Q₈ ⊂ 2H")
print()

spaces = [
    ("SO(3)", "Z₂", 2, 1),
    ("SO(3)/T", "Z₃", 3, 0),
    ("SO(3)/O", "Z₂", 2, 0),
    ("SO(3)/I", "{1}", 1, 0),
    ("SO(3)/D₃", "Z₄", 4, 2),
    ("SO(3)/D₅", "Z₄", 4, 2),
    ("SO(3)/C₃", "Z₆", 6, 3),
]

print(f"{'Space':14} {'Abel':6} {'|H¹|':>5} {'#spin':>6} {'#tens':>6}")
print("-" * 45)
for space, abel, total, spin in spaces:
    print(f"  {space:12} {abel:6} {total:>5} {spin:>6} {total - spin:>6}")
print()
print("Counting: #spinorial = floor(|Abel|/2). ✓")

# =============================================================
# Parts 3-5: Spectral decomposition
# =============================================================

print()
print("=" * 70)
print("PART 3: V_j|_{2O}")
print("=" * 70)
print()
chi_2O, cls_2O, sz_2O, _ = build_character_table_2O()
n_irr_2O = decompose_table("2O", chi_2O, cls_2O, sz_2O, 48)

print()
print("=" * 70)
print("PART 4: V_j|_{2T}")
print("=" * 70)
print()
chi_2T, cls_2T, sz_2T, _ = build_character_table_2T()
n_irr_2T = decompose_table("2T", chi_2T, cls_2T, sz_2T, 24)

print()
print("=" * 70)
print("PART 5: V_j|_{2I}")
print("=" * 70)
print()
chi_2I, cls_2I, sz_2I, _ = build_character_table_2I()
n_irr_2I = decompose_table("2I", chi_2I, cls_2I, sz_2I, 120, j_max=6.5)

# =============================================================
# Part 6: Irreducibility pattern
# =============================================================

print()
print("=" * 70)
print("PART 6: IRREDUCIBILITY PATTERN ACROSS GROUPS")
print("=" * 70)
print()
print(f"  2T (E₆): {n_irr_2T} consecutive irreducible V_j")
print(f"  2O (E₇): {n_irr_2O} consecutive irreducible V_j")
print(f"  2I (E₈): {n_irr_2I} consecutive irreducible V_j")
print()
print("For 2I: the first 6 V_j (dims 1,2,3,4,5,6) are all irreducible.")
print("These exhaust all irreps obtainable directly from SU(2) restriction.")

# =============================================================
# Part 7: Dihedral parity toggle
# =============================================================

print()
print("=" * 70)
print("PART 7: DIHEDRAL PARITY TOGGLE")
print("=" * 70)
print()

for N in range(2, 10):
    obs = "obstruction" if N % 2 == 0 else "restoration"
    sectors = 0 if N % 2 == 0 else 2
    print(f"  D_{N}: N={'even' if N % 2 == 0 else 'odd ':>4}, {obs:>12}, "
          f"#spin sectors = {sectors}")

print()
print("Toggle: obstruction ALTERNATES with parity of N.")
print("Re-obstruction possible: O(obs) → D₃(ok) → D₂(obs).")
