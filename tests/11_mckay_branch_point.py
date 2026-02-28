#!/usr/bin/env python3
"""
McKay graph and branch-point theorem for consecutive irreducibility.

THEOREM: For a binary polyhedral group 2G with McKay graph (= affine Dynkin
diagram Ê_n), the number of consecutive irreducible V_j|_{2G} equals
  distance(extending node, branch point) + 1
where the extending node is the trivial irrep Γ₀.

  2T (Ê₆): dist = 2, consecutive = 3  ✓
  2O (Ê₇): dist = 3, consecutive = 4  ✓
  2I (Ê₈): dist = 5, consecutive = 6  ✓

MECHANISM: The SU(2) recursion V_{j+1/2} = V_{1/2} ⊗ V_j - V_{j-1/2}
translates to a walk in the McKay graph. At each node with degree 2,
tensoring with V_{1/2} gives 2 neighbors; subtracting V_{j-1/2} (the
node we came from) leaves 1 new irrep → still irreducible.

At the branch point (degree ≥ 3), tensoring gives ≥ 3 neighbors;
subtracting the predecessor leaves ≥ 2 new irreps → reducible.

RAW OUTPUT:
===========

======================================================================
McKay graph: 2T (Ê₆)
======================================================================

  7 irreps, sum dim² = 24 (should be 24): ✓

  McKay graph (V_{1/2} ⊗ Γ_i decomposition):
      A (dim 1, deg 1): → F₁
     F₁ (dim 2, deg 2): → A + T
      T (dim 3, deg 3): → F₁ + F₂ + F₃
     Eω (dim 1, deg 1): → F₂
    Eω² (dim 1, deg 1): → F₃
     F₂ (dim 2, deg 2): → T + Eω
     F₃ (dim 2, deg 2): → T + Eω²

  Adjacency symmetric: True  ✓
  Branch point: T (degree 3)
  Distance A → T: 2
  Consecutive irreducible = dist + 1 = 3
  Expected: 3
  Match: ✓

  Path (V_j walk through McKay graph):
    j=0.0: A (dim 1) *
    j=0.5: F₁ (dim 2) *
    j=1.0: T (dim 3) * ← BRANCH
    j=1.5: F₂ (dim 2)

======================================================================
McKay graph: 2O (Ê₇)
======================================================================

  8 irreps, sum dim² = 48 (should be 48): ✓

  McKay graph (V_{1/2} ⊗ Γ_i decomposition):
     A₀ (dim 1, deg 1): → H₁
     H₁ (dim 2, deg 2): → A₀ + T₂
     T₂ (dim 3, deg 2): → H₁ + L
      L (dim 4, deg 3): → T₂ + T₁ + E
     A₁ (dim 1, deg 1): → H₂
     T₁ (dim 3, deg 2): → L + H₂
     H₂ (dim 2, deg 2): → A₁ + T₁
      E (dim 2, deg 1): → L

  Adjacency symmetric: True  ✓
  Branch point: L (degree 3)
  Distance A₀ → L: 3
  Consecutive irreducible = dist + 1 = 4
  Expected: 4
  Match: ✓

  Path (V_j walk through McKay graph):
    j=0.0: A₀ (dim 1) *
    j=0.5: H₁ (dim 2) *
    j=1.0: T₂ (dim 3) *
    j=1.5: L (dim 4) * ← BRANCH
    j=2.0: T₁ (dim 3)

======================================================================
McKay graph: 2I (Ê₈)
======================================================================

  9 irreps, sum dim² = 120 (should be 120): ✓

  McKay graph (V_{1/2} ⊗ Γ_i decomposition):
      A (dim 1, deg 1): → F₁
     F₁ (dim 2, deg 2): → A + T
      T (dim 3, deg 2): → F₁ + G'
     G' (dim 4, deg 2): → T + H
      H (dim 5, deg 2): → G' + D
      D (dim 6, deg 3): → H + T' + G
     T' (dim 3, deg 1): → D
     F₂ (dim 2, deg 1): → G
      G (dim 4, deg 2): → D + F₂

  Adjacency symmetric: True  ✓
  Branch point: D (degree 3)
  Distance A → D: 5
  Consecutive irreducible = dist + 1 = 6
  Expected: 6
  Match: ✓

  Path (V_j walk through McKay graph):
    j=0.0: A (dim 1) *
    j=0.5: F₁ (dim 2) *
    j=1.0: T (dim 3) *
    j=1.5: G' (dim 4) *
    j=2.0: H (dim 5) *
    j=2.5: D (dim 6) * ← BRANCH
    j=3.0: T' (dim 3)

======================================================================
SUMMARY
======================================================================

  Group  Dynkin  dist(Γ₀,branch)  consec_irred  match
  ------------------------------------------------------
  2T     Ê₆     2                 3             ✓
  2O     Ê₇     3                 4             ✓
  2I     Ê₈     5                 6             ✓

Formula: consecutive_irreducible = dist(extending_node, branch_point) + 1

The V_j sequence traces the unique path from the extending node
(trivial irrep) to the branch point in the affine Dynkin diagram.
At the branch point, the path bifurcates and V_j becomes reducible.
"""

import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
from collections import deque
from src.quaternion import qmul, qinv, qkey
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                        build_binary_icosahedral, compute_conjugacy_classes)


def chi_su2(j, trace):
    """SU(2) spin-j character at given trace."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


# =============================================================
# Character table builders (from script 9, proven correct)
# =============================================================

def build_character_table_2T():
    """Build complete character table for 2T (7 irreps)."""
    elements = build_binary_tetrahedral()
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    nc = len(classes)

    chi = {}
    names = []
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
    names = ['A', 'F₁', 'T', 'Eω', 'Eω²', 'F₂', 'F₃']
    return names, chi, classes, sizes, elements


def build_character_table_2O():
    """Build complete character table for 2O (8 irreps)."""
    elements = build_binary_octahedral()
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    nc = len(classes)

    chi = {}
    chi['A₀'] = np.array([chi_su2(0, c['trace']) for c in classes])
    chi['H₁'] = np.array([chi_su2(0.5, c['trace']) for c in classes])
    chi['T₂'] = np.array([chi_su2(1, c['trace']) for c in classes])

    # A₁: sign character (2O/[2O,2O] = Z₂)
    elements_2T = build_binary_tetrahedral()
    keys_2T = {qkey(e) for e in elements_2T}
    chi_A1 = np.zeros(nc)
    for i, c in enumerate(classes):
        rep_key = list(c['keys'])[0]
        chi_A1[i] = 1.0 if rep_key in keys_2T else -1.0
    chi['A₁'] = chi_A1

    chi['T₁'] = chi_A1 * chi['T₂']
    chi['H₂'] = chi_A1 * chi['H₁']

    chi_V2 = np.array([chi_su2(2, c['trace']) for c in classes])
    chi['E'] = chi_V2 - chi['T₁']

    chi['L'] = np.array([chi_su2(1.5, c['trace']) for c in classes])

    names = ['A₀', 'H₁', 'T₂', 'L', 'A₁', 'T₁', 'H₂', 'E']
    return names, chi, classes, sizes, elements


def build_character_table_2I():
    """Build complete character table for 2I (9 irreps)."""
    elements = build_binary_icosahedral()
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    nc = len(classes)

    phi = (1 + np.sqrt(5)) / 2
    inv_phi = 1 / phi

    chi = {}
    for name, j in [('A', 0), ('F₁', 0.5), ('T', 1), ("G'", 1.5),
                     ('H', 2), ('D', 2.5)]:
        chi[name] = np.array([chi_su2(j, c['trace']) for c in classes])

    sorted_traces = [c['trace'] for c in classes]

    def find_idx(t):
        return int(np.argmin([abs(st - t) for st in sorted_traces]))

    sigma = list(range(nc))
    sigma[find_idx(phi)] = find_idx(-inv_phi)
    sigma[find_idx(-inv_phi)] = find_idx(phi)
    sigma[find_idx(-phi)] = find_idx(inv_phi)
    sigma[find_idx(inv_phi)] = find_idx(-phi)

    chi["T'"] = np.array([chi['T'][sigma[i]] for i in range(nc)])
    chi['F₂'] = np.array([chi['F₁'][sigma[i]] for i in range(nc)])

    chi_V3 = np.array([chi_su2(3, c['trace']) for c in classes])
    chi_remainder = chi_V3.copy()
    for name_check in chi:
        m = sum(sizes[i] * np.conj(chi[name_check][i]) * chi_V3[i]
                for i in range(nc)) / 120
        m_int = int(round(m.real))
        if m_int > 0:
            chi_remainder = chi_remainder - m_int * chi[name_check]
    chi['G'] = chi_remainder

    names = ['A', 'F₁', 'T', "G'", 'H', 'D', "T'", 'F₂', 'G']
    return names, chi, classes, sizes, elements


# =============================================================
# McKay graph computation
# =============================================================

def compute_mckay_adjacency(names, chi, classes, sizes, order):
    """Compute McKay graph: A_{ik} = ⟨V_{1/2}⊗Γ_i, Γ_k⟩."""
    nc = len(classes)
    chi_fund = np.array([chi_su2(0.5, c['trace']) for c in classes])
    n = len(names)
    adj = np.zeros((n, n), dtype=int)

    for i, ni in enumerate(names):
        product = chi_fund * chi[ni]
        for k, nk in enumerate(names):
            m = sum(sizes[ii] * np.conj(chi[nk][ii]) * product[ii]
                    for ii in range(nc)) / order
            adj[i][k] = int(round(m.real))

    return adj


def bfs_to_branch(adj, start=0):
    """BFS from start to first node with degree >= 3."""
    n = adj.shape[0]
    degrees = [sum(adj[i]) for i in range(n)]
    visited = {start}
    queue = deque([(start, 0)])

    while queue:
        node, dist = queue.popleft()
        if degrees[node] >= 3:
            return dist, node, degrees[node]
        for j in range(n):
            if adj[node][j] > 0 and j not in visited:
                visited.add(j)
                queue.append((j, dist + 1))

    return -1, -1, -1


# =============================================================
# Main computation
# =============================================================

builders = [
    ("2T", build_character_table_2T, 24, "Ê₆", 3),
    ("2O", build_character_table_2O, 48, "Ê₇", 4),
    ("2I", build_character_table_2I, 120, "Ê₈", 6),
]

for group_name, builder, order, dynkin, expected in builders:
    print("=" * 70)
    print(f"McKay graph: {group_name} ({dynkin})")
    print("=" * 70)
    print()

    names, chi, classes, sizes, elements = builder()
    nc = len(classes)
    dims = {n: int(round(abs(chi[n][0]))) for n in names}

    # Verify: sum dim² = order
    sum_d2 = sum(d**2 for d in dims.values())
    print(f"  {len(names)} irreps, sum dim² = {sum_d2} (should be {order}): "
          f"{'✓' if sum_d2 == order else '✗'}")
    print()

    adj = compute_mckay_adjacency(names, chi, classes, sizes, order)

    # Print McKay graph
    print(f"  McKay graph (V_{{1/2}} ⊗ Γ_i decomposition):")
    for i, nm in enumerate(names):
        deg = sum(adj[i])
        neighbors = [names[j] for j in range(len(names)) if adj[i][j] > 0]
        print(f"    {nm:>3} (dim {dims[nm]}, deg {deg}): "
              f"→ {' + '.join(neighbors)}")
    print()

    # Symmetry check
    sym = np.allclose(adj, adj.T)
    print(f"  Adjacency symmetric: {sym}  {'✓' if sym else '✗'}")

    # BFS to branch point
    dist, bp, bp_deg = bfs_to_branch(adj)
    bp_name = names[bp] if bp >= 0 else "none"
    consec = dist + 1

    print(f"  Branch point: {bp_name} (degree {bp_deg})")
    print(f"  Distance {names[0]} → {bp_name}: {dist}")
    print(f"  Consecutive irreducible = dist + 1 = {consec}")
    print(f"  Expected: {expected}")
    print(f"  Match: {'✓' if consec == expected else '✗'}")

    # Path trace
    print()
    print(f"  Path (V_j walk through McKay graph):")
    node = 0
    prev = -1
    for step in range(dist + 2):
        j = step / 2.0
        nm_node = names[node]
        deg_node = sum(adj[node])
        tag_bp = " ← BRANCH" if deg_node >= 3 else ""
        tag_ir = " *" if step <= dist else ""
        print(f"    j={j:.1f}: {nm_node} (dim {dims[nm_node]}){tag_ir}{tag_bp}")

        if step <= dist:
            for j2 in range(len(names)):
                if adj[node][j2] > 0 and j2 != prev:
                    prev = node
                    node = j2
                    break

    print()

# =============================================================
# Summary
# =============================================================

print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("  Group  Dynkin  dist(Γ₀,branch)  consec_irred  match")
print("  " + "-" * 54)

for group_name, builder, order, dynkin, expected in builders:
    names, chi, classes, sizes, elements = builder()
    adj = compute_mckay_adjacency(names, chi, classes, sizes, order)
    dist, bp, bp_deg = bfs_to_branch(adj)
    consec = dist + 1
    ok = "✓" if consec == expected else "✗"
    print(f"  {group_name:5}  {dynkin:5}  {dist:<17} {consec:<14}{ok}")

print()
print("Formula: consecutive_irreducible = dist(extending_node, branch_point) + 1")
print()
print("The V_j sequence traces the unique path from the extending node")
print("(trivial irrep) to the branch point in the affine Dynkin diagram.")
print("At the branch point, the path bifurcates and V_j becomes reducible.")
