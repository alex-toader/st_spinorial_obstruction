#!/usr/bin/env python3
"""
H₁ bundle spectrum under symmetry breaking O → D₃.

QUESTION: A spinorial particle on SO(3)/O lives in the H₁ vector bundle
(rank 2). When symmetry reduces O → D₃ (trigonal distortion), what happens
to its spectrum?

FROM SCRIPT 3 (branching):
  H₁|_{2D₃} = ψ₁ (rank-2, remains irreducible)
  H₂|_{2D₃} = ψ₁ (same)
  L|_{2D₃}  = χ₁ ⊕ χ₃ ⊕ ψ₁  (rank-4 splits: two scalar spinorial + rank-2)

So H₁ stays rank-2 under O→D₃. The particle cannot become scalar.
But its SPECTRUM changes: multiplicities m(H₁, j) on 2O vs m(ψ₁, j) on 2D₃.

INVESTIGATION:
  Part 1: Group setup (2O, 2D₃, character tables)
  Part 2: H₁ spectrum on SO(3)/O — multiplicities, gap, density
  Part 3: ψ₁ spectrum on SO(3)/D₃ — same quantities
  Part 4: New scalar spinorial sectors χ₁, χ₃ — what L contributes
  Part 5: Comparison table and Weyl densities
  Part 6: Partition functions and Schottky peaks
  Part 7: Complete spinorial picture (H₁, H₂, L) + branching consistency
  Part 8: Physical interpretation for ST_ Sector 2

Asserts: 296

RAW OUTPUT:
===========

2O: 8 classes, |G| = 48
  Spinorial irreps: H₁(dim 2), H₂(dim 2), L(dim 4)
2D₃: 6 classes, |G| = 12
  Spinorial: χ₁(dim 1), χ₃(dim 1), ψ₁(dim 2)

H₁ SPECTRUM ON SO(3)/O:
    j   m(H₁)     Spectral gap: j=0.5 → j=3.5, ΔE = 15B
   0.5     1       Weyl density: 4/48 = 1/12
   3.5     1
   4.5     1

ψ₁ SPECTRUM ON SO(3)/D₃ (= H₁ after O→D₃):
    j   m(ψ₁)     NO spectral gap: j=0.5, 1.5, 2.5 consecutive
   0.5     1       Weyl density: 4/12 = 1/3
   1.5     1       Densification: 4× vs H₁ on O
   2.5     2

NEW SCALAR SPINORIAL SECTORS χ₁, χ₃ ON SO(3)/D₃:
    j   m(χ₁)     First at j=3/2, E=3.75B
   1.5     1       Topologically forbidden on SO(3)/O
   2.5     1       χ₁, χ₃ conjugate: identical spectra
   3.5     1       Weyl density per sector: 1/12

COMPARISON TABLE (j ≤ 14.5):
      j   H₁(O)   ψ₁(D₃)   χ₁(D₃)   χ₃(D₃)  total(D₃)
    0.5       1        1        0        0          2
    1.5       0        1        1        1          4
    3.5       1        3        1        1          8
  Total: H₁(O)=40 states, D₃=240 states, ratio 6.0×

WEYL DENSITIES:
  H₁ on O:          4/48  = 0.0833
  ψ₁ on D₃:         4/12  = 0.3333
  χ₁+χ₃ on D₃:      2/12  = 0.1667
  Total spin(D₃):    6/12  = 0.5000
  Ratio: 6.0× increase in spinorial state density

PARTITION FUNCTIONS (Z_spin(D₃)/Z_H₁(O)):
  T→∞: ratio → 6.0 (Weyl limit)
  T=1K: ratio = 5.75
  T=0.1K: ratio = 1.08 (ground state dominated)

SCHOTTKY PEAKS:
  H₁ on O:   T = 0.81 K (gap 15B)
  ψ₁ on D₃:  T = 0.44 K (gap 3B)

COMPLETE SPINORIAL PICTURE (2O):
  First appearances: H₁ at j=0.5, L at j=1.5, H₂ at j=2.5
  Branching consistency (100 checks, j ≤ 50):
    m(ψ₁; D₃) = m(H₁; O) + m(H₂; O) + m(L; O)  ✓
    m(χ₁; D₃) = m(L; O)                           ✓

KEY INSIGHT: New scalar spinorial sectors on D₃ come EXCLUSIVELY from L
(rank-4 spinorial irrep of 2O). H₁, H₂ remain rank-2 bundles.
A particle in H₁ stays in a bundle sector at the grain boundary —
its gap closes (15B→3B) but it does not become scalar.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.quaternion import qkey, qmul, qconj
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
                       compute_conjugacy_classes, compute_commutator_subgroup)

# ==========================================================
# Constants
# ==========================================================
B_CM = 0.09107      # SF₆ rotational constant [cm⁻¹]
CM_K = 1.4388       # 1 cm⁻¹ = 1.4388 K
B_K  = B_CM * CM_K  # B in Kelvin

SQ2 = np.sqrt(2)
IDENTITY = np.array([1., 0., 0., 0.])

n_asserts = 0

# ==========================================================
# Infrastructure
# ==========================================================

def chi_su2(j, trace):
    """SU(2) spin-j character at element with given trace."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def generate_group(generators, max_order=300):
    """Generate finite group from generators by closure."""
    elements = [IDENTITY.copy()]
    seen = {qkey(elements[0])}
    queue = list(elements)
    idx = 0
    while idx < len(queue):
        g = queue[idx]; idx += 1
        for gen in generators:
            for h in [qmul(g, gen), qmul(gen, g)]:
                k = qkey(h)
                if k not in seen:
                    seen.add(k)
                    queue.append(h)
                    elements.append(h)
        if len(elements) > max_order:
            raise RuntimeError(f"Group too large: {len(elements)}")
    return elements


def multiplicity_irrep(j, chi_irrep, classes, sizes, G):
    """Multiplicity of irrep (given by class character values) in V_j|_G."""
    chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
    m = np.sum(np.conj(chi_irrep) * chi_Vj * sizes) / G
    assert abs(m.imag) < 1e-8, f"Non-real multiplicity at j={j}: imag={m.imag}"
    return int(round(m.real))


# ==========================================================
# Part 1: Build groups and character tables
# ==========================================================

print("=" * 70)
print("PART 1: GROUP SETUP")
print("=" * 70)

# --- 2O ---
elements_2O = build_binary_octahedral()
classes_2O = compute_conjugacy_classes(elements_2O)
classes_2O.sort(key=lambda c: -c['trace'])
sizes_2O = np.array([c['size'] for c in classes_2O], dtype=float)
G_2O = 48
assert int(sum(sizes_2O)) == G_2O; n_asserts += 1

# 2O character table (from script 3)
# Class signature: (size, trace)
# 2O character table indexed by (size, round(trace, 2))
CHAR_TABLE_2O = {
    'A0': {(1,2.0):1, (6,1.41):1,  (8,1.0):1,  (12,0.0):1,
           (6,0.0):1, (8,-1.0):1, (6,-1.41):1, (1,-2.0):1},
    'A1': {(1,2.0):1, (6,1.41):-1, (8,1.0):1,  (12,0.0):-1,
           (6,0.0):1, (8,-1.0):1, (6,-1.41):-1,(1,-2.0):1},
    'E':  {(1,2.0):2, (6,1.41):0,  (8,1.0):-1, (12,0.0):0,
           (6,0.0):2, (8,-1.0):-1,(6,-1.41):0, (1,-2.0):2},
    'T1': {(1,2.0):3, (6,1.41):-1, (8,1.0):0,  (12,0.0):1,
           (6,0.0):-1,(8,-1.0):0, (6,-1.41):-1,(1,-2.0):3},  # = V₁ restricted
    'T2': {(1,2.0):3, (6,1.41):1,  (8,1.0):0,  (12,0.0):-1,
           (6,0.0):-1,(8,-1.0):0, (6,-1.41):1, (1,-2.0):3},
    'H1': {(1,2.0):2, (6,1.41):SQ2,(8,1.0):1,  (12,0.0):0,
           (6,0.0):0, (8,-1.0):-1,(6,-1.41):-SQ2,(1,-2.0):-2},
    'H2': {(1,2.0):2, (6,1.41):-SQ2,(8,1.0):1, (12,0.0):0,
           (6,0.0):0, (8,-1.0):-1,(6,-1.41):SQ2,(1,-2.0):-2},
    'L':  {(1,2.0):4, (6,1.41):0,  (8,1.0):-1, (12,0.0):0,
           (6,0.0):0, (8,-1.0):1, (6,-1.41):0, (1,-2.0):-4},
}


def get_chi_on_classes(irrep_name, classes):
    """Build character array for irrep on given classes."""
    table = CHAR_TABLE_2O[irrep_name]
    vals = []
    for c in classes:
        s, t = c['size'], round(c['trace'], 2)
        matched = False
        for (ks, kt), v in table.items():
            if ks == s and abs(kt - t) < 0.05:
                vals.append(v)
                matched = True
                break
        assert matched, f"No match for {irrep_name}: size={s}, trace={t}"
    return np.array(vals, dtype=float)


chi_H1_2O = get_chi_on_classes('H1', classes_2O)
chi_H2_2O = get_chi_on_classes('H2', classes_2O)
chi_L_2O  = get_chi_on_classes('L', classes_2O)
chi_A0_2O = get_chi_on_classes('A0', classes_2O)
chi_A1_2O = get_chi_on_classes('A1', classes_2O)

# Verify: H1 is spinorial (trace at -1 = -2)
idx_m1_2O = [i for i, c in enumerate(classes_2O) if abs(c['trace'] + 2) < 0.01][0]
assert chi_H1_2O[idx_m1_2O] == -2; n_asserts += 1  # spinorial
assert chi_H2_2O[idx_m1_2O] == -2; n_asserts += 1  # spinorial
assert chi_L_2O[idx_m1_2O] == -4; n_asserts += 1   # spinorial
assert chi_A0_2O[idx_m1_2O] == 1; n_asserts += 1   # tensorial

# Verify dimensions
assert chi_H1_2O[0] == 2; n_asserts += 1
assert chi_L_2O[0] == 4; n_asserts += 1

print(f"\n  2O: {len(classes_2O)} classes, |G| = {G_2O}")
print(f"  Spinorial irreps: H₁(dim 2), H₂(dim 2), L(dim 4)")
print(f"  H₁(-1) = {chi_H1_2O[idx_m1_2O]:.0f}, L(-1) = {chi_L_2O[idx_m1_2O]:.0f}")

# --- 2D₃ ---
a_gen = np.array([0.5, 0.5, 0.5, 0.5])       # (1+i+j+k)/2, order 6
b_gen = np.array([0.0, 1/SQ2, -1/SQ2, 0.0])  # (i-j)/√2, order 4
elements_2D3 = generate_group([a_gen, b_gen, qconj(a_gen)])
assert len(elements_2D3) == 12; n_asserts += 1
G_2D3 = 12

# Verify 2D₃ ⊂ 2O (embedding check)
keys_2O = {qkey(g) for g in elements_2O}
assert all(qkey(g) in keys_2O for g in elements_2D3), "2D₃ not embedded in 2O"
n_asserts += 1

# Conjugacy classes
def conjugacy_classes_local(elements):
    remaining = list(range(len(elements)))
    classes = []
    while remaining:
        rep_idx = remaining[0]
        rep = elements[rep_idx]
        cls_indices, new_remaining = [], []
        for idx in remaining:
            g = elements[idx]
            is_conj = any(np.allclose(qmul(qmul(h, g), qconj(h)), rep, atol=1e-10)
                          for h in elements)
            (cls_indices if is_conj else new_remaining).append(idx)
        classes.append({
            'size': len(cls_indices),
            'trace': round(2 * elements[cls_indices[0]][0], 6),
            'rep': elements[cls_indices[0]].copy(),
            'keys': {qkey(elements[i]) for i in cls_indices},
        })
        remaining = new_remaining
    return classes

classes_2D3 = conjugacy_classes_local(elements_2D3)
classes_2D3.sort(key=lambda c: -c['trace'])
sizes_2D3 = np.array([c['size'] for c in classes_2D3], dtype=float)
assert len(classes_2D3) == 6; n_asserts += 1

# Identify classes by (size, trace) and distinguish the two size-3
# classes by which one contains b. Both are outside [G,G] = {1, a², a⁴},
# so commutator membership cannot distinguish them.
comm_2D3 = compute_commutator_subgroup(elements_2D3)
comm_2D3_keys = set(comm_2D3.keys())

# Verify: -1 ∉ [2D₃, 2D₃] — non-obstructed (Theorem 2.3)
minus_one_key = qkey(np.array([-1.0, 0, 0, 0]))
assert minus_one_key not in comm_2D3_keys, "-1 should NOT be in [2D₃,2D₃]"
n_asserts += 1

# Standard 2D₃ (= Dic₃) character table
# Classes: {1}(1,+2), {a,a⁵}(2,+1), {a²,a⁴}(2,-1), {-1}(1,-2),
#          {b,a²b,a⁴b}(3,0), {ab,a³b,a⁵b}(3,0)

# Find which size-3 class contains b and which contains ab
b_key = qkey(b_gen)
ab_key = qkey(qmul(a_gen, b_gen))
c_b_idx, c_ab_idx = None, None
for i, c in enumerate(classes_2D3):
    if c['size'] == 3:
        if b_key in c['keys']:
            c_b_idx = i
        if ab_key in c['keys']:
            c_ab_idx = i
assert c_b_idx is not None and c_ab_idx is not None
assert c_b_idx != c_ab_idx

def map_classes_to_standard(classes, c_b, c_ab):
    """Map 2D₃ classes to standard Dic₃ order.
    Standard: 0={1}, 1={a,a⁵}, 2={a²,a⁴}, 3={-1}, 4={b-class}, 5={ab-class}."""
    mapping = {}
    for i, c in enumerate(classes):
        s, t = c['size'], c['trace']
        if s == 1 and t > 0:
            mapping[i] = 0
        elif s == 2 and t > 0:
            mapping[i] = 1
        elif s == 2 and t < 0:
            mapping[i] = 2
        elif s == 1 and t < 0:
            mapping[i] = 3
        elif i == c_b:
            mapping[i] = 4
        elif i == c_ab:
            mapping[i] = 5
        else:
            raise ValueError(f"Unmatched class: size={s}, trace={t}")
    return mapping

cls_map = map_classes_to_standard(classes_2D3, c_b_idx, c_ab_idx)
assert len(cls_map) == 6; n_asserts += 1

# Standard character table for Dic₃ (6 irreps, 6 classes)
# Order: {1}, {a,a⁵}, {a²,a⁴}, {-1}, {b-type}, {ab-type}
STD_CHARS = {
    'chi_0': np.array([1, 1, 1, 1, 1, 1], dtype=complex),         # trivial, T
    'chi_1': np.array([1, -1, 1, -1, 1j, -1j], dtype=complex),    # spinorial
    'chi_2': np.array([1, 1, 1, 1, -1, -1], dtype=complex),       # sign, T
    'chi_3': np.array([1, -1, 1, -1, -1j, 1j], dtype=complex),    # spinorial
    'psi_1': np.array([2, 1, -1, -2, 0, 0], dtype=complex),       # defining, S
    'psi_2': np.array([2, -1, -1, 2, 0, 0], dtype=complex),       # T, rank-2
}

# Convert to our class ordering
CHAR_2D3 = {}
for name, std_vals in STD_CHARS.items():
    our_vals = np.array([std_vals[cls_map[i]] for i in range(6)], dtype=complex)
    CHAR_2D3[name] = our_vals

# Verify orthogonality
for n1 in CHAR_2D3:
    for n2 in CHAR_2D3:
        ip = np.sum(sizes_2D3 * np.conj(CHAR_2D3[n1]) * CHAR_2D3[n2]) / G_2D3
        expected = 1.0 if n1 == n2 else 0.0
        assert abs(ip - expected) < 1e-10, f"Orthogonality fail: <{n1},{n2}> = {ip}"
n_asserts += 1  # bulk orthogonality

# Spinorial status
idx_m1_2D3 = [i for i in range(6) if cls_map[i] == 3][0]
assert CHAR_2D3['chi_1'][idx_m1_2D3].real == -1; n_asserts += 1  # spinorial
assert CHAR_2D3['chi_3'][idx_m1_2D3].real == -1; n_asserts += 1  # spinorial
assert CHAR_2D3['psi_1'][idx_m1_2D3].real == -2; n_asserts += 1  # spinorial rank-2
assert CHAR_2D3['chi_0'][idx_m1_2D3].real == 1; n_asserts += 1   # tensorial

print(f"  2D₃: {len(classes_2D3)} classes, |G| = {G_2D3}")
print(f"  Spinorial: χ₁(dim 1), χ₃(dim 1), ψ₁(dim 2)")
print(f"  Tensorial: χ₀(dim 1), χ₂(dim 1), ψ₂(dim 2)")

# ==========================================================
# Part 2: H₁ spectrum on SO(3)/O
# ==========================================================

print("\n" + "=" * 70)
print("PART 2: H₁ SPECTRUM ON SO(3)/O")
print("=" * 70)

J_MAX = 50  # half-integer j up to this

# H₁ multiplicities: only half-integer j (spinorial)
m_H1 = {}
for two_j in range(1, 2 * J_MAX + 1, 2):
    j = two_j / 2.0
    m = multiplicity_irrep(j, chi_H1_2O, classes_2O, sizes_2O, G_2O)
    if m > 0:
        m_H1[j] = m

print(f"\n  m(H₁, j) for j ≤ 15 (half-integer only):")
print(f"  {'j':>5}  {'E/B':>7}  {'m(H₁)':>6}  {'states':>6}")
print("  " + "-" * 32)
for j in sorted(m_H1.keys()):
    if j > 15:
        break
    E_B = j * (j + 1)
    states = m_H1[j] * 2  # rank-2 bundle → 2 states per multiplicity
    print(f"  {j:5.1f}  {E_B:7.2f}  {m_H1[j]:6d}  {states:6d}")

# Gap analysis
j_vals_H1 = sorted(m_H1.keys())
j_first_H1 = j_vals_H1[0]
j_second_H1 = j_vals_H1[1]
gap_H1 = j_second_H1 * (j_second_H1 + 1) - j_first_H1 * (j_first_H1 + 1)
assert j_first_H1 == 0.5; n_asserts += 1
assert j_second_H1 == 3.5; n_asserts += 1
print(f"\n  First level: j = {j_first_H1} (E = {j_first_H1*(j_first_H1+1):.2f}B)")
print(f"  Second level: j = {j_second_H1} (E = {j_second_H1*(j_second_H1+1):.2f}B)")
print(f"  Spectral gap: ΔE = {gap_H1:.2f}B = {gap_H1*B_K:.3f} K")

# Weyl density for H₁ on 2O
# Asymptotic: m(H₁, j) ~ dim(H₁)² × (2j+1) / |G| = 4(2j+1)/48 = (2j+1)/12
# But only half-integer j, so effective density = 2/12 = 1/6 per unit j
weyl_H1 = 4.0 / G_2O  # dim²/|G|
print(f"\n  Weyl density: dim(H₁)²/|2O| = 4/48 = {weyl_H1:.4f}")

# ==========================================================
# Part 3: ψ₁ spectrum on SO(3)/D₃ (H₁ descendant)
# ==========================================================

print("\n" + "=" * 70)
print("PART 3: ψ₁ SPECTRUM ON SO(3)/D₃ (= H₁ after O→D₃)")
print("=" * 70)

chi_psi1 = CHAR_2D3['psi_1']

m_psi1 = {}
for two_j in range(1, 2 * J_MAX + 1, 2):
    j = two_j / 2.0
    m = multiplicity_irrep(j, chi_psi1.real, classes_2D3, sizes_2D3, G_2D3)
    if m > 0:
        m_psi1[j] = m

print(f"\n  m(ψ₁, j) for j ≤ 15 (half-integer only):")
print(f"  {'j':>5}  {'E/B':>7}  {'m(ψ₁)':>6}  {'states':>6}")
print("  " + "-" * 32)
for j in sorted(m_psi1.keys()):
    if j > 15:
        break
    E_B = j * (j + 1)
    states = m_psi1[j] * 2
    print(f"  {j:5.1f}  {E_B:7.2f}  {m_psi1[j]:6d}  {states:6d}")

# Gap analysis
j_vals_psi1 = sorted(m_psi1.keys())
assert j_vals_psi1[0] == 0.5; n_asserts += 1
assert j_vals_psi1[1] == 1.5; n_asserts += 1  # NO gap! consecutive
print(f"\n  First level: j = 0.5, Second: j = 1.5 — NO spectral gap")

weyl_psi1 = 4.0 / G_2D3
print(f"  Weyl density: dim(ψ₁)²/|2D₃| = 4/12 = {weyl_psi1:.4f}")
print(f"  Densification: {weyl_psi1/weyl_H1:.1f}× vs H₁ on O")

# ==========================================================
# Part 4: Scalar spinorial sectors χ₁, χ₃ on SO(3)/D₃
# ==========================================================

print("\n" + "=" * 70)
print("PART 4: NEW SCALAR SPINORIAL SECTORS χ₁, χ₃ ON SO(3)/D₃")
print("=" * 70)

chi_chi1 = CHAR_2D3['chi_1']
chi_chi3 = CHAR_2D3['chi_3']

m_chi1 = {}
m_chi3 = {}
for two_j in range(1, 2 * J_MAX + 1, 2):
    j = two_j / 2.0
    m1 = multiplicity_irrep(j, chi_chi1, classes_2D3, sizes_2D3, G_2D3)
    m3 = multiplicity_irrep(j, chi_chi3, classes_2D3, sizes_2D3, G_2D3)
    if m1 > 0:
        m_chi1[j] = m1
    if m3 > 0:
        m_chi3[j] = m3

# χ₁ and χ₃ are complex conjugate, so their spectra must be identical
for j in set(list(m_chi1.keys()) + list(m_chi3.keys())):
    assert m_chi1.get(j, 0) == m_chi3.get(j, 0), f"χ₁ ≠ χ₃ at j={j}"
n_asserts += 1

print(f"\n  χ₁ and χ₃ are conjugate: identical spectra ✓")
print(f"\n  m(χ₁, j) for j ≤ 15:")
print(f"  {'j':>5}  {'E/B':>7}  {'m(χ₁)':>6}")
print("  " + "-" * 22)
for j in sorted(m_chi1.keys()):
    if j > 15:
        break
    print(f"  {j:5.1f}  {j*(j+1):7.2f}  {m_chi1[j]:6d}")

j_first_chi1 = min(m_chi1.keys())
assert j_first_chi1 == 1.5; n_asserts += 1
print(f"\n  First scalar spinorial level: j = 3/2, E = 3.75B")
print(f"  These sectors DO NOT EXIST on SO(3)/O — topologically forbidden")

weyl_chi1 = 1.0 / G_2D3
print(f"  Weyl density per χ: dim(χ)²/|2D₃| = 1/12 = {weyl_chi1:.4f}")
print(f"  Two conjugate sectors: total scalar spinorial density = {2*weyl_chi1:.4f}")

# ==========================================================
# Part 5: Comparison table
# ==========================================================

print("\n" + "=" * 70)
print("PART 5: SPECTRUM COMPARISON TABLE")
print("=" * 70)

print(f"\n  All spinorial sectors, j ≤ 14.5:")
print(f"  {'j':>5}  {'H₁(O)':>6}  {'ψ₁(D₃)':>7}  {'χ₁(D₃)':>7}  {'χ₃(D₃)':>7}  {'total(D₃)':>9}")
print("  " + "-" * 50)

total_H1_O = 0
total_psi1_D3 = 0
total_chi_D3 = 0
for two_j in range(1, 30, 2):
    j = two_j / 2.0
    mH1 = m_H1.get(j, 0)
    mpsi1 = m_psi1.get(j, 0)
    mc1 = m_chi1.get(j, 0)
    mc3 = m_chi3.get(j, 0)
    # Total states on D₃: ψ₁ contributes dim=2 states, χ₁,χ₃ each dim=1
    states_D3 = mpsi1 * 2 + mc1 + mc3
    states_O = mH1 * 2
    total_H1_O += states_O
    total_psi1_D3 += mpsi1 * 2
    total_chi_D3 += mc1 + mc3
    if mH1 > 0 or mpsi1 > 0 or mc1 > 0:
        print(f"  {j:5.1f}  {mH1:6d}  {mpsi1:7d}  {mc1:7d}  {mc3:7d}  {states_D3:9d}")

total_D3 = total_psi1_D3 + total_chi_D3
print(f"\n  Total states (j ≤ 14.5):")
print(f"    H₁ on O:       {total_H1_O} (bundle only)")
print(f"    ψ₁ on D₃:      {total_psi1_D3} (bundle)")
print(f"    χ₁+χ₃ on D₃:   {total_chi_D3} (NEW scalar)")
print(f"    Total on D₃:   {total_D3}")
print(f"    Ratio D₃/O:    {total_D3/total_H1_O:.2f}×")

# Asymptotic Weyl densities
print(f"\n  Asymptotic Weyl densities (states per unit (2j+1)):")
print(f"    H₁ on O:          dim²/|G| = 4/48  = {4/48:.4f}")
print(f"    ψ₁ on D₃:         dim²/|G| = 4/12  = {4/12:.4f}")
print(f"    χ₁ on D₃:         dim²/|G| = 1/12  = {1/12:.4f}")
print(f"    χ₃ on D₃:         dim²/|G| = 1/12  = {1/12:.4f}")
print(f"    Total spin(D₃):   4/12+2/12         = {6/12:.4f}")
print(f"    Ratio:            (6/12)/(4/48)     = {(6/12)/(4/48):.1f}×")

ratio_weyl = (6.0/12) / (4.0/48)
assert abs(ratio_weyl - 6.0) < 0.01; n_asserts += 1
print(f"\n  Factor 6 increase in spinorial state density at O→D₃")

# ==========================================================
# Part 6: Partition functions and Schottky peaks
# ==========================================================

print("\n" + "=" * 70)
print("PART 6: THERMODYNAMICS")
print("=" * 70)

def Z_from_dict(m_dict, dim_irrep, T):
    """Partition function Z = Σ m(j) × dim × (2j+1) × exp(-Bj(j+1)/T)."""
    Z = 0.0
    for j, m in m_dict.items():
        Z += m * dim_irrep * (2*j + 1) * np.exp(-B_K * j * (j + 1) / T)
    return Z

def C_from_dict(m_dict, dim_irrep, T):
    """Specific heat C/k_B from multiplicity dictionary."""
    Z = 0.0; E1 = 0.0; E2 = 0.0
    for j, m in m_dict.items():
        Ej = B_K * j * (j + 1)
        w = m * dim_irrep * (2*j + 1) * np.exp(-Ej / T)
        Z += w; E1 += w * Ej; E2 += w * Ej**2
    if Z < 1e-30:
        return 0.0
    return (E2/Z - (E1/Z)**2) / T**2

print(f"\n  B = {B_K:.4f} K (SF₆ energy scale)")

# 6a. Partition functions comparison
print(f"\n  {'T [K]':>7}  {'Z_H₁(O)':>9}  {'Z_ψ₁(D₃)':>10}  {'Z_χ₁(D₃)':>10}  "
      f"{'Z_spin(D₃)':>11}  {'ratio':>6}")
print("  " + "-" * 60)

for T in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]:
    Z_H1 = Z_from_dict(m_H1, 2, T)       # H₁ dim=2
    Z_psi1 = Z_from_dict(m_psi1, 2, T)    # ψ₁ dim=2
    Z_chi1 = Z_from_dict(m_chi1, 1, T)    # χ₁ dim=1
    Z_spin_D3 = Z_psi1 + 2 * Z_chi1       # total spinorial on D₃
    ratio = Z_spin_D3 / Z_H1 if Z_H1 > 1e-30 else float('inf')
    print(f"  {T:7.2f}  {Z_H1:9.3f}  {Z_psi1:10.3f}  {Z_chi1:10.3f}  "
          f"{Z_spin_D3:11.3f}  {ratio:6.2f}")

# High-T ratio should → |2O|/|2D₃| × (relevant dim ratio)
# Z_spin(D₃)/Z_H₁(O) → (4+1+1)/12 ÷ (4/48) = (6/12)/(4/48) = 0.5/0.0833 = 6
# But Z counts states weighted by degeneracy, so need to be careful
# Actually Z ~ T × density, ratio → density ratio = 6

# 6b. Specific heat per sector — find Schottky peaks
print(f"\n  Schottky anomaly (specific heat peaks):")

T_grid = np.linspace(0.02, 10.0, 2000)
C_H1_grid = [C_from_dict(m_H1, 2, T) for T in T_grid]
C_psi1_grid = [C_from_dict(m_psi1, 2, T) for T in T_grid]
C_chi1_grid = [C_from_dict(m_chi1, 1, T) for T in T_grid]

idx_H1 = np.argmax(C_H1_grid)
idx_psi1 = np.argmax(C_psi1_grid)
idx_chi1 = np.argmax(C_chi1_grid)

T_peak_H1 = T_grid[idx_H1]
T_peak_psi1 = T_grid[idx_psi1]
T_peak_chi1 = T_grid[idx_chi1]

print(f"    H₁ on O:  T_peak = {T_peak_H1:.3f} K, C_max/k_B = {C_H1_grid[idx_H1]:.4f}")
print(f"    ψ₁ on D₃: T_peak = {T_peak_psi1:.3f} K, C_max/k_B = {C_psi1_grid[idx_psi1]:.4f}")
print(f"    χ₁ on D₃: T_peak = {T_peak_chi1:.3f} K, C_max/k_B = {C_chi1_grid[idx_chi1]:.4f}")

# Schottky peak position scales with spectral gap
# H₁ gap = 15B = 1.96 K → peak ~ 0.4 × gap ~ 0.8 K
# ψ₁ gap = 2.25B (j=0.5→1.5) → peak ~ 0.3 K
# χ₁ has no j=0.5, starts at j=1.5
print(f"\n    H₁ gap (j=0.5→3.5): {15*B_K:.3f} K → Schottky at {T_peak_H1:.3f} K")
print(f"    ψ₁ gap (j=0.5→1.5): {(3.75-0.75)*B_K:.3f} K → Schottky at {T_peak_psi1:.3f} K")

# ==========================================================
# Part 7: H₂ and L spectra (complete picture)
# ==========================================================

print("\n" + "=" * 70)
print("PART 7: COMPLETE SPINORIAL PICTURE")
print("=" * 70)

# H₂ on 2O
m_H2 = {}
for two_j in range(1, 2 * J_MAX + 1, 2):
    j = two_j / 2.0
    m = multiplicity_irrep(j, chi_H2_2O, classes_2O, sizes_2O, G_2O)
    if m > 0:
        m_H2[j] = m

# L on 2O
m_L = {}
for two_j in range(1, 2 * J_MAX + 1, 2):
    j = two_j / 2.0
    m = multiplicity_irrep(j, chi_L_2O, classes_2O, sizes_2O, G_2O)
    if m > 0:
        m_L[j] = m

print(f"\n  All spinorial irreps of 2O, multiplicities for j ≤ 14.5:")
print(f"  {'j':>5}  {'m(H₁)':>6}  {'m(H₂)':>6}  {'m(L)':>5}  {'total states':>12}")
print("  " + "-" * 40)
for two_j in range(1, 30, 2):
    j = two_j / 2.0
    mH1 = m_H1.get(j, 0)
    mH2 = m_H2.get(j, 0)
    mL = m_L.get(j, 0)
    # States: H₁ × 2 + H₂ × 2 + L × 4
    total = mH1 * 2 + mH2 * 2 + mL * 4
    if total > 0:
        print(f"  {j:5.1f}  {mH1:6d}  {mH2:6d}  {mL:5d}  {total:12d}")

# Verify: at j=0.5, only H₁ contributes (V_{1/2} = H₁ as 2O rep)
assert m_H1.get(0.5, 0) == 1; n_asserts += 1
assert m_H2.get(0.5, 0) == 0; n_asserts += 1
assert m_L.get(0.5, 0) == 0; n_asserts += 1

# H₂ first appearance
j_first_H2 = min(m_H2.keys())
j_first_L = min(m_L.keys())
print(f"\n  First appearances:")
print(f"    H₁: j = {min(m_H1.keys())} (E = {min(m_H1.keys())*(min(m_H1.keys())+1):.2f}B)")
print(f"    H₂: j = {j_first_H2} (E = {j_first_H2*(j_first_H2+1):.2f}B)")
print(f"    L:  j = {j_first_L} (E = {j_first_L*(j_first_L+1):.2f}B)")

# Branching verification: H₁|_{2D₃} = ψ₁, H₂|_{2D₃} = ψ₁, L|_{2D₃} = χ₁+χ₃+ψ₁
# So: m(ψ₁, j; D₃) should equal m(H₁,j;O) + m(H₂,j;O) + m(L,j;O)
#     m(χ₁, j; D₃) should equal m(L,j;O)
print(f"\n  Branching consistency check:")
print(f"  H₁→ψ₁, H₂→ψ₁, L→χ₁+χ₃+ψ₁")
print(f"  So m(ψ₁; D₃) = m(H₁;O) + m(H₂;O) + m(L;O)")
print(f"     m(χ₁; D₃) = m(L;O)")

n_checks = 0
for two_j in range(1, 2 * J_MAX + 1, 2):
    j = two_j / 2.0
    mH1 = m_H1.get(j, 0)
    mH2 = m_H2.get(j, 0)
    mL = m_L.get(j, 0)
    mpsi1 = m_psi1.get(j, 0)
    mc1 = m_chi1.get(j, 0)

    # ψ₁ branching
    expected_psi1 = mH1 + mH2 + mL
    assert mpsi1 == expected_psi1, \
        f"j={j}: m(ψ₁)={mpsi1} ≠ m(H₁)+m(H₂)+m(L)={expected_psi1}"
    n_checks += 1

    # χ₁ branching
    assert mc1 == mL, f"j={j}: m(χ₁)={mc1} ≠ m(L)={mL}"
    n_checks += 1

n_asserts += 1  # bulk branching check
print(f"  Verified: {n_checks} checks (2 per j, j ≤ {J_MAX}), ALL pass ✓")

# ==========================================================
# Part 8: Physical interpretation
# ==========================================================

print("\n" + "=" * 70)
print("PART 8: PHYSICAL INTERPRETATION")
print("=" * 70)

gap_O_K = 15 * B_K
gap_D3_K = 3 * B_K
scalar_E_K = 3.75 * B_K
print(f"""
  A spinorial particle on SO(3)/O lives in the H1 bundle (rank 2).
  Its spectrum has a LARGE gap: j=1/2 to j=7/2 (Delta E = 15B).

  At a trigonal grain boundary (O to D3):
  1. H1 -> psi1 (remains rank-2, stays a bundle sector)
  2. The spectral gap CLOSES: j=1/2, 3/2, 5/2 all present
  3. Multiplicities INCREASE: factor 4 asymptotically
  4. NEW scalar spinorial sectors (chi1, chi3) appear at j=3/2
     These come from L|_2D3 = chi1 + chi3 + psi1
     L does not contribute to scalar sectors on O (rank-4),
     but DOES on D3 (chi1, chi3 are rank-1)

  The particle sees the grain boundary through:
  - Closing of its spectral gap (low-energy states become available)
  - Appearance of NEW scalar sectors (bundle to scalar transition)
  - 6x overall increase in spinorial state density

  Key numbers:
    H1 gap on O:    15B = {gap_O_K:.3f} K
    psi1 gap on D3:  3B = {gap_D3_K:.3f} K  (5x smaller)
    New scalar at: 3.75B = {scalar_E_K:.3f} K
    Density ratio:   6x
""")

# ==========================================================
# Part 9: Vector bundle duality at half-integer j
# V_j + V_{j*-j}|_G = rho_reg^- for half-integer j
# => m(rho,j) + m(rho,j*-j) = dim(rho) for spinorial rho
# ==========================================================

print("\n" + "=" * 70)
print("PART 9: VECTOR BUNDLE DUALITY (half-integer j)")
print("=" * 70)

j_star = G_2O // 4 - 1  # = 48/4 - 1 = 11
print(f"\n  j* = |2O|/4 - 1 = {j_star}")

# Build character arrays for tensorial irreps
chi_E_2O  = get_chi_on_classes('E', classes_2O)
chi_T1_2O = get_chi_on_classes('T1', classes_2O)
chi_T2_2O = get_chi_on_classes('T2', classes_2O)

# Part 9a: Spinorial irreps — m(rho,j)+m(rho,j*-j) = dim(rho)
print(f"\n  9a. Spinorial irreps: m(rho,j)+m(rho,j*-j) = dim(rho)")
print(f"  {'j':>5} {'j*-j':>5}  {'H1+H1':>6}  {'H2+H2':>6}  {'L+L':>5}")
print("  " + "-" * 40)

n_duality = 0
# Only iterate j <= j*/2 to avoid checking each pair twice
for two_j in range(1, j_star + 1, 2):  # half-integer j from 0.5 to 5.5
    j = two_j / 2.0
    jp = j_star - j  # j* - j, also half-integer

    s_H1 = m_H1.get(j, 0) + m_H1.get(jp, 0)
    s_H2 = m_H2.get(j, 0) + m_H2.get(jp, 0)
    s_L = m_L.get(j, 0) + m_L.get(jp, 0)

    print(f"  {j:5.1f} {jp:5.1f}  {s_H1:6d}  {s_H2:6d}  {s_L:5d}")

    assert s_H1 == 2, f"H1 duality fails at j={j}: {s_H1} != 2"
    assert s_H2 == 2, f"H2 duality fails at j={j}: {s_H2} != 2"
    assert s_L == 4, f"L duality fails at j={j}: {s_L} != 4"
    n_duality += 3

n_asserts += n_duality
print(f"  {n_duality} checks pass")

# Part 9b: Tensorial irreps — must have ZERO multiplicity at half-integer j
# This is required for V_j+V_{j*-j}|_G = rho_reg^- (fermionic half)
print(f"\n  9b. Tensorial irreps: m(rho,j) = 0 at half-integer j")
n_tens = 0
for two_j in range(1, 2 * J_MAX + 1, 2):
    j = two_j / 2.0
    m_A0 = multiplicity_irrep(j, chi_A0_2O, classes_2O, sizes_2O, G_2O)
    m_A1 = multiplicity_irrep(j, chi_A1_2O, classes_2O, sizes_2O, G_2O)
    m_E  = multiplicity_irrep(j, chi_E_2O,  classes_2O, sizes_2O, G_2O)
    m_T1 = multiplicity_irrep(j, chi_T1_2O, classes_2O, sizes_2O, G_2O)
    m_T2 = multiplicity_irrep(j, chi_T2_2O, classes_2O, sizes_2O, G_2O)
    assert m_A0 == 0, f"A0 nonzero at half-integer j={j}"
    assert m_A1 == 0, f"A1 nonzero at half-integer j={j}"
    assert m_E  == 0, f"E nonzero at half-integer j={j}"
    assert m_T1 == 0, f"T1 nonzero at half-integer j={j}"
    assert m_T2 == 0, f"T2 nonzero at half-integer j={j}"
    n_tens += 5

n_asserts += n_tens
print(f"  {n_tens} checks pass (5 irreps x {2*J_MAX} half-integer j values)")
print(f"\n  CONCLUSION: V_j + V_{{j*-j}}|_G = rho_reg^- (fermionic half)")
print(f"  at half-integer j. All {n_duality + n_tens} checks pass.")

# ==========================================================
# Summary
# ==========================================================
print("=" * 70)
print(f"TOTAL ASSERTS: {n_asserts}")
print("=" * 70)
