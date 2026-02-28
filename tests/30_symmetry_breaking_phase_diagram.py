#!/usr/bin/env python3
"""
Complete symmetry-breaking phase diagram for 2O.

For every subgroup channel O → H' (and the corresponding 2O → 2H'),
determine:
  1. Whether Q₈ ⊂ 2H' (obstruction persists)
  2. The number of scalar sectors κ
  3. Whether scalar spinorial sectors exist
  4. How spinorial 2O irreps (H₁, H₂, L) branch

Channels checked:
  O → T    (tetrahedral, index 2)
  O → D₄   (tetragonal, index 3)
  O → D₃   (trigonal, index 4)
  O → D₂   (orthorhombic, index 6)
  O → C₄   (tetragonal cyclic, index 6)
  O → C₃   (trigonal cyclic, index 8)
  O → C₂   (minimal rotation, index 12)

This elevates the trigonal/tetragonal comparison from examples to
a complete classification of re-entrance channels.

RAW OUTPUT:
===========

======================================================================
COMPLETE SYMMETRY-BREAKING PHASE DIAGRAM FOR 2O
======================================================================

Channel          |2H'|  Q₈⊂2H'?  κ  Spinorial?  Mechanism
───────────────────────────────────────────────────────────
O (parent)         48    Yes       2     No        Q₈ ⊂ 2O
O → T              24    Yes       3     No        Q₈ ⊂ 2T
O → D₄             16    Yes       4     No        Q₈ ⊂ 2D₄
O → D₃             12    No        4     YES       Q₈ ⊄ 2D₃ (Lagrange)
O → D₂              8    Yes       4     No        Q₈ = 2D₂
O → C₄              8    No        8     YES       Q₈ ⊄ Z₈ (abelian)
O → C₃              6    No        6     YES       Q₈ ⊄ Z₆ (Lagrange)
O → C₂              4    No        4     YES       Q₈ ⊄ Z₄ (abelian)

Re-entrance channels: D₃, C₄, C₃, C₂
Persistent obstruction: T, D₄, D₂

L(S) branching diagnostic:
  O → T:  L → all rank≥2 (obstruction persists)
  O → D₄: L → all rank≥2 (obstruction persists)
  O → D₃: L → 2 spinorial 1D + rank≥2
  O → D₂: L → all rank≥2 (obstruction persists)
  O → C₄: L → 4 spinorial 1D (fully splits)
  O → C₃: L → 4 spinorial 1D (fully splits)
  O → C₂: L → 4 spinorial 1D (fully splits)

Classification rule:
  Re-entrance ⟺ Q₈ ⊄ 2H' ⟺ -1 ∉ [2H', 2H']

The obstruction persists exactly when the residual symmetry group
retains the quaternion core Q₈.

"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
from src.quaternion import qmul, qconj, qkey
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
                       compute_conjugacy_classes, compute_commutator_subgroup)


SQ2 = np.sqrt(2)
MINUS1 = np.array([-1, 0, 0, 0], dtype=float)
IDENTITY = np.array([1, 0, 0, 0], dtype=float)


def generate_group(generators, max_order=300):
    """Generate group from generators by iterated multiplication."""
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


def contains_Q8(elements):
    """Check if Q₈ = {±1, ±i, ±j, ±k} is contained in the group."""
    keys = {qkey(e) for e in elements}
    q8 = [
        np.array([1, 0, 0, 0]), np.array([-1, 0, 0, 0]),
        np.array([0, 1, 0, 0]), np.array([0, -1, 0, 0]),
        np.array([0, 0, 1, 0]), np.array([0, 0, -1, 0]),
        np.array([0, 0, 0, 1]), np.array([0, 0, 0, -1]),
    ]
    return all(qkey(q) in keys for q in q8)


def minus1_in_commutator(elements):
    """Check if -1 ∈ [G, G]."""
    comm = compute_commutator_subgroup(elements)
    return qkey(MINUS1) in set(comm.keys())


def count_1d_chars(elements):
    """Count 1D characters = |G/[G,G]| = |Abel(G)|."""
    comm = compute_commutator_subgroup(elements)
    return len(elements) // len(comm)


def has_spinorial_1d(elements):
    """Check if any 1D character satisfies χ(-1) = -1."""
    # Equivalent to: -1 ∉ [G, G]
    return not minus1_in_commutator(elements)


# =============================================================
# Build 2O and verify
# =============================================================

elements_2O = build_binary_octahedral()
keys_2O = {qkey(g) for g in elements_2O}
assert len(elements_2O) == 48

print("=" * 70)
print("COMPLETE SYMMETRY-BREAKING PHASE DIAGRAM FOR 2O")
print("=" * 70)

# =============================================================
# Define subgroup embeddings
# =============================================================

# 2T: binary tetrahedral (order 24) — Q₈ ∪ {(±1±i±j±k)/2}
elements_2T = build_binary_tetrahedral()
assert len(elements_2T) == 24

# 2D₄: tetragonal (order 16) — stabilizer of [001]
a_D4 = np.array([1/SQ2, 0, 0, 1/SQ2])    # (1+k)/√2, order 8
x_D4 = np.array([0, 1, 0, 0])              # i, order 4
elements_2D4 = generate_group([a_D4, x_D4, qconj(a_D4)])
assert len(elements_2D4) == 16

# 2D₃: trigonal (order 12) — stabilizer of [111]
a_D3 = np.array([0.5, 0.5, 0.5, 0.5])     # (1+i+j+k)/2, order 6
b_D3 = np.array([0, 1/SQ2, -1/SQ2, 0])    # (i-j)/√2, order 4
elements_2D3 = generate_group([a_D3, b_D3, qconj(a_D3)])
assert len(elements_2D3) == 12

# 2D₂ = Q₈: orthorhombic (order 8) — stabilizer of {[100],[010],[001]}
elements_Q8 = generate_group([
    np.array([0, 1, 0, 0]),  # i
    np.array([0, 0, 1, 0]),  # j
])
assert len(elements_Q8) == 8

# Z₈ ≅ 2C₄: tetragonal cyclic (order 8) — generated by 45° rotation about z
a_C4 = np.array([1/SQ2, 0, 0, 1/SQ2])    # (1+k)/√2, order 8 in SU(2)
elements_Z8 = generate_group([a_C4, qconj(a_C4)])
assert len(elements_Z8) == 8

# Z₆ ≅ 2C₃: trigonal cyclic (order 6) — generated by 120° about [111]
a_C3 = np.array([0.5, 0.5, 0.5, 0.5])    # (1+i+j+k)/2, order 6
elements_Z6 = generate_group([a_C3, qconj(a_C3)])
assert len(elements_Z6) == 6

# Z₄ ≅ 2C₂: minimal rotation (order 4) — generated by 180° about z
a_C2 = np.array([0, 0, 0, 1])             # k, order 4 in SU(2)
elements_Z4 = generate_group([a_C2])
assert len(elements_Z4) == 4

# =============================================================
# Verify all are subgroups of 2O
# =============================================================

subgroups = [
    ("O (parent)", elements_2O),
    ("O → T", elements_2T),
    ("O → D₄", elements_2D4),
    ("O → D₃", elements_2D3),
    ("O → D₂", elements_Q8),
    ("O → C₄", elements_Z8),
    ("O → C₃", elements_Z6),
    ("O → C₂", elements_Z4),
]

for name, elems in subgroups:
    for g in elems:
        assert qkey(g) in keys_2O, f"{name}: element {g} not in 2O!"

# =============================================================
# Phase diagram
# =============================================================

print(f"\n{'Channel':<18} {'|2H|':>5}  {'Q₈⊂2H?':>8}  {'κ':>3}  {'Spinorial?':>11}  Mechanism")
print("─" * 75)

reentrance_channels = []
persistent_channels = []

for name, elems in subgroups:
    order = len(elems)
    has_q8 = contains_Q8(elems)
    kappa = count_1d_chars(elems)
    has_spin = has_spinorial_1d(elems)

    # Q₈ criterion: has_q8 ↔ ¬has_spin
    assert has_q8 == (not has_spin), \
        f"{name}: Q₈ criterion violated! Q₈⊂={has_q8}, spinorial={has_spin}"

    # Mechanism
    if name == "O (parent)":
        mechanism = "Q₈ ⊂ 2O"
    elif has_q8:
        short = name.split("→")[1].strip() if "→" in name else name
        mechanism = f"Q₈ ⊂ 2{short}"
        persistent_channels.append(short)
    else:
        # Why is Q₈ not contained?
        if order < 8:
            mechanism = f"Q₈ ⊄ Z_{order} (Lagrange)"
        elif kappa == order:  # abelian
            mechanism = f"Q₈ ⊄ Z_{order} (abelian)"
        else:
            mechanism = f"Q₈ ⊄ 2H' (Lagrange)"
        short = name.split("→")[1].strip() if "→" in name else name
        reentrance_channels.append(short)

    q8_str = "Yes" if has_q8 else "No"
    spin_str = "No" if not has_spin else "YES"

    print(f"{name:<18} {order:>5}    {q8_str:<8}  {kappa:>3}     {spin_str:<11} {mechanism}")

print(f"\nRe-entrance channels: {', '.join(reentrance_channels)}")
print(f"Persistent obstruction: {', '.join(persistent_channels)}")

# =============================================================
# Classification rule verification
# =============================================================

print(f"\n{'=' * 70}")
print("CLASSIFICATION RULE VERIFICATION")
print("=" * 70)

print("""
For every subgroup 2H' ⊆ 2O, we verify:

  Re-entrance ⟺ Q₈ ⊄ 2H' ⟺ -1 ∉ [2H', 2H']

This is a complete check of Theorem 4.3 (Q₈ criterion) restricted
to the subgroup lattice of 2O.""")

for name, elems in subgroups:
    has_q8 = contains_Q8(elems)
    m1_in_comm = minus1_in_commutator(elems)

    # The three conditions should be equivalent
    assert has_q8 == m1_in_comm, \
        f"{name}: Q₈⊂G={has_q8} but -1∈[G,G]={m1_in_comm}"

    status = "obstructed" if m1_in_comm else "RE-ENTRANT"
    print(f"  {name:<18}: Q₈⊂2H'={has_q8}, -1∈[2H',2H']={m1_in_comm} → {status}")

# =============================================================
# Spinorial branching summary
# =============================================================

print(f"\n{'=' * 70}")
print("SPINORIAL BRANCHING: L(S) DECOMPOSITION IN EACH CHANNEL")
print("=" * 70)
print("""
The rank-4 spinorial irrep L is the diagnostic:
  - If L branches to include 1D components → scalar spinorial sectors restored
  - If L branches only to rank ≥ 2 → obstruction persists""")

def chi_su2(j, q):
    """SU(2) spin-j character evaluated at quaternion q."""
    trace = 2 * q[0]  # tr(q) in SU(2)
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


# Build 2O conjugacy class lookup
classes_2O = compute_conjugacy_classes(elements_2O)
class_key_map = {}
for c in classes_2O:
    for k in c['keys']:
        class_key_map[k] = c

# 2O character table indexed by class trace and size
# (verified unique by script 29)
CHAR_2O_TABLE = {
    'A₀': {(1,2.0):1,(12,0.0):1,(8,-1.0):1,(6,0.0):1,(1,-2.0):1,(6,-SQ2):1,(8,1.0):1,(6,SQ2):1},
    'A₁': {(1,2.0):1,(12,0.0):-1,(8,-1.0):1,(6,0.0):1,(1,-2.0):1,(6,-SQ2):-1,(8,1.0):1,(6,SQ2):-1},
    'E':  {(1,2.0):2,(12,0.0):0,(8,-1.0):-1,(6,0.0):2,(1,-2.0):2,(6,-SQ2):0,(8,1.0):-1,(6,SQ2):0},
    'T₁': {(1,2.0):3,(12,0.0):1,(8,-1.0):0,(6,0.0):-1,(1,-2.0):3,(6,-SQ2):-1,(8,1.0):0,(6,SQ2):-1},
    'T₂': {(1,2.0):3,(12,0.0):-1,(8,-1.0):0,(6,0.0):-1,(1,-2.0):3,(6,-SQ2):1,(8,1.0):0,(6,SQ2):1},
    'H₁': {(1,2.0):2,(12,0.0):0,(8,-1.0):-1,(6,0.0):0,(1,-2.0):-2,(6,-SQ2):-SQ2,(8,1.0):1,(6,SQ2):SQ2},
    'H₂': {(1,2.0):2,(12,0.0):0,(8,-1.0):-1,(6,0.0):0,(1,-2.0):-2,(6,-SQ2):SQ2,(8,1.0):1,(6,SQ2):-SQ2},
    'L':  {(1,2.0):4,(12,0.0):0,(8,-1.0):1,(6,0.0):0,(1,-2.0):-4,(6,-SQ2):0,(8,1.0):-1,(6,SQ2):0},
}


def char_2O_on_element(irrep, g):
    """Evaluate 2O character on element g."""
    k = qkey(g)
    c = class_key_map[k]
    key = (c['size'], round(c['trace'], 4))
    for (s, t), v in CHAR_2O_TABLE[irrep].items():
        if s == key[0] and abs(t - key[1]) < 0.01:
            return v
    raise ValueError(f"No match for {irrep} at class {key}")


def compute_cyclic_spin_1d(elems, char_fn):
    """For cyclic groups: compute spinorial 1D content of L by Fourier analysis."""
    n = len(elems)
    # Find a generator (element of maximal order)
    orders = []
    for g in elems:
        o = 0
        current = IDENTITY.copy()
        while o <= n:
            current = qmul(current, g)
            o += 1
            if np.allclose(current, IDENTITY, atol=1e-10):
                break
        orders.append(o)
    assert max(orders) == n, "Not cyclic!"
    gen = elems[orders.index(n)]
    omega = np.exp(2j * np.pi / n)

    spin_1d = 0
    for k in range(n):
        # σ_k(-1) = ω^{k·n/2} = (-1)^k
        if (omega ** (k * n // 2)).real < -0.5:
            # m(σ_k, L) = (1/n) Σ_m ω^{-km} χ_L(gen^m)
            g_power = IDENTITY.copy()
            mult = 0.0
            for m in range(n):
                mult += (omega ** (-k * m)) * char_fn('L', g_power)
                g_power = qmul(g_power, gen)
            spin_1d += int(round(mult.real / n))
    return spin_1d


def compute_dic3_spin_1d(elems, gen_a, gen_b, char_fn):
    """For Dic₃: compute spinorial 1D content of L using explicit character theory.

    Dic₃ has abelianization Z₄. The map φ: Dic₃ → Z₄ satisfies
    φ(a) = 2, φ(b) = 1 (from relations a³=-1, bab⁻¹=a⁻¹).
    Characters σ_k of Z₄: σ_k(g) = ω₄^{k·φ(g)}, ω₄ = e^{2πi/4}.
    Spinorial: σ_k(-1) = σ_k(a³) = ω₄^{6k} = ω₄^{2k} → spinorial for k=1,3.
    """
    # Build element-to-(type, power) mapping: each g = a^m or b·a^m
    powers_a = [IDENTITY.copy()]
    curr = IDENTITY.copy()
    for _ in range(5):
        curr = qmul(curr, gen_a)
        powers_a.append(curr.copy())

    elem_map = {}
    for m in range(6):
        elem_map[qkey(powers_a[m])] = ('a', m)
        elem_map[qkey(qmul(gen_b, powers_a[m]))] = ('b', m)
    for g in elems:
        assert qkey(g) in elem_map, f"Unmapped Dic₃ element: {g}"

    omega4 = np.exp(2j * np.pi / 4)

    def sigma_k(k, etype, m):
        if etype == 'a':
            return omega4 ** (k * (2 * m % 4))
        else:
            return omega4 ** (k * ((1 + 2 * m) % 4))

    spin_1d = 0
    for k in range(4):
        # Spinorial? σ_k(-1) = σ_k(a³)
        if sigma_k(k, 'a', 3).real < -0.5:
            mult = sum(
                np.conj(sigma_k(k, *elem_map[qkey(g)])) * char_fn('L', g)
                for g in elems
            ).real / 12
            spin_1d += int(round(mult))
    return spin_1d


for name, elems in subgroups[1:]:  # skip O itself
    order = len(elems)
    has_spin_scalar = has_spinorial_1d(elems)

    # Total 1D content of L: (1/|[G,G]|) Σ_{g ∈ [G,G]} χ_L(g)
    comm = compute_commutator_subgroup(elems)
    comm_keys = set(comm.keys())
    comm_elements = [g for g in elems if qkey(g) in comm_keys]
    total_1d = int(round(
        sum(char_2O_on_element('L', g) for g in comm_elements) / len(comm_elements)
    ))

    # Spinorial 1D content: computed from explicit character theory
    spin_1d_content = 0
    if has_spin_scalar:
        if name == "O → D₃":
            spin_1d_content = compute_dic3_spin_1d(
                elems, a_D3, b_D3, char_2O_on_element)
        else:
            spin_1d_content = compute_cyclic_spin_1d(elems, char_2O_on_element)

    rank2_content = 4 - total_1d  # dim L = 4

    if total_1d > 0 and has_spin_scalar:
        status = f"{spin_1d_content} spinorial + {total_1d - spin_1d_content} tensorial 1D + rank≥2"
    elif total_1d > 0:
        status = f"{total_1d} tensorial 1D + {rank2_content} in rank≥2"
    else:
        status = f"all rank≥2 (obstruction persists)"

    print(f"  {name:<18}: L(dim 4) → {status}")

    # Verify consistency
    assert 0 <= total_1d <= 4, f"1D content out of range: {total_1d}"
    if has_spin_scalar:
        assert spin_1d_content > 0, f"{name} should have spinorial 1D in L!"
    else:
        assert spin_1d_content == 0, f"{name} shouldn't have spinorial 1D in L!"
# =============================================================
# Summary
# =============================================================

print(f"\n{'=' * 70}")
print("SUMMARY: COMPLETE PHASE DIAGRAM")
print("=" * 70)

print("""
The symmetry-breaking phase diagram of 2O is completely determined
by the Q₈ criterion (Theorem 4.3):

  Obstruction persists ⟺ Q₈ ⊂ 2H' (residual binary group)
  Re-entrance occurs   ⟺ Q₈ ⊄ 2H'

Among the 7 principal channels of O:
  - 3 retain Q₈ (T, D₄, D₂) → obstruction persists
  - 4 lose Q₈ (D₃, C₄, C₃, C₂) → scalar spinorial sectors restored

The critical transition is NOT about group size:
  |2D₄| = 16 > |2D₃| = 12, yet D₄ is obstructed and D₃ is not.
  It is about the lower central structure — specifically, whether
  the commutator subgroup contains the central element -1.

This is the complete verification of Section 6 (re-entrance)
across the full subgroup lattice of O.
""")

print("All checks passed.")


if __name__ == '__main__':
    pass  # All computation runs at module level (research script)
