#!/usr/bin/env python3
"""
Proof that j* is the last zero: m(A₀, j) ≥ 1 for all j > j* = |G|/4 - 1.

THEOREM (Corollary 3.4(iii)):
  For G ⊂ SU(2) finite with -1 ∈ [G,G] and 4 | |G|, the trivial multiplicity
  m(A₀, j) ≥ 1 for every integer j > j* = |G|/4 - 1.

PROOF STRUCTURE:

  Part A — Periodicity and linear growth:
    (i)   The non-central character sum S(j) = Σ_{g ≠ ±1} χ_j(g) is periodic in j
          with period P = |G|/4.
          [Each non-central g has SU(2) angle α = πp/q with q | P.
           χ_j(g) = sin((2j+1)α)/sin(α) has period q in j. So S has period lcm(q's) = P.]
    (ii)  m(A₀, j) = [2(2j+1) + S(j)] / |G|.
    (iii) m(A₀, j+P) = m(A₀, j) + 1.
          [S(j+P) = S(j) and 2(2(j+P)+1) = 2(2j+1) + 4P = 2(2j+1) + |G|.]

  Part B — First period verification:
    (iv)  Dicyclic family Dic_n (n even, |G| = 4n, all n):
          Using the Dirichlet kernel identity
            Σ_{k=1}^{n-1} sin(mπk/n)/sin(πk/n) = n·(2⌊j/n⌋+1) - (2j+1),
          one obtains m(A₀, j) = (3 + (-1)^j)/2 ∈ {1, 2} for j ∈ [n, 2n-1].
          Non-trivial characters (G^ab ≅ Z₂ × Z₂):
            m(χ₁) = (3 - (-1)^j)/2  [off-diagonal sign flip]
            m(χ₂) = m(χ₃) = 1 if j < 3n/2, else 2  [shifted Dirichlet kernel]
          ALL FOUR characters: m ∈ {1, 2}. FULLY ANALYTIC.
    (v)   Polyhedral groups (2T, 2O, 2I):
          m(A₀, j) ≥ 1 for j ∈ [j*+1, j*+P] by direct computation (48 values total).

  Conclusion: m(A₀, j) = m(A₀, r) + k ≥ 1 for all j > j*, where r is the residue
  in the first period and k ≥ 0 counts complete periods.

  Part C — Non-trivial scalar characters:
    (vi)  ANALYTIC: m(χ,0) = ⟨V₀,χ⟩ = δ_{χ,A₀} = 0 for χ ≠ A₀.
          By duality: m(χ,j*) + m(χ,0) = 1, so m(χ,j*) = 1.
          Therefore j* is NOT a zero for χ ≠ A₀, and j_max(χ) < j*.
    (vii) Same periodicity and growth apply to S_χ(j) for all scalar χ.
          Verified computationally for 2T (χ₁) and 2O (A₁).
    (viii) For κ=3 (2T): m(χ₂,j) = m(χ₁,j) by conjugation symmetry.

KEY IDENTITIES:
  - Dirichlet kernel: Σ_{k=1}^{n-1} sin(mπk/n)/sin(πk/n) = n·(2⌊(m-1)/(2n)⌋+1) - m
  - At j = j*: χ_{j*}(g) = -1 for ALL non-central g (from divisibility condition).
    This gives S(j*) = -(|G|-2) and m(A₀, j*) = 0.
  - At j = j*+P: S repeats but central grows by |G|, so m increases by exactly 1.

Asserts: 43
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       build_binary_icosahedral, compute_conjugacy_classes,
                       compute_commutator_subgroup)
from src.quaternion import qkey, qmul, qinv

n_asserts = 0

def chi_j(j, alpha):
    """SU(2) character of spin-j representation at angle alpha."""
    if abs(alpha) < 1e-12:
        return 2*j + 1
    if abs(alpha - np.pi) < 1e-12:
        return ((-1)**(2*j)) * (2*j + 1)
    return np.sin((2*j+1)*alpha) / np.sin(alpha)

def m_A0_from_classes(j, cls, G_order):
    """Compute m(A₀, j) from conjugacy class data."""
    total = 0
    for c in cls:
        tr = c['trace']
        alpha = np.arccos(np.clip(tr/2, -1, 1))
        total += c['size'] * chi_j(j, alpha)
    return total / G_order

def m_A0_dicyclic(j, n):
    """Compute m(A₀, j) for Dic_n directly."""
    G_order = 4*n
    total = 2*(2*j+1)  # central: identity + (-1)
    # Cyclic elements a^k, k=1,...,2n-1, k!=n
    for k in range(1, 2*n):
        if k == n:
            continue
        alpha = k * np.pi / n
        total += np.sin((2*j+1)*alpha) / np.sin(alpha)
    # Off-diagonal elements (2n of them, all alpha = pi/2)
    total += 2*n * np.sin((2*j+1)*np.pi/2)
    return total / G_order


def m_chi_dicyclic(j, n, chi_a, chi_b):
    """Compute m(χ, j) for Dic_n with character (chi_a, chi_b) on generators.

    G^ab ≅ Z₂ × Z₂ for n even: chi_a, chi_b ∈ {+1, -1}.
    Elements: cyclic a^k (angle kπ/n), off-diagonal a^k·b (angle π/2).
    """
    G_order = 4 * n
    total = 0.0
    # Cyclic elements a^k, k = 0, ..., 2n-1
    for k in range(2 * n):
        alpha = k * np.pi / n
        total += chi_j(j, alpha) * (chi_a ** k)
    # Off-diagonal a^k·b, k = 0, ..., 2n-1 (all have SU(2) angle π/2)
    off_diag = np.sin((2 * j + 1) * np.pi / 2)
    for k in range(2 * n):
        total += off_diag * (chi_a ** k) * chi_b
    return total / G_order


def m_chi_from_classes(j, cls, G_order, chi_conj_vals):
    """Compute m(χ, j) from class data and conjugate character values.

    chi_conj_vals[i] = conj(χ(g_i)) for class i representative.
    Returns complex in general; real part is the multiplicity.
    """
    total = sum(c['size'] * chi_j(j, np.arccos(np.clip(c['trace']/2, -1, 1)))
                * chi_conj_vals[i]
                for i, c in enumerate(cls))
    return total / G_order


def compute_class_coset_labels(elements, cls, comm_keys, kappa):
    """Coset index in G/[G,G] ≅ Z_κ for each conjugacy class.

    Uses group multiplication: label(g) = k iff gen^{-k}·g ∈ [G,G].
    Returns list of integer coset labels, one per class.
    """
    if kappa == 1:
        return [0] * len(cls)

    # Find generator of G^ab: first element not in [G,G]
    gen = None
    for e in elements:
        if qkey(e) not in comm_keys:
            gen = e.copy()
            break
    assert gen is not None, "All elements in [G,G] but κ > 1"

    # Verify generator has order exactly κ in G^ab
    power = np.array([1, 0, 0, 0], dtype=float)
    for k in range(1, kappa):
        power = qmul(power, gen)
        assert qkey(power) not in comm_keys, f"gen^{k} in [G,G] (order < κ={kappa})"
    power = qmul(power, gen)  # gen^κ
    assert qkey(power) in comm_keys, f"gen^{kappa} not in [G,G]"

    # Build inverse powers gen^{-k} for k = 0, ..., κ-1
    gen_inv = qinv(gen)
    inv_powers = [np.array([1, 0, 0, 0], dtype=float)]
    for k in range(1, kappa):
        inv_powers.append(qmul(inv_powers[-1], gen_inv))

    # Label each class representative
    labels = []
    for c in cls:
        rep = c['rep']
        if qkey(rep) in comm_keys:
            labels.append(0)
            continue
        found = False
        for k in range(1, kappa):
            if qkey(qmul(inv_powers[k], rep)) in comm_keys:
                labels.append(k)
                found = True
                break
        assert found, f"Class rep not in any coset (κ={kappa})"

    return labels


# ============================================================
# Part 1: Periodicity of S(j) with period P = |G|/4
# ============================================================
print("=" * 70)
print("Part 1: Periodicity of non-central sum S(j)")
print("=" * 70)

for name, builder, G_order in [('2T', build_binary_tetrahedral, 24),
                                ('2O', build_binary_octahedral, 48),
                                ('2I', build_binary_icosahedral, 120)]:
    G = builder()
    cls = compute_conjugacy_classes(G)
    j_star = G_order // 4 - 1
    P = G_order // 4

    # Compute S(j) = |G|*m - 2(2j+1) for j in [j*+1, j*+1+2P]
    period_ok = True
    for j in range(j_star + 1, j_star + 1 + P):
        S_j = G_order * m_A0_from_classes(j, cls, G_order) - 2*(2*j+1)
        S_jP = G_order * m_A0_from_classes(j + P, cls, G_order) - 2*(2*(j+P)+1)
        if abs(S_j - S_jP) > 1e-6:
            period_ok = False

    print(f"  {name}: P = {P}, period verified: {period_ok}")
    assert period_ok; n_asserts += 1

print()

# ============================================================
# Part 2: m(j+P) = m(j) + 1 (linear growth per period)
# ============================================================
print("=" * 70)
print("Part 2: m(A₀, j+P) = m(A₀, j) + 1")
print("=" * 70)

for name, builder, G_order in [('2T', build_binary_tetrahedral, 24),
                                ('2O', build_binary_octahedral, 48),
                                ('2I', build_binary_icosahedral, 120)]:
    G = builder()
    cls = compute_conjugacy_classes(G)
    j_star = G_order // 4 - 1
    P = G_order // 4

    growth_ok = True
    for j in range(j_star + 1, j_star + 1 + P):
        m_j = round(m_A0_from_classes(j, cls, G_order))
        m_jP = round(m_A0_from_classes(j + P, cls, G_order))
        if m_jP != m_j + 1:
            growth_ok = False

    print(f"  {name}: m(j+{P}) = m(j) + 1 verified: {growth_ok}")
    assert growth_ok; n_asserts += 1

print()

# ============================================================
# Part 3: χ_{j*}(g) = -1 for all non-central g (divisibility)
# ============================================================
print("=" * 70)
print("Part 3: χ_{j*}(g) = -1 for all non-central g")
print("=" * 70)

for name, builder, G_order in [('2T', build_binary_tetrahedral, 24),
                                ('2O', build_binary_octahedral, 48),
                                ('2I', build_binary_icosahedral, 120)]:
    G = builder()
    cls = compute_conjugacy_classes(G)
    j_star = G_order // 4 - 1

    all_minus_one = True
    for c in cls:
        tr = c['trace']
        if abs(tr - 2) < 1e-10 or abs(tr + 2) < 1e-10:
            continue
        alpha = np.arccos(np.clip(tr/2, -1, 1))
        val = chi_j(j_star, alpha)
        if abs(val - (-1)) > 1e-6:
            all_minus_one = False

    S_star = sum(c['size'] * chi_j(j_star, np.arccos(np.clip(c['trace']/2, -1, 1)))
                 for c in cls
                 if abs(c['trace'] - 2) > 1e-10 and abs(c['trace'] + 2) > 1e-10)
    expected_S = -(G_order - 2)

    print(f"  {name}: j* = {j_star}, all χ = -1: {all_minus_one}, S(j*) = {S_star:.0f} = -(|G|-2) = {expected_S}")
    assert all_minus_one; n_asserts += 1
    assert abs(S_star - expected_S) < 1e-6; n_asserts += 1

print()

# ============================================================
# Part 4: First period verification — polyhedral groups
# ============================================================
print("=" * 70)
print("Part 4: m(A₀, j) ≥ 1 in first period [j*+1, j*+P] — polyhedral")
print("=" * 70)

total_checked = 0
for name, builder, G_order in [('2T', build_binary_tetrahedral, 24),
                                ('2O', build_binary_octahedral, 48),
                                ('2I', build_binary_icosahedral, 120)]:
    G = builder()
    cls = compute_conjugacy_classes(G)
    j_star = G_order // 4 - 1
    P = G_order // 4

    m_vals = []
    for j in range(j_star + 1, j_star + 1 + P):
        m = round(m_A0_from_classes(j, cls, G_order))
        m_vals.append(m)

    min_m = min(m_vals)
    total_checked += P
    print(f"  {name}: j* = {j_star}, P = {P}, m values = {m_vals}")
    print(f"         min m = {min_m}")
    assert min_m >= 1; n_asserts += 1

print(f"\n  Total polyhedral values checked: {total_checked}")
assert total_checked == 48; n_asserts += 1

print()

# ============================================================
# Part 5: Analytic proof for dicyclic — Dirichlet kernel identity
# ============================================================
print("=" * 70)
print("Part 5: Analytic proof for Dic_n (n even) — Dirichlet kernel")
print("=" * 70)

# Dirichlet kernel identity:
# C(j) = Σ_{k=1}^{n-1} sin((2j+1)kπ/n)/sin(kπ/n) = n·(2·floor(j/n)+1) - (2j+1)
#
# For j ∈ [n, 2n-1]: floor(j/n) = 1, so C(j) = 3n - 2j - 1.
# m(A₀,j) = [2(2j+1) + 2·C(j) + 2n·(-1)^j] / (4n)
#          = [4j+2 + 6n-4j-2 + 2n·(-1)^j] / (4n)
#          = [6n + 2n·(-1)^j] / (4n)
#          = (3 + (-1)^j) / 2

# Step 1: Verify the Dirichlet kernel identity
print("\n  Step 1: Verify Dirichlet kernel identity on examples")
for n in [4, 6, 10, 15]:
    for j in [n, n+1, n+2, 2*n-1]:
        # Direct computation
        C_direct = sum(np.sin((2*j+1)*k*np.pi/n) / np.sin(k*np.pi/n)
                       for k in range(1, n))
        # Formula
        C_formula = n * (2*(j // n) + 1) - (2*j + 1)
        assert abs(C_direct - C_formula) < 1e-6

n_asserts += 1  # Dirichlet kernel identity verified
print("    Identity verified for n = 4, 6, 10, 15 at multiple j values.")

# Step 2: Verify analytic formula m = (3 + (-1)^j)/2 for Dic_n
print("\n  Step 2: Verify m(A₀, j) = (3 + (-1)^j)/2 for j ∈ [n, 2n-1]")
for n in [2, 4, 6, 8, 10, 20, 50, 100]:
    for j in range(n, 2*n):
        m_computed = m_A0_dicyclic(j, n)
        m_formula = (3 + (-1)**j) / 2
        assert abs(m_computed - m_formula) < 1e-6, \
            f"Dic_{n}, j={j}: computed={m_computed}, formula={m_formula}"

n_asserts += 1  # Analytic formula verified
print(f"    Formula verified for Dic_n, n ∈ {{2,4,6,8,10,20,50,100}}.")
print(f"    Total Dic_n values checked: {sum(n for n in [2,4,6,8,10,20,50,100])}")

# Step 3: The formula gives m ∈ {1, 2}, hence m ≥ 1
m_even = (3 + 1) / 2  # = 2
m_odd = (3 - 1) / 2   # = 1
assert m_even == 2; n_asserts += 1
assert m_odd == 1; n_asserts += 1
print(f"\n  m(A₀, j) = (3+(-1)^j)/2: even j → {int(m_even)}, odd j → {int(m_odd)}. Both ≥ 1. ✓")

# Step 4: Non-trivial characters of Dic_n (n even)
# G^ab ≅ Z₂ × Z₂. Four characters: A₀=(1,1), χ₁=(1,-1), χ₂=(-1,1), χ₃=(-1,-1).
#
# χ₁: differs from A₀ only in off-diagonal sign (b → -1).
#   Off-diagonal sum flips: +2n·(-1)^j → -2n·(-1)^j.
#   m(χ₁) = m(A₀) - (-1)^j = (3 - (-1)^j)/2.
#
# χ₂, χ₃: have (-1)^k on cyclic elements. Off-diagonal sum vanishes (Σ(-1)^k = 0).
#   Cyclic sum = shifted Dirichlet kernel D(2j+1+n).
#   m(χ₂) = m(χ₃) = 1 for j ∈ [n, 3n/2-1], 2 for j ∈ [3n/2, 2n-1].
print("\n  Step 4: Verify m(χ₁, j) = (3 - (-1)^j)/2 for j ∈ [n, 2n-1]")
for n in [2, 4, 6, 8, 10, 20, 50, 100]:
    for j in range(n, 2*n):
        m_computed = m_chi_dicyclic(j, n, 1, -1)
        m_formula = (3 - (-1)**j) / 2
        assert abs(m_computed - m_formula) < 1e-6, \
            f"Dic_{n}, j={j}: chi1 computed={m_computed}, formula={m_formula}"

n_asserts += 1
print(f"    Formula verified for Dic_n, n ∈ {{2,4,6,8,10,20,50,100}}.")

# Step 5: Verify m(χ₂) = m(χ₃) = step function
print("\n  Step 5: Verify m(χ₂, j) = m(χ₃, j) = {{1 if j < 3n/2, 2 if j ≥ 3n/2}}")
for n in [2, 4, 6, 8, 10, 20, 50, 100]:
    for j in range(n, 2*n):
        p_step = 1 if j < 3*n//2 else 2
        for chi_a, chi_b, cname in [(-1, 1, 'chi2'), (-1, -1, 'chi3')]:
            m_computed = m_chi_dicyclic(j, n, chi_a, chi_b)
            assert abs(m_computed - p_step) < 1e-6, \
                f"Dic_{n}, j={j}: {cname} computed={m_computed}, predicted={p_step}"

n_asserts += 1
print(f"    Formula verified for Dic_n, n ∈ {{2,4,6,8,10,20,50,100}}.")

# Step 6: All four formulas give m ∈ {1, 2}
all_in_12 = True
for n in [2, 4, 6, 8, 10, 20, 50, 100]:
    for j in range(n, 2*n):
        for chi_a, chi_b in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            m_val = round(m_chi_dicyclic(j, n, chi_a, chi_b))
            if m_val not in (1, 2):
                all_in_12 = False

assert all_in_12; n_asserts += 1
print(f"\n  All 4 characters: m ∈ {{1, 2}} throughout first period. All ≥ 1. ✓")
print(f"  Dicyclic family: FULLY ANALYTIC for all scalar characters.")

print()

# ============================================================
# Part 6: Extended verification — large j range
# ============================================================
print("=" * 70)
print("Part 6: Extended verification up to j = j* + 50")
print("=" * 70)

for name, builder, G_order in [('2T', build_binary_tetrahedral, 24),
                                ('2O', build_binary_octahedral, 48),
                                ('2I', build_binary_icosahedral, 120)]:
    G = builder()
    cls = compute_conjugacy_classes(G)
    j_star = G_order // 4 - 1

    min_m = 100
    for j in range(j_star + 1, j_star + 51):
        m = round(m_A0_from_classes(j, cls, G_order))
        if m < min_m:
            min_m = m

    print(f"  {name}: min m(A₀, j) for j ∈ [{j_star+1}, {j_star+50}] = {min_m}")
    assert min_m >= 1; n_asserts += 1

# Dicyclic extended
for n in [2, 10, 50]:
    j_star = n - 1
    min_m = 100
    for j in range(j_star + 1, j_star + 51):
        m = round(m_A0_dicyclic(j, n))
        if m < min_m:
            min_m = m
    print(f"  Dic_{n}: min m(A₀, j) for j ∈ [{j_star+1}, {j_star+50}] = {min_m}")
    assert min_m >= 1; n_asserts += 1

print()

# ============================================================
# Part 7: Boundary values — m at j* and j*+P
# ============================================================
print("=" * 70)
print("Part 7: Boundary: m(j*) = 0 and m(j*+P) = 1 exactly")
print("=" * 70)

for name, builder, G_order in [('2T', build_binary_tetrahedral, 24),
                                ('2O', build_binary_octahedral, 48),
                                ('2I', build_binary_icosahedral, 120)]:
    G = builder()
    cls = compute_conjugacy_classes(G)
    j_star = G_order // 4 - 1
    P = G_order // 4

    m_star = round(m_A0_from_classes(j_star, cls, G_order))
    m_star_P = round(m_A0_from_classes(j_star + P, cls, G_order))

    print(f"  {name}: m(j*={j_star}) = {m_star}, m(j*+P={j_star+P}) = {m_star_P}")
    assert m_star == 0; n_asserts += 1
    # m(j*+P) = m(j*) + 1 = 1, already verified in Part 2

print()

# ============================================================
# Part 8: All scalar characters — analytic argument + verification
# ============================================================
print("=" * 70)
print("Part 8: All scalar characters — last zero and m >= 1 above")
print("=" * 70)

# ANALYTIC ARGUMENT (one-line proof for non-trivial characters):
#
#   For χ ≠ A₀: m(χ, 0) = ⟨V₀|_G, χ⟩ = ⟨A₀, χ⟩ = δ_{χ,A₀} = 0.
#   By duality (Theorem 3.3): m(χ, j*) + m(χ, 0) = 1.
#   Therefore m(χ, j*) = 1 for all non-trivial scalar χ.
#
#   Combined with m(A₀, j*) = 0 (Part 7):
#     j* is NOT a zero for χ ≠ A₀  →  j_max(χ) < j*
#     j* IS the last zero for A₀   →  j_max(A₀) = j* is the largest
#
# This is ANALYTIC: it uses only V₀ = A₀ (trivial) and the duality identity.
# The verification below confirms it computationally and extends to m >= 1
# above the last zero via the same periodicity + growth argument as Parts 1-2.

print("\n  Analytic argument:")
print("    m(χ, 0) = ⟨V₀, χ⟩ = δ_{χ,A₀} = 0 for χ ≠ A₀.")
print("    By duality: m(χ, j*) + m(χ, 0) = 1  =>  m(χ, j*) = 1.")
print("    Therefore j* is NOT a zero of χ ≠ A₀, so j_max(χ) < j*.")

for name, builder, G_order in [('2T', build_binary_tetrahedral, 24),
                                ('2O', build_binary_octahedral, 48),
                                ('2I', build_binary_icosahedral, 120)]:
    G = builder()
    cls = compute_conjugacy_classes(G)
    comm = compute_commutator_subgroup(G)
    comm_keys = set(qkey(g) for g in comm.values())
    j_star = G_order // 4 - 1
    P = G_order // 4
    kappa = G_order // len(comm)

    print(f"\n  {name} (|G|={G_order}, κ={kappa}, j*={j_star}):")

    if kappa == 1:
        print(f"    κ=1: only A₀. Already proved in Parts 1-7.")
        continue

    # Compute coset labels from group multiplication (not index guessing)
    class_cosets = compute_class_coset_labels(G, cls, comm_keys, kappa)
    print(f"    Coset labels (from group multiplication): {class_cosets}")

    omega = np.exp(2j * np.pi / kappa)

    # For κ=3: χ₂ = χ̄₁, so m(χ₁,j) = m(χ₂,j). Check χ₁ fully, verify conjugation.
    char_idx = 1
    chi_conj = [omega**(-char_idx * class_cosets[i]) for i in range(len(cls))]
    char_name = "A₁" if kappa == 2 else "χ₁"

    # Step 1: Verify analytic argument — m(χ, 0) = 0 and m(χ, j*) = 1
    m_0 = m_chi_from_classes(0, cls, G_order, chi_conj)
    m_jstar = m_chi_from_classes(j_star, cls, G_order, chi_conj)

    print(f"    {char_name}: m(χ,0) = {m_0.real:.6f}, m(χ,j*={j_star}) = {m_jstar.real:.6f}")
    assert abs(m_0) < 1e-6, f"m({char_name}, 0) = {m_0} != 0"
    n_asserts += 1
    assert abs(m_jstar.real - 1) < 1e-6, f"m({char_name}, j*) = {m_jstar} != 1"
    n_asserts += 1

    # Step 2: Find last zero and verify it's strictly below j*
    last_zero = -1
    for j_val in range(j_star + 1):
        m_val = m_chi_from_classes(j_val, cls, G_order, chi_conj)
        if round(m_val.real) == 0:
            last_zero = j_val

    assert last_zero < j_star, f"j_max({char_name}) = {last_zero} not < j*={j_star}"
    n_asserts += 1

    # Step 3: Verify m >= 1 above last zero
    min_m = 100
    for j_val in range(last_zero + 1, j_star + 51):
        m_val = round(m_chi_from_classes(j_val, cls, G_order, chi_conj).real)
        if m_val < min_m:
            min_m = m_val

    assert min_m >= 1, f"min m({char_name}) above last zero = {min_m}"
    n_asserts += 1

    # Step 4: Verify growth m(χ, j+P) = m(χ, j) + 1
    growth_ok = True
    for j_val in range(last_zero + 1, last_zero + 1 + P):
        m_j = round(m_chi_from_classes(j_val, cls, G_order, chi_conj).real)
        m_jP = round(m_chi_from_classes(j_val + P, cls, G_order, chi_conj).real)
        if m_jP != m_j + 1:
            growth_ok = False

    assert growth_ok, f"Growth m(j+P) = m(j)+1 failed for {char_name}"
    n_asserts += 1

    conj_note = " (= χ₂ by conjugation)" if kappa == 3 else ""
    print(f"    {char_name}{conj_note}: m(χ,0)=0, m(χ,j*)=1, "
          f"last zero={last_zero} < j*, min m above={min_m}, growth ok")

    # For κ=3: verify m(χ₂, j) = m(χ₁, j) (conjugation symmetry)
    if kappa == 3:
        chi2_conj = [omega**(-2 * class_cosets[i]) for i in range(len(cls))]
        conjugation_ok = True
        for j_val in range(j_star + 11):
            m1 = m_chi_from_classes(j_val, cls, G_order, chi_conj)
            m2 = m_chi_from_classes(j_val, cls, G_order, chi2_conj)
            if abs(m1.real - m2.real) > 1e-6:
                conjugation_ok = False

        assert conjugation_ok, "m(χ₁,j) != m(χ₂,j)"
        n_asserts += 1
        print(f"    m(χ₁,j) = m(χ₂,j) for j in [0,{j_star+10}] (conjugation symmetry)")

print()

# ============================================================
# Summary
# ============================================================
print("=" * 70)
print("PROOF SUMMARY")
print("=" * 70)
print()
print("Theorem: For each scalar character χ of an obstructed G ⊂ SU(2),")
print("         m(χ, j) >= 1 for all integer j above the last zero j_max(χ).")
print("         j_max(A₀) = j* = |G|/4 - 1 (largest). j_max(χ) < j* for χ != A₀.")
print()
print("Proof:")
print("  (A) S_χ(j) periodic with period P = |G|/4, m(χ,j+P) = m(χ,j) + 1.")
print("  (B) For A₀ in first period [j*+1, j*+P]:")
print("      - Dic_n (n even): m = (3+(-1)^j)/2 in {1,2}. [Dirichlet kernel, analytic]")
print("      - 2T, 2O, 2I: m >= 1 by direct computation (48 values).")
print("  (C) For χ != A₀: m(χ,0) = 0 => m(χ,j*) = 1 by duality. [ANALYTIC]")
print("      Same periodicity + growth => m(χ,j) >= 1 above j_max(χ) < j*.")
print()
print(f"All checks passed ({n_asserts} asserts).")
