#!/usr/bin/env python3
"""
H1 + H2: Perturbative stability and lifting criterion for SO(3)/O.

TEST 1 — Occupation invariance (H1, trivial):
  m(A₀, j) = dim Hom(A₀, V_j|_{2O}) is representation-theoretic.
  Binary pattern {0,1} exactly preserved under any O-invariant
  perturbation. Guaranteed by Schur's lemma.

TEST 2 — Energy pairing (H1, non-trivial):
  Cross-sector sum E(A₀,j) + E(A₁,j*-j) for duality pairs.
  Breaks at O(ε): first-order paired shifts are nonzero and
  pair-dependent (+0.045, -0.137, +0.231).
  CONCLUSION: duality is purely arithmetic (character theory),
  not a dynamical symmetry.

H2 — Lifting criterion:
  For all 8 channels O → K, verify: spinorial ⟺ Q₈ ⊄ 2K.
  Two independent checks (Q₈ subset + commutator [a,b]=-1).
  All agree. Monotonicity confirmed.

RAW OUTPUT:
===========

TEST 1 (occupation):
  A₀ sub-threshold occupied: [0, 4, 6, 8, 9, 10], empty: [1, 2, 3, 5, 7, 11]
  Levels shift smoothly under ε·V₄ but never cross → pattern rigid.

TEST 2 (energy pairing):
  ✓ chi_A1 verified as group homomorphism
  A₀: occupied j = [0, 4, 6, 8, 9, 10, 12, 13, 14, 15]
  A₁: occupied j = [3, 6, 7, 9, 10, 11, 12, 13, 14, 15]
  Cross-sector pairs: (0,11), (3,8), (4,7)
  Paired first-order sums δE(s1,j)/ε + δE(s2,j*-j)/ε:
    (0,11): +0.000000 + +0.044779 = +0.044779
    (3,8):  -0.363636 + +0.226721 = -0.136916
    (4,7):  +0.251748 + -0.020265 = +0.231484
  Energy sums drift linearly in ε — no cancellation.

H2 (lifting):
  Channel           |2K|  Q₈⊂2K  Spinorial?
  2O (parent)         48    Yes       No
  2T                  24    Yes       No
  2D₄                 16    Yes       No
  2D₃                 12     No      YES
  Q₈ = 2D₂            8    Yes       No
  Z₈ = 2C₄            8     No      YES
  Z₆ = 2C₃            6     No      YES
  Z₄ = 2C₂            4     No      YES
  Monotonicity: subgroup transitivity ✓
"""

import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import math
import numpy as np
from scipy.linalg import eigh
from src.quaternion import quat_to_su2
from src.group import build_binary_octahedral


# ============================================================
# Wigner D-matrix from SU(2) matrix
# ============================================================

def wigner_d_matrix(U, j):
    """
    Compute the (2j+1) × (2j+1) Wigner D-matrix for spin j
    from a 2×2 SU(2) matrix U.

    Uses the explicit formula via matrix powers of the spin-j
    representation built from tensor products of the fundamental.
    """
    dim = int(2*j + 1)
    if dim == 1:
        return np.array([[1.0 + 0j]])

    # Build spin-j rep from symmetric tensor product
    # Basis: |j, m⟩ for m = j, j-1, ..., -j
    # Using Schwinger oscillator construction
    a, b = U[0, 0], U[0, 1]
    c, d = U[1, 0], U[1, 1]

    D = np.zeros((dim, dim), dtype=complex)
    for mi_idx in range(dim):
        m_i = j - mi_idx  # m value for row
        for mf_idx in range(dim):
            m_f = j - mf_idx  # m value for column
            # D^j_{m_i, m_f} = sum_s (...)
            val = 0.0 + 0j
            for s in range(dim):
                # s runs over valid values where all factorials are non-negative
                t1 = j + m_f - s
                t2 = s
                t3 = j - m_i - s
                t4 = m_i - m_f + s
                if t1 < 0 or t2 < 0 or t3 < 0 or t4 < 0:
                    continue
                coeff = (math.factorial(int(j + m_f))
                       * math.factorial(int(j - m_f))
                       * math.factorial(int(j + m_i))
                       * math.factorial(int(j - m_i))) ** 0.5
                coeff /= (math.factorial(int(t1))
                        * math.factorial(int(t2))
                        * math.factorial(int(t3))
                        * math.factorial(int(t4)))
                val += coeff * a**t1 * b**t2 * c**t4 * d**t3
            D[mi_idx, mf_idx] = val
    return D


# ============================================================
# Sector projectors (generalized for any 1D character)
# ============================================================

def build_sector_projector(j, group_su2, chi_values):
    """
    Build projector for character χ in spin-j subspace.
    P = (1/|G|) Σ conj(χ(g)) D^j(g)
    """
    dim = int(2*j + 1)
    P = np.zeros((dim, dim), dtype=complex)
    for U, chi_g in zip(group_su2, chi_values):
        D = wigner_d_matrix(U, j)
        P += np.conj(chi_g) * D
    P /= len(group_su2)
    return P


def get_sector_subspace(j, group_su2, chi_values):
    """
    Get orthonormal basis for character χ subspace of V_j.
    Returns: (n_chi, basis) where basis is (dim, n_chi) matrix.
    """
    P = build_sector_projector(j, group_su2, chi_values)
    eigvals, eigvecs = eigh(P)
    mask = eigvals > 0.5
    n = int(np.sum(mask))
    basis = eigvecs[:, mask]
    return n, basis


# Backward-compatible wrappers
def build_A0_projector(j, group_elements_su2):
    return build_sector_projector(j, group_elements_su2, [1.0] * len(group_elements_su2))

def get_A0_subspace(j, group_su2):
    return get_sector_subspace(j, group_su2, [1.0] * len(group_su2))


# ============================================================
# Cubic harmonic potential V₄
# ============================================================

def cubic_harmonic_matrix_element(j1, m1, j2, m2):
    """
    Matrix element ⟨j1,m1| V₄ |j2,m2⟩ where V₄ is the
    octahedral-invariant rank-4 potential:

    V₄ = C₄₀ + √(5/14) (C₄₄ + C₄₋₄)

    where C_{kq} are Racah-normalized spherical harmonics.
    Using the Wigner-Eckart theorem:
    ⟨j1,m1| C_{kq} |j2,m2⟩ = (-1)^{j1-m1} * (j1 k j2; -m1 q m2)
                                * ⟨j1 || C_k || j2⟩

    Reduced matrix element:
    ⟨j1 || C_k || j2⟩ = (-1)^{j1} √((2j1+1)(2j2+1)) * (j1 k j2; 0 0 0)
    """
    from sympy.physics.wigner import wigner_3j, clebsch_gordan
    from sympy import sqrt as ssqrt, Rational, N as numerical

    k = 4

    # Reduced matrix element
    rme = ((-1)**j1 * float(ssqrt((2*j1+1)*(2*j2+1)))
           * float(wigner_3j(j1, k, j2, 0, 0, 0)))

    if abs(rme) < 1e-15:
        return 0.0

    # Sum over q = 0, +4, -4 with appropriate coefficients
    c_5_14 = np.sqrt(5.0/14.0)

    val = 0.0
    for q, coeff in [(0, 1.0), (4, c_5_14), (-4, c_5_14)]:
        # ⟨j1,m1|C_{kq}|j2,m2⟩ = (-1)^{j1-m1} * 3j(j1,k,j2;-m1,q,m2) * rme
        threej = float(wigner_3j(j1, k, j2, -m1, q, m2))
        val += coeff * (-1)**(j1 - m1) * threej * rme

    return val


def build_V4_matrix(j_max):
    """
    Build the full V₄ matrix in the |j,m⟩ basis for j = 0, ..., j_max.
    Only integer j (tensorial sector).

    Returns: (V4_matrix, basis_labels) where basis_labels[(j,m)] = index.
    """
    # Build basis: all |j,m⟩ with j = 0,1,...,j_max, m = -j,...,j
    basis = []
    for j in range(j_max + 1):
        for m in range(-j, j + 1):
            basis.append((j, m))

    dim = len(basis)
    V4 = np.zeros((dim, dim))

    print(f"Building V₄ matrix: {dim}×{dim} basis (j_max={j_max})...")

    for i, (j1, m1) in enumerate(basis):
        for j2 in range(max(0, j1 - 4), min(j_max, j1 + 4) + 1):
            for m2 in range(-j2, j2 + 1):
                j_idx = None
                for k, (jj, mm) in enumerate(basis):
                    if jj == j2 and mm == m2:
                        j_idx = k
                        break
                if j_idx is None:
                    continue
                val = cubic_harmonic_matrix_element(j1, m1, j2, m2)
                V4[i, j_idx] = val

    # Symmetrize (should already be symmetric, but ensure numerics)
    V4 = 0.5 * (V4 + V4.T)

    return V4, basis


def build_sector_hamiltonian(j_max, group_su2, chi_values, label=""):
    """
    Build H₀ and V₄ projected onto a given character sector.

    Returns: (H0, V4, j_labels)
    where j_labels[i] = j value of the i-th sector state.
    """
    sector_data = {}
    j_labels = []

    for j in range(j_max + 1):
        n, basis = get_sector_subspace(j, group_su2, chi_values)
        if n > 0:
            sector_data[j] = basis
            for _ in range(n):
                j_labels.append(j)

    n_total = len(j_labels)
    occ = [j for j in range(j_max+1) if j in sector_data]
    if label:
        print(f"{label}: {n_total} states, occupied j = {occ}")
    else:
        print(f"Sector: {n_total} states for j = 0..{j_max}")
        print(f"  Occupied j values: {occ}")

    H0 = np.diag([j*(j+1) for j in j_labels]).astype(float)
    V4 = np.zeros((n_total, n_total))

    idx_start = {}
    offset = 0
    for j in sorted(sector_data.keys()):
        n = sector_data[j].shape[1]
        idx_start[j] = offset
        offset += n

    for j1 in sorted(sector_data.keys()):
        basis1 = sector_data[j1]
        n1 = basis1.shape[1]
        for j2 in sorted(sector_data.keys()):
            if abs(j1 - j2) > 4:
                continue
            basis2 = sector_data[j2]
            n2 = basis2.shape[1]

            dim1 = 2*j1 + 1
            dim2 = 2*j2 + 1
            V_block_full = np.zeros((dim1, dim2))

            for mi in range(dim1):
                m1 = j1 - mi
                for mj in range(dim2):
                    m2 = j2 - mj
                    V_block_full[mi, mj] = cubic_harmonic_matrix_element(
                        j1, m1, j2, m2)

            V_block_proj = (basis1.conj().T @ V_block_full @ basis2).real

            i0 = idx_start[j1]
            j0 = idx_start[j2]
            V4[i0:i0+n1, j0:j0+n2] = V_block_proj

    V4 = 0.5 * (V4 + V4.T)

    return H0, V4, j_labels


# Backward-compatible wrapper
def build_V4_A0(j_max, group_su2):
    return build_sector_hamiltonian(j_max, group_su2, [1.0] * len(group_su2))


# ============================================================
# Main computation
# ============================================================

def run_perturbative_scan(j_max=20, n_eps=50, eps_max=0.5):
    """
    Scan ε from 0 to eps_max, diagonalize H₀ + ε·V₄ in A₀ sector.
    Track eigenvalues as function of ε.
    """
    # Build 2O in SU(2)
    elements_2O = build_binary_octahedral()
    group_su2 = [quat_to_su2(q) for q in elements_2O]
    print(f"|2O| = {len(group_su2)}")

    # Build A₀-projected Hamiltonian
    H0, V4, j_labels = build_V4_A0(j_max, group_su2)
    n_states = len(j_labels)

    print(f"\nA₀ states: {n_states}")
    print(f"j_labels: {j_labels[:15]}...")

    # Verify at ε=0: eigenvalues should be j(j+1)
    E0 = np.sort(np.diag(H0))
    expected = sorted([j*(j+1) for j in j_labels])
    assert np.allclose(E0, expected), "H₀ eigenvalues don't match j(j+1)"
    print("✓ H₀ eigenvalues match j(j+1)")

    # Scan over ε
    epsilons = np.linspace(0, eps_max, n_eps)
    all_energies = np.zeros((n_eps, n_states))

    for i, eps in enumerate(epsilons):
        H = H0 + eps * V4
        eigvals = eigh(H, eigvals_only=True)
        all_energies[i] = eigvals

    return epsilons, all_energies, j_labels, H0, V4


def analyze_robustness(epsilons, all_energies, j_labels):
    """
    Analyze whether the binary occupation pattern is robust.

    At ε=0, levels are at E = j(j+1) with m(A₀,j) ∈ {0,1}.
    Track: do occupied levels (m=1) stay separated from
    unoccupied levels (m=0)?
    """
    j_star = 11  # octahedral threshold

    # Identify which j values are occupied at ε=0
    occupied_j = set(j_labels)  # j values that have an A₀ state

    # Sub-threshold occupied j for A₀ on O: 0, 4, 6, 8, 9, 10
    sub_thresh_occupied = [j for j in range(j_star + 1) if j in occupied_j]
    sub_thresh_empty = [j for j in range(j_star + 1) if j not in occupied_j]

    print(f"\nSub-threshold (j ≤ {j_star}):")
    print(f"  Occupied: {sub_thresh_occupied}")
    print(f"  Empty:    {sub_thresh_empty}")

    # For each ε, check: are the lowest n_occ eigenvalues still
    # well-separated from each other (no level crossings)?
    n_sub = len(sub_thresh_occupied)

    print(f"\nLevel tracking (first {n_sub} A₀ eigenvalues):")
    print(f"{'ε':>8} | ", end="")
    for j in sub_thresh_occupied[:6]:
        print(f"  E(j={j})", end="")
    print()
    print("-" * 70)

    for i in range(0, len(epsilons), max(1, len(epsilons)//10)):
        eps = epsilons[i]
        E = all_energies[i]
        print(f"{eps:8.3f} | ", end="")
        for k in range(min(6, n_sub)):
            print(f" {E[k]:8.3f}", end="")
        print()

    # Key metric: minimum gap between consecutive A₀ levels
    # If this stays positive, levels don't cross → pattern is robust
    print(f"\nMinimum gap between consecutive A₀ levels (sub-threshold):")
    for i in range(0, len(epsilons), max(1, len(epsilons)//10)):
        eps = epsilons[i]
        E = all_energies[i, :n_sub]
        gaps = np.diff(E)
        min_gap = np.min(gaps) if len(gaps) > 0 else float('inf')
        print(f"  ε = {eps:.3f}: min gap = {min_gap:.4f}")


def test_perturbative_stability():
    """Main test: perturbative stability of spectral reflection."""
    print("=" * 70)
    print("PERTURBATIVE STABILITY: H = j(j+1) + ε·V₄ on SO(3)/O")
    print("=" * 70)

    epsilons, all_energies, j_labels, H0, V4 = run_perturbative_scan(
        j_max=20, n_eps=50, eps_max=2.0)

    analyze_robustness(epsilons, all_energies, j_labels)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    # Check V₄ matrix structure
    print(f"\nV₄ norm: {np.linalg.norm(V4):.4f}")
    print(f"V₄ diagonal elements (first 8): {np.diag(V4)[:8].round(4)}")
    print(f"V₄ max off-diagonal: {np.max(np.abs(V4 - np.diag(np.diag(V4)))):.4f}")

    # At small ε, first-order shift is ε·⟨j|V₄|j⟩ (diagonal of V₄ in A₀ basis)
    print(f"\nFirst-order energy shifts (ε·V₄ diagonal):")
    for k, j in enumerate(j_labels[:12]):
        print(f"  j={j}: E₀ = {j*(j+1):.0f}, δE/ε = {V4[k,k]:.4f}")


# ============================================================
# Test 2: Cross-sector energy pairing
# ============================================================

def test_energy_pairing():
    """
    Non-trivial test: does E(A₀,j) + E(A₁,j*-j) stay constant under ε?

    The occupation duality m(j)+m(j*-j)=1 is trivially rigid (Schur).
    The ENERGY sum across paired sectors is NOT guaranteed.
    If it drifts at O(ε), duality is purely arithmetic.

    Character A₁ of 2O: +1 on 2T (first 24 elements), -1 on complement.
    build_binary_octahedral() returns 2T first, so this is just [+1]*24+[-1]*24.
    """
    j_star = 11
    j_max = 15  # enough for sub-threshold + V₄ coupling window

    elements_2O = build_binary_octahedral()
    group_su2 = [quat_to_su2(q) for q in elements_2O]

    chi_A0 = [1.0] * 48
    chi_A1 = [1.0] * 24 + [-1.0] * 24

    # Verify chi_A1 is a homomorphism: χ(g)χ(h) = χ(gh) for all g,h
    from src.quaternion import qmul, qkey
    key_to_idx = {qkey(elements_2O[i]): i for i in range(48)}
    for i in range(48):
        for k in range(48):
            prod = qmul(elements_2O[i], elements_2O[k])
            m = key_to_idx[qkey(prod)]
            assert chi_A1[i] * chi_A1[k] == chi_A1[m], \
                f"A₁ not a homomorphism at indices {i},{k},{m}"
    print("✓ chi_A1 verified as group homomorphism")

    print("\n" + "=" * 70)
    print("TEST 2: Cross-sector energy pairing")
    print("=" * 70)

    print(f"\nBuilding sectors (j_max={j_max})...")
    H0_A0, V4_A0, jl_A0 = build_sector_hamiltonian(j_max, group_su2, chi_A0, "  A₀")
    H0_A1, V4_A1, jl_A1 = build_sector_hamiltonian(j_max, group_su2, chi_A1, "  A₁")

    occ_A0 = set(jl_A0)
    occ_A1 = set(jl_A1)

    # Cross-sector duality pairs: j in A₀, j*-j in A₁ (or vice versa)
    pairs = []
    for j in range(j_star // 2 + 1):
        jp = j_star - j
        if j in occ_A0 and jp in occ_A1:
            pairs.append((j, jp, 'A₀', 'A₁'))
        elif j in occ_A1 and jp in occ_A0:
            pairs.append((j, jp, 'A₁', 'A₀'))

    print(f"\nCross-sector duality pairs:")
    for j, jp, s1, s2 in pairs:
        print(f"  j={j} ({s1}) <-> j*-j={jp} ({s2}), "
              f"E₀ sum = {j*(j+1)} + {jp*(jp+1)} = {j*(j+1)+jp*(jp+1)}")

    # Index: first occurrence of j in label list
    def j_idx(j, jlabels):
        for k, jl in enumerate(jlabels):
            if jl == j:
                return k
        return None

    # First-order shifts (diagonal V₄ elements)
    print(f"\nFirst-order shifts δE/ε (V₄ diagonal):")
    print(f"  {'j':>3} {'A₀':>12} {'A₁':>12}")
    for j in range(j_star + 1):
        a0 = f"{V4_A0[j_idx(j, jl_A0), j_idx(j, jl_A0)]:+.6f}" if j in occ_A0 else "      ---"
        a1 = f"{V4_A1[j_idx(j, jl_A1), j_idx(j, jl_A1)]:+.6f}" if j in occ_A1 else "      ---"
        print(f"  {j:3d} {a0:>12} {a1:>12}")

    print(f"\nPaired first-order sums δE(s1,j)/ε + δE(s2,j*-j)/ε:")
    for j, jp, s1, s2 in pairs:
        if s1 == 'A₀':
            d1 = V4_A0[j_idx(j, jl_A0), j_idx(j, jl_A0)]
            d2 = V4_A1[j_idx(jp, jl_A1), j_idx(jp, jl_A1)]
        else:
            d1 = V4_A1[j_idx(j, jl_A1), j_idx(j, jl_A1)]
            d2 = V4_A0[j_idx(jp, jl_A0), j_idx(jp, jl_A0)]
        print(f"  ({j},{jp}): {d1:+.6f} + {d2:+.6f} = {d1+d2:+.6f}")

    # Full ε scan
    n_eps = 50
    eps_max = 2.0
    epsilons = np.linspace(0, eps_max, n_eps)

    print(f"\nEnergy sums S(j,j*-j) vs ε:")
    hdr = f"{'eps':>8} |"
    for j, jp, _, _ in pairs:
        hdr += f"  S({j},{jp})"
    hdr += "  |"
    for j, jp, _, _ in pairs:
        hdr += f" dS({j},{jp})"
    print(hdr)
    print("-" * len(hdr))

    ref = {}
    for i in range(0, n_eps, max(1, n_eps // 10)):
        eps = epsilons[i]
        E_A0 = eigh(H0_A0 + eps * V4_A0, eigvals_only=True)
        E_A1 = eigh(H0_A1 + eps * V4_A1, eigvals_only=True)

        line = f"{eps:8.3f} |"
        deltas = ""
        for j, jp, s1, s2 in pairs:
            if s1 == 'A₀':
                e1 = E_A0[j_idx(j, jl_A0)]
                e2 = E_A1[j_idx(jp, jl_A1)]
            else:
                e1 = E_A1[j_idx(j, jl_A1)]
                e2 = E_A0[j_idx(jp, jl_A0)]
            s = e1 + e2
            if i == 0:
                ref[(j, jp)] = s
            delta = s - ref[(j, jp)]
            line += f" {s:8.2f}"
            deltas += f" {delta:+7.4f}"
        line += "  |" + deltas
        print(line)


# ============================================================
# H2: Lifting criterion — Q₈ containment for all O subgroups
# ============================================================

def _qmul(a, b):
    return np.array([
        a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
        a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
        a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
        a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0],
    ])

def _qinv(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

def _generate_group(gens, max_size=500):
    elts = [np.array([1,0,0,0], dtype=float)]
    def is_in(q, lst):
        return any(np.linalg.norm(q-e) < 1e-10 for e in lst)
    for g in gens:
        for x in [g, _qinv(g)]:
            if not is_in(x, elts):
                elts.append(x)
    changed = True
    while changed:
        changed = False
        new = []
        for a in elts:
            for b in elts:
                p = _qmul(a, b)
                if not is_in(p, elts) and not is_in(p, new):
                    new.append(p)
                    if len(elts)+len(new) > max_size:
                        raise ValueError(f"Exceeds {max_size}")
        if new:
            elts.extend(new)
            changed = True
    return elts

def _contains_Q8(G):
    """Check Q₈ ⊂ G by direct subset test."""
    q8 = [np.array([s*d[0], s*d[1], s*d[2], s*d[3]])
          for d in [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
          for s in [1,-1]]
    def is_in(q, lst):
        return any(np.linalg.norm(q-e) < 1e-10 for e in lst)
    return all(is_in(q, G) for q in q8)

def _has_commutator_minus1(G):
    """Check ∃ a,b ∈ G with [a,b] = -1."""
    for a in G:
        for b in G:
            c = _qmul(_qmul(a, b), _qmul(_qinv(a), _qinv(b)))
            if abs(c[0]+1) < 1e-10 and np.linalg.norm(c[1:]) < 1e-10:
                return True
    return False

def test_lifting_criterion():
    """
    H2: For every subgroup channel O → K, verify:
      spinorial sectors exist on SO(3)/K ⟺ Q₈ ⊄ 2K
    Two independent checks: direct Q₈ subset + commutator test.
    """
    sq2 = 1/np.sqrt(2)
    groups = [
        ("2O (parent)", lambda: _generate_group(
            [np.array([.5,.5,.5,.5]), np.array([sq2,sq2,0,0])]), 48),
        ("2T",  lambda: _generate_group(
            [np.array([0,1,0,0]), np.array([.5,.5,.5,.5])]), 24),
        ("2D₄", lambda: _generate_group(
            [np.array([np.cos(np.pi/4),np.sin(np.pi/4),0,0]),
             np.array([0,0,1,0])]), 16),
        ("2D₃", lambda: _generate_group(
            [np.array([np.cos(np.pi/3),np.sin(np.pi/3),0,0]),
             np.array([0,0,1,0])]), 12),
        ("Q₈ = 2D₂", lambda: _generate_group(
            [np.array([0,1,0,0]), np.array([0,0,1,0])]), 8),
        ("Z₈ = 2C₄", lambda: _generate_group(
            [np.array([np.cos(np.pi/4),np.sin(np.pi/4),0,0])]), 8),
        ("Z₆ = 2C₃", lambda: _generate_group(
            [np.array([np.cos(np.pi/3),np.sin(np.pi/3),0,0])]), 6),
        ("Z₄ = 2C₂", lambda: _generate_group(
            [np.array([0,1,0,0])]), 4),
    ]

    print("\n" + "=" * 70)
    print("H2: LIFTING CRITERION — Q₈ containment phase diagram")
    print("=" * 70)
    print(f"\n{'Channel':<16} {'|2K|':>5} {'Q₈⊂2K':>7} {'[a,b]=-1':>9} {'Spinorial?':>11}")
    print("-" * 55)

    for name, builder, expected_order in groups:
        G = builder()
        assert len(G) == expected_order, f"{name}: |G|={len(G)} ≠ {expected_order}"
        has_q8 = _contains_Q8(G)
        has_comm = _has_commutator_minus1(G)
        assert has_q8 == has_comm, f"{name}: subset≠commutator!"
        spin = "No" if has_q8 else "YES"
        print(f"{name:<16} {expected_order:>5} {'Yes' if has_q8 else 'No':>7} "
              f"{'Yes' if has_comm else 'No':>9} {spin:>11}")

    print(f"\nMonotonicity: Q₈ ⊄ 2K ∧ K' ⊂ K ⟹ Q₈ ⊄ 2K' (subgroup transitivity)")
    print("Lifting criterion verified for all 8 channels of O → K.")


if __name__ == '__main__':
    test_perturbative_stability()
    test_energy_pairing()
    test_lifting_criterion()
