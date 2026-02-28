#!/usr/bin/env python3
"""
I7 — Thermodynamic consequence of spectral reflection.

INVESTIGATION: The duality m(χ,j) + m(χ,j*-j) = 1 for j ≤ j*
implies a constraint on the sector partition function:
  Z_χ(T) + Z_χ^refl(T) = Z_trunc(T)
where Z_trunc sums over j ≤ j* (all integers, matching parity of χ).

This is a direct restatement of m + m' = 1 at the partition function
level. Verified numerically but not a T-duality (Z_refl mixes
multiplicities at j' = j*-j with Boltzmann factors at j).

CONCLUSION: correct but not transformative — the identity adds
no content beyond the spectral duality itself. Could appear as
a remark in the paper.

RAW OUTPUT:
  Reflection identity Z_A₀ + Z_A₀^refl = Z_trunc: ✓ (5 temperatures)
  Reflection identity Z_A₁ + Z_A₁^refl = Z_trunc: ✓ (5 temperatures)
  At T=10B:  Z_A₀ = 2.4283,  Z_A₀^refl = 7.9118
    Sum = 10.3401 = Z_trunc → not a simple T-duality
  Low-T:  Z_A₀ → 1
  High-T: Z_A₀/Z_trunc = 80/144 = 0.5556 (≠ 0.5, reflection ≠ (2j+1))
  5 passed in 0.37s

STATUS: COMPLETE — verified, mild interest
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest


# ==========================================================
# § 0. Infrastructure: compute m(χ, j) for SO(3)/H
# ==========================================================

def chi_su2(j, g_trace):
    """Character of V_j at element with tr(g) = g_trace.
    χ_j(g) = sin((2j+1)α/2) / sin(α/2) where tr(g) = 2cos(α/2)."""
    if abs(g_trace - 2.0) < 1e-10:
        return 2*j + 1
    if abs(g_trace + 2.0) < 1e-10:
        return (2*j + 1) * ((-1)**(int(2*j)))
    alpha_half = np.arccos(g_trace / 2.0)
    return np.sin((2*j + 1) * alpha_half) / np.sin(alpha_half)


def multiplicity(j, chi_sector, group_elements, group_order):
    """m(χ, j) = (1/|G|) Σ_g conj(χ(g)) · χ_j(g).
    chi_sector: dict mapping element_key → character value.
    group_elements: list of (key, trace) pairs."""
    s = 0.0
    for key, tr in group_elements:
        s += np.conj(chi_sector[key]) * chi_su2(j, tr)
    return s.real / group_order


# ==========================================================
# § 1. Build 2O and its scalar characters
# ==========================================================

def build_2O_data():
    """Build binary octahedral group data for partition function computation.
    Returns elements, characters, and group parameters."""
    from src.group import build_binary_octahedral
    from src.quaternion import qkey

    G = build_binary_octahedral()
    order = len(G)  # 48
    j_star = order // 4 - 1  # |H|/2 - 1 = 24/2 - 1 = 11

    # Elements as (key, trace) pairs
    elements = []
    for g in G:
        key = qkey(g)
        tr = 2 * g[0]  # trace of SU(2) matrix = 2·Re(q)
        elements.append((key, tr))

    # Scalar characters of 2O: A₀ (trivial) and A₁ (sign)
    # A₀: χ(g) = 1 for all g
    # A₁: χ(g) = det(g as SO(3)) lifted: +1 for even, -1 for odd rotations
    #     Equivalently: A₁(g) = +1 if g ∈ 2T, -1 otherwise
    from src.group import compute_commutator_subgroup
    comm = compute_commutator_subgroup(G)
    comm_keys = {qkey(g) for g in comm}

    chi_A0 = {qkey(g): 1.0 for g in G}
    # A₁: trivial on [G,G] = 2T, sign on coset
    chi_A1 = {}
    for g in G:
        k = qkey(g)
        chi_A1[k] = 1.0 if k in comm_keys else -1.0

    return {
        'elements': elements,
        'order': order,
        'j_star': j_star,
        'chi_A0': chi_A0,
        'chi_A1': chi_A1,
    }


# ==========================================================
# § 2. Test the reflection identity
#
# For sector χ: m(χ, j) + m(χ, j*-j) = 1 for all integer j ≤ j*.
#
# Define Z_refl(T) = Σ_{j=0}^{j*} m(χ, j*-j) · (2j+1) · exp(-j(j+1)/T)
#
# Then: Z_χ(T) + Z_refl(T) = Σ_{j=0}^{j*} (2j+1) · exp(-j(j+1)/T)
#                            = Z_trunc(T)
#
# Z_refl is NOT Z_χ at a different temperature — it mixes multiplicities
# at j' = j*-j with Boltzmann factors at j. Not a T-duality.
# ==========================================================

class TestReflectionIdentity:
    """Test Z_χ + Z_χ^refl = Z_trunc for sub-threshold levels."""

    def test_reflection_identity_A0(self):
        """For A₀ of 2O: Z_A₀(T) + Z_A₀^refl(T) = Z_trunc(T)
        at several temperatures. Sum over ALL j from 0 to j*."""
        data = build_2O_data()
        j_star = data['j_star']  # 11

        for T in [1.0, 5.0, 10.0, 50.0, 100.0]:
            Z_chi = 0.0
            Z_refl = 0.0
            Z_trunc = 0.0
            for j in range(0, j_star + 1):
                m_j = multiplicity(j, data['chi_A0'], data['elements'],
                                   data['order'])
                m_dual = multiplicity(j_star - j, data['chi_A0'],
                                      data['elements'], data['order'])
                boltz = (2*j + 1) * np.exp(-j*(j+1) / T)
                Z_chi += m_j * boltz
                Z_refl += m_dual * boltz
                Z_trunc += boltz

            assert abs(Z_chi + Z_refl - Z_trunc) < 1e-10, \
                f"Identity fails at T={T}: {Z_chi} + {Z_refl} ≠ {Z_trunc}"

        print(f"\nReflection identity Z_A₀ + Z_A₀^refl = Z_trunc: ✓ (5 temperatures)")

    def test_reflection_identity_A1(self):
        """Same identity for A₁ of 2O. Sum over ALL j from 0 to j*."""
        data = build_2O_data()
        j_star = data['j_star']  # 11

        for T in [1.0, 5.0, 10.0, 50.0, 100.0]:
            Z_chi = 0.0
            Z_refl = 0.0
            Z_trunc = 0.0
            for j in range(0, j_star + 1):
                m_j = multiplicity(j, data['chi_A1'], data['elements'],
                                   data['order'])
                m_dual = multiplicity(j_star - j, data['chi_A1'],
                                      data['elements'], data['order'])
                boltz = (2*j + 1) * np.exp(-j*(j+1) / T)
                Z_chi += m_j * boltz
                Z_refl += m_dual * boltz
                Z_trunc += boltz

            assert abs(Z_chi + Z_refl - Z_trunc) < 1e-10, \
                f"Identity fails at T={T}"

        print(f"Reflection identity Z_A₁ + Z_A₁^refl = Z_trunc: ✓ (5 temperatures)")

    def test_identity_content(self):
        """The identity is Z_χ + Z_χ^refl = Z_trunc.
        Z_refl has DIFFERENT Boltzmann weights than Z_χ
        (m(χ,j*-j) weighted by exp(-j(j+1)/T), not exp(-(j*-j)(j*-j+1)/T)).
        So this is NOT a simple duality Z(T) = Z(T') — it's a constraint
        that mixes multiplicities at j with Boltzmann factors at j."""
        data = build_2O_data()
        j_star = data['j_star']
        T = 10.0

        # Compute both Z and Z^refl to show they're different
        Z_chi = 0.0
        Z_refl = 0.0
        for j in range(0, j_star + 1):
            m_j = multiplicity(j, data['chi_A0'], data['elements'],
                               data['order'])
            m_dual = multiplicity(j_star - j, data['chi_A0'],
                                  data['elements'], data['order'])
            boltz = (2*j + 1) * np.exp(-j*(j+1) / T)
            Z_chi += m_j * boltz
            Z_refl += m_dual * boltz

        # They should be different (not symmetric)
        assert abs(Z_chi - Z_refl) > 0.01, \
            "Z and Z_refl should differ (different j weights)"
        print(f"\nAt T=10B:  Z_A₀ = {Z_chi:.4f},  Z_A₀^refl = {Z_refl:.4f}")
        print(f"  Sum = {Z_chi + Z_refl:.4f} = Z_trunc")
        print(f"  Z and Z^refl are different → not a simple T-duality")


# ==========================================================
# § 4. Limiting behaviour
#
# Low T:  Z_A₀ → 1 (only j=0 contributes, m(A₀,0)=1).
# High T: Z_A₀/Z_trunc → Σ_{j:m=1} (2j+1) / (j*+1)².
#   NOT 1/2 because reflection j ↔ j*-j does not preserve (2j+1).
#   (It preserves the sum: (2j+1) + (2(j*-j)+1) = 2j*+2 = |H|.)
#
# Remark: Cor 3.5 (regular-rep tiling) gives a stronger identity
# summing over ALL irreps ρ (not just scalar):
#   Σ_ρ (dim ρ) · Z_ρ(T) + Σ_ρ (dim ρ) · Z_ρ^refl(T) = |H| · Z_trunc(T)
# This constrains the full partition function (including non-scalar
# sectors) but is not tested here.
# ==========================================================

class TestThermodynamicConsequences:
    """Physical implications of the reflection identity."""

    def test_low_T_limit(self):
        """At T→0, Z_A₀ → 1 (only j=0 contributes, m(A₀,0)=1).
        So Z_refl → Z_trunc - 1."""
        data = build_2O_data()
        j_star = data['j_star']
        T = 0.01  # very low

        Z_chi = 0.0
        Z_trunc = 0.0
        for j in range(0, j_star + 1):
            m_j = multiplicity(j, data['chi_A0'], data['elements'],
                               data['order'])
            boltz = (2*j + 1) * np.exp(-j*(j+1) / T)
            Z_chi += m_j * boltz
            Z_trunc += boltz

        assert abs(Z_chi - 1.0) < 1e-10, "Z_A₀ → 1 at low T"
        print(f"\nLow-T limit (T=0.01B):")
        print(f"  Z_A₀ = {Z_chi:.10f} → 1")

    def test_high_T_ratio(self):
        """At T→∞, exp(-j(j+1)/T) → 1, so:
          Z_A₀ → Σ_{j: m=1} (2j+1),  Z_trunc → Σ_j (2j+1) = (j*+1)²
        Ratio → Σ_{j: m=1} (2j+1) / (j*+1)².
        NOT 1/2 because occupied levels (j=0,4,6,8,9,10 for A₀ of O)
        have different (2j+1) than unoccupied ones."""
        data = build_2O_data()
        j_star = data['j_star']

        # Exact high-T ratio from (2j+1) weights
        num = 0.0   # Σ_{j: m=1} (2j+1)
        den = 0.0   # Σ_j (2j+1) = (j*+1)²
        for j in range(0, j_star + 1):
            m_j = multiplicity(j, data['chi_A0'], data['elements'],
                               data['order'])
            deg = 2*j + 1
            den += deg
            if round(m_j) == 1:
                num += deg
        exact_ratio = num / den
        assert abs(den - (j_star + 1)**2) < 1e-10  # Σ(2j+1) = (j*+1)²

        # Verify at very high T
        T = 1e6
        Z_chi = 0.0
        Z_trunc = 0.0
        for j in range(0, j_star + 1):
            m_j = multiplicity(j, data['chi_A0'], data['elements'],
                               data['order'])
            boltz = (2*j + 1) * np.exp(-j*(j+1) / T)
            Z_chi += m_j * boltz
            Z_trunc += boltz

        ratio = Z_chi / Z_trunc
        assert abs(ratio - exact_ratio) < 1e-4, \
            f"High-T ratio {ratio} ≠ exact {exact_ratio}"
        print(f"\nHigh-T limit:")
        print(f"  Z_A₀/Z_trunc = {ratio:.6f}")
        print(f"  Exact: Σ(2j+1)[m=1] / (j*+1)² = {num:.0f}/{den:.0f}"
              f" = {exact_ratio:.6f}")
        print(f"  ≠ 0.5 because reflection doesn't preserve (2j+1)")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
