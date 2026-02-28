#!/usr/bin/env python3
"""
I31 — Quadratic rigidity on S³/Γ.

QUESTION: Does the S² reflection identity lift to S³/Γ, and does the
quadratic structure of m_S³ discriminate condition (B)?

RAW OUTPUT:
  §0: CG cross-check passes for 2T, 2O, 2I (ℓ=0..2|H|):
      cumulative Σ m_S²(k) ≡ (1/|G|) Σ_g [χ_{ℓ/2}(g)]².
  §1: Δ²(ℓ) = m_S³(ℓ+2P)−2m_S³(ℓ+P)+m_S³(ℓ):
      (B) true  → Δ² = P constant (2T:6, 2O:12, 2I:30, D₄:4, D₆:6).
      (B) false → Δ² oscillates (D₃:{2,4}, D₅:{4,6}).
  §2: Exact integer identity m_S³(nP−1) = Pn²/2 verified for n=1..5:
      2T: m_S³(5,11,17,23,29) = (3,12,27,48,75) = 6·n²/2.
      2O: m_S³(11,23,...) = 12·n²/2. 2I: 30·n²/2.
  §3: D_n odd oscillation: Δ² takes exactly 2 values = {P−1, P+1}.
      D₃:{2,4}, D₅:{4,6}, D₇:{6,8}.

ANSWER: POSITIVE. The S² staircase lifts to S³/Γ as quadratic rigidity:
  m_S³(nP−1) = Pn²/2  (exact integer identity at period boundaries).
  The constant second difference Δ² = P is equivalent to condition (B).
  Without (B), Δ² oscillates between {P−1, P+1}.
  Reflection itself does NOT lift (cumulative sums are not constant-sum),
  but the periodicity does, as quadratic polynomial growth.

PLAN:
  §0. CG cross-check: cumulative sum = squared-character formula [3 tests]
  §1. Constant second difference ⟺ (B) [7 tests]
  §2. Exact quadratic formula at period boundaries [3 tests]
  §3. Non-(B) oscillation pattern [3 tests]

STATUS: COMPLETE — 16/16 tests pass
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
                       build_binary_icosahedral, compute_conjugacy_classes)
from src.quaternion import qmul


# ---- helpers ----

def build_binary_dihedral(n):
    """2D_n = Dic_n in SU(2). Order 4n."""
    j_quat = np.array([0.0, 0.0, 1.0, 0.0])
    result = []
    for k in range(2 * n):
        theta = k * np.pi / n
        q = np.array([np.cos(theta), 0, 0, np.sin(theta)])
        result.append(q)
        result.append(qmul(j_quat, q))
    G = np.array(result)
    x2 = qmul(j_quat, j_quat)
    a_n = np.array([np.cos(np.pi), 0, 0, np.sin(np.pi)])
    assert np.allclose(x2, a_n, atol=1e-12), "Dicyclic relation x²=aⁿ failed"
    return G


def chi_su2(j, trace):
    """SU(2) spin-j character at given trace = 2cos(α)."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def prepare_group(G):
    """Prepare conjugacy class data from element array."""
    classes = compute_conjugacy_classes(G)
    classes.sort(key=lambda c: -c['trace'])
    return classes, len(G)


def m_S2(j, classes, order):
    """m(A₀, j) on S² — multiplicity of trivial character in V_j."""
    s = sum(chi_su2(j, c['trace']) * c['size'] for c in classes)
    val = (s / order).real
    assert abs(val - round(val)) < 1e-6, f"j={j}: non-integer {val}"
    return int(round(val))


def m_S3_cumulative(ell, classes, order):
    """m_S³(ℓ) via CG: Σ_{k=0}^{ℓ} m_S²(A₀, k).
    V_{ℓ/2} ⊗ V_{ℓ/2} = V₀ ⊕ V₁ ⊕ ... ⊕ V_ℓ (integer spins)."""
    return sum(m_S2(k, classes, order) for k in range(ell + 1))


def m_S3_direct(ell, classes, order):
    """m_S³(ℓ) via squared character: (1/|G|) Σ_g [χ_{ℓ/2}(g)]²."""
    s = sum(chi_su2(ell / 2, c['trace']) ** 2 * c['size'] for c in classes)
    val = (s / order).real
    assert abs(val - round(val)) < 1e-6, f"ell={ell}: non-integer {val}"
    return int(round(val))


def s3_mults(classes, order, ell_max):
    """Compute m_S³(ℓ) for ℓ=0..ell_max via cumulative sum."""
    s2 = [m_S2(k, classes, order) for k in range(ell_max + 1)]
    s3 = []
    running = 0
    for k in range(ell_max + 1):
        running += s2[k]
        s3.append(running)
    return s3


# ==========================================================
# § 0. CG cross-check
#
# Two independent computations of m_S³(ℓ):
#   Method A: Σ_{k=0}^{ℓ} m_S²(k)          (Clebsch-Gordan)
#   Method B: (1/|G|) Σ_g [χ_{ℓ/2}(g)]²    (squared character)
# Must agree exactly (both integer-valued).
# ==========================================================

class TestCGCrossCheck:
    """Clebsch-Gordan cumulative sum = squared-character formula."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_polyhedral(self, name, builder):
        """Cumulative S² mults = squared-character S³ mults."""
        G = builder()
        classes, order = prepare_group(G)
        H_order = order // 2
        for ell in range(2 * H_order):
            cum = m_S3_cumulative(ell, classes, order)
            direct = m_S3_direct(ell, classes, order)
            assert cum == direct, \
                f"{name} ell={ell}: cumulative={cum}, direct={direct}"


# ==========================================================
# § 1. Constant second difference ⟺ (B)
#
# Define Δ²(ℓ) = m_S³(ℓ+2P) − 2·m_S³(ℓ+P) + m_S³(ℓ).
#
# Claim: Δ²(ℓ) = P (constant) for all ℓ  ⟺  condition (B).
#
# Proof sketch (⇐): Under (B), m_S²(j+P) = m_S²(j)+1 (staircase).
# Cumulative sum telescopes: Δ² = P.
# Without (B), the staircase fails: increment per period is
# non-uniform, so Δ² oscillates.
# ==========================================================

class TestConstantSecondDiff:
    """Δ²(ℓ) = P iff (B)."""

    @pytest.mark.parametrize("name,builder,H_order", [
        ("2T", build_binary_tetrahedral, 12),
        ("2O", build_binary_octahedral, 24),
        ("2I", build_binary_icosahedral, 60),
    ])
    def test_polyhedral_constant(self, name, builder, H_order):
        """Polyhedral (B=True): Δ² = P for all ℓ in [0, 2P]."""
        G = builder()
        classes, order = prepare_group(G)
        P = H_order // 2
        s3 = s3_mults(classes, order, 4 * P)
        for ell in range(2 * P + 1):
            d2 = s3[ell + 2 * P] - 2 * s3[ell + P] + s3[ell]
            assert d2 == P, \
                f"{name} ell={ell}: Δ²={d2}, expected P={P}"

    @pytest.mark.parametrize("n", [4, 6])
    def test_even_dihedral_constant(self, n):
        """D_n even (B=True): Δ² = n for all ℓ."""
        G = build_binary_dihedral(n)
        classes, order = prepare_group(G)
        P = n  # |H|/2 = 2n/2 = n
        s3 = s3_mults(classes, order, 4 * P)
        for ell in range(2 * P + 1):
            d2 = s3[ell + 2 * P] - 2 * s3[ell + P] + s3[ell]
            assert d2 == P, \
                f"D_{n} ell={ell}: Δ²={d2}, expected P={P}"

    @pytest.mark.parametrize("n", [3, 5])
    def test_odd_dihedral_not_constant(self, n):
        """D_n odd (B=False): Δ² is NOT constant."""
        G = build_binary_dihedral(n)
        classes, order = prepare_group(G)
        P = n
        s3 = s3_mults(classes, order, 4 * P)
        d2_vals = [s3[ell + 2 * P] - 2 * s3[ell + P] + s3[ell]
                   for ell in range(2 * P + 1)]
        assert len(set(d2_vals)) > 1, \
            f"D_{n}: Δ² is constant={d2_vals[0]} but (B) fails"


# ==========================================================
# § 2. Exact quadratic formula at period boundaries
#
# Under (B), m_S²(k) = ⌊k/P⌋ + ε(k mod P) with Σε = P/2.
# Summing over k = 0..nP-1:
#   m_S³(nP-1) = Σ_{k=0}^{nP-1} ⌊k/P⌋ + n·(P/2)
#              = P·n(n-1)/2 + nP/2
#              = Pn²/2.
# This is EXACT (integer identity), not asymptotic.
# ==========================================================

class TestExactQuadratic:
    """m_S³(nP-1) = Pn²/2 at period boundaries."""

    @pytest.mark.parametrize("name,builder,H_order", [
        ("2T", build_binary_tetrahedral, 12),
        ("2O", build_binary_octahedral, 24),
        ("2I", build_binary_icosahedral, 60),
    ])
    def test_period_boundary(self, name, builder, H_order):
        """m_S³(nP-1) = Pn²/2 for n = 1..5."""
        G = builder()
        classes, order = prepare_group(G)
        P = H_order // 2
        s3 = s3_mults(classes, order, 5 * P)
        for n in range(1, 6):
            ell = n * P - 1
            expected = P * n * n // 2
            assert s3[ell] == expected, \
                f"{name} n={n}: m_S3({ell})={s3[ell]}, expected Pn²/2={expected}"


# ==========================================================
# § 3. Non-(B) oscillation structure
#
# For D_n odd (B fails), Δ² oscillates between two values.
# The S² staircase fails: m_S²(j+P) - m_S²(j) alternates
# between 0 and 2 (because some modes have period n, others 2n).
# The second difference Δ² inherits a 2-valued oscillation
# with values {P-1, P+1}.
# ==========================================================

class TestNonBOscillation:
    """D_n odd: Δ² oscillates between exactly two values."""

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_two_valued(self, n):
        """D_n odd: Δ² takes exactly two distinct values."""
        G = build_binary_dihedral(n)
        classes, order = prepare_group(G)
        P = n
        s3 = s3_mults(classes, order, 4 * P)
        d2_vals = set(s3[ell + 2 * P] - 2 * s3[ell + P] + s3[ell]
                      for ell in range(2 * P + 1))
        assert len(d2_vals) == 2, \
            f"D_{n}: expected 2 Δ² values, got {len(d2_vals)}: {sorted(d2_vals)}"
        # Values should be P-1 and P+1
        assert d2_vals == {P - 1, P + 1}, \
            f"D_{n}: Δ² values = {sorted(d2_vals)}, expected {{{P-1}, {P+1}}}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
