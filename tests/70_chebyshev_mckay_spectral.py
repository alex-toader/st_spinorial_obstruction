#!/usr/bin/env python3
"""
I31/R6 — Chebyshev–McKay spectral interpretation of reflection.

QUESTION: Is reflection a Chebyshev polynomial identity on the McKay spectrum?

RAW OUTPUT:
  §0: Factorization χ_j + χ_{j*-j} = 2·T_{|2j-j*|}·U_{P-1} verified
      for 2T, 2O, 2I at all j ≤ j* and all non-trivial angles.
  §1: U_{P-1}(cos α) = 0 at all non-central classes ⟺ (B):
      2T(✓), 2O(✓), 2I(✓), D₄(✓), D₆(✓) — all zeros.
      D₃(✗), D₅(✗), D₇(✗) — U_{P-1} ≠ 0 at half-turn α=π/2.
  §2: Selective duality: U_{2k+1}(0) = 0 for all k ≥ 0.
      Half-integer modes vanish at half-turns → spinorial sector immune
      to (B)-violating classes. Verified D₃, D₅, D₇.
  §3: μ(ℓ) = m_S³(ℓ) − ℓ is palindromic under ℓ ↔ 2j*−ℓ.
      Universal (holds for ALL groups, regardless of (B)).
      Verified 2T, 2O, 2I, D₃–D₇.

ANSWER: POSITIVE. Reflection factors as mode × spectral:
  χ_j(g) + χ_{j*-j}(g) = 2·T_{|2j-j*|}(cos α) · U_{P-1}(cos α).
  The spectral factor U_{P-1} vanishes at all non-central group angles
  iff (B) holds. This is a Chebyshev vanishing condition on the McKay
  eigenvalues. Selective duality arises because half-integer modes kill
  the violating factor independently.

PLAN:
  §0. Chebyshev factorization identity [3 tests]
  §1. U_{P-1} vanishing ⟺ (B) [5+3 tests]
  §2. Selective duality via mode vanishing [3 tests]
  §3. S³ palindrome μ(ℓ) = m_S³(ℓ)−ℓ [5+3 tests]

STATUS: COMPLETE — 22/22 tests pass
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
    return np.array(result)


def chi_su2(j, trace):
    """SU(2) spin-j character at trace = 2cos(α)."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def chebyshev_U(n, x):
    """Chebyshev U_n(x) = sin((n+1)·arccos(x))/sin(arccos(x))."""
    if abs(x - 1) < 1e-14:
        return float(n + 1)
    if abs(x + 1) < 1e-14:
        return float((n + 1) * (-1) ** n)
    a = np.arccos(np.clip(x, -1, 1))
    return np.sin((n + 1) * a) / np.sin(a)


def chebyshev_T(n, x):
    """Chebyshev T_n(x) = cos(n·arccos(x))."""
    if n == 0:
        return 1.0
    a = np.arccos(np.clip(x, -1, 1))
    return np.cos(n * a)


def prepare_group(G):
    """Return (classes, order, P, jstar)."""
    classes = compute_conjugacy_classes(G)
    classes.sort(key=lambda c: -c['trace'])
    order = len(G)
    P = order // 4
    return classes, order, P, P - 1


def m_S2(j, classes, order):
    """m(A₀, j) on S²."""
    s = sum(chi_su2(j, c['trace']) * c['size'] for c in classes)
    return int(round(s / order))


def m_S3_seq(classes, order, ell_max):
    """m_S³(ℓ) = Σ_{k=0}^{ℓ} m_S²(k) for ℓ = 0..ell_max."""
    result = []
    running = 0
    for k in range(ell_max + 1):
        running += m_S2(k, classes, order)
        result.append(running)
    return result


# ==========================================================
# § 0. Chebyshev factorization
#
# Identity: χ_j(g) + χ_{j*-j}(g) = 2·T_{|2j-j*|}(cos α)·U_{P-1}(cos α)
#
# Proof: sin A + sin B = 2cos((A-B)/2)sin((A+B)/2) with
#   A = (2j+1)α, B = (2(j*-j)+1)α.
#   A+B = 2Pα, (A-B)/2 = (2j-j*)α.
# Divide by sin(α): U_{2j} + U_{2(j*-j)} = 2·T_{|2j-j*|}·U_{P-1}.
# ==========================================================

class TestChebyshevFactorization:
    """χ_j + χ_{j*-j} = 2·T·U_{P-1} at all angles."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_polyhedral(self, name, builder):
        """Factorization holds at all group element angles."""
        G = builder()
        classes, order, P, jstar = prepare_group(G)
        for j in range(jstar + 1):
            m = abs(2 * j - jstar)
            for c in classes:
                x = c['trace'] / 2  # cos(alpha)
                direct = chi_su2(j, c['trace']) + chi_su2(jstar - j, c['trace'])
                factored = 2 * chebyshev_T(m, x) * chebyshev_U(P - 1, x)
                assert abs(direct - factored) < 1e-8, \
                    f"{name} j={j} trace={c['trace']:.4f}: {direct:.6f} vs {factored:.6f}"


# ==========================================================
# § 1. U_{P-1} vanishing ⟺ (B)
#
# Under (B): every non-central element has ord(g) | P,
# so α = mπ/ord(g) with Pα = multiple of π → sin(Pα) = 0
# → U_{P-1}(cos α) = 0.
#
# Without (B): ∃ g with ord(g) ∤ P → sin(Pα) ≠ 0.
# ==========================================================

class TestUVanishing:
    """U_{P-1} vanishes at all non-central McKay eigenvalues iff (B)."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_polyhedral_all_zero(self, name, builder):
        """(B) true: U_{P-1}(cos α) = 0 at every non-central class."""
        G = builder()
        classes, order, P, jstar = prepare_group(G)
        for c in classes:
            if abs(c['trace'] - 2) < 1e-10 or abs(c['trace'] + 2) < 1e-10:
                continue  # skip central
            x = c['trace'] / 2
            val = chebyshev_U(P - 1, x)
            assert abs(val) < 1e-8, \
                f"{name} trace={c['trace']:.4f}: U_{P-1}({x:.4f}) = {val:.6f}"

    @pytest.mark.parametrize("n", [4, 6])
    def test_even_dihedral_all_zero(self, n):
        """D_n even (B true): U_{P-1} = 0 at non-central classes."""
        G = build_binary_dihedral(n)
        classes, order, P, jstar = prepare_group(G)
        for c in classes:
            if abs(c['trace'] - 2) < 1e-10 or abs(c['trace'] + 2) < 1e-10:
                continue
            x = c['trace'] / 2
            val = chebyshev_U(P - 1, x)
            assert abs(val) < 1e-8, \
                f"D_{n} trace={c['trace']:.4f}: U_{P-1}({x:.4f}) = {val:.6f}"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_odd_dihedral_nonzero(self, n):
        """D_n odd (B false): U_{P-1} ≠ 0 at half-turn class (trace=0).

        For D_n odd, half-turns (order 4, α=π/2) are the ONLY violating
        class: all other non-central elements have ord | P = n.
        So checking U_{P-1}(0) ≠ 0 suffices for the iff."""
        G = build_binary_dihedral(n)
        classes, order, P, jstar = prepare_group(G)
        # Half-turns have trace = 0, alpha = pi/2
        half_turn_classes = [c for c in classes
                            if abs(c['trace']) < 0.01
                            and c['size'] > 1]  # exclude identity
        assert len(half_turn_classes) > 0, f"D_{n}: no half-turn class found"
        for c in half_turn_classes:
            val = chebyshev_U(P - 1, 0.0)
            assert abs(val) > 0.5, \
                f"D_{n}: U_{P-1}(0) = {val:.6f}, expected nonzero"


# ==========================================================
# § 2. Selective duality via mode vanishing
#
# For half-integer j: χ_j(g) = U_{2j}(cos α). At α = π/2:
#   U_{2j}(0) = sin((2j+1)π/2) / sin(π/2) = sin(kπ) = 0
# because 2j+1 is even when j is half-integer.
#
# So even when U_{P-1}(0) ≠ 0 (D_n odd), the CHARACTER
# itself vanishes at the violating class → spinorial sector
# is immune to (B)-violation.
# ==========================================================

class TestSelectiveDuality:
    """Half-integer modes vanish at half-turns."""

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_spinorial_mode_vanishing(self, n):
        """U_{2j}(0) = 0 for all half-integer j up to 2P."""
        P = n
        for m in range(4 * P):
            j = (2 * m + 1) / 2  # half-integer
            val = chebyshev_U(int(2 * j), 0.0)
            assert abs(val) < 1e-10, \
                f"D_{n} j={j}: U_{int(2*j)}(0) = {val:.6f}, expected 0"


# ==========================================================
# § 3. S³ palindrome: μ(ℓ) = m_S³(ℓ) − ℓ
#
# Define μ(ℓ) = m_S³(ℓ) − ℓ for ℓ = 0, ..., 2j*.
# Claim: μ(ℓ) = μ(2j* − ℓ) (palindromic).
#
# Proof: Under (B), χ_j = −χ_{j*−j} at non-central g,
# so χ_j² = χ_{j*−j}² there. At central g = ±1:
# (2j+1)² − (2(j*−j)+1)² = (sum)(diff) = 2P·2(ℓ−j*).
# Weighted by 2/|Γ| = 1/(2P): difference = ℓ − ℓ'.
# Hence m_S³(ℓ) − m_S³(ℓ') = ℓ − ℓ', i.e. μ(ℓ) = μ(ℓ').
#
# For D_n odd (B false): palindrome still holds. Proof:
# - Regular classes (ord | P): χ_j² = χ_{j*-j}² (same as (B) case).
# - Half-turns (trace=0, α=π/2): χ_j(0) = sin((2j+1)π/2) = (−1)^j,
#   so χ_j(0)² = 1 for ALL j. Hence χ_j² = χ_{j*-j}² = 1 at half-turns.
# - Central (±1): (2j+1)² ≠ (2(j*-j)+1)², difference = 2P·(2ℓ−2j*).
#   Weighted by 2/|Γ| = 1/(2P): gives ℓ − ℓ'.
# So m_S³(ℓ) − m_S³(ℓ') = ℓ − ℓ' universally, i.e. μ(ℓ) = μ(ℓ').
# ==========================================================

class TestS3Palindrome:
    """μ(ℓ) = m_S³(ℓ) − ℓ is palindromic under ℓ ↔ 2j*−ℓ."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_polyhedral(self, name, builder):
        """Polyhedral (B true): μ palindrome."""
        G = builder()
        classes, order, P, jstar = prepare_group(G)
        L = 2 * jstar
        ms3 = m_S3_seq(classes, order, L)
        mu = [ms3[ell] - ell for ell in range(L + 1)]
        for ell in range(L + 1):
            assert mu[ell] == mu[L - ell], \
                f"{name} ℓ={ell}: μ({ell})={mu[ell]} ≠ μ({L-ell})={mu[L-ell]}"

    @pytest.mark.parametrize("n", [4, 6])
    def test_even_dihedral(self, n):
        """D_n even (B true): μ palindrome."""
        G = build_binary_dihedral(n)
        classes, order, P, jstar = prepare_group(G)
        L = 2 * jstar
        ms3 = m_S3_seq(classes, order, L)
        mu = [ms3[ell] - ell for ell in range(L + 1)]
        for ell in range(L + 1):
            assert mu[ell] == mu[L - ell], \
                f"D_{n} ℓ={ell}: μ({ell})={mu[ell]} ≠ μ({L-ell})={mu[L-ell]}"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_odd_dihedral(self, n):
        """D_n odd (B false): μ palindrome still holds."""
        G = build_binary_dihedral(n)
        classes, order, P, jstar = prepare_group(G)
        L = 2 * jstar
        ms3 = m_S3_seq(classes, order, L)
        mu = [ms3[ell] - ell for ell in range(L + 1)]
        for ell in range(L + 1):
            assert mu[ell] == mu[L - ell], \
                f"D_{n} ℓ={ell}: μ({ell})={mu[ell]} ≠ μ({L-ell})={mu[L-ell]}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
