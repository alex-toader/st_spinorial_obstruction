#!/usr/bin/env python3
"""
I26 — Molien series and Gorenstein property.

QUESTION: Is m(A₀,j) the j-th coefficient of the Molien series
of C[x,y]^{2H}, and does reflection (C) follow from the
Gorenstein property?

ANSWER: YES (consolidation). m(A₀,j) = [t^{2j}] M(t) where
M(t) = (1/|G|) Σ 1/det(I - t·g) is the Molien series of C[x,y]^{2H}.
The reflection identity m(j) + m(j*-j) = 1 is the spectral shadow
of the Gorenstein property (Watanabe 1974): M(t) = (1+t^s)/((1-t^{d₁})(1-t^{d₂}))
with palindromic numerator and s = d₁+d₂-2.  Rational forms verified:
  2T: (1+t^{12})/((1-t^6)(1-t^8))     [hsop 6,8;   s=12]
  2O: (1+t^{18})/((1-t^8)(1-t^{12}))  [hsop 8,12;  s=18]
  2I: (1+t^{30})/((1-t^{12})(1-t^{20})) [hsop 12,20; s=30]

RAW OUTPUT:
  2T (|G|=24, j*=5):
    m(A0,j) j=0..10: [1, 0, 0, 1, 1, 0, 2, 1, 1, 2, 2]
    reflection sums:  [1, 1, 1, 1, 1, 1]
    series vs rational at t=0.3: |diff| = 4.44e-16

  2O (|G|=48, j*=11):
    m(A0,j) j=0..16: [1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 2, 1, 1, 1, 2]
    reflection sums:  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    series vs rational at t=0.3: |diff| = 0.00e+00

  2I (|G|=120, j*=29):
    m(A0,j) j=0..34: [1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,2,...]
    reflection sums:  all 1
    series vs rational at t=0.3: |diff| = 0.00e+00

PLAN:
  §0. Cross-check: character-sum coefficients vs matrix-det rational function [3 tests]
  §1. Reflection identity m(j)+m(j*-j)=1 with j* computed from |G| [6 tests]
  §2. Gorenstein rational form: lattice-point counting vs coefficients [3 tests]
  §3. Gorenstein structure: binary range + first m=2 at j*+1 [3 tests]

STATUS: COMPLETE — 15/15 tests pass
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
                       build_binary_icosahedral)
from src.quaternion import quat_to_su2


# ==========================================================
# § 0. Molien series coefficients vs m(A₀, j)
# ==========================================================

def molien_coefficients(G, max_j):
    """Compute Molien series M(t) = (1/|G|) Σ_g 1/det(I - t·g)
    as a power series up to t^{max_j}.

    For g ∈ SU(2) with eigenvalues e^{±iα}, det(I - t·g) =
    1 - t·e^{iα} - t·e^{-iα} + t² = 1 - 2t·cos(α) + t².

    M(t) = (1/|G|) Σ_g 1/(1 - 2t·cos(α_g) + t²)

    Expand in power series using 1/(1 - 2tx + t²) = Σ_j U_j(x) t^j
    where U_j is the Chebyshev polynomial of the second kind.
    But U_j(cos α) = sin((j+1)α)/sin(α) = χ_j(g).

    So coefficient of t^j in M(t) = (1/|G|) Σ_g χ_j(g) = m(A₀, j).

    Note: only integer j (via range); half-integer j not supported.
    """
    order = len(G)
    coeffs = []
    for j in range(max_j + 1):
        s = 0.0
        for g in G:
            tr = 2.0 * g[0]  # SU(2) trace
            if abs(tr - 2.0) < 1e-10:
                s += 2*j + 1
            elif abs(tr + 2.0) < 1e-10:
                s += (2*j + 1) * ((-1)**(2*j))
            else:
                alpha_half = np.arccos(np.clip(tr / 2.0, -1, 1))
                s += np.sin((2*j + 1) * alpha_half) / np.sin(alpha_half)
        val = s / order
        assert abs(val - round(val)) < 1e-6, \
            f"j={j}: non-integer multiplicity {val} (group may be broken)"
        coeffs.append(round(val))
    return coeffs


def molien_rational(G, t_val):
    """Evaluate M(t) directly as rational function.
    M(t) = (1/|G|) Σ_g 1/det(I - t·U_g)
    where U_g is the 2×2 SU(2) matrix.

    Caution: diverges as t → 1 (det → 0 for g = identity).
    Use |t| < 1 well away from roots of unity."""
    order = len(G)
    s = 0.0
    for g in G:
        U = quat_to_su2(g)
        det_val = np.linalg.det(np.eye(2) - t_val * U)
        s += 1.0 / det_val
    return (s / order).real


class TestMolienCrossCheck:
    """Cross-check: character-sum coefficients vs matrix-det rational function.

    Method A (molien_coefficients): (1/|G|) Σ_g χ_j(g) via sin formula.
    Method B (molien_rational):     (1/|G|) Σ_g 1/det(I - t·U_g) via 2×2 matrices.

    These are algebraically equivalent but numerically independent.
    """

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_series_vs_rational(self, name, builder):
        """Σ_j m(A₀,j) t^{2j} from character sum must equal
        (1/|G|) Σ_g 1/det(I-tg) from matrix determinants."""
        G = builder()
        # 80 terms: truncation error ~ t^{160}. At t=0.5: 7e-49. Safe.
        coeffs = molien_coefficients(G, 80)
        for t in [0.2, 0.3, 0.5]:
            from_series = sum(c * t**(2*j) for j, c in enumerate(coeffs))
            from_rational = molien_rational(G, t)
            assert abs(from_series - from_rational) < 1e-10, \
                f"{name} at t={t}: series={from_series}, rational={from_rational}"


# ==========================================================
# § 1. Functional equation for truncated polynomial
#
# P(t) = Σ_{j=0}^{j*} m(A₀,j) t^j
#
# Claim: P(t) + t^{j*} P(1/t) = (1 - t^{j*+1})/(1-t)
#                                = 1 + t + t² + ... + t^{j*}
#
# Proof: coefficient of t^j on LHS = m(j) + m(j*-j) = 1.
# ==========================================================

class TestReflectionIdentity:
    """Verify m(A₀,j) + m(A₀,j*-j) = 1 for j=0..j*, j* = |G|/4 - 1."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_functional_equation(self, name, builder):
        """P(t) + t^{j*} P(1/t) = 1 + t + ... + t^{j*}
        evaluated at several t values."""
        G = builder()
        j_star = len(G) // 4 - 1
        coeffs = molien_coefficients(G, j_star)

        for t in [0.3, 0.5, 0.7, 0.9]:
            P_t = sum(coeffs[j] * t**j for j in range(j_star + 1))
            P_inv = sum(coeffs[j] * t**(-j) for j in range(j_star + 1))
            lhs = P_t + t**j_star * P_inv
            rhs = sum(t**j for j in range(j_star + 1))
            assert abs(lhs - rhs) < 1e-10, \
                f"{name} at t={t}: LHS={lhs}, RHS={rhs}"

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_coefficient_level(self, name, builder):
        """Direct coefficient check: m(j) + m(j*-j) = 1."""
        G = builder()
        j_star = len(G) // 4 - 1
        coeffs = molien_coefficients(G, j_star)
        for j in range(j_star + 1):
            assert coeffs[j] + coeffs[j_star - j] == 1, \
                f"{name}: m({j}) + m({j_star - j}) ≠ 1"


# ==========================================================
# § 2. Gorenstein form of M(t)
#
# For G ⊂ SL(2,C), C[x,y]^G is Gorenstein (Watanabe 1974).
# Cohen-Macaulay (2 variables) gives
# M(t) = (1+t^s) / ((1-t^{d1})(1-t^{d2}))
# with s = d₁+d₂-2 and hsop degrees d₁, d₂:
#   2T: (1+t^{12})/((1-t^6)(1-t^8))
#   2O: (1+t^{18})/((1-t^8)(1-t^{12}))
#   2I: (1+t^{30})/((1-t^{12})(1-t^{20}))
#
# Cross-check: [t^{2j}] of rational form (lattice-point counting)
# must equal molien_coefficients(G, j).
# ==========================================================

class TestGorensteinForm:
    """Verify Gorenstein property of the Molien series."""

    def test_2O_rational_form(self):
        """2O: M(t) = (1 + t^{18}) / ((1-t^8)(1-t^{12})).
        Klein Grundformen: Φ (deg 8), Ψ (deg 12), χ (deg 18).
        Key: m(A₀, j) = [t^{2j}] M(t), since V_j = Sym^{2j}(C²)."""
        G = build_binary_octahedral()
        coeffs = molien_coefficients(G, 30)

        for j in range(len(coeffs)):
            d = 2 * j  # polynomial degree
            count = 0
            for a in range(d // 8 + 1):
                for b in range((d - 8*a) // 12 + 1):
                    if 8*a + 12*b == d:
                        count += 1
                    if 8*a + 12*b == d - 18:
                        count += 1
            assert coeffs[j] == count, \
                f"j={j}: Molien={coeffs[j]}, rational={count}"

    def test_2T_rational_form(self):
        """2T: M(t) = (1 + t^{12}) / ((1-t^6)(1-t^8)).
        hsop degrees 6, 8; module generator degree 12 (= 6+8-2)."""
        G = build_binary_tetrahedral()
        coeffs = molien_coefficients(G, 20)

        for j in range(len(coeffs)):
            d = 2 * j
            count = 0
            for a in range(d // 6 + 1):
                for b in range((d - 6*a) // 8 + 1):
                    if 6*a + 8*b == d:
                        count += 1
                    if 6*a + 8*b == d - 12:
                        count += 1
            assert coeffs[j] == count, \
                f"j={j}: Molien={coeffs[j]}, rational={count}"

    def test_2I_rational_form(self):
        """2I: M(t) = (1 + t^{30}) / ((1-t^{12})(1-t^{20})).
        Klein Grundformen: Φ (deg 12), Ψ (deg 20), χ (deg 30)."""
        G = build_binary_icosahedral()
        coeffs = molien_coefficients(G, 40)

        for j in range(len(coeffs)):
            d = 2 * j
            count = 0
            for a in range(d // 12 + 1):
                for b in range((d - 12*a) // 20 + 1):
                    if 12*a + 20*b == d:
                        count += 1
                    if 12*a + 20*b == d - 30:
                        count += 1
            assert coeffs[j] == count, \
                f"j={j}: Molien={coeffs[j]}, rational={count}"

    @pytest.mark.parametrize("name,builder,d1,d2", [
        ("2T", build_binary_tetrahedral, 6, 8),
        ("2O", build_binary_octahedral, 8, 12),
        ("2I", build_binary_icosahedral, 12, 20),
    ])
    def test_gorenstein_structure(self, name, builder, d1, d2):
        """Verify:
        (a) binary range: m ∈ {0,1} for j ≤ j* = |G|/4-1,
        (b) first repeated invariant: m(j*+1) = 2."""
        G = builder()
        j_star = len(G) // 4 - 1
        coeffs = molien_coefficients(G, j_star + 5)
        # (a) binary below threshold
        for j in range(j_star + 1):
            assert coeffs[j] in (0, 1), \
                f"{name}: m(A₀, {j}) = {coeffs[j]}, expected 0 or 1"
        # (b) first repeated invariant at j*+1
        assert coeffs[j_star + 1] == 2, \
            f"{name}: m(A₀, {j_star+1}) = {coeffs[j_star+1]}, expected 2"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
