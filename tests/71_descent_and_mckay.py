#!/usr/bin/env python3
"""
I32 — Descent theorem, McKay eigenvalues, and Coxeter check.

QUESTION: What structural connections exist between S³ palindrome,
S² reflection, McKay graph spectrum, and Coxeter exponents?

RAW OUTPUT:
  §0: Universal identity m(j+P)+m(P-1-j)=2 from S³ palindrome.
      Verified 2T, 2O, 2I, D₃, D₄, D₅, D₆, D₇.
  §1: Staircase m(j+P)=m(j)+1 iff (B). Increments uniform {1} under (B),
      alternating {0,2} without (B).
  §2: Descent: (I)+(II)⟹reflection. Verified explicitly.
  §3: McKay graph eigenvalues = conjugacy class traces for Ê₆, Ê₇, Ê₈.
  §4: Coxeter eigenvalues 2cos(πm_k/h) ≠ affine eigenvalues. No match.
  §5: C2 negative: palindrome = trig identity from Lagrange, not Poincaré duality.

ANSWER:
  C1 POSITIVE: Reflection = universal palindrome + arithmetic staircase.
    (B) is the descent condition for S³ symmetry to project to S² reflection.
  C2 NEGATIVE: Palindrome is sin²((ℓ+1)α)=sin²((L-ℓ+1)α) from sin(|H|α)=0.
    Group theory (Lagrange), not topology (Poincaré duality).
  C4 PARTIAL NEGATIVE: McKay eigenvalues = affine Dynkin eigenvalues (confirmed),
    but these ≠ Coxeter eigenvalues 2cos(πm_k/h). No clean P-vs-h formula.

PLAN:
  §0. Universal identity from palindrome [8 tests]
  §1. Staircase iff (B) [8 tests]
  §2. Descent: universal + staircase ⟹ reflection [5 tests]
  §3. McKay = affine Dynkin eigenvalues [3 tests]
  §4. Coxeter ≠ affine [3 tests]
  §5. Palindrome mechanism: trig not topology [11 tests]

STATUS: COMPLETE — 38/38 tests pass
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
    assert len(G) == 4 * n, f"Expected |2D_{n}|={4*n}, got {len(G)}"
    return G


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


def m_s2(j, classes, order):
    """m(A₀, j) on S²."""
    s = sum(chi_su2(j, c['trace']) * c['size'] for c in classes)
    return int(round(s / order))


def prepare(G):
    """Return (classes, order, P)."""
    classes = compute_conjugacy_classes(G)
    classes.sort(key=lambda c: -c['trace'])
    order = len(G)
    P = order // 4
    return classes, order, P


def graph_eigenvalues(edges, n):
    """Eigenvalues of adjacency matrix of a graph."""
    A = np.zeros((n, n))
    for i, j in edges:
        A[i, j] = A[j, i] = 1
    return sorted(np.linalg.eigvalsh(A), reverse=True)


# ==========================================================
# § 0. Universal identity: m(j+P) + m(P-1-j) = 2
#
# From S³ palindrome μ(ℓ) = μ(L-ℓ) where L = 2P-2:
#   Differencing gives m(ℓ) + m(2P-1-ℓ) = 2 for all ℓ.
#   Setting ℓ = j+P: m(j+P) + m(P-1-j) = 2.
# This is UNIVERSAL — holds for ALL groups regardless of (B).
# ==========================================================

class TestUniversalIdentity:
    """m(j+P) + m(P-1-j) = 2 for all groups."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_polyhedral(self, name, builder):
        """(B) true: universal identity holds."""
        G = builder()
        classes, order, P = prepare(G)
        s2 = [m_s2(j, classes, order) for j in range(2 * P)]
        for j in range(P):
            assert s2[j + P] + s2[P - 1 - j] == 2, \
                f"{name} j={j}: m({j+P})+m({P-1-j}) = {s2[j+P]}+{s2[P-1-j]}"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_odd_dihedral(self, n):
        """(B) false: universal identity STILL holds."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        s2 = [m_s2(j, classes, order) for j in range(2 * P)]
        for j in range(P):
            assert s2[j + P] + s2[P - 1 - j] == 2, \
                f"D_{n} j={j}: m({j+P})+m({P-1-j}) = {s2[j+P]}+{s2[P-1-j]}"

    @pytest.mark.parametrize("n", [4, 6])
    def test_even_dihedral(self, n):
        """(B) true: universal identity holds."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        s2 = [m_s2(j, classes, order) for j in range(2 * P)]
        for j in range(P):
            assert s2[j + P] + s2[P - 1 - j] == 2, \
                f"D_{n} j={j}: m({j+P})+m({P-1-j}) = {s2[j+P]}+{s2[P-1-j]}"


# ==========================================================
# § 1. Staircase m(j+P) = m(j) + 1 iff (B)
#
# Under (B): every non-central mode has period P, central
# contribution increases by 1 per period → uniform increment.
# Without (B): some modes have period 2P → increments {0,2}.
# ==========================================================

class TestStaircase:
    """Staircase holds iff (B)."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_polyhedral_uniform(self, name, builder):
        """(B) true: all increments = 1."""
        G = builder()
        classes, order, P = prepare(G)
        s2 = [m_s2(j, classes, order) for j in range(2 * P)]
        increments = [s2[j + P] - s2[j] for j in range(P)]
        assert all(inc == 1 for inc in increments), \
            f"{name}: increments = {increments}, expected all 1"

    @pytest.mark.parametrize("n", [4, 6])
    def test_even_dihedral_uniform(self, n):
        """D_n even (B true): all increments = 1."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        s2 = [m_s2(j, classes, order) for j in range(2 * P)]
        increments = [s2[j + P] - s2[j] for j in range(P)]
        assert all(inc == 1 for inc in increments), \
            f"D_{n}: increments = {increments}"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_odd_dihedral_nonuniform(self, n):
        """D_n odd (B false): increments alternate {0, 2}."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        s2 = [m_s2(j, classes, order) for j in range(2 * P)]
        increments = [s2[j + P] - s2[j] for j in range(P)]
        assert set(increments) == {0, 2}, \
            f"D_{n}: increments = {set(increments)}, expected {{0,2}}"


# ==========================================================
# § 2. Descent: universal identity + staircase ⟹ reflection
#
# From §0: m(j+P) + m(P-1-j) = 2.
# From §1: m(j+P) = m(j) + 1 under (B).
# Substituting: m(j) + 1 + m(P-1-j) = 2 ⟹ m(j) + m(j*-j) = 1.
# ==========================================================

class TestDescent:
    """Reflection follows from universal identity + staircase."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_polyhedral_reflection(self, name, builder):
        """(B) true: descent gives reflection."""
        G = builder()
        classes, order, P = prepare(G)
        jstar = P - 1
        s2 = [m_s2(j, classes, order) for j in range(2 * P)]
        # Check all three parts
        for j in range(P):
            universal = s2[j + P] + s2[P - 1 - j]
            staircase = s2[j + P] - s2[j]
            reflection = s2[j] + s2[jstar - j]
            assert universal == 2, f"{name} j={j}: universal={universal}"
            assert staircase == 1, f"{name} j={j}: staircase={staircase}"
            assert reflection == 1, f"{name} j={j}: reflection={reflection}"

    @pytest.mark.parametrize("n", [3, 5])
    def test_odd_dihedral_no_descent(self, n):
        """(B) false: universal holds, staircase fails, reflection fails."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        jstar = P - 1
        s2 = [m_s2(j, classes, order) for j in range(2 * P)]
        universal_ok = all(s2[j + P] + s2[P - 1 - j] == 2 for j in range(P))
        staircase_ok = all(s2[j + P] - s2[j] == 1 for j in range(P))
        reflection_ok = all(s2[j] + s2[jstar - j] == 1 for j in range(P))
        assert universal_ok, f"D_{n}: universal should hold"
        assert not staircase_ok, f"D_{n}: staircase should fail"
        assert not reflection_ok, f"D_{n}: reflection should fail"


# ==========================================================
# § 3. McKay graph eigenvalues = affine Dynkin eigenvalues
#
# The McKay adjacency matrix of 2H has eigenvalues equal to
# the conjugacy class traces {2cos(α_g)}.
# These are also the eigenvalues of the affine Dynkin diagram.
# ==========================================================

class TestMcKayAffine:
    """Affine Dynkin diagram eigenvalues = class traces."""

    def test_E6_hat(self):
        """Ê₆ eigenvalues = 2T class traces."""
        # Ê₆: three arms of length 2 from central node
        edges = [(1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (6, 0)]
        eigs = sorted(graph_eigenvalues(edges, 7), reverse=True)
        G = build_binary_tetrahedral()
        classes, order, P = prepare(G)
        traces = sorted([c['trace'] for c in classes], reverse=True)
        for e, t in zip(eigs, traces):
            assert abs(e - t) < 1e-10, \
                f"Ê₆: graph eig {e:.4f} != trace {t:.4f}"

    def test_E7_hat(self):
        """Ê₇ eigenvalues = 2O class traces."""
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (3, 7)]
        eigs = sorted(graph_eigenvalues(edges, 8), reverse=True)
        G = build_binary_octahedral()
        classes, order, P = prepare(G)
        traces = sorted([c['trace'] for c in classes], reverse=True)
        for e, t in zip(eigs, traces):
            assert abs(e - t) < 1e-10, \
                f"Ê₇: graph eig {e:.4f} != trace {t:.4f}"

    def test_E8_hat(self):
        """Ê₈ eigenvalues = 2I class traces."""
        # Ê₈: backbone 0-1-2-3-4-5-6-7, branch at 2 to node 8
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (2, 8)]
        eigs = sorted(graph_eigenvalues(edges, 9), reverse=True)
        G = build_binary_icosahedral()
        classes, order, P = prepare(G)
        traces = sorted([c['trace'] for c in classes], reverse=True)
        for e, t in zip(eigs, traces):
            assert abs(e - t) < 1e-10, \
                f"Ê₈: graph eig {e:.4f} != trace {t:.4f}"


# ==========================================================
# § 4. Coxeter eigenvalues ≠ affine eigenvalues
#
# The Coxeter eigenvalues 2cos(πm_k/h) (from the FINITE
# Dynkin diagram) are NOT the same as the affine eigenvalues
# (= class traces). The connection to Coxeter theory is indirect.
# ==========================================================

class TestCoxeterMismatch:
    """Coxeter eigenvalues ≠ affine eigenvalues."""

    @pytest.mark.parametrize("name,builder,h,exponents", [
        ("E6", build_binary_tetrahedral, 12, [1, 4, 5, 7, 8, 11]),
        ("E7", build_binary_octahedral, 18, [1, 5, 7, 9, 11, 13, 17]),
        ("E8", build_binary_icosahedral, 30, [1, 7, 11, 13, 17, 19, 23, 29]),
    ])
    def test_coxeter_ne_affine(self, name, builder, h, exponents):
        """2cos(πm_k/h) at Coxeter exponents ≠ class traces."""
        G = builder()
        classes, _, _ = prepare(G)
        traces = sorted([c['trace'] for c in classes], reverse=True)
        # Coxeter eigenvalues (with affine m=0)
        coxeter = sorted([2 * np.cos(np.pi * m / h) for m in [0] + exponents],
                         reverse=True)
        # They should NOT match: same count, but different values
        assert len(coxeter) == len(traces), \
            f"{name}: unexpected count mismatch {len(coxeter)} vs {len(traces)}"
        mismatches = sum(abs(c - t) > 0.01 for c, t in zip(coxeter, traces))
        assert mismatches >= 1, \
            f"{name}: Coxeter eigenvalues unexpectedly equal affine eigenvalues"


# ==========================================================
# § 5. Palindrome mechanism: trig, not topology
#
# sin²((ℓ+1)α) = sin²((L-ℓ+1)α) when sin(|H|α) = 0.
# This follows from Lagrange's theorem (ord(g) | |2H|),
# not from Poincaré duality of S³/Γ.
# ==========================================================

class TestPalindromeMechanism:
    """Palindrome comes from sin(|H|α)=0, not topology."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_sin_H_alpha_zero(self, name, builder):
        """sin(|H|·α) = 0 at ALL non-central group elements."""
        G = builder()
        classes, order, _ = prepare(G)
        H_order = order // 2
        for c in classes:
            if abs(c['trace'] - 2) < 0.01 or abs(c['trace'] + 2) < 0.01:
                continue
            alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
            val = np.sin(H_order * alpha)
            assert abs(val) < 1e-10, \
                f"{name} trace={c['trace']:.4f}: sin({H_order}·α) = {val:.6f}"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_sin_H_alpha_zero_dihedral(self, n):
        """sin(|H|·α) = 0 for D_n odd too — universal."""
        G = build_binary_dihedral(n)
        classes, order, _ = prepare(G)
        H_order = order // 2  # = 2n
        for c in classes:
            if abs(c['trace'] - 2) < 0.01 or abs(c['trace'] + 2) < 0.01:
                continue
            alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
            val = np.sin(H_order * alpha)
            assert abs(val) < 1e-10, \
                f"D_{n} trace={c['trace']:.4f}: sin({H_order}·α) = {val:.6f}"

    @pytest.mark.parametrize("name,builder", [
        ("2O", build_binary_octahedral),
        ("2T", build_binary_tetrahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_squared_character_symmetry_B_true(self, name, builder):
        """χ_j(g)² = χ_{j*-j}(g)² at non-central g, (B) true groups.

        Under (B): χ_{j*-j} = -χ_j so squares match trivially.
        Pairing is j ↔ j*-j (= P-1-j)."""
        G = builder()
        classes, order, P = prepare(G)
        jstar = P - 1
        for c in classes:
            if abs(c['trace'] - 2) < 0.01 or abs(c['trace'] + 2) < 0.01:
                continue
            for j in range(jstar + 1):
                chi_j = chi_su2(j, c['trace'])
                chi_partner = chi_su2(jstar - j, c['trace'])
                assert abs(chi_j ** 2 - chi_partner ** 2) < 1e-8, \
                    f"{name} j={j}: χ²={chi_j**2:.4f} vs χ²_partner={chi_partner**2:.4f}"

    @pytest.mark.parametrize("n", [3, 5])
    def test_squared_character_symmetry_B_false(self, n):
        """χ_j(g)² = χ_{j*-j}(g)² at non-central g, (B) false (D_n odd).

        Without (B): at violating classes χ_{j*-j} = +χ_j (not -χ_j),
        so squares still match but for a different reason. Universality."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        jstar = P - 1
        for c in classes:
            if abs(c['trace'] - 2) < 0.01 or abs(c['trace'] + 2) < 0.01:
                continue
            for j in range(jstar + 1):
                chi_j = chi_su2(j, c['trace'])
                chi_partner = chi_su2(jstar - j, c['trace'])
                assert abs(chi_j ** 2 - chi_partner ** 2) < 1e-8, \
                    f"D_{n} j={j}: χ²={chi_j**2:.4f} vs χ²_partner={chi_partner**2:.4f}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
