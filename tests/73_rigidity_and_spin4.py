#!/usr/bin/env python3
"""
I34 — Rank-1 rigidity, McKay involution, and Spin(4) test.

QUESTION: Is reflection a rank-1 phenomenon? Can it be reformulated
as a McKay graph property? Does it survive in Spin(4)?

RAW OUTPUT:
  §0 McKay bipartiteness: D̂_n, Ê₆, Ê₇, Ê₈ are trees (|E|=|V|-1).
     Trees always bipartite → test trivial. Â_n (cycle) bipartite iff
     n+1 even = |H| even, not (B). D_n odd: D̂ bipartite but (B) fails.
  §1 Rank-1 rigidity: m(j+P)-m(j)=1 verified on 2T(j<12), 2O(j<24),
     2I(j<60) — real data, constant slope 1. SU(3): dim growth = k+1
     (non-constant), staircase impossible.
  §2 Spin(4) product: 2T×2T, 2O×2O, 2I×2I, 2T×2O all fail.
     Failure rate = exactly 50% (algebraic: zero fraction = 1/2 per
     factor → m1≠m2 iff failure). Mixed 2T×2O also 50%.
  §3 Spin(4) diagonal: Clebsch-Gordan sums m_diag(j1,j2) > 1 for
     all of 2T, 2O, 2I. Binary occupation fails → no reflection.

ANSWER:
  §0 NEGATIVE — bipartiteness ≠ (B). All affine Dynkin except Â are
     trees (trivially bipartite). Cannot discriminate.
  §1 POSITIVE — rank-1 rigidity confirmed on real SU(2) data (staircase
     slope 1 exact) + SU(3) dim argument (non-constant growth).
     See also file 55 (#87) for SU(3)/Δ(27) actual negative.
  §2 NEGATIVE — Spin(4) product: reflection fails, exactly 50% rate,
     universal for all (B)-satisfying groups.
  §3 NEGATIVE — Spin(4) diagonal: values > 1, binary fails.
  Overall: reflection is an SU(2)-specific rank-1 phenomenon.

PLAN:
  §0. McKay bipartiteness vs (B) — all ADE families are trees except Â
  §1. Rank-1 rigidity: staircase on real SU(2) data + SU(3) dim argument
  §2. Spin(4) product subgroups: reflection fails (50% rate)
  §3. Spin(4) diagonal subgroups: binary occupation fails (values > 1)

STATUS: COMPLETE — 20/20 tests pass
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
                       build_binary_icosahedral, compute_conjugacy_classes)
from src.quaternion import qmul, qkey


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
    assert len(G) == 4 * n
    assert len(set(qkey(q) for q in G)) == 4 * n, \
        f"2D_{n}: expected {4*n} distinct elements"
    return G


def chi_su2(j, trace):
    """SU(2) spin-j character at trace = 2cos(α).

    Special cases for trace ≈ ±2 (α ≈ 0 or π) avoid sin(dα)/sin(α)
    instability. Safety margin: |trace ∓ 2| < 1e-10 catches all
    near-degenerate cases; the sin(ha) < 1e-14 fallback is redundant
    (never reached for physical quaternion data)."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def m_A0(j, classes, order):
    """m(A₀, j) for SU(2) subgroup."""
    s = sum(chi_su2(j, c['trace']) * c['size'] for c in classes)
    return int(round(s / order))


def prepare(G):
    """Return (classes, order, P).

    Classes sorted by trace descending for deterministic ordering."""
    classes = compute_conjugacy_classes(G)
    classes.sort(key=lambda c: -c['trace'])
    order = len(G)
    P = order // 4
    return classes, order, P


def is_bipartite(edges, n):
    """Check if graph with n nodes and given edges is bipartite."""
    adj = [[] for _ in range(n)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)
    color = [-1] * n
    color[0] = 0
    queue = [0]
    while queue:
        u = queue.pop(0)
        for v in adj[u]:
            if color[v] == -1:
                color[v] = 1 - color[u]
                queue.append(v)
            elif color[v] == color[u]:
                return False
    return True


def A_hat(n):
    """Affine Â_n: cycle of n+1 nodes. Only affine Dynkin with cycles."""
    edges = [(i, (i + 1) % (n + 1)) for i in range(n + 1)]
    return edges, n + 1


def D_hat(n):
    """Affine D̂_n (n≥4): tree with branches at both ends. n+1 nodes, n edges."""
    assert n >= 4
    edges = [(0, 2), (1, 2)]
    for i in range(2, n - 2):
        edges.append((i, i + 1))
    edges.append((n - 2, n - 1))
    edges.append((n - 2, n))
    return edges, n + 1


# ==========================================================
# § 0. McKay bipartiteness vs (B)
#
# Hypothesis: (B) ⟺ McKay graph bipartite.
# Result: NEGATIVE. All affine Dynkin diagrams except Â_n are
# trees (|E| = |V|-1), hence trivially bipartite.
# Â_n (cycle) is bipartite iff n+1 even (= |H| even), which
# discriminates |H| parity, not (B).
# ==========================================================

class TestMcKayBipartite:
    """McKay bipartiteness does NOT characterize (B).

    Root cause: D̂_n, Ê₆, Ê₇, Ê₈ are trees (no cycles).
    Trees are always bipartite. Only Â_n has a cycle."""

    def test_A_hat_bipartite_iff_even(self):
        """Â_n (cyclic Z_{n+1}): bipartite iff n+1 even.
        This is the only non-trivial case — Â has a cycle."""
        for n in range(2, 10):
            e, nn = A_hat(n)
            assert len(e) == nn, f"Â_{n}: should have {nn} edges (cycle)"
            expected = ((n + 1) % 2 == 0)
            assert is_bipartite(e, nn) == expected, \
                f"Â_{n}: bipartite={not expected}, expected={expected}"

    def test_D_hat_is_tree_hence_bipartite(self):
        """D̂_n (n≥4): tree (|E|=|V|-1), hence ALWAYS bipartite."""
        for n in range(4, 12):
            e, nn = D_hat(n)
            assert len(e) == nn - 1, \
                f"D̂_{n}: {len(e)} edges, {nn} nodes — not a tree"
            assert is_bipartite(e, nn), f"D̂_{n} should be bipartite (tree)"

    def test_E_hat_is_tree_hence_bipartite(self):
        """Ê₆, Ê₇, Ê₈: trees (|E|=|V|-1), hence always bipartite."""
        # Ê₆: 7 nodes, 6 edges
        e6 = [(1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (6, 0)]
        assert len(e6) == 7 - 1, "Ê₆ should be a tree"
        assert is_bipartite(e6, 7)
        # Ê₇: 8 nodes, 7 edges
        e7 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (3, 7)]
        assert len(e7) == 8 - 1, "Ê₇ should be a tree"
        assert is_bipartite(e7, 8)
        # Ê₈: 9 nodes, 8 edges
        e8 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (2, 8)]
        assert len(e8) == 9 - 1, "Ê₈ should be a tree"
        assert is_bipartite(e8, 9)

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_dihedral_odd_bipartite_but_B_fails(self, n):
        """D_n odd: McKay graph D̂_{n+2} is bipartite (tree) but (B) fails."""
        e, nn = D_hat(n + 2)
        assert len(e) == nn - 1, "D̂ is a tree"
        assert is_bipartite(e, nn), "D̂ should be bipartite (tree)"
        # (B) fails for D_n odd: verify staircase breaks
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        s2 = [m_A0(j, classes, order) for j in range(2 * P)]
        increments = [s2[j + P] - s2[j] for j in range(P)]
        assert set(increments) != {1}, \
            f"D_{n}: staircase should fail (increments={set(increments)})"


# ==========================================================
# § 1. Rank-1 rigidity: dim V linear ⟺ rank 1
#
# For SU(N), dim V grows as polynomial of degree N-1.
# Constant-slope staircase requires degree 1 (linear).
# So descent only works for SU(2).
#
# SU(2) side: verified on actual multiplicities (2T, 2O, 2I).
# SU(3) side: dimensional formula shows non-constant growth.
# (See also file 55 for actual SU(3)/Δ(27) negative result #87.)
# ==========================================================

class TestRank1Rigidity:
    """Descent requires dim(V) linear in j, i.e., rank 1."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_su2_staircase_on_data(self, name, builder):
        """SU(2): m(A₀, j+P) - m(A₀, j) = 1 for all j, on real data.

        This is the core rank-1 property: the staircase has constant
        slope 1 because dim(V_j) = 2j+1 grows linearly and |2H| = 4P."""
        G = builder()
        classes, order, P = prepare(G)
        mults = [m_A0(j, classes, order) for j in range(3 * P)]
        for j in range(2 * P):
            assert mults[j + P] - mults[j] == 1, \
                f"{name}: m({j+P}) - m({j}) = {mults[j+P] - mults[j]}, expected 1"

    def test_su3_dim_growth_non_constant(self):
        """SU(3): dim V_{(k,0)} = (k+1)(k+2)/2 has growth rate k+1 (non-constant).

        Dimensional formula check (no group computation).
        A constant-slope staircase is impossible when dim growth varies.
        For actual SU(3)/Δ(27) negative result, see file 55 (#87)."""
        deltas = []
        for k in range(1, 15):
            d = (k + 1) * (k + 2) // 2
            d_prev = k * (k + 1) // 2
            delta = d - d_prev
            assert delta == k + 1, f"k={k}: delta={delta}, expected {k+1}"
            deltas.append(delta)
        assert len(set(deltas)) > 1, "SU(3) growth should vary"
        assert min(deltas) < max(deltas), \
            f"SU(3): growth range [{min(deltas)}, {max(deltas)}] should be non-trivial"


# ==========================================================
# § 2. Spin(4) product subgroups: reflection fails
#
# For G1 × G2 ⊂ SU(2)×SU(2), m factorizes:
# m(j1,j2) = m1(j1)·m2(j2).
# Reflection sum: m1·m2 + (1-m1)·(1-m2) = 2m1m2 - m1 - m2 + 1.
# Equals 1 only when m1 = m2. Under (B), zero fraction = 1/2,
# so exactly 50% of pairs fail.
# ==========================================================

class TestSpin4Product:
    """Reflection fails for product subgroups of Spin(4)."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_product_counterexample(self, name, builder):
        """G × G: find specific (j1,j2) where reflection fails."""
        G = builder()
        classes, order, P = prepare(G)
        jstar = P - 1
        m_vals = [m_A0(j, classes, order) for j in range(jstar + 1)]
        found_failure = False
        for j1 in range(jstar + 1):
            for j2 in range(jstar + 1):
                m = m_vals[j1] * m_vals[j2]
                m_paired = m_vals[jstar - j1] * m_vals[jstar - j2]
                if m + m_paired != 1:
                    found_failure = True
                    break
            if found_failure:
                break
        assert found_failure, f"{name}×{name}: should find reflection failure"

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_product_failure_rate_exactly_half(self, name, builder):
        """G × G: exactly 50% of (j1,j2) pairs fail reflection.

        Algebraic: under (B), zero fraction = 1/2.
        Failures iff m1 ≠ m2. Count: 2·(P/2)² = P²/2 out of P²."""
        G = builder()
        classes, order, P = prepare(G)
        jstar = P - 1
        m_vals = [m_A0(j, classes, order) for j in range(jstar + 1)]
        failures = 0
        total = 0
        for j1 in range(jstar + 1):
            for j2 in range(jstar + 1):
                m = m_vals[j1] * m_vals[j2]
                m_p = m_vals[jstar - j1] * m_vals[jstar - j2]
                total += 1
                if m + m_p != 1:
                    failures += 1
        assert failures > 0, f"{name}: no failures found"
        assert failures * 2 == total, \
            f"{name}: expected exactly 50% failures, got {failures}/{total}"

    def test_mixed_product_2T_x_2O(self):
        """2T × 2O: mixed product also fails with 50% rate.

        Each factor has zero fraction = 1/2 independently."""
        G1 = build_binary_tetrahedral()
        c1, o1, P1 = prepare(G1)
        G2 = build_binary_octahedral()
        c2, o2, P2 = prepare(G2)
        js1 = P1 - 1
        js2 = P2 - 1
        m1_vals = [m_A0(j, c1, o1) for j in range(js1 + 1)]
        m2_vals = [m_A0(j, c2, o2) for j in range(js2 + 1)]
        failures = 0
        total = 0
        for j1 in range(js1 + 1):
            for j2 in range(js2 + 1):
                m = m1_vals[j1] * m2_vals[j2]
                m_p = m1_vals[js1 - j1] * m2_vals[js2 - j2]
                total += 1
                if m + m_p != 1:
                    failures += 1
        assert failures * 2 == total, \
            f"2T×2O: expected 50% failures, got {failures}/{total}"


# ==========================================================
# § 3. Spin(4) diagonal subgroups: binary occupation fails
#
# For diagonal Γ = {(g,g) : g ∈ G},
# m_diag(j1,j2) = Σ_{j=|j1-j2|}^{j1+j2} m(j) (Clebsch-Gordan).
# Values > 1 ⟹ binary occupation fails ⟹ no reflection.
# ==========================================================

class TestSpin4Diagonal:
    """Binary occupation fails for diagonal subgroups of Spin(4)."""

    @pytest.mark.parametrize("name,builder", [
        ("2T", build_binary_tetrahedral),
        ("2O", build_binary_octahedral),
        ("2I", build_binary_icosahedral),
    ])
    def test_diagonal_exceeds_binary(self, name, builder):
        """Diagonal embedding: m_diag(j1,j2) > 1 for some (j1,j2).

        Checks all sub-threshold (j1,j2), reports max value."""
        G = builder()
        classes, order, P = prepare(G)
        jstar = P - 1

        def m_diag(j1, j2):
            return sum(m_A0(j, classes, order)
                       for j in range(abs(j1 - j2), j1 + j2 + 1))

        max_val = 0
        for j1 in range(jstar + 1):
            for j2 in range(jstar + 1):
                val = m_diag(j1, j2)
                if val > max_val:
                    max_val = val
        assert max_val > 1, \
            f"{name} diagonal: max m_diag = {max_val}, expected > 1 (binary fails)"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
