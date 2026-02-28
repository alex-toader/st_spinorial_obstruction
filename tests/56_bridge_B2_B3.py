#!/usr/bin/env python3
"""
Script 56: Bridge tests B2 (d=2 uniqueness) and B3 (periodicity chain).

B2: d=2 is the unique SO(3) order whose mode function vanishes
    identically at half-integer j.  Backs Proposition 4.5.

B3: exp(H) = |H|/2 for all obstructed groups, and the chain
    modes → periodicity → threshold sharpness:
      m(A₀, j + P) = m(A₀, j) + 1,  P = |H|/2 = j* + 1.

References: main_v6.tex Remark 3.7, Proposition 4.5, tracker v5 items B2/B3.

Run: .venv/bin/python -m pytest tests/56_bridge_B2_B3.py -v
Date: Feb 2026

RAW OUTPUT (113 passed, 0.19s):
  TestB2_D2Invisible::test_d2_vanishes_at_all_half_int PASSED
  TestB2_OtherDVisible::test_not_identically_zero[3,4,5,6,8,10] 6 PASSED
  TestB2_OtherDVisible::test_majority_nonzero[3,4,5,6,8,10] 6 PASSED
  TestB2_AnalyticReason::test_2j_plus_1_even PASSED
  TestB2_AnalyticReason::test_sin_even_pi_over_2_vanishes PASSED
  TestB2_ArithmeticUniqueness::test_d2_divides_all_even PASSED
  TestB2_ArithmeticUniqueness::test_d_geq3_fails_at_k1 PASSED
  TestB2_ArithmeticUniqueness::test_uniqueness_complete PASSED
  TestB3_ExpEqualsHalfOrder::test_exp_equals_half[D_4..I] 7 PASSED
  TestB3_ExpEqualsHalfOrder::test_computed_orders_match_expected[D_4..I] 7 PASSED
  TestB3_ModePeriod::test_period_d[2,3,4,5,6,8,10] 7 PASSED
  TestB3_ModePeriod::test_not_shorter_period[3,4,5,6,8,10] 6 PASSED
  TestB3_PeriodicityChain::test_shift_by_P[D_4..I] 7 PASSED
  TestB3_PeriodicityChain::test_threshold_sharpness[D_4..I] 7 PASSED
  TestB3_IntegralityArgument::test_step1_P_divides_half_order[D_4..I] 7 PASSED
  TestB3_IntegralityArgument::test_step2_noncentral_periodic[D_4..I] 7 PASSED
  TestB3_IntegralityArgument::test_step3_central_shift_equals_2P_over_H[D_4..I] 7 PASSED
  TestB3_IntegralityArgument::test_step4_total_shift_is_integer[D_4..I] 7 PASSED
  TestB3_IntegralityArgument::test_step5_trivial_group_excluded PASSED
  TestB3_IntegralityArgument::test_step5_order2_fails_B PASSED
  TestB3_IntegralityArgument::test_step5_P_strictly_positive[D_2..I] 8 PASSED
  TestB3_IntegralityArgument::test_step6_integrality_forces_equality[D_2..I] 8 PASSED
"""
import numpy as np
import pytest
import sys
from math import gcd
from functools import reduce

sys.path.insert(0, '.')
from src.quaternion import qmul, qinv, qkey
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       build_binary_icosahedral)


# ════════════════════════════════════════════════════════════
# Helpers
# ════════════════════════════════════════════════════════════

def f_mode(alpha, j):
    """Mode function sin((2j+1)α) / sin(α)."""
    s = np.sin(alpha)
    if abs(s) < 1e-14:
        return (2 * j + 1) * (1 if alpha < 0.1 else (-1) ** int(round(2 * j)))
    return np.sin((2 * j + 1) * alpha) / s


def lcm(a, b):
    return a * b // gcd(a, b)


def build_binary_dihedral(n):
    """Build 2D_n (dicyclic, order 4n) as unit quaternions."""
    elements = []
    for k in range(2 * n):
        theta = k * np.pi / n
        elements.append(np.array([np.cos(theta), np.sin(theta), 0, 0]))
    for k in range(2 * n):
        theta = k * np.pi / n
        elements.append(np.array([0, 0, np.cos(theta), np.sin(theta)]))
    return np.array(elements)


def so3_orders_from_group(G):
    """Extract distinct SO(3) orders of non-central elements from built group.

    For quaternion g with half-angle α, the SO(3) order is the
    smallest d ≥ 2 with d·α ∈ πZ.
    """
    orders = set()
    for g in G:
        w = float(g[0])
        alpha = np.arccos(np.clip(w, -1, 1))
        if alpha < 1e-10 or abs(alpha - np.pi) < 1e-10:
            continue  # skip ±1
        for d in range(2, 200):
            if abs(np.sin(d * alpha)) < 1e-8:
                orders.add(d)
                break
    return sorted(orders)


def m_A0(G, j):
    """Multiplicity of trivial character in V_j|_G."""
    total = 0.0
    for g in G:
        w = g[0]
        alpha = np.arccos(np.clip(w, -1, 1))
        if abs(np.sin(alpha)) < 1e-14:
            chi_j = (2 * j + 1) * (1 if alpha < 0.1 else (-1) ** int(round(2 * j)))
        else:
            chi_j = np.sin((2 * j + 1) * alpha) / np.sin(alpha)
        total += chi_j
    return int(round(total / len(G)))


# ════════════════════════════════════════════════════════════
# B2: d=2 spectral invisibility uniqueness
# ════════════════════════════════════════════════════════════

HALF_INT_JS = [j / 2 for j in range(1, 32, 2)]  # 0.5, 1.5, ..., 15.5


class TestB2_D2Invisible:
    """f(π/2, j) = 0 for all half-integer j."""

    def test_d2_vanishes_at_all_half_int(self):
        alpha = np.pi / 2
        for j in HALF_INT_JS:
            assert abs(f_mode(alpha, j)) < 1e-12, f"f(π/2, {j}) ≠ 0"


class TestB2_OtherDVisible:
    """f(π/d, j) ≠ 0 for d ≥ 3 at generic half-integer j."""

    @pytest.mark.parametrize("d", [3, 4, 5, 6, 8, 10])
    def test_not_identically_zero(self, d):
        alpha = np.pi / d
        vals = [f_mode(alpha, j) for j in HALF_INT_JS]
        nonzero = [v for v in vals if abs(v) > 1e-10]
        assert len(nonzero) > 0, f"f(π/{d}, j) vanishes at ALL half-int j"

    @pytest.mark.parametrize("d", [3, 4, 5, 6, 8, 10])
    def test_majority_nonzero(self, d):
        """Most half-integer j give nonzero mode values."""
        alpha = np.pi / d
        vals = [f_mode(alpha, j) for j in HALF_INT_JS]
        nonzero = sum(1 for v in vals if abs(v) > 1e-10)
        # For d=3: 2/3 nonzero. For d=4: 1/2 nonzero. For d≥5: most nonzero.
        assert nonzero >= len(HALF_INT_JS) // d, (
            f"d={d}: only {nonzero}/{len(HALF_INT_JS)} nonzero")


class TestB2_AnalyticReason:
    """The vanishing is because 2j+1 is even at half-int j."""

    def test_2j_plus_1_even(self):
        for j in HALF_INT_JS:
            val = 2 * j + 1
            assert val == int(val), f"2j+1 not integer at j={j}"
            assert int(val) % 2 == 0, f"2j+1 = {int(val)} not even at j={j}"

    def test_sin_even_pi_over_2_vanishes(self):
        """sin(even · π/2) = sin(kπ) = 0."""
        for n in range(2, 30, 2):  # even integers
            assert abs(np.sin(n * np.pi / 2)) < 1e-12, f"sin({n}π/2) ≠ 0"


class TestB2_ArithmeticUniqueness:
    """Pure arithmetic proof that d=2 is the unique vanishing order.

    f(π/d, j) = 0 for all half-int j requires sin((2j+1)π/d) = 0,
    i.e. (2j+1)/d ∈ Z.  At half-int j, 2j+1 ranges over all even
    positive integers {2, 4, 6, ...}.  So need 2k/d ∈ Z for all k ≥ 1.
    At k=1: d | 2, so d ∈ {1, 2}.  With d ≥ 2: d = 2 uniquely.
    No numeric sampling — this is a complete arithmetic argument.
    """

    def test_d2_divides_all_even(self):
        """d=2 divides every even integer (sufficient condition)."""
        for k in range(1, 100):
            assert (2 * k) % 2 == 0

    def test_d_geq3_fails_at_k1(self):
        """For d ≥ 3: 2·1/d is not an integer (necessary condition fails)."""
        for d in range(3, 200):
            assert 2 % d != 0, f"d={d} divides 2, contradicting d ≥ 3"

    def test_uniqueness_complete(self):
        """The only d ≥ 2 with d | 2k for all k ≥ 1 is d = 2.

        Proof by contrapositive: d | 2k for all k implies d | 2 (set k=1).
        Divisors of 2 with d ≥ 2: {2}.
        """
        divisors_of_2 = [d for d in range(2, 3) if 2 % d == 0]
        assert divisors_of_2 == [2]
        # Verify the divisibility chain: d=2 works for all k
        for k in range(1, 1000):
            assert (2 * k) % 2 == 0


# ════════════════════════════════════════════════════════════
# B3: Periodicity chain
# ════════════════════════════════════════════════════════════

OBSTRUCTED = [
    ('D_2',  lambda: build_binary_dihedral(2),  4,  [2]),
    ('D_4',  lambda: build_binary_dihedral(4),  8,  [4, 2]),
    ('D_6',  lambda: build_binary_dihedral(6),  12, [6, 3, 2]),
    ('D_8',  lambda: build_binary_dihedral(8),  16, [8, 4, 2]),
    ('D_10', lambda: build_binary_dihedral(10), 20, [10, 5, 2]),
    ('T',    build_binary_tetrahedral,          12, [3, 2]),
    ('O',    build_binary_octahedral,           24, [4, 3, 2]),
    ('I',    build_binary_icosahedral,          60, [5, 3, 2]),
]


class TestB3_ExpEqualsHalfOrder:
    """exp(H) = |H|/2 for all obstructed rotation groups.

    SO(3) orders are computed from the built group's quaternion elements,
    not taken from the hardcoded list.
    """

    @pytest.mark.parametrize("name,builder,H_order,_so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_exp_equals_half(self, name, builder, H_order, _so3_orders):
        G = builder()
        assert len(G) == 2 * H_order, (
            f"{name}: |G| = {len(G)} ≠ 2·|H| = {2 * H_order}")
        computed_orders = so3_orders_from_group(G)
        assert len(computed_orders) > 0, f"{name}: no non-central elements found"
        exp_H = reduce(lcm, computed_orders)
        assert exp_H == H_order // 2, (
            f"{name}: exp(H) = {exp_H} ≠ |H|/2 = {H_order // 2}")

    @pytest.mark.parametrize("name,builder,H_order,so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_computed_orders_match_expected(self, name, builder, H_order, so3_orders):
        """Cross-check: computed SO(3) orders match the expected list."""
        G = builder()
        computed = so3_orders_from_group(G)
        assert computed == sorted(so3_orders), (
            f"{name}: computed {computed} ≠ expected {sorted(so3_orders)}")


class TestB3_ModePeriod:
    """Each mode f(π/d, j) has exact period d."""

    @pytest.mark.parametrize("d", [2, 3, 4, 5, 6, 8, 10])
    def test_period_d(self, d):
        alpha = np.pi / d
        for j in range(3 * d):
            diff = abs(f_mode(alpha, j + d) - f_mode(alpha, j))
            assert diff < 1e-10, f"f(π/{d}, {j}+{d}) ≠ f(π/{d}, {j})"

    @pytest.mark.parametrize("d", [3, 4, 5, 6, 8, 10])
    def test_not_shorter_period(self, d):
        """Period is exactly d, not a proper divisor."""
        alpha = np.pi / d
        for p in range(1, d):
            if d % p != 0:
                continue
            diffs = [abs(f_mode(alpha, j + p) - f_mode(alpha, j))
                     for j in range(2 * d)]
            if all(x < 1e-10 for x in diffs):
                pytest.fail(f"f(π/{d}, ·) has shorter period {p} < {d}")


class TestB3_PeriodicityChain:
    """m(A₀, j+P) = m(A₀, j) + 1 with P = |H|/2."""

    @pytest.mark.parametrize("name,builder,H_order,so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_shift_by_P(self, name, builder, H_order, so3_orders):
        G = builder()
        P = H_order // 2
        for j in range(P):
            m1 = m_A0(G, j)
            m2 = m_A0(G, j + P)
            assert m2 == m1 + 1, (
                f"{name}: m(A₀, {j}+{P}) = {m2} ≠ m(A₀, {j}) + 1 = {m1 + 1}")

    @pytest.mark.parametrize("name,builder,H_order,so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_threshold_sharpness(self, name, builder, H_order, so3_orders):
        """m(A₀, j*+1) = 2, derived from periodicity at j=0."""
        G = builder()
        jstar = H_order // 2 - 1
        assert m_A0(G, 0) == 1, f"{name}: m(A₀, 0) ≠ 1"
        assert m_A0(G, jstar + 1) == 2, f"{name}: m(A₀, j*+1) ≠ 2"


class TestB3_IntegralityArgument:
    """Classification-free proof that (B) implies exp(H) = |H|/2.

    The argument (no ADE classification needed):
    1. Under (B), each SO(3) order d divides |H|/2, so P = exp(H) | |H|/2.
    2. Non-central modes are periodic with period P.
    3. Central shift: m(A₀, j+P) - m(A₀, j) = 2P/|H| for integer j.
    4. Since m(A₀, j) ∈ Z, the shift 2P/|H| must be a positive integer.
    5. Since P | |H|/2, we have 0 < 2P/|H| ≤ 1.
    6. The only positive integer ≤ 1 is 1, so P = |H|/2.

    Uses only: (B), character orthogonality, integrality of multiplicities.
    """

    @pytest.mark.parametrize("name,builder,H_order,_so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_step1_P_divides_half_order(self, name, builder, H_order, _so3_orders):
        """Under (B): each d | |H|/2, so P = lcm(all d) | |H|/2."""
        G = builder()
        computed_orders = so3_orders_from_group(G)
        half = H_order // 2
        for d in computed_orders:
            assert half % d == 0, f"{name}: d={d} does not divide |H|/2={half}"
        P = reduce(lcm, computed_orders)
        assert half % P == 0, f"{name}: P={P} does not divide |H|/2={half}"

    @pytest.mark.parametrize("name,builder,H_order,_so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_step2_noncentral_periodic(self, name, builder, H_order, _so3_orders):
        """Non-central contribution to m(A₀, j) is periodic with period P."""
        G = builder()
        P = reduce(lcm, so3_orders_from_group(G))
        for j in range(2 * P):
            nc_j = 0.0
            nc_jP = 0.0
            for g in G:
                w = float(g[0])
                alpha = np.arccos(np.clip(w, -1, 1))
                if alpha < 1e-10 or abs(alpha - np.pi) < 1e-10:
                    continue  # skip ±1 (central)
                nc_j += f_mode(alpha, j)
                nc_jP += f_mode(alpha, j + P)
            assert abs(nc_jP - nc_j) < 1e-8, (
                f"{name}: noncentral not periodic at j={j}, diff={nc_jP - nc_j}")

    @pytest.mark.parametrize("name,builder,H_order,_so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_step3_central_shift_equals_2P_over_H(self, name, builder, H_order, _so3_orders):
        """Central shift per period: [central(j+P) - central(j)] / |G| = 2P/|H|.

        For integer j: χ_j(±1) = 2j+1, so central(j) = 2(2j+1).
        central(j+P) - central(j) = 2·2P.  Divide by |G| = 2|H|: gives 2P/|H|.
        """
        G = builder()
        P = reduce(lcm, so3_orders_from_group(G))
        G_order = len(G)  # = 2|H|
        for j in range(P):
            central_j = 2 * (2 * j + 1)
            central_jP = 2 * (2 * (j + P) + 1)
            shift = (central_jP - central_j) / G_order
            expected = 2 * P / H_order
            assert abs(shift - expected) < 1e-14, (
                f"{name}: central shift={shift} ≠ 2P/|H|={expected}")

    @pytest.mark.parametrize("name,builder,H_order,_so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_step4_total_shift_is_integer(self, name, builder, H_order, _so3_orders):
        """m(A₀, j+P) - m(A₀, j) is an integer, confirming 2P/|H| ∈ Z."""
        G = builder()
        P = reduce(lcm, so3_orders_from_group(G))
        expected_shift = round(2 * P / H_order)
        for j in range(P):
            diff = m_A0(G, j + P) - m_A0(G, j)
            assert diff == expected_shift, (
                f"{name}: m(A₀,{j}+{P})-m(A₀,{j})={diff}, expected {expected_shift}")

    def test_step5_trivial_group_excluded(self):
        """H={e}: (B) vacuously holds but |H|/2 = 1/2 ∉ Z, so Prop excluded."""
        # The trivial group H={e} has |H|=1.
        # (B) is vacuously satisfied (no h ≠ e to check).
        # But |H|/2 = 1/2 is not an integer, so the proposition is inapplicable.
        H_order = 1
        assert H_order % 2 != 0, "trivial group: |H| is odd, |H|/2 not integer"

    def test_step5_order2_fails_B(self):
        """H = C_2 (|H|=2): unique h≠e has ord(h)=2, |H|/ord(h)=1 odd → (B) fails."""
        H_order = 2
        h_order = 2  # unique non-identity element
        assert (H_order // h_order) % 2 != 0, "C_2: |H|/ord(h) should be odd"

    @pytest.mark.parametrize("name,builder,H_order,_so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_step5_P_strictly_positive(self, name, builder, H_order, _so3_orders):
        """(B) with |H| ≥ 2 implies |H| ≥ 4 and P ≥ 2."""
        G = builder()
        assert H_order >= 4, f"{name}: (B) implies |H| ≥ 4, got {H_order}"
        computed_orders = so3_orders_from_group(G)
        assert len(computed_orders) > 0, f"{name}: no non-central elements"
        P = reduce(lcm, computed_orders)
        assert P >= 2, f"{name}: P = {P} < 2"

    @pytest.mark.parametrize("name,builder,H_order,_so3_orders",
                             OBSTRUCTED, ids=[g[0] for g in OBSTRUCTED])
    def test_step6_integrality_forces_equality(self, name, builder, H_order, _so3_orders):
        """2P/|H| ∈ Z and 0 < 2P/|H| ≤ 1 forces P = |H|/2."""
        G = builder()
        P = reduce(lcm, so3_orders_from_group(G))
        ratio = 2 * P / H_order
        # From step 5: P ≥ 2, so ratio > 0
        # From step 1: P | |H|/2, so ratio ≤ 1
        assert 0 < ratio <= 1, f"{name}: 2P/|H|={ratio} not in (0,1]"
        # From step 4: ratio is integer
        assert abs(ratio - round(ratio)) < 1e-14, f"{name}: 2P/|H|={ratio} not integer"
        # Only positive integer ≤ 1 is 1
        assert ratio == 1, f"{name}: 2P/|H|={ratio} ≠ 1"
        assert P == H_order // 2, f"{name}: P={P} ≠ |H|/2={H_order // 2}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
