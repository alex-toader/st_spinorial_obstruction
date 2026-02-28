#!/usr/bin/env python3
"""
I3 — Gap universality: is the first spinorial level always j=1/2?

INVESTIGATION: For every K with Q₈ ⊄ 2K (non-obstructed),
does there exist a spinorial character χ with m(χ, 1/2) ≥ 1?

ANSWER: NO. Gap universality fails for dihedral groups.

RESULT:
  Cyclic C_n:         j_min^spin = 1/2,   gap = 3B/4       (universal)
  Dihedral D_n odd:   j_min^spin = n/2,   gap = n(n+2)B/4  (grows with n)

PROOF SKETCH (for D_n):
  V_j|_{Z_{2n}} has weights {2m mod 2n : m = -j,...,j}.
  The 1D spinorial chars of Dic_n restrict to χ_n on Z_{2n}.
  Need weight n in V_j, i.e., 2m = n ⟹ m = n/2.
  Requires j ≥ n/2. First half-integer: j = n/2.

RAW OUTPUT (42 passed in 0.10s):
  TestCyclicGap::test_j_min_half[1..12]         10 PASSED
  TestDihedralGap::test_m_half_zero[3..11]        5 PASSED
  TestDihedralGap::test_j_min_equals_n_over_2[3..11]  5 PASSED
  TestDihedralGap::test_gap_formula[3..11]        5 PASSED
  TestDihedralGap::test_below_j_min_all_zero[3..11]   5 PASSED
  TestProofVerification::test_coset_trace_zero[3..7]   3 PASSED
  TestProofVerification::test_spinorial_char_restricts_to_chi_n[3..7]  3 PASSED
  TestProofVerification::test_V_half_weights_miss_n[3..7]  3 PASSED
  TestProofVerification::test_V_nhalf_weights_hit_n[3..7]  3 PASSED

STATUS: COMPLETE — clean negative result with exact formula
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest


# ==========================================================
# § 0. SU(2) characters
# ==========================================================

def chi_su2(j, g_trace):
    """Character of V_j at element with tr(g) = g_trace.
    χ_j(g) = sin((2j+1)α/2) / sin(α/2) where tr(g) = 2cos(α/2)."""
    if abs(g_trace - 2.0) < 1e-10:
        return 2*j + 1
    if abs(g_trace + 2.0) < 1e-10:
        return (2*j + 1) * ((-1)**(int(2*j)))
    alpha_half = np.arccos(np.clip(g_trace / 2.0, -1, 1))
    return np.sin((2*j + 1) * alpha_half) / np.sin(alpha_half)


# ==========================================================
# § 1. Non-obstructed groups
#
# K ⊂ SO(3) with Q₈ ⊄ 2K:
#   Cyclic: C_n for all n ≥ 1
#   Dihedral: D_n for n odd (n = 3, 5, 7, ...)
#
# Obstructed (skip): D_n for n even, T, O, I
# ==========================================================

def build_2Cn(n):
    """Build binary cyclic group 2C_n = Z_{2n} ⊂ SU(2).
    Elements: e^{kπi/n} for k = 0,...,2n-1.
    As quaternions: (cos(πk/n), 0, 0, sin(πk/n))."""
    elements = []
    for k in range(2*n):
        angle = np.pi * k / n
        elements.append((np.cos(angle), 0.0, 0.0, np.sin(angle)))
    return elements


def build_2Dn_odd(n):
    """Build binary dihedral group 2D_n for n odd (= Dic_n).
    Order 4n. Presentation: a^{2n} = 1, b^2 = a^n, bab^{-1} = a^{-1}.
    Elements: a^k (k=0..2n-1) and a^k·j (k=0..2n-1).
    Here a = (cos(π/n), 0, 0, sin(π/n)), j = (0, 0, 1, 0)."""
    elements = []
    # First half: a^k
    for k in range(2*n):
        angle = np.pi * k / n
        elements.append((np.cos(angle), 0.0, 0.0, np.sin(angle)))
    # Second half: a^k · j via Hamilton product
    # (c, 0, 0, s) · (0, 0, 1, 0) = (0, -s, c, 0)
    for k in range(2*n):
        angle = np.pi * k / n
        c, s = np.cos(angle), np.sin(angle)
        elements.append((0.0, -s, c, 0.0))
    return elements


def get_scalar_chars_2Cn(n):
    """Scalar (1D) characters of 2C_n = Z_{2n}.
    χ_k(a^m) = e^{2πi·k·m/(2n)}, spinorial iff k odd."""
    chars = {}
    elements = build_2Cn(n)
    order = 2 * n
    for k in range(order):
        char_vals = {}
        for m, q in enumerate(elements):
            key = tuple(round(x, 10) for x in q)
            char_vals[key] = np.exp(2.0j * np.pi * k * m / order)
        chars[f'chi_{k}'] = {
            'values': char_vals,
            'spinorial': (k % 2 == 1),
        }
    return chars, elements


def get_scalar_chars_2Dn_odd(n):
    """Scalar (1D) characters of 2D_n (n odd) = Dic_n.
    Abelianization: Z_4 with π(a) = g², π(b) = g.
    χ_k(a^m) = i^{2km}, χ_k(a^m·b) = i^{k(2m+1)}.
    Spinorial (k odd): χ_k(-1) = (-1)^k = -1."""
    elements = build_2Dn_odd(n)
    order = 4 * n
    zeta = 1j  # fourth root of unity
    chars = {}
    for k in range(4):
        char_vals = {}
        for m in range(2*n):
            q = elements[m]
            key = tuple(round(x, 10) for x in q)
            char_vals[key] = zeta ** (2 * k * m)
        for m in range(2*n):
            q = elements[2*n + m]
            key = tuple(round(x, 10) for x in q)
            char_vals[key] = zeta ** (k * (2 * m + 1))
        chars[f'chi_{k}'] = {
            'values': char_vals,
            'spinorial': (k % 2 == 1),
        }
    return chars, elements


# ==========================================================
# § 2. Multiplicity computation
# ==========================================================

def multiplicity_spin(chi_vals, elements, order, j):
    """m(χ, j) = (1/|G|) Σ_g conj(χ(g)) · χ_j(g)."""
    s = 0.0
    for q in elements:
        key = tuple(round(x, 10) for x in q)
        s += np.conj(chi_vals.get(key, 0.0)) * chi_su2(j, 2.0 * q[0])
    return (s / order).real


def find_j_min_spin(chars, elements, order, max_j=30):
    """Find smallest half-integer j with m(χ_spin, j) ≥ 1."""
    for j_twice in range(1, 2*max_j + 1, 2):
        j = j_twice / 2.0
        for name, data in chars.items():
            if data['spinorial']:
                m = multiplicity_spin(data['values'], elements, order, j)
                if round(m) >= 1:
                    return j
    return None


# ==========================================================
# § 3. Tests: cyclic groups (gap universal at j=1/2)
# ==========================================================

class TestCyclicGap:
    """C_n: gap universality holds. j_min^spin = 1/2 always."""

    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 7, 8, 10, 12])
    def test_j_min_half(self, n):
        """V_{1/2}|_{Z_{2n}} = χ_1 ⊕ χ_{2n-1}, both spinorial."""
        chars, elements = get_scalar_chars_2Cn(n)
        order = 2 * n
        j_min = find_j_min_spin(chars, elements, order)
        assert j_min == 0.5, f"C_{n}: j_min = {j_min}, expected 0.5"


# ==========================================================
# § 4. Tests: dihedral groups (gap universality FAILS)
# ==========================================================

class TestDihedralGap:
    """D_n (n odd): gap universality fails. j_min^spin = n/2."""

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11])
    def test_m_half_zero(self, n):
        """V_{1/2}|_{Dic_n} is an irreducible 2D rep → no 1D components.
        So m(χ_spin, 1/2) = 0 for all spinorial χ."""
        chars, elements = get_scalar_chars_2Dn_odd(n)
        order = 4 * n
        for name, data in chars.items():
            if data['spinorial']:
                m = multiplicity_spin(data['values'], elements, order, 0.5)
                assert abs(m) < 1e-10, \
                    f"D_{n}, {name}: m(χ, 1/2) = {m}, expected 0"

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11])
    def test_j_min_equals_n_over_2(self, n):
        """j_min^spin(D_n) = n/2. Proof: need weight n in V_j|_{Z_{2n}},
        i.e., 2m = n, so m = n/2, requiring j ≥ n/2."""
        chars, elements = get_scalar_chars_2Dn_odd(n)
        order = 4 * n
        j_min = find_j_min_spin(chars, elements, order)
        expected = n / 2.0
        assert j_min == expected, \
            f"D_{n}: j_min = {j_min}, expected {expected}"

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11])
    def test_gap_formula(self, n):
        """Gap = j_min(j_min+1)B = n(n+2)/4 · B."""
        j_min = n / 2.0
        gap = j_min * (j_min + 1)
        expected = n * (n + 2) / 4.0
        assert abs(gap - expected) < 1e-10, \
            f"D_{n}: gap = {gap}, expected n(n+2)/4 = {expected}"

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11])
    def test_below_j_min_all_zero(self, n):
        """For j < n/2 (half-integer), m(χ_spin, j) = 0."""
        chars, elements = get_scalar_chars_2Dn_odd(n)
        order = 4 * n
        j_min = n / 2.0
        for j_twice in range(1, n, 2):  # j = 1/2, 3/2, ..., < n/2
            j = j_twice / 2.0
            for name, data in chars.items():
                if data['spinorial']:
                    m = multiplicity_spin(data['values'], elements, order, j)
                    assert abs(m) < 1e-10, \
                        f"D_{n}: m(χ_spin, {j}) = {m} ≠ 0 below j_min"


# ==========================================================
# § 5. Proof verification: why V_{1/2}|_{Dic_n} has no 1D part
#
# V_j|_{Z_{2n}} has Z_{2n}-weights {2m mod 2n : m = -j,...,j}.
# The 1D spinorial chars of Dic_n restrict to χ_n on Z_{2n}.
# For V_{1/2}: weights are {-1, +1}. Neither equals n (for n ≥ 3).
# For V_{n/2}: weights include ±n, which equals n mod 2n. ✓
# ==========================================================

class TestProofVerification:
    """Verify the representation-theoretic argument."""

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_coset_trace_zero(self, n):
        """All coset elements a^k·b have trace 0.
        So they never contribute to m(χ, j) for half-integer j.
        (trace 0 → α/2 = π/2 → χ_j = sin((2j+1)π/2).
        At j half-integer, 2j+1 is even, so sin(even·π/2) = 0.
        This is Lemma 4.4 (spectral invisibility, d=2) from the paper.)"""
        elements = build_2Dn_odd(n)
        for k in range(2*n):
            q = elements[2*n + k]
            assert abs(q[0]) < 1e-10, \
                f"Coset element a^{k}·b has nonzero trace {2*q[0]}"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_spinorial_char_restricts_to_chi_n(self, n):
        """χ_1(Dic_n) restricted to Z_{2n} = <a> gives χ_n(Z_{2n}).
        I.e., χ_1(a^m) = e^{iπnm/n} = e^{iπm} = (-1)^m."""
        chars, _ = get_scalar_chars_2Dn_odd(n)
        elements = build_2Dn_odd(n)
        chi1 = chars['chi_1']['values']
        for m in range(2*n):
            q = elements[m]
            key = tuple(round(x, 10) for x in q)
            expected = (-1)**m
            actual = chi1[key]
            assert abs(actual - expected) < 1e-10, \
                f"χ_1(a^{m}) = {actual}, expected (-1)^{m} = {expected}"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_V_half_weights_miss_n(self, n):
        """V_{1/2} on Z_{2n} has weights {-1, +1} mod 2n.
        Neither equals n for n ≥ 3, so no overlap with spinorial chars."""
        weights = {(-1) % (2*n), 1 % (2*n)}  # {2n-1, 1}
        assert n not in weights, \
            f"Weight n={n} found in V_{{1/2}} weights {weights}"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_V_nhalf_weights_hit_n(self, n):
        """V_{n/2} on Z_{2n} has weight 2·(n/2) = n, hitting the target."""
        j = n / 2.0
        weights = set()
        for m_twice in range(-n, n + 1):  # 2m from -n to +n
            weights.add(m_twice % (2*n))
        assert n in weights, \
            f"Weight n={n} not found in V_{{{j}}} weights"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
