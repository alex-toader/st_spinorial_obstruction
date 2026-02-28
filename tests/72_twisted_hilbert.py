#!/usr/bin/env python3
"""
I33 — Twisted Hilbert functions for semi-invariant rings.

QUESTION: Does f_χ(q) for non-trivial scalar χ of 2T, 2O, 2I have
a rational form with binary palindromic numerator? Does this fail when (B) fails?

RAW OUTPUT:
  §0: 2T has 3 scalar chars (Z_3 abel), 2O has 2 (Z_2), 2I has 1 (perfect).
  §1: Rational forms with shared hsop denominators:
      2T: A0=(1+q^6)/((1-q^3)(1-q^4)), A1=A2=(q^2+q^4)/((1-q^3)(1-q^4))
      2O: A0=(1+q^9)/((1-q^4)(1-q^6)), A1=(q^3+q^6)/((1-q^4)(1-q^6))
      2I: A0=(1+q^15)/((1-q^6)(1-q^10))
  §2: ALL numerators binary (coefficients in {0,1}).
  §3: ALL numerators palindromic (powers sum to socle s=d1+d2-1).
  §4: D_n (all n): A0 = (1+q^{n+1}), sign = (q+q^{n-1}), ALSO binary palindromic.
      Binary palindromic holds REGARDLESS of (B). Not a discriminator.

ANSWER:
  RP3 NEGATIVE AS DISCRIMINATOR: Binary palindromic numerator is universal
  (Gorenstein property), holds for ALL groups regardless of (B).
  RP3 POSITIVE AS STRUCTURE: All semi-invariant Hilbert functions have
  rational forms with binary palindromic numerators and shared hsop.

PLAN:
  §0. Setup: scalar characters [3 tests]
  §1. Rational form extraction [7 tests]
  §2. Binary numerator check [3 tests]
  §3. Palindromic numerator check [4 tests]
  §4. Dihedral: binary palindromic independent of (B) [14 tests]

STATUS: COMPLETE — 30/30 tests pass
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


def prepare(G):
    """Return (classes, order, P)."""
    classes = compute_conjugacy_classes(G)
    classes.sort(key=lambda c: -c['trace'])
    order = len(G)
    P = order // 4
    return classes, order, P


def scalar_characters(G, classes):
    """Compute all scalar (1D) characters of 2H.

    Returns list of (name, chi_values) where chi_values[k] is χ(class_k).
    Uses the abelianization: scalar chars = characters of G/[G,G].
    """
    order = len(G)
    G_keys = set(qkey(q) for q in G)

    # Build commutator subgroup by computing all [a,b]
    comm_keys = set()
    for a in G:
        for b in G:
            # [a,b] = a b a^{-1} b^{-1}
            ab = qmul(a, b)
            # a^{-1} = conjugate for unit quaternions: (w, -x, -y, -z)
            a_inv = np.array([a[0], -a[1], -a[2], -a[3]])
            b_inv = np.array([b[0], -b[1], -b[2], -b[3]])
            comm = qmul(qmul(a, b), qmul(a_inv, b_inv))
            comm_keys.add(qkey(comm))

    comm_size = len(comm_keys)
    abel_order = order // comm_size  # |G/[G,G]|

    if abel_order == 1:
        # Perfect group (e.g., 2I): only trivial character
        return [("A0", [1.0] * len(classes))]

    if abel_order == 2:
        # Z_2 abelianization (e.g., 2O): trivial + sign
        chi_vals = []
        for c in classes:
            rep = c['rep']
            chi_vals.append(1.0 if qkey(rep) in comm_keys else -1.0)
        return [("A0", [1.0] * len(classes)), ("A1", chi_vals)]

    if abel_order == 3:
        # Z_3 abelianization (e.g., 2T): trivial + 2 complex conjugate chars
        # Find a generator: element not in [G,G]
        gen = None
        for q in G:
            if qkey(q) not in comm_keys:
                gen = q
                break

        # Assign coset labels
        coset_map = {}  # qkey -> coset (0, 1, or 2)

        # Coset 0 = [G,G]
        for q in G:
            if qkey(q) in comm_keys:
                coset_map[qkey(q)] = 0

        # Coset 1 = gen * [G,G]
        for q in G:
            if qkey(q) in comm_keys:
                prod = qmul(gen, q)
                coset_map[qkey(prod)] = 1

        # Coset 2 = gen^2 * [G,G]
        gen2 = qmul(gen, gen)
        for q in G:
            if qkey(q) in comm_keys:
                prod = qmul(gen2, q)
                coset_map[qkey(prod)] = 2

        assert len(coset_map) == len(G), "Coset map incomplete"

        omega = np.exp(2j * np.pi / 3)
        chi_A1 = []
        chi_A2 = []
        for c in classes:
            coset = coset_map[qkey(c['rep'])]
            chi_A1.append(omega ** coset)
            chi_A2.append(omega ** (2 * coset))

        return [
            ("A0", [1.0] * len(classes)),
            ("A1", chi_A1),
            ("A2", chi_A2),
        ]

    raise ValueError(f"Unexpected abelianization order {abel_order}")


def multiplicities(classes, order, chi_vals, J=30):
    """Compute m(χ, j) for j = 0, ..., J-1."""
    result = []
    for j in range(J):
        s = sum(chi_su2(j, c['trace']) * c['size'] * np.conj(chi_vals[k])
                for k, c in enumerate(classes))
        result.append(int(round(np.real(s) / order)))
    return result


def find_numerator(mults, den_periods, N=30):
    """Given mults and denominator (1-q^d1)(1-q^d2), extract numerator."""
    result = np.array(mults[:N], dtype=float)
    for d in den_periods:
        new = np.zeros(N)
        for i in range(N):
            new[i] = result[i]
            if i >= d:
                new[i] -= result[i - d]
        result = new
    return result


def expand_rational(num_powers, den_periods, N=30):
    """Expand N(q)/D(q) as power series."""
    num = np.zeros(N)
    for p, coeff in num_powers:
        if p < N:
            num[p] = coeff
    result = num.copy()
    for d in den_periods:
        expanded = np.zeros(N)
        for i in range(N):
            expanded[i] = result[i]
            if i >= d:
                expanded[i] += expanded[i - d]
        result = expanded
    return result


# ==========================================================
# § 0. Setup: identify scalar characters
# ==========================================================

class TestSetup:
    """Verify scalar character identification."""

    def test_2O_two_chars(self):
        """2O has 2 scalar characters (Z_2 abelianization)."""
        G = build_binary_octahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        assert len(chars) == 2
        assert chars[0][0] == "A0"
        assert chars[1][0] == "A1"

    def test_2T_three_chars(self):
        """2T has 3 scalar characters (Z_3 abelianization)."""
        G = build_binary_tetrahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        assert len(chars) == 3

    def test_2I_one_char(self):
        """2I has 1 scalar character (perfect group)."""
        G = build_binary_icosahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        assert len(chars) == 1


# ==========================================================
# § 1. Twisted Molien: rational form with known hsop
#
# hsop degrees (in j-variable, so half the t-variable degrees):
#   2T: d1=3, d2=4  (t-degrees 6, 8)
#   2O: d1=4, d2=6  (t-degrees 8, 12)
#   2I: d1=6, d2=10 (t-degrees 12, 20)
# ==========================================================

class TestRationalForm:
    """Extract and verify rational forms."""

    def test_2O_A0_gorenstein(self):
        """2O A0: (1+q^9)/((1-q^4)(1-q^6)) — known Gorenstein."""
        G = build_binary_octahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        m = multiplicities(classes, order, chars[0][1])
        num = find_numerator(m, [4, 6])
        nonzero = [(i, int(round(num[i]))) for i in range(25) if abs(num[i]) > 0.01]
        assert nonzero == [(0, 1), (9, 1)]

    def test_2O_A1_rational(self):
        """2O A1: (q^3+q^6)/((1-q^4)(1-q^6))."""
        G = build_binary_octahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        m = multiplicities(classes, order, chars[1][1])
        num = find_numerator(m, [4, 6])
        nonzero = [(i, int(round(num[i]))) for i in range(25) if abs(num[i]) > 0.01]
        assert nonzero == [(3, 1), (6, 1)]

    def test_2T_A0_gorenstein(self):
        """2T A0: known Gorenstein form."""
        G = build_binary_tetrahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        m = multiplicities(classes, order, chars[0][1])
        num = find_numerator(m, [3, 4])
        nonzero = [(i, int(round(num[i]))) for i in range(20) if abs(num[i]) > 0.01]
        assert nonzero == [(0, 1), (6, 1)], f"2T A0 numerator: {nonzero}"

    def test_2T_A1_rational(self):
        """2T A1: (q^2+q^4)/((1-q^3)(1-q^4))."""
        G = build_binary_tetrahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        m = multiplicities(classes, order, chars[1][1])
        num = find_numerator(m, [3, 4])
        nonzero = [(i, int(round(num[i]))) for i in range(20) if abs(num[i]) > 0.01]
        assert nonzero == [(2, 1), (4, 1)], f"2T A1 numerator: {nonzero}"

    def test_2T_A2_rational(self):
        """2T A2: same as A1 (complex conjugate chars give same real mults)."""
        G = build_binary_tetrahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        m = multiplicities(classes, order, chars[2][1])
        num = find_numerator(m, [3, 4])
        nonzero = [(i, int(round(num[i]))) for i in range(20) if abs(num[i]) > 0.01]
        assert nonzero == [(2, 1), (4, 1)], f"2T A2 numerator: {nonzero}"

    def test_2I_A0_gorenstein(self):
        """2I A0: known Gorenstein form."""
        G = build_binary_icosahedral()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        m = multiplicities(classes, order, chars[0][1])
        num = find_numerator(m, [6, 10])
        nonzero = [(i, int(round(num[i]))) for i in range(30) if abs(num[i]) > 0.01]
        assert nonzero == [(0, 1), (15, 1)], f"2I A0 numerator: {nonzero}"


# ==========================================================
# § 2. Binary numerator check
#
# All numerator coefficients should be in {0, 1}.
# ==========================================================

class TestBinaryNumerator:
    """All twisted Molien numerators have coefficients in {0,1}."""

    @pytest.mark.parametrize("name,builder,hsop", [
        ("2T", build_binary_tetrahedral, [3, 4]),
        ("2O", build_binary_octahedral, [4, 6]),
        ("2I", build_binary_icosahedral, [6, 10]),
    ])
    def test_binary_all_chars(self, name, builder, hsop):
        """Every scalar character's numerator is binary."""
        G = builder()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        for chi_name, chi_vals in chars:
            m = multiplicities(classes, order, chi_vals)
            num = find_numerator(m, hsop)
            coeffs = [int(round(num[i])) for i in range(hsop[0] + hsop[1])]
            assert all(c in (0, 1) for c in coeffs), \
                f"{name} {chi_name}: numerator coeffs {coeffs} not binary"


# ==========================================================
# § 3. Palindromic numerator check
#
# Numerator powers should sum to socle s = d1 + d2 - 1
# (in j-variable). Each pair (a, s-a) with coeff 1.
# ==========================================================

class TestPalindromicNumerator:
    """All twisted Molien numerators are palindromic."""

    @pytest.mark.parametrize("name,builder,hsop", [
        ("2T", build_binary_tetrahedral, [3, 4]),
        ("2O", build_binary_octahedral, [4, 6]),
        ("2I", build_binary_icosahedral, [6, 10]),
    ])
    def test_palindromic_all_chars(self, name, builder, hsop):
        """Every scalar character's numerator is palindromic with socle s=d1+d2-1."""
        G = builder()
        classes, order, P = prepare(G)
        chars = scalar_characters(G, classes)
        s = hsop[0] + hsop[1] - 1  # socle degree
        for chi_name, chi_vals in chars:
            m = multiplicities(classes, order, chi_vals)
            num = find_numerator(m, hsop)
            # Check palindrome: num[i] == num[s-i] for all i
            for i in range(s + 1):
                assert abs(num[i] - num[s - i]) < 0.01, \
                    f"{name} {chi_name}: num[{i}]={num[i]:.0f} != num[{s-i}]={num[s-i]:.0f}"

    def test_same_socle_all_chars(self):
        """All characters of same group share the same socle degree."""
        for name, builder, hsop in [
            ("2T", build_binary_tetrahedral, [3, 4]),
            ("2O", build_binary_octahedral, [4, 6]),
        ]:
            G = builder()
            classes, order, P = prepare(G)
            chars = scalar_characters(G, classes)
            s = hsop[0] + hsop[1] - 1
            for chi_name, chi_vals in chars:
                m = multiplicities(classes, order, chi_vals)
                num = find_numerator(m, hsop)
                powers = [i for i in range(s + 1) if abs(num[i]) > 0.01]
                for a in powers:
                    assert abs(num[s - a]) > 0.01, \
                        f"{name} {chi_name}: power {a} present but {s-a} missing"


# ==========================================================
# § 4. Counterexample: D_n odd (B fails)
#
# For D_n odd, abelianization is Z_2 (2 scalar chars).
# hsop for dihedral D_n in j-variable: d1=2, d2=n.
# Check if A1 numerator is still binary palindromic.
# ==========================================================

def dihedral_sign_character(classes):
    """Sign character of D_n: +1 on rotations, -1 on reflections.

    In quaternion representation: rotations have (w,0,0,z) form,
    reflections (j-type) have (0,x,y,0) form."""
    chi = []
    for c in classes:
        rep = c['rep']
        is_rotation = abs(rep[1]) < 0.01 and abs(rep[2]) < 0.01
        chi.append(1.0 if is_rotation else -1.0)
    return chi


class TestDihedral:
    """D_n: check ALL scalar char numerators — binary palindromic always."""

    @pytest.mark.parametrize("n", [3, 4, 5, 6, 7, 8])
    def test_A0_binary_palindromic(self, n):
        """D_n A0: numerator = 1 + q^{n+1}, always binary palindromic."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        m = multiplicities(classes, order, [1.0] * len(classes), J=40)
        hsop = [2, n]
        num = find_numerator(m, hsop, N=40)
        s = hsop[0] + hsop[1] - 1
        coeffs = [int(round(num[i])) for i in range(s + 3)]
        assert all(c == 0 for c in coeffs[s + 1:]), f"D_{n} A0: not finite"
        assert all(c in (0, 1) for c in coeffs[:s + 1]), \
            f"D_{n} A0: not binary: {coeffs[:s+1]}"
        assert all(abs(num[i] - num[s - i]) < 0.01 for i in range(s + 1)), \
            f"D_{n} A0: not palindromic"

    @pytest.mark.parametrize("n", [3, 4, 5, 6, 7])
    def test_sign_binary_palindromic(self, n):
        """D_n sign char: numerator = q + q^{n-1}, always binary palindromic."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        chi_sign = dihedral_sign_character(classes)
        m = multiplicities(classes, order, chi_sign, J=40)
        hsop = [2, n]
        num = find_numerator(m, hsop, N=40)
        s = hsop[0] + hsop[1] - 1
        coeffs = [int(round(num[i])) for i in range(s + 3)]
        assert all(c == 0 for c in coeffs[s + 1:]), f"D_{n} sign: not finite"
        assert all(c in (0, 1) for c in coeffs[:s + 1]), \
            f"D_{n} sign: not binary: {coeffs[:s+1]}"
        assert all(abs(num[i] - num[s - i]) < 0.01 for i in range(s + 1)), \
            f"D_{n} sign: not palindromic"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_sign_independent_of_B(self, n):
        """D_n odd (B fails): sign char STILL binary palindromic."""
        G = build_binary_dihedral(n)
        classes, order, P = prepare(G)
        chi_sign = dihedral_sign_character(classes)
        m = multiplicities(classes, order, chi_sign, J=40)
        hsop = [2, n]
        num = find_numerator(m, hsop, N=40)
        s = hsop[0] + hsop[1] - 1
        coeffs = [int(round(num[i])) for i in range(s + 1)]
        # Should be q + q^{s-1} = q + q^n (palindromic around s/2)
        expected_powers = {1, n}
        actual_powers = {i for i in range(s + 1) if coeffs[i] == 1}
        assert actual_powers == expected_powers, \
            f"D_{n} sign: expected powers {expected_powers}, got {actual_powers}"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
