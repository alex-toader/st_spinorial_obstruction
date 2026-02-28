#!/usr/bin/env python3
"""
Script 52: Converse formula verification.

For each non-obstructed group with |H| even (C_n n=2..20 even, D_n n=3..15 odd):
  1. Compute S_odd from group elements directly (not from analytic formula)
  2. Assert |S_odd| = |H|/2
  3. Compute m(A_0, j*) from the multiplicity formula
  4. Assert 2*|S_odd| = |H| * m(A_0, j*) (exact, no integer division)
  5. Assert m(A_0, 0) + m(A_0, j*) = 2 (duality fails at j=0)

For each obstructed group (D_n n=2..10 even, T, O, I):
  - Compute S_odd from group elements, assert S_odd = 0
  - Verify condition (B) directly: |H|/d even for all non-identity SO(3) orders d
  - Assert m(A_0, j*) = 0

S_odd is computed from the actual group elements by projecting to SO(3)
and checking |H|/ord(h) parity, not from analytic formulas.

Closes the converse direction of the duality theorem computationally.
References: main_v5.tex eqs (chi-jstar), (m-jstar).

Run: /usr/bin/python3 -m pytest tests/52_converse_formula.py -v
Date: Feb 2026

RAW OUTPUT (112 passed, 0.45s):
  TestConverseFormulaCyclic::test_S_odd_equals_H_half[2,4,...,20] 10 PASSED
  TestConverseFormulaCyclic::test_non_identity_count[2,4,...,20] 10 PASSED
  TestConverseFormulaCyclic::test_m_jstar_equals_1[2,4,...,20] 10 PASSED
  TestConverseFormulaCyclic::test_formula_exact[2,4,...,20] 10 PASSED
  TestConverseFormulaCyclic::test_duality_fails_at_j0[2,4,...,20] 10 PASSED
  TestConverseFormulaDihedralOdd::test_S_odd_equals_H_half[3,5,...,15] 7 PASSED
  TestConverseFormulaDihedralOdd::test_non_identity_count[3,5,...,15] 7 PASSED
  TestConverseFormulaDihedralOdd::test_m_jstar_equals_1[3,5,...,15] 7 PASSED
  TestConverseFormulaDihedralOdd::test_formula_exact[3,5,...,15] 7 PASSED
  TestConverseFormulaDihedralOdd::test_duality_fails_at_j0[3,5,...,15] 7 PASSED
  TestObstructedDihedralEven::test_S_odd_zero[2,4,6,8,10] 5 PASSED
  TestObstructedDihedralEven::test_condition_B[2,4,6,8,10] 5 PASSED
  TestObstructedDihedralEven::test_m_jstar_zero[2,4,6,8,10] 5 PASSED
  TestObstructedPolyhedral::test_S_odd_zero[T,O,I] 3 PASSED
  TestObstructedPolyhedral::test_condition_B[T,O,I] 3 PASSED
  TestObstructedPolyhedral::test_m_jstar_zero[T,O,I] 3 PASSED
  TestObstructedPolyhedral::test_non_identity_count[T,O,I] 3 PASSED
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
import pytest
from src.quaternion import qkey, qmul
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       build_binary_icosahedral, compute_conjugacy_classes)


def chi_su2(j, trace):
    """SU(2) spin-j character at given trace = 2*cos(alpha)."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def m_A0(j, classes, G_order):
    """Multiplicity of trivial character in spin-j rep."""
    total = 0.0
    for c in classes:
        total += c['size'] * chi_su2(j, c['trace'])
    return int(round(total / G_order))


def so3_order(q):
    """SO(3) order of the rotation corresponding to quaternion q.
    q and -q give the same SO(3) element. Uses |q[0]| = cos(half-angle).
    """
    cos_ha = abs(q[0])
    if cos_ha > 1 - 1e-12:
        return 1
    ha = np.arccos(np.clip(cos_ha, -1, 1))
    for d in range(1, 200):
        if abs(d * ha / np.pi - round(d * ha / np.pi)) < 1e-8:
            return d
    raise ValueError(f"Could not determine SO(3) order for q={q}")


def compute_S_odd(elements, H_order):
    """Compute |S_odd| directly from SU(2) group elements.

    Projects to SO(3) by picking one from each {q, -q} pair.
    For each non-identity h in H = G/{+/-1}, computes ord(h) and
    checks whether |H|/ord(h) is odd.

    Returns (s_odd_count, non_identity_count).
    """
    seen = set()
    s_odd = 0
    non_id = 0
    for i in range(len(elements)):
        q = elements[i]
        k = qkey(q)
        neg_k = qkey(-q)
        if neg_k in seen:
            continue
        seen.add(k)
        d = so3_order(q)
        if d == 1:
            continue
        non_id += 1
        K = H_order // d
        assert H_order % d == 0, f"|H|={H_order} not divisible by ord={d}"
        if K % 2 == 1:
            s_odd += 1
    return s_odd, non_id


def condition_B_holds(elements, H_order):
    """Check condition (B) directly: |H|/ord(h) is even for all h != e.
    Returns True iff (B) holds (all quotients even = obstructed).
    """
    s_odd, _ = compute_S_odd(elements, H_order)
    return s_odd == 0


# ============================================================
# Group builders
# ============================================================

def build_binary_cyclic(n):
    """2C_n = C_{2n} in SU(2). Order 2n."""
    elements = []
    for k in range(2 * n):
        theta = k * np.pi / n
        elements.append(np.array([np.cos(theta), 0, 0, np.sin(theta)]))
    return np.array(elements)


def build_binary_dihedral(n):
    """2D_n (= Dic_n) in SU(2). Order 4n."""
    j_quat = np.array([0.0, 0.0, 1.0, 0.0])
    result = []
    for k in range(2 * n):
        theta = k * np.pi / n
        q = np.array([np.cos(theta), 0, 0, np.sin(theta)])
        result.append(q)
        result.append(qmul(j_quat, q))
    return np.array(result)


def prepare_group(name, elements, H_order):
    """Prepare group data dict from elements."""
    G_order = len(elements)
    j_star = H_order // 2 - 1
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    s_odd, non_id = compute_S_odd(elements, H_order)
    return {
        'name': name, 'G': G_order, 'H': H_order,
        'j_star': j_star, 'classes': classes, 'elements': elements,
        'S_odd': s_odd, 'non_id': non_id,
    }


# ============================================================
# Non-obstructed groups: C_n (n even), D_n (n odd)
# C_n with n odd excluded: |H| odd => j* not integer.
# ============================================================

class TestConverseFormulaCyclic:
    """Converse formula for cyclic groups C_n (n=2,4,...,20)."""

    @pytest.fixture(params=[2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
                    scope='class')
    def grp(self, request):
        n = request.param
        return prepare_group(f'C_{n}', build_binary_cyclic(n), n)

    def test_S_odd_equals_H_half(self, grp):
        """Assert |S_odd| = |H|/2, computed from group elements."""
        assert grp['S_odd'] == grp['H'] // 2

    def test_non_identity_count(self, grp):
        """Sanity: |H|-1 non-identity elements in SO(3) projection."""
        assert grp['non_id'] == grp['H'] - 1

    def test_m_jstar_equals_1(self, grp):
        """Assert m(A_0, j*) = 1 via character formula."""
        m = m_A0(grp['j_star'], grp['classes'], grp['G'])
        assert m == 1

    def test_formula_exact(self, grp):
        """Assert 2*|S_odd| = |H| * m(A_0, j*) (exact integer identity)."""
        m = m_A0(grp['j_star'], grp['classes'], grp['G'])
        assert 2 * grp['S_odd'] == grp['H'] * m

    def test_duality_fails_at_j0(self, grp):
        """m(A_0, 0) + m(A_0, j*) = 2 != 1."""
        m0 = m_A0(0, grp['classes'], grp['G'])
        mj = m_A0(grp['j_star'], grp['classes'], grp['G'])
        assert m0 + mj == 2


class TestConverseFormulaDihedralOdd:
    """Converse formula for dihedral groups D_n (n=3,5,...,15)."""

    @pytest.fixture(params=[3, 5, 7, 9, 11, 13, 15], scope='class')
    def grp(self, request):
        n = request.param
        return prepare_group(f'D_{n}', build_binary_dihedral(n), 2 * n)

    def test_S_odd_equals_H_half(self, grp):
        assert grp['S_odd'] == grp['H'] // 2

    def test_non_identity_count(self, grp):
        assert grp['non_id'] == grp['H'] - 1

    def test_m_jstar_equals_1(self, grp):
        m = m_A0(grp['j_star'], grp['classes'], grp['G'])
        assert m == 1

    def test_formula_exact(self, grp):
        m = m_A0(grp['j_star'], grp['classes'], grp['G'])
        assert 2 * grp['S_odd'] == grp['H'] * m

    def test_duality_fails_at_j0(self, grp):
        m0 = m_A0(0, grp['classes'], grp['G'])
        mj = m_A0(grp['j_star'], grp['classes'], grp['G'])
        assert m0 + mj == 2


# ============================================================
# Obstructed groups: S_odd = 0, condition (B), m(A_0, j*) = 0
# ============================================================

class TestObstructedDihedralEven:
    """For D_n (n even), S_odd = 0 and m(A_0, j*) = 0."""

    @pytest.fixture(params=[2, 4, 6, 8, 10], scope='class')
    def grp(self, request):
        n = request.param
        return prepare_group(f'D_{n}', build_binary_dihedral(n), 2 * n)

    def test_S_odd_zero(self, grp):
        """S_odd = 0 from group elements (condition (B) holds)."""
        assert grp['S_odd'] == 0

    def test_condition_B(self, grp):
        """Condition (B): |H|/ord(h) even for every non-identity h."""
        assert condition_B_holds(grp['elements'], grp['H'])

    def test_m_jstar_zero(self, grp):
        m = m_A0(grp['j_star'], grp['classes'], grp['G'])
        assert m == 0


class TestObstructedPolyhedral:
    """For T, O, I: S_odd = 0, condition (B), m(A_0, j*) = 0."""

    @pytest.fixture(params=[
        ('T', build_binary_tetrahedral, 12),
        ('O', build_binary_octahedral, 24),
        ('I', build_binary_icosahedral, 60),
    ], scope='class')
    def grp(self, request):
        name, build_fn, H_order = request.param
        return prepare_group(name, build_fn(), H_order)

    def test_S_odd_zero(self, grp):
        """S_odd = 0 from group elements."""
        assert grp['S_odd'] == 0

    def test_condition_B(self, grp):
        """Condition (B): |H|/d even for every non-identity SO(3) order d."""
        assert condition_B_holds(grp['elements'], grp['H'])

    def test_m_jstar_zero(self, grp):
        m = m_A0(grp['j_star'], grp['classes'], grp['G'])
        assert m == 0

    def test_non_identity_count(self, grp):
        """Sanity: correct number of non-identity H elements."""
        assert grp['non_id'] == grp['H'] - 1


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
