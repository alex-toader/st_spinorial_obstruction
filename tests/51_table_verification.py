#!/usr/bin/env python3
"""
Script 51: Entry-level verification of paper tables.

Computes m(chi, j) for every scalar character and every j <= j*,
asserts each value against the corrected paper table, and
checks cross-consistency: the dual column m(chi, j*-j) as printed
in the paper must agree with m(chi, j*-j) from the character formula.

This is the test that would have caught the 8 wrong entries
in the octahedral table (main_v5.tex, Table 2).

Groups: 2O (j*=11, kappa=2), 2T (j*=5, kappa=3), 2I (j*=29, kappa=1).

Run: /usr/bin/python3 -m pytest tests/51_table_verification.py -v
Date: Feb 2026

RAW OUTPUT (197 passed, 1.46s):
  TestOctahedralTable::test_group_order PASSED
  TestOctahedralTable::test_j_star PASSED
  TestOctahedralTable::test_kappa PASSED
  TestOctahedralTable::test_mA0_value[0..11] 12 PASSED
  TestOctahedralTable::test_mA1_value[0..11] 12 PASSED
  TestOctahedralTable::test_A0_dual_column_consistency[0..11] 12 PASSED
  TestOctahedralTable::test_A1_dual_column_consistency[0..11] 12 PASSED
  TestOctahedralTable::test_A0_duality_sum[0..11] 12 PASSED
  TestOctahedralTable::test_A1_duality_sum[0..11] 12 PASSED
  TestOctahedralTable::test_m_A0_jstar_zero PASSED
  TestOctahedralTable::test_m_A1_jstar_one PASSED
  TestOctahedralTable::test_first_nonzero PASSED
  TestOctahedralTable::test_zero_count PASSED
  TestTetrahedralTable::test_group_order PASSED
  TestTetrahedralTable::test_j_star PASSED
  TestTetrahedralTable::test_kappa PASSED
  TestTetrahedralTable::test_A0_binary[0..5] 6 PASSED
  TestTetrahedralTable::test_A1_binary[0..5] 6 PASSED
  TestTetrahedralTable::test_A2_binary[0..5] 6 PASSED
  TestTetrahedralTable::test_A0_duality[0..5] 6 PASSED
  TestTetrahedralTable::test_A1_duality[0..5] 6 PASSED
  TestTetrahedralTable::test_A2_duality[0..5] 6 PASSED
  TestTetrahedralTable::test_A1_A2_conjugate[0..5] 6 PASSED
  TestTetrahedralTable::test_A0_jstar_zero PASSED
  TestTetrahedralTable::test_nonA0_jstar_one PASSED
  TestTetrahedralTable::test_print_table PASSED
  TestIcosahedralTable::test_group_order PASSED
  TestIcosahedralTable::test_j_star PASSED
  TestIcosahedralTable::test_kappa_1 PASSED
  TestIcosahedralTable::test_binary_multiplicity[0..29] 30 PASSED
  TestIcosahedralTable::test_duality[0..29] 30 PASSED
  TestIcosahedralTable::test_first_nonzero PASSED
  TestIcosahedralTable::test_zero_count PASSED
  TestIcosahedralTable::test_m_A0_jstar_zero PASSED
  TestIcosahedralTable::test_print_table PASSED
  TestObstructionKappaConsistency::test_minus1_in_commutator[2T,2O,2I] 3 PASSED
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
import pytest
from src.quaternion import qmul, qkey
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       build_binary_icosahedral,
                       compute_conjugacy_classes, compute_commutator_subgroup)

IDENTITY_KEY = qkey(np.array([1, 0, 0, 0], dtype=float))
MINUS1_KEY = qkey(np.array([-1, 0, 0, 0], dtype=float))


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


def multiplicity(j, classes, G_order, chi_values):
    """Multiplicity of character chi in spin-j rep.
    chi_values[i] = chi(class_rep_i) for each conjugacy class.
    """
    total = 0.0
    for i, c in enumerate(classes):
        total += c['size'] * np.conj(chi_values[i]) * chi_su2(j, c['trace'])
    return int(round(total.real / G_order))


# ============================================================
# Build groups
# ============================================================

def build_group_data(name, build_fn):
    elements = build_fn()
    G = len(elements)
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    comm = compute_commutator_subgroup(elements)
    comm_keys = set(comm.keys())
    return {
        'name': name,
        'G': G,
        'H': G // 2,
        'j_star': G // 4 - 1,
        'classes': classes,
        'elements': elements,
        'comm_keys': comm_keys,
    }


# ============================================================
# 2O: octahedral (j*=11, kappa=2, characters A_0 and A_1)
# ============================================================

class TestOctahedralTable:
    """Verify every entry of the octahedral table (paper Table 2)."""

    @pytest.fixture(scope='class')
    def grp(self):
        return build_group_data('2O', build_binary_octahedral)

    @pytest.fixture(scope='class')
    def multiplicities(self, grp):
        """Compute all m(A_0,j) and m(A_1,j) for j=0..11."""
        classes = grp['classes']
        G = grp['G']
        comm_keys = grp['comm_keys']

        # A_0: trivial character (all +1)
        chi_A0 = [1.0] * len(classes)

        # A_1: sign character (+1 on [G,G]=2T, -1 off)
        chi_A1 = [1.0 if qkey(c['rep']) in comm_keys else -1.0
                  for c in classes]

        mA0 = {}
        mA1 = {}
        for j in range(12):
            mA0[j] = multiplicity(j, classes, G, chi_A0)
            mA1[j] = multiplicity(j, classes, G, chi_A1)

        return mA0, mA1

    # ----------------------------------------------------------
    # Group structure
    # ----------------------------------------------------------

    def test_group_order(self, grp):
        assert grp['G'] == 48

    def test_j_star(self, grp):
        assert grp['j_star'] == 11

    def test_kappa(self, grp):
        """[2O,2O] = 2T (order 24) => kappa = 48/24 = 2."""
        assert len(grp['comm_keys']) == 24

    # ----------------------------------------------------------
    # Paper table values (corrected, main_v5.tex Table 2)
    # 7 columns: j, m(A0,j), m(A1,j), m(A0,11-j), Sigma_A0, m(A1,11-j), Sigma_A1
    # ----------------------------------------------------------

    # Column 2: m(A_0, j) for j=0..11
    PAPER_A0 = [1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0]
    # Column 3: m(A_1, j) for j=0..11
    PAPER_A1 = [0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1]
    # Column 4: m(A_0, 11-j) as printed in row j
    PAPER_A0_DUAL = [0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1]
    # Column 6: m(A_1, 11-j) as printed in row j
    PAPER_A1_DUAL = [1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0]

    # ----------------------------------------------------------
    # Value-level tests: computed m(chi,j) vs paper entries
    # ----------------------------------------------------------

    @pytest.mark.parametrize("j", range(12))
    def test_mA0_value(self, multiplicities, j):
        """Assert m(A_0, j) matches paper table column 2."""
        mA0, _ = multiplicities
        assert mA0[j] == self.PAPER_A0[j], (
            f"m(A_0, {j}) = {mA0[j]}, paper says {self.PAPER_A0[j]}")

    @pytest.mark.parametrize("j", range(12))
    def test_mA1_value(self, multiplicities, j):
        """Assert m(A_1, j) matches paper table column 3."""
        _, mA1 = multiplicities
        assert mA1[j] == self.PAPER_A1[j], (
            f"m(A_1, {j}) = {mA1[j]}, paper says {self.PAPER_A1[j]}")

    # ----------------------------------------------------------
    # Cross-consistency: dual columns match primary columns
    # This is the test that catches the original 8-entry bug.
    # The paper prints m(chi, 11-j) in a separate column;
    # that column must equal PAPER_chi[11-j] from column 2/3.
    # We encode both columns independently and check agreement.
    # ----------------------------------------------------------

    @pytest.mark.parametrize("j", range(12))
    def test_A0_dual_column_consistency(self, j):
        """Paper column 4 (m(A_0, 11-j) in row j) must equal PAPER_A0[11-j]."""
        assert self.PAPER_A0_DUAL[j] == self.PAPER_A0[11 - j], (
            f"Row {j}: paper dual col shows m(A_0,{11-j})={self.PAPER_A0_DUAL[j]}, "
            f"but primary col of row {11-j} has {self.PAPER_A0[11 - j]}")

    @pytest.mark.parametrize("j", range(12))
    def test_A1_dual_column_consistency(self, j):
        """Paper column 6 (m(A_1, 11-j) in row j) must equal PAPER_A1[11-j]."""
        assert self.PAPER_A1_DUAL[j] == self.PAPER_A1[11 - j], (
            f"Row {j}: paper dual col shows m(A_1,{11-j})={self.PAPER_A1_DUAL[j]}, "
            f"but primary col of row {11-j} has {self.PAPER_A1[11 - j]}")

    # ----------------------------------------------------------
    # Duality: m(chi, j) + m(chi, 11-j) = 1
    # ----------------------------------------------------------

    @pytest.mark.parametrize("j", range(12))
    def test_A0_duality_sum(self, multiplicities, j):
        """m(A_0, j) + m(A_0, 11-j) = 1."""
        mA0, _ = multiplicities
        assert mA0[j] + mA0[11 - j] == 1

    @pytest.mark.parametrize("j", range(12))
    def test_A1_duality_sum(self, multiplicities, j):
        """m(A_1, j) + m(A_1, 11-j) = 1."""
        _, mA1 = multiplicities
        assert mA1[j] + mA1[11 - j] == 1

    # ----------------------------------------------------------
    # Named results from Corollary 3.5
    # ----------------------------------------------------------

    def test_m_A0_jstar_zero(self, multiplicities):
        """Corollary 3.5(iii): m(A_0, j*) = 0."""
        mA0, _ = multiplicities
        assert mA0[11] == 0

    def test_m_A1_jstar_one(self, multiplicities):
        """Corollary 3.5(iii): m(chi != A_0, j*) = 1."""
        _, mA1 = multiplicities
        assert mA1[11] == 1

    def test_first_nonzero(self, multiplicities):
        """j_1 = 3 for SO(3)/O: first j > 0 with any scalar m > 0."""
        mA0, mA1 = multiplicities
        for j in range(1, 3):
            assert mA0[j] == 0 and mA1[j] == 0, (
                f"j={j}: m(A_0)={mA0[j]}, m(A_1)={mA1[j]}, both should be 0")
        assert mA1[3] == 1, f"m(A_1, 3) = {mA1[3]}, should be 1"

    def test_zero_count(self, multiplicities):
        """Exactly half the levels have m = 0 (for each character)."""
        mA0, mA1 = multiplicities
        assert sum(1 for j in range(12) if mA0[j] == 0) == 6
        assert sum(1 for j in range(12) if mA1[j] == 0) == 6


# ============================================================
# 2T: tetrahedral (j*=5, kappa=3, characters A_0, A_1, A_2)
# ============================================================

class TestTetrahedralTable:
    """Verify multiplicities for 2T (not in paper, but validates pipeline)."""

    @pytest.fixture(scope='class')
    def grp(self):
        return build_group_data('2T', build_binary_tetrahedral)

    @pytest.fixture(scope='class')
    def multiplicities(self, grp):
        """Compute m(chi, j) for all three scalar characters, j=0..5.

        A_1 and A_2 are complex conjugate (omega vs omega^2 on cosets).
        Coset labeling depends on iteration order, so A_1 <-> A_2 may
        be swapped. Tests only use symmetric properties (duality holds
        for both) or the unordered multiset {m(A_1,j), m(A_2,j)}.
        """
        classes = grp['classes']
        G = grp['G']
        comm_keys = grp['comm_keys']
        elements = grp['elements']

        # A_0: trivial
        chi_A0 = [1.0] * len(classes)

        # [2T, 2T] = Q_8 (order 8), abelianization Z_3
        comm_sub = compute_commutator_subgroup(elements)

        # Coset map: element -> coset index 0,1,2
        coset_assigned = {}
        coset_idx = 0
        for e in elements:
            k = qkey(e)
            if k in coset_assigned:
                continue
            coset = set()
            for q8_key, q8_elem in comm_sub.items():
                prod = qmul(e, q8_elem)
                coset.add(qkey(prod))
            for ck in coset:
                coset_assigned[ck] = coset_idx
            coset_idx += 1

        omega = np.exp(2j * np.pi / 3)

        # A_1 (omega on coset 1) and A_2 (omega^2 on coset 1)
        chi_A1 = [omega ** coset_assigned[qkey(c['rep'])]
                  for c in classes]
        chi_A2 = [omega ** (2 * coset_assigned[qkey(c['rep'])])
                  for c in classes]

        result = {}
        for label, chi in [('A0', chi_A0), ('A1', chi_A1), ('A2', chi_A2)]:
            m = {}
            for j in range(6):
                m[j] = multiplicity(j, classes, G, chi)
            result[label] = m

        return result

    def test_group_order(self, grp):
        assert grp['G'] == 24

    def test_j_star(self, grp):
        assert grp['j_star'] == 5

    def test_kappa(self, grp):
        """Abelianization has order 3: |2T| / |[2T,2T]| = 24 / 8 = 3."""
        assert len(grp['comm_keys']) == 8

    @pytest.mark.parametrize("j", range(6))
    def test_A0_binary(self, multiplicities, j):
        """m(A_0, j) in {0, 1} below threshold."""
        assert multiplicities['A0'][j] in (0, 1)

    @pytest.mark.parametrize("j", range(6))
    def test_A1_binary(self, multiplicities, j):
        assert multiplicities['A1'][j] in (0, 1)

    @pytest.mark.parametrize("j", range(6))
    def test_A2_binary(self, multiplicities, j):
        assert multiplicities['A2'][j] in (0, 1)

    @pytest.mark.parametrize("j", range(6))
    def test_A0_duality(self, multiplicities, j):
        m = multiplicities['A0']
        assert m[j] + m[5 - j] == 1

    @pytest.mark.parametrize("j", range(6))
    def test_A1_duality(self, multiplicities, j):
        m = multiplicities['A1']
        assert m[j] + m[5 - j] == 1

    @pytest.mark.parametrize("j", range(6))
    def test_A2_duality(self, multiplicities, j):
        m = multiplicities['A2']
        assert m[j] + m[5 - j] == 1

    @pytest.mark.parametrize("j", range(6))
    def test_A1_A2_conjugate(self, multiplicities, j):
        """A_1 and A_2 are complex conjugate: m(A_1,j) = m(A_2,j)."""
        assert multiplicities['A1'][j] == multiplicities['A2'][j]

    def test_A0_jstar_zero(self, multiplicities):
        """m(A_0, j*) = 0."""
        assert multiplicities['A0'][5] == 0

    def test_nonA0_jstar_one(self, multiplicities):
        """m(chi != A_0, j*) = 1."""
        assert multiplicities['A1'][5] == 1
        assert multiplicities['A2'][5] == 1

    def test_print_table(self, multiplicities):
        """Print the tetrahedral table for visual inspection."""
        print("\n2T tetrahedral table (j*=5):")
        print(f"  j   A0   A1   A2")
        for j in range(6):
            a0 = multiplicities['A0'][j]
            a1 = multiplicities['A1'][j]
            a2 = multiplicities['A2'][j]
            print(f"  {j}    {a0}    {a1}    {a2}")


# ============================================================
# 2I: icosahedral (j*=29, kappa=1, only A_0)
# ============================================================

class TestIcosahedralTable:
    """Verify multiplicities for 2I -- 30-row stress test."""

    @pytest.fixture(scope='class')
    def grp(self):
        return build_group_data('2I', build_binary_icosahedral)

    @pytest.fixture(scope='class')
    def mA0(self, grp):
        """Compute m(A_0, j) for j=0..29."""
        classes = grp['classes']
        G = grp['G']
        chi_A0 = [1.0] * len(classes)
        return {j: multiplicity(j, classes, G, chi_A0) for j in range(30)}

    def test_group_order(self, grp):
        assert grp['G'] == 120

    def test_j_star(self, grp):
        assert grp['j_star'] == 29

    def test_kappa_1(self, grp):
        """2I is perfect: [2I, 2I] = 2I, so kappa=1."""
        assert len(grp['comm_keys']) == 120

    @pytest.mark.parametrize("j", range(30))
    def test_binary_multiplicity(self, mA0, j):
        """m(A_0, j) in {0, 1} for j <= 29."""
        assert mA0[j] in (0, 1), f"m(A_0, {j}) = {mA0[j]}"

    @pytest.mark.parametrize("j", range(30))
    def test_duality(self, mA0, j):
        """m(A_0, j) + m(A_0, 29-j) = 1."""
        assert mA0[j] + mA0[29 - j] == 1, (
            f"m(A_0, {j}) + m(A_0, {29-j}) = {mA0[j]} + {mA0[29-j]}")

    def test_first_nonzero(self, mA0):
        """First nonzero level after j=0 is j_1=6 (from paper Table 1)."""
        assert mA0[0] == 1, "m(A_0, 0) = 1 always (trivial rep)"
        for j in range(1, 6):
            assert mA0[j] == 0, f"m(A_0, {j}) = {mA0[j]} should be 0"
        assert mA0[6] == 1, f"m(A_0, 6) = {mA0[6]} should be 1"

    def test_zero_count(self, mA0):
        """Exactly half the levels are zero."""
        n_zero = sum(1 for j in range(30) if mA0[j] == 0)
        assert n_zero == 15

    def test_m_A0_jstar_zero(self, mA0):
        """m(A_0, j*=29) = 0."""
        assert mA0[29] == 0

    def test_print_table(self, mA0):
        """Print the icosahedral table for visual inspection."""
        print("\n2I icosahedral table (j*=29):")
        print("  j  m(A_0,j)  j  m(A_0,j)")
        for j in range(15):
            jd = 29 - j
            print(f"  {j:2d}    {mA0[j]}     {jd:2d}    {mA0[jd]}")


# ============================================================
# Global obstruction-kappa consistency
# ============================================================

class TestObstructionKappaConsistency:
    """Verify: -1 in [G,G] iff kappa = |G|/|[G,G]| has no spinorial sector.

    For obstructed groups (T, O, I): -1 in [G,G], all scalar chi have chi(-1)=+1.
    For non-obstructed: -1 not in [G,G], exists chi with chi(-1)=-1.
    """

    @pytest.fixture(scope='class', params=[
        ('2T', build_binary_tetrahedral, True),
        ('2O', build_binary_octahedral, True),
        ('2I', build_binary_icosahedral, True),
    ])
    def group_data(self, request):
        name, build_fn, expected_obstructed = request.param
        elements = build_fn()
        G = len(elements)
        comm = compute_commutator_subgroup(elements)
        comm_keys = set(comm.keys())
        return {
            'name': name,
            'G': G,
            'comm_keys': comm_keys,
            'expected_obstructed': expected_obstructed,
        }

    def test_minus1_in_commutator(self, group_data):
        """-1 in [G,G] iff group is obstructed."""
        assert (MINUS1_KEY in group_data['comm_keys']) == \
               group_data['expected_obstructed']


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
