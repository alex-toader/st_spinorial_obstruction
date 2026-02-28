#!/usr/bin/env python3
"""
Script 53: Above-threshold breakdown verification.

For each obstructed polyhedral group (2T, 2O, 2I):
  1. Compute m(A₀, j) for j = j*+1 .. j*+P (P = |G|/4, one full period)
  2. Assert m(A₀, j*+1) = 2 (binary occupation fails immediately)
  3. Assert all m ≥ 1 above threshold (j* is truly the last zero)
  4. Assert periodicity: m(j+P) = m(j) + 1 for j = 0 .. j*+P
  5. Pin all above-threshold values against hardcoded table

For 2O: also verify A₁ (sign character) above threshold.

The periodicity m(j+P) = m(j) + 1 follows from S(j+P) = S(j) (non-central
sum is periodic with period P = |G|/4) and the central term growing by |G|.
This was proved analytically in file 38; here we verify all values.

Key result: the threshold j* is maximally sharp. Since j*+1 = P = |G|/4,
periodicity gives m(A₀, j*+1) = m(A₀, 0) + 1 = 2. The binary occupation
law m ∈ {0,1} breaks at the very first level above threshold.

Above threshold, the duality sum shifts: for j in [j*+1, 2j*+1],
m(A₀, j) + m(A₀, 2j*+1-j) = 2 (instead of 1 below threshold).
Proof: j = j*+1+k maps to dual j' = j*-k, and m(j*+1+k) = m(k)+1
by periodicity, while m(j*-k) = 1-m(k) by sub-threshold duality.
Sum = (m(k)+1) + (1-m(k)) = 2.

References: main_v5.tex Corollary 3.5, file 38 (periodicity proof).

Run: /usr/bin/python3 -m pytest tests/53_above_threshold.py -v
Date: Feb 2026

RAW OUTPUT (266 passed, 0.42s):
  TestAboveThreshold2O::test_commutator_subgroup_order PASSED
  TestAboveThreshold2O::test_j_star PASSED
  TestAboveThreshold2O::test_period PASSED
  TestAboveThreshold2O::test_jstar_plus_1_equals_P PASSED
  TestAboveThreshold2O::test_m_A0_jstar_plus_1_equals_2 PASSED
  TestAboveThreshold2O::test_m_A0_jstar_zero PASSED
  TestAboveThreshold2O::test_hardcoded_table_length PASSED
  TestAboveThreshold2O::test_m_A0_above_pinned[12..23] 12 PASSED
  TestAboveThreshold2O::test_m_A1_above_pinned[12..23] 12 PASSED
  TestAboveThreshold2O::test_m_A0_geq_1[12..23] 12 PASSED
  TestAboveThreshold2O::test_periodicity_A0[0..23] 24 PASSED
  TestAboveThreshold2O::test_count_geq2_A0 PASSED (6/12)
  TestAboveThreshold2O::test_count_geq2_A1 PASSED (6/12)
  TestAboveThreshold2O::test_hardcoded_self_consistency_A0 PASSED
  TestAboveThreshold2O::test_hardcoded_self_consistency_A1 PASSED
  TestAboveThreshold2O::test_extended_duality_sum_A0[12..23] 12 PASSED (sum=2)
  TestAboveThreshold2O::test_extended_duality_sum_A1[12..23] 12 PASSED (sum=2)
  TestAboveThreshold2T::test_j_star PASSED
  TestAboveThreshold2T::test_period PASSED
  TestAboveThreshold2T::test_jstar_plus_1_equals_P PASSED
  TestAboveThreshold2T::test_m_A0_jstar_plus_1_equals_2 PASSED
  TestAboveThreshold2T::test_m_A0_jstar_zero PASSED
  TestAboveThreshold2T::test_hardcoded_table_length PASSED
  TestAboveThreshold2T::test_m_A0_above_pinned[6..11] 6 PASSED
  TestAboveThreshold2T::test_m_A0_geq_1[6..11] 6 PASSED
  TestAboveThreshold2T::test_periodicity_A0[0..11] 12 PASSED
  TestAboveThreshold2T::test_count_geq2_A0 PASSED (3/6)
  TestAboveThreshold2T::test_hardcoded_self_consistency_A0 PASSED
  TestAboveThreshold2T::test_extended_duality_sum_A0[6..11] 6 PASSED (sum=2)
  TestAboveThreshold2I::test_j_star PASSED
  TestAboveThreshold2I::test_period PASSED
  TestAboveThreshold2I::test_jstar_plus_1_equals_P PASSED
  TestAboveThreshold2I::test_m_A0_jstar_plus_1_equals_2 PASSED
  TestAboveThreshold2I::test_m_A0_jstar_zero PASSED
  TestAboveThreshold2I::test_hardcoded_table_length PASSED
  TestAboveThreshold2I::test_m_A0_above_pinned[30..59] 30 PASSED
  TestAboveThreshold2I::test_m_A0_geq_1[30..59] 30 PASSED
  TestAboveThreshold2I::test_periodicity_A0[0..29] 30 PASSED
  TestAboveThreshold2I::test_count_geq2_A0 PASSED (15/30)
  TestAboveThreshold2I::test_extended_duality_sum_A0[30..59] 30 PASSED (sum=2)
  TestUniversalThresholdSharpness::test_jstar_plus_1_equals_P[2T,2O,2I] 3 PASSED
  TestUniversalThresholdSharpness::test_m_jstar_plus_1_equals_2[2T,2O,2I] 3 PASSED
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
import pytest
from src.quaternion import qkey
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       build_binary_icosahedral, compute_conjugacy_classes,
                       compute_commutator_subgroup)


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


def m_chi(j, classes, G_order, chi_bar_fn=None):
    """Multiplicity of character chi in spin-j rep.
    chi_bar_fn: function(class_dict) -> chi_bar value. None = trivial (A₀).
    """
    total = 0.0
    for c in classes:
        w = 1.0 if chi_bar_fn is None else chi_bar_fn(c)
        total += c['size'] * w * chi_su2(j, c['trace'])
    return int(round(total / G_order))


def prepare_polyhedral(name, builder, G_order):
    """Build group data for a polyhedral group."""
    elems = builder()
    classes = compute_conjugacy_classes(elems)
    classes.sort(key=lambda c: -c['trace'])
    j_star = G_order // 4 - 1
    P = G_order // 4
    return {
        'name': name, 'G': G_order, 'j_star': j_star, 'P': P,
        'classes': classes, 'elements': elems,
    }


# ============================================================
# Hardcoded above-threshold values (one full period past j*)
# Computed from character formula; pinned here for regression.
# ============================================================

# 2T: j*=5, P=6. Above-threshold j=6..11.
ABOVE_2T_A0 = {6: 2, 7: 1, 8: 1, 9: 2, 10: 2, 11: 1}

# 2O: j*=11, P=12. Above-threshold j=12..23.
ABOVE_2O_A0 = {12: 2, 13: 1, 14: 1, 15: 1, 16: 2, 17: 1,
               18: 2, 19: 1, 20: 2, 21: 2, 22: 2, 23: 1}
ABOVE_2O_A1 = {12: 1, 13: 1, 14: 1, 15: 2, 16: 1, 17: 1,
               18: 2, 19: 2, 20: 1, 21: 2, 22: 2, 23: 2}

# 2I: j*=29, P=30. Above-threshold j=30..59.
ABOVE_2I_A0 = {30: 2, 31: 1, 32: 1, 33: 1, 34: 1, 35: 1,
               36: 2, 37: 1, 38: 1, 39: 1, 40: 2, 41: 1,
               42: 2, 43: 1, 44: 1, 45: 2, 46: 2, 47: 1,
               48: 2, 49: 1, 50: 2, 51: 2, 52: 2, 53: 1,
               54: 2, 55: 2, 56: 2, 57: 2, 58: 2, 59: 1}


# Sub-threshold values (from paper tables, independently encoded).
# Used for self-consistency: ABOVE[j] == SUB[j-P] + 1.
SUB_2T_A0 = {0: 1, 1: 0, 2: 0, 3: 1, 4: 1, 5: 0}
SUB_2O_A0 = {0: 1, 1: 0, 2: 0, 3: 0, 4: 1, 5: 0,
             6: 1, 7: 0, 8: 1, 9: 1, 10: 1, 11: 0}
SUB_2O_A1 = {0: 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0,
             6: 1, 7: 1, 8: 0, 9: 1, 10: 1, 11: 1}

# Count of levels with m ≥ 2 in first above-threshold period.
COUNT_GEQ2_2T_A0 = 3   # j=6,9,10
COUNT_GEQ2_2O_A0 = 6   # j=12,16,18,20,21,22
COUNT_GEQ2_2O_A1 = 6   # j=15,18,19,21,22,23
COUNT_GEQ2_2I_A0 = 15  # 30,36,40,42,45,46,48,50,51,52,54,55,56,57,58


# ============================================================
# 2O tests — A₀ and A₁
# ============================================================

class TestAboveThreshold2O:
    """Above-threshold breakdown for the octahedral group."""

    @pytest.fixture(scope='class')
    def grp(self):
        g = prepare_polyhedral('2O', build_binary_octahedral, 48)
        comm = compute_commutator_subgroup(g['elements'])
        comm_keys = set(comm.keys())
        g['comm_keys'] = comm_keys
        return g

    def _chi_bar_A1(self, c, comm_keys):
        return 1.0 if qkey(c['rep']) in comm_keys else -1.0

    def test_commutator_subgroup_order(self, grp):
        """[2O,2O] = 2T has order 24. Guards A₁ construction."""
        assert len(grp['comm_keys']) == 24

    def test_j_star(self, grp):
        assert grp['j_star'] == 11

    def test_period(self, grp):
        assert grp['P'] == 12

    def test_jstar_plus_1_equals_P(self, grp):
        """j*+1 = P, so m(j*+1) = m(0) + 1 = 2 by periodicity."""
        assert grp['j_star'] + 1 == grp['P']

    def test_m_A0_jstar_plus_1_equals_2(self, grp):
        """Binary occupation breaks immediately: m(A₀, 12) = 2."""
        m = m_chi(12, grp['classes'], 48)
        assert m == 2

    def test_m_A0_jstar_zero(self, grp):
        """m(A₀, j*) = 0 (last zero)."""
        m = m_chi(11, grp['classes'], 48)
        assert m == 0

    def test_hardcoded_table_length(self, grp):
        assert len(ABOVE_2O_A0) == grp['P']
        assert len(ABOVE_2O_A1) == grp['P']

    @pytest.mark.parametrize('j', range(12, 24))
    def test_m_A0_above_pinned(self, grp, j):
        """Pin each above-threshold A₀ value."""
        m = m_chi(j, grp['classes'], 48)
        assert m == ABOVE_2O_A0[j]

    @pytest.mark.parametrize('j', range(12, 24))
    def test_m_A1_above_pinned(self, grp, j):
        """Pin each above-threshold A₁ value."""
        chi_bar = lambda c: self._chi_bar_A1(c, grp['comm_keys'])
        m = m_chi(j, grp['classes'], 48, chi_bar)
        assert m == ABOVE_2O_A1[j]

    @pytest.mark.parametrize('j', range(12, 24))
    def test_m_A0_geq_1(self, grp, j):
        """All m(A₀, j) ≥ 1 above threshold."""
        m = m_chi(j, grp['classes'], 48)
        assert m >= 1

    @pytest.mark.parametrize('j', range(0, 24))
    def test_periodicity_A0(self, grp, j):
        """m(A₀, j+P) = m(A₀, j) + 1."""
        m = m_chi(j, grp['classes'], 48)
        m_shifted = m_chi(j + grp['P'], grp['classes'], 48)
        assert m_shifted == m + 1

    def test_count_geq2_A0(self, grp):
        """Exactly 6/12 above-threshold A₀ levels have m ≥ 2."""
        count = sum(1 for j in range(12, 24)
                    if m_chi(j, grp['classes'], 48) >= 2)
        assert count == COUNT_GEQ2_2O_A0

    def test_count_geq2_A1(self, grp):
        """Exactly 6/12 above-threshold A₁ levels have m ≥ 2."""
        chi_bar = lambda c: self._chi_bar_A1(c, grp['comm_keys'])
        count = sum(1 for j in range(12, 24)
                    if m_chi(j, grp['classes'], 48, chi_bar) >= 2)
        assert count == COUNT_GEQ2_2O_A1

    def test_hardcoded_self_consistency_A0(self):
        """Hardcoded above = hardcoded sub + 1 (periodicity cross-check)."""
        for j in range(12, 24):
            assert ABOVE_2O_A0[j] == SUB_2O_A0[j - 12] + 1, f"j={j}"

    def test_hardcoded_self_consistency_A1(self):
        """Hardcoded above A₁ = hardcoded sub A₁ + 1."""
        for j in range(12, 24):
            assert ABOVE_2O_A1[j] == SUB_2O_A1[j - 12] + 1, f"j={j}"

    @pytest.mark.parametrize('j', range(12, 24))
    def test_extended_duality_sum_A0(self, grp, j):
        """Above threshold: m(j) + m(2j*+1-j) = 2."""
        dual = 2 * grp['j_star'] + 1 - j
        m1 = m_chi(j, grp['classes'], 48)
        m2 = m_chi(dual, grp['classes'], 48)
        assert m1 + m2 == 2

    @pytest.mark.parametrize('j', range(12, 24))
    def test_extended_duality_sum_A1(self, grp, j):
        """Above threshold: m_A1(j) + m_A1(2j*+1-j) = 2."""
        chi_bar = lambda c: self._chi_bar_A1(c, grp['comm_keys'])
        dual = 2 * grp['j_star'] + 1 - j
        m1 = m_chi(j, grp['classes'], 48, chi_bar)
        m2 = m_chi(dual, grp['classes'], 48, chi_bar)
        assert m1 + m2 == 2


# ============================================================
# 2T tests — A₀ only (κ=3, but A₁=A₂ by conjugacy)
# ============================================================

class TestAboveThreshold2T:
    """Above-threshold breakdown for the tetrahedral group."""

    @pytest.fixture(scope='class')
    def grp(self):
        return prepare_polyhedral('2T', build_binary_tetrahedral, 24)

    def test_j_star(self, grp):
        assert grp['j_star'] == 5

    def test_period(self, grp):
        assert grp['P'] == 6

    def test_jstar_plus_1_equals_P(self, grp):
        """j*+1 = P, so m(j*+1) = m(0) + 1 = 2 by periodicity."""
        assert grp['j_star'] + 1 == grp['P']

    def test_m_A0_jstar_plus_1_equals_2(self, grp):
        """Binary occupation breaks immediately: m(A₀, 6) = 2."""
        m = m_chi(6, grp['classes'], 24)
        assert m == 2

    def test_m_A0_jstar_zero(self, grp):
        m = m_chi(5, grp['classes'], 24)
        assert m == 0

    def test_hardcoded_table_length(self, grp):
        assert len(ABOVE_2T_A0) == grp['P']

    @pytest.mark.parametrize('j', range(6, 12))
    def test_m_A0_above_pinned(self, grp, j):
        """Pin each above-threshold A₀ value."""
        m = m_chi(j, grp['classes'], 24)
        assert m == ABOVE_2T_A0[j]

    @pytest.mark.parametrize('j', range(6, 12))
    def test_m_A0_geq_1(self, grp, j):
        m = m_chi(j, grp['classes'], 24)
        assert m >= 1

    @pytest.mark.parametrize('j', range(0, 12))
    def test_periodicity_A0(self, grp, j):
        m = m_chi(j, grp['classes'], 24)
        m_shifted = m_chi(j + grp['P'], grp['classes'], 24)
        assert m_shifted == m + 1

    def test_count_geq2_A0(self, grp):
        """Exactly 3/6 above-threshold A₀ levels have m ≥ 2."""
        count = sum(1 for j in range(6, 12)
                    if m_chi(j, grp['classes'], 24) >= 2)
        assert count == COUNT_GEQ2_2T_A0

    def test_hardcoded_self_consistency_A0(self):
        """Hardcoded above = hardcoded sub + 1."""
        for j in range(6, 12):
            assert ABOVE_2T_A0[j] == SUB_2T_A0[j - 6] + 1, f"j={j}"

    @pytest.mark.parametrize('j', range(6, 12))
    def test_extended_duality_sum_A0(self, grp, j):
        """Above threshold: m(j) + m(2j*+1-j) = 2."""
        dual = 2 * grp['j_star'] + 1 - j
        m1 = m_chi(j, grp['classes'], 24)
        m2 = m_chi(dual, grp['classes'], 24)
        assert m1 + m2 == 2


# ============================================================
# 2I tests — A₀ only (κ=1)
# ============================================================

class TestAboveThreshold2I:
    """Above-threshold breakdown for the icosahedral group."""

    @pytest.fixture(scope='class')
    def grp(self):
        return prepare_polyhedral('2I', build_binary_icosahedral, 120)

    def test_j_star(self, grp):
        assert grp['j_star'] == 29

    def test_period(self, grp):
        assert grp['P'] == 30

    def test_jstar_plus_1_equals_P(self, grp):
        """j*+1 = P, so m(j*+1) = m(0) + 1 = 2 by periodicity."""
        assert grp['j_star'] + 1 == grp['P']

    def test_m_A0_jstar_plus_1_equals_2(self, grp):
        """Binary occupation breaks immediately: m(A₀, 30) = 2."""
        m = m_chi(30, grp['classes'], 120)
        assert m == 2

    def test_m_A0_jstar_zero(self, grp):
        m = m_chi(29, grp['classes'], 120)
        assert m == 0

    def test_hardcoded_table_length(self, grp):
        assert len(ABOVE_2I_A0) == grp['P']

    @pytest.mark.parametrize('j', range(30, 60))
    def test_m_A0_above_pinned(self, grp, j):
        """Pin each above-threshold A₀ value."""
        m = m_chi(j, grp['classes'], 120)
        assert m == ABOVE_2I_A0[j]

    @pytest.mark.parametrize('j', range(30, 60))
    def test_m_A0_geq_1(self, grp, j):
        m = m_chi(j, grp['classes'], 120)
        assert m >= 1

    @pytest.mark.parametrize('j', range(0, 30))
    def test_periodicity_A0(self, grp, j):
        m = m_chi(j, grp['classes'], 120)
        m_shifted = m_chi(j + grp['P'], grp['classes'], 120)
        assert m_shifted == m + 1

    def test_count_geq2_A0(self, grp):
        """Exactly 15/30 above-threshold A₀ levels have m ≥ 2."""
        count = sum(1 for j in range(30, 60)
                    if m_chi(j, grp['classes'], 120) >= 2)
        assert count == COUNT_GEQ2_2I_A0

    @pytest.mark.parametrize('j', range(30, 60))
    def test_extended_duality_sum_A0(self, grp, j):
        """Above threshold: m(j) + m(2j*+1-j) = 2."""
        dual = 2 * grp['j_star'] + 1 - j
        m1 = m_chi(j, grp['classes'], 120)
        m2 = m_chi(dual, grp['classes'], 120)
        assert m1 + m2 == 2


# ============================================================
# Cross-group: universal pattern m(A₀, j*+1) = 2
# ============================================================

class TestUniversalThresholdSharpness:
    """m(A₀, j*+1) = 2 for all three polyhedral groups."""

    @pytest.fixture(params=[
        ('2T', build_binary_tetrahedral, 24),
        ('2O', build_binary_octahedral, 48),
        ('2I', build_binary_icosahedral, 120),
    ], scope='class')
    def grp(self, request):
        name, builder, G_order = request.param
        return prepare_polyhedral(name, builder, G_order)

    def test_jstar_plus_1_equals_P(self, grp):
        """j*+1 = P by definition: j* = |H|/2-1, P = |H|/2."""
        assert grp['j_star'] + 1 == grp['P']

    def test_m_jstar_plus_1_equals_2(self, grp):
        """m(A₀, j*+1) = m(A₀, 0)+1 = 2: threshold is maximally sharp."""
        m = m_chi(grp['j_star'] + 1, grp['classes'], grp['G'])
        assert m == 2


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
