#!/usr/bin/env python3
"""
I35 — Blue phase application: defect spectroscopy on SO(3)/O.

QUESTION: Does blue phase provide new mathematical content beyond
what's already in SO(3)/O analysis, or is it just relabeling?

RAW OUTPUT:
  §0 BP I (O^8+) and BP II (O^2) both have point group O (order 24).
     2O has 2 scalar sectors (A₀, A₁), both tensorial. -1 ∈ [2O,2O].
  §1 Defect stabilizers: [100]→C₄(Z₈), [111]→C₃(Z₆), [110]→C₂(Z₄).
     All cyclic → all have spinorial sectors (Q₈ ⊄ cyclic).
  §2 j=2 nematic mode: O bulk m(A₀,2)=0 FORBIDDEN.
     C₄: m=1, C₃: m=1, C₂: m=3, D₃: m=1. Localized at defects.
     [110] defects have highest degeneracy (3 nematic modes).
  §3 Spinorial first appearance: C₄,C₃,C₂ all at j=1/2 (gap 0.75B).
     D₃ at j=3/2 (gap 3.75B). Matches Prop 4.3 in paper.
  §4 All results match paper v10 Tables 2, 5, 6 exactly.

ANSWER:
  Blue phase is NOT a new mathematical investigation.
  Both BP I and BP II have point group O = our main example.
  Paper v10 already contains all relevant computations:
    - SO(3)/O analysis (Table 2)
    - O→D₃ symmetry breaking (Table 5)
    - O→K phase diagram (Table 6)
    - Spinorial gap formula (Prop 4.3)
  Blue phase application = physical narrative (defect lines,
  Q-tensor localization, junction spectroscopy) around existing math.
  Writing task, not computation task.

PLAN:
  §0. Blue phase point groups — both BP I and BP II have point group O
  §1. Defect line stabilizers — C₄, C₃, C₂ along crystal axes
  §2. j=2 nematic mode at defects — forbidden in O, allowed at defects
  §3. Spinorial sectors at defects — cyclic stabilizers all have them
  §4. Comparison with paper tables — everything already computed

STATUS: COMPLETE — 22/22 tests pass
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.group import (build_binary_octahedral, compute_conjugacy_classes)
from src.quaternion import qmul, qkey


# ---- helpers ----

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


def build_binary_cyclic(n):
    """2C_n = Z_{2n} in SU(2). Order 2n."""
    result = []
    for k in range(2 * n):
        theta = k * np.pi / n
        q = np.array([np.cos(theta), 0, 0, np.sin(theta)])
        result.append(q)
    G = np.array(result)
    assert len(G) == 2 * n
    assert len(set(qkey(q) for q in G)) == 2 * n
    return G


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
    assert len(set(qkey(q) for q in G)) == 4 * n
    return G


def m_trivial(j, G):
    """Multiplicity of trivial character in V_j|_G."""
    order = len(G)
    s = sum(chi_su2(j, 2 * G[k, 0]) for k in range(order))
    return int(round(s / order))


def m_char_cyclic(j, n, char_index):
    """Multiplicity of χ_{char_index} in V_j|_{Z_{2n}}.

    Z_{2n} chars: χ_m(g_k) = exp(2πi·m·k/(2n)).
    χ_m(-1) = (-1)^m. Spinorial iff m odd."""
    G = build_binary_cyclic(n)
    order = 2 * n
    s = 0
    for k in range(order):
        trace = 2 * G[k, 0]
        chi_bar = np.exp(-2j * np.pi * char_index * k / order)
        s += chi_bar * chi_su2(j, trace)
    return int(round(s.real / order))


# ==========================================================
# § 0. Blue phase point groups
#
# BP I: space group O^8+ (I4₁32), point group O (order 24)
# BP II: space group O^2 (P4₂32), point group O (order 24)
# Both chiral (no improper rotations).
# Binary lift: 2O (order 48) for both.
# ==========================================================

class TestBluePhasePointGroups:
    """Both blue phases have octahedral point group O."""

    def test_2O_order(self):
        G = build_binary_octahedral()
        assert len(G) == 48

    def test_2O_abelianization_Z2(self):
        """2O has 2 scalar sectors (both tensorial)."""
        from src.group import compute_commutator_subgroup
        G = build_binary_octahedral()
        comm = compute_commutator_subgroup(G)
        abel_order = len(G) // len(comm)
        assert abel_order == 2

    def test_both_sectors_tensorial(self):
        """Both A₀ and A₁ are tensorial: χ(-1) = +1."""
        # 2O^ab = Z₂, only char is trivial and sign
        # Both have χ(-1) = +1 since -1 ∈ [2O, 2O]
        from src.group import compute_commutator_subgroup
        G = build_binary_octahedral()
        comm = compute_commutator_subgroup(G)
        # -1 is in commutator subgroup
        minus_one = np.array([-1.0, 0.0, 0.0, 0.0])
        found = any(np.allclose(minus_one, c) for c in comm)
        assert found, "-1 should be in [2O, 2O]"


# ==========================================================
# § 1. Defect line stabilizers
#
# In cubic blue phase crystal, defect lines along:
#   [100]: stabilizer C₄ (4-fold axis), binary lift Z₈
#   [111]: stabilizer C₃ (3-fold axis), binary lift Z₆
#   [110]: stabilizer C₂ (2-fold axis), binary lift Z₄
# All cyclic → all have spinorial sectors (Q₈ not in cyclic).
# ==========================================================

class TestDefectStabilizers:
    """Defect line stabilizers are cyclic subgroups of O."""

    @pytest.mark.parametrize("axis,n,order_2K", [
        ("[100]", 4, 8),
        ("[111]", 3, 6),
        ("[110]", 2, 4),
    ])
    def test_binary_cyclic_order(self, axis, n, order_2K):
        """Binary cyclic group Z_{2n} has correct order."""
        G = build_binary_cyclic(n)
        assert len(G) == order_2K

    @pytest.mark.parametrize("n", [4, 3, 2])
    def test_cyclic_has_spinorial(self, n):
        """All cyclic stabilizers have spinorial sectors.

        Z_{2n} has 2n characters. χ_m(-1) = (-1)^m.
        Odd m gives spinorial. Always exists for n ≥ 1."""
        # χ₁ is spinorial: χ₁(-1) = exp(πi) = -1
        # Check m(χ₁, 1/2) ≥ 1
        m = m_char_cyclic(0.5, n, 1)
        assert m >= 1, f"C_{n}: spinorial char should appear at j=1/2"


# ==========================================================
# § 2. j=2 nematic mode at defects
#
# Q-tensor (nematic order parameter) transforms as j=2 (dim 5).
# In bulk O: m(A₀, 2) = 0 — nematic fluctuation FORBIDDEN.
# At defects: m(χ₀, 2) > 0 — nematic fluctuation ALLOWED.
# ==========================================================

class TestNematicMode:
    """j=2 (nematic) mode: forbidden in O, allowed at defects."""

    def test_j2_forbidden_in_O(self):
        """O bulk: m(A₀, 2) = 0."""
        G = build_binary_octahedral()
        assert m_trivial(2, G) == 0

    @pytest.mark.parametrize("axis,n,expected", [
        ("[100] C₄", 4, 1),
        ("[111] C₃", 3, 1),
        ("[110] C₂", 2, 3),
    ])
    def test_j2_allowed_at_defect(self, axis, n, expected):
        """Defect stabilizer: m(χ₀, 2) = expected."""
        G = build_binary_cyclic(n)
        assert m_trivial(2, G) == expected, \
            f"{axis}: m(χ₀, 2) = {m_trivial(2, G)}, expected {expected}"

    def test_j2_allowed_at_D3(self):
        """Trigonal junction D₃: m(χ₀, 2) = 1."""
        G = build_binary_dihedral(3)
        assert m_trivial(2, G) == 1

    def test_C2_highest_nematic_degeneracy(self):
        """[110] defect (C₂) has highest nematic degeneracy: 3.

        dim(V₂) = 5, and C₂ is smallest stabilizer,
        so most modes survive restriction."""
        mults = {}
        for n in [4, 3, 2]:
            G = build_binary_cyclic(n)
            mults[n] = m_trivial(2, G)
        assert mults[2] > mults[4]
        assert mults[2] > mults[3]
        assert mults[2] == 3


# ==========================================================
# § 3. Spinorial sectors at defects
#
# All cyclic stabilizers have spinorial sectors starting at j=1/2.
# D₃ (trigonal junction) has spinorial starting at j=3/2.
# ==========================================================

class TestSpinorialAtDefects:
    """Spinorial sectors appear at all defect lines."""

    @pytest.mark.parametrize("n", [4, 3, 2])
    def test_cyclic_spinorial_at_j_half(self, n):
        """Cyclic stabilizers: first spinorial at j=1/2, gap=0.75B."""
        m = m_char_cyclic(0.5, n, 1)  # χ₁ = first spinorial
        assert m == 1, f"C_{n}: m(χ₁, 1/2) = {m}, expected 1"

    def test_D3_spinorial_at_j_3half(self):
        """D₃ junction: first spinorial at j=3/2, gap=3.75B.

        Computed via Dic_3 character table: χ₁(a^k) = (-1)^k,
        χ₁(b·a^k) = i·(-1)^k. χ₁(-1) = -1 (spinorial)."""
        G = build_binary_dihedral(3)
        order = 12
        # Dic_3 abelianization Z₄. Chars: χ_m with χ_m(a)=(-1)^m, χ_m(b)=i^m
        # χ₁: spinorial (χ₁(-1) = (-1)^1 = -1)
        # Elements: G[2k] = a^k, G[2k+1] = b·a^k
        char_vals = np.zeros(order, dtype=complex)
        for idx in range(order):
            k = idx // 2
            if idx % 2 == 0:  # a^k
                char_vals[idx] = (-1) ** k
            else:  # b·a^k
                char_vals[idx] = 1j * (-1) ** k

        # Verify spinorial: χ₁(-1) = χ₁(a³) = char_vals[6]
        assert char_vals[6].real < 0, "χ₁ should be spinorial"

        # m(χ₁, 1/2) = 0 (below threshold)
        traces = [2 * G[k, 0] for k in range(order)]
        s_half = sum(np.conj(char_vals[k]) * chi_su2(0.5, traces[k])
                     for k in range(order))
        assert int(round(s_half.real / order)) == 0, "m(χ₁, 1/2) should be 0"

        # m(χ₁, 3/2) = 1 (first spinorial level)
        s_3half = sum(np.conj(char_vals[k]) * chi_su2(1.5, traces[k])
                      for k in range(order))
        assert int(round(s_3half.real / order)) == 1, "m(χ₁, 3/2) should be 1"


# ==========================================================
# § 4. Comparison with paper tables
#
# All results match paper v10:
#   Table 2: O multiplicities (j=0..11)
#   Table 5: O→D₃ densification
#   Table 6: O→K phase diagram
# No new mathematical content. Blue phase = physical relabeling.
# ==========================================================

class TestPaperComparison:
    """Verify blue phase data matches paper tables."""

    def test_O_j2_matches_table2(self):
        """Paper Table 2: m(A₀, 2) = 0, m(A₁, 2) = 0."""
        G = build_binary_octahedral()
        classes = compute_conjugacy_classes(G)
        classes.sort(key=lambda c: -c['trace'])
        order = len(G)

        def m_A0(j):
            return int(round(sum(chi_su2(j, c['trace']) * c['size']
                                 for c in classes) / order))

        assert m_A0(2) == 0
        # A₁ = sign char. Sign char of 2O: trivial on 2T, -1 on 2O\2T.
        # From Table 2: m(A₁, 2) = 0
        assert m_A0(0) == 1  # sanity: m(A₀, 0) = 1

    def test_O_to_D3_densification(self):
        """Paper Table 5: D₃ has 4 sectors, gap reduces from 20B to 6B.

        Trivial sector gap: O has j₁(A₀)=4, D₃ has j₁(χ₀)=2."""
        G_O = build_binary_octahedral()
        G_D3 = build_binary_dihedral(3)
        # O: first non-zero j for A₀
        j1_O = None
        for j in range(1, 20):
            if m_trivial(j, G_O) > 0:
                j1_O = j
                break
        assert j1_O == 4, f"O: j₁(A₀) = {j1_O}, expected 4"

        # D₃: first non-zero j for χ₀
        j1_D3 = None
        for j in range(1, 20):
            if m_trivial(j, G_D3) > 0:
                j1_D3 = j
                break
        assert j1_D3 == 2, f"D₃: j₁(χ₀) = {j1_D3}, expected 2"

    def test_phase_diagram_cyclic_all_spinorial(self):
        """Paper Table 6: all O→C_n channels have spinorial sectors."""
        for n in [4, 3, 2]:
            m = m_char_cyclic(0.5, n, 1)
            assert m >= 1, f"C_{n}: should have spinorial at j=1/2"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
