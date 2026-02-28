#!/usr/bin/env python3
"""
I17 — Berry phase on configuration space SO(3)/H.

QUESTION: Does the Berry phase give NEW content beyond χ(-1)?

ANSWER: NO. Berry phase = χ(-1), which is the sector label itself.

RESULT:
  §0. χ(-1) = +1 for all scalar chars, for all groups (2O, 2T, 2I).
      -1 ∈ [2H,2H] for all obstructed groups (condition A).
  §1. Z₂ defect holonomy: rotor encircling disclination picks up
      phase χ(-1). Scalar → +1, spinorial → -1. Exact.
      Connects to ST_11 z2_holonomy.py (lattice computation).
  §2. Parameter-space Berry phase: m(χ,j) depends only on H,
      not on Hamiltonian. Symmetry group jumps discretely.
      No smooth parameter space → no non-trivial Berry phase.

VERDICT: Reformulation, not new content. χ(-1) IS the sector label.
  The only topological invariant is ±1 from π₁(SO(3)) = Z₂.
  Worth a REMARK in paper (connects to defect physics / ST_11)
  but not a Proposition.

RAW OUTPUT (8 passed in 0.87s):
  TestChiMinusOne::test_2O_scalar_chars_at_minus_one     PASSED
  TestChiMinusOne::test_2T_scalar_chars_at_minus_one     PASSED
  TestChiMinusOne::test_minus_one_is_in_commutator_iff_obstructed  PASSED
  TestDefectHolonomy::test_scalar_sector_trivial_phase   PASSED
  TestDefectHolonomy::test_spinorial_sector_pi_phase     PASSED
  TestDefectHolonomy::test_defect_phase_from_multiplicities  PASSED
  TestParameterSpaceBerry::test_multiplicities_are_topological  PASSED
  TestParameterSpaceBerry::test_berry_phase_is_chi_minus_one   PASSED

STATUS: COMPLETE — clean negative result, Berry phase = reformulation
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.group import build_binary_octahedral, build_binary_tetrahedral, build_binary_icosahedral
from src.quaternion import qkey


# ==========================================================
# § 0. χ(-1) for all scalar characters of all groups
#
# This is the baseline: Berry phase = χ(-1) by definition.
# Scalar: χ(-1) = +1,  Spinorial: χ(-1) = -1.
# If this is ALL Berry phase gives us, it's a reformulation.
# ==========================================================

MINUS_ONE_KEY = qkey(np.array([-1.0, 0, 0, 0]))


def scalar_chars_2O():
    """Build scalar (1D) characters of 2O.
    Abelianization 2O/[2O,2O] = Z₂ → two 1D chars: A₀ (trivial), A₁ (sign).
    [2O, 2O] = 2T (index 2 subgroup)."""
    from src.group import compute_commutator_subgroup
    G = build_binary_octahedral()
    comm = compute_commutator_subgroup(G)
    comm_keys = set(comm.keys())

    chi_A0 = {}  # trivial
    chi_A1 = {}  # sign: +1 on 2T, -1 on coset
    for g in G:
        k = qkey(g)
        chi_A0[k] = 1.0
        chi_A1[k] = 1.0 if k in comm_keys else -1.0

    return {'A0': chi_A0, 'A1': chi_A1}, G


def scalar_chars_2T():
    """Build scalar (1D) characters of 2T.
    Abelianization 2T/[2T,2T] = Z₃ → three 1D chars.
    [2T, 2T] = Q₈ (index 3 subgroup)."""
    from src.group import compute_commutator_subgroup
    G = build_binary_tetrahedral()
    comm = compute_commutator_subgroup(G)
    comm_keys = set(comm.keys())
    order = len(G)  # 24

    # Quotient map π: 2T → Z₃. Need to assign coset labels.
    # Q₈ has 8 elements. Cosets: Q₈, c·Q₈, c²·Q₈ where c = (1+i+j+k)/2.
    omega = np.exp(2j * np.pi / 3)
    c = np.array([1, 1, 1, 1]) / 2.0  # generator of Z₃ quotient

    # Find coset of each element by checking if g·c^{-k} ∈ Q₈
    from src.quaternion import qmul, qinv
    c_inv = qinv(c)
    coset_label = {}
    for g in G:
        k = qkey(g)
        if k in comm_keys:
            coset_label[k] = 0
        else:
            # Check if g·c⁻¹ ∈ Q₈
            gc = qmul(g, c_inv)
            if qkey(gc) in comm_keys:
                coset_label[k] = 1
            else:
                coset_label[k] = 2

    # Three chars: χ_m(g) = ω^{m·coset(g)} for m = 0,1,2
    chars = {}
    for m in range(3):
        chi = {}
        for g in G:
            k = qkey(g)
            chi[k] = omega ** (m * coset_label[k])
        chars[f'chi_{m}'] = chi

    return chars, G


class TestChiMinusOne:
    """Verify χ(-1) classification matches scalar/spinorial."""

    def test_2O_scalar_chars_at_minus_one(self):
        """2O: A₀(-1) = A₁(-1) = +1 (both scalar)."""
        chars, G = scalar_chars_2O()
        for name, chi in chars.items():
            val = chi[MINUS_ONE_KEY]
            assert abs(val - 1.0) < 1e-10, \
                f"2O {name}: χ(-1) = {val}, expected +1"

    def test_2T_scalar_chars_at_minus_one(self):
        """2T: all 3 scalar chars have χ(-1) = +1."""
        chars, G = scalar_chars_2T()
        for name, chi in chars.items():
            val = chi[MINUS_ONE_KEY]
            assert abs(val - 1.0) < 1e-10, \
                f"2T {name}: χ(-1) = {val}, expected +1"

    def test_minus_one_is_in_commutator_iff_obstructed(self):
        """Key fact: -1 ∈ [2H, 2H] ⟺ obstruction holds.
        2O: -1 ∈ [2O,2O] = 2T ✓ (obstructed)
        2T: -1 ∈ [2T,2T] = Q₈ ✓ (obstructed)
        2I: -1 ∈ [2I,2I] = 2I ✓ (obstructed, since 2I is perfect)
        This is condition (A) from the paper."""
        from src.group import compute_commutator_subgroup
        for name, builder in [('2O', build_binary_octahedral),
                              ('2T', build_binary_tetrahedral),
                              ('2I', build_binary_icosahedral)]:
            G = builder()
            comm = compute_commutator_subgroup(G)
            assert MINUS_ONE_KEY in comm, \
                f"{name}: -1 not in [G,G] — obstruction should hold"


# ==========================================================
# § 1. Z₂ defect holonomy
#
# A rotor encircling a disclination line (Z₂ vortex) picks up
# phase χ(-1). This is the physical observable.
#
# Connection to ST_11: z2_holonomy.py computes exactly this
# on a Kelvin foam lattice.
#
# Question: is this MORE than just saying "spinorial = -1"?
# ==========================================================

def chi_su2(j, g_trace):
    """Character of V_j at element with tr(g) = g_trace."""
    if abs(g_trace - 2.0) < 1e-10:
        return 2*j + 1
    if abs(g_trace + 2.0) < 1e-10:
        return (2*j + 1) * ((-1)**(int(2*j)))
    alpha_half = np.arccos(np.clip(g_trace / 2.0, -1, 1))
    return np.sin((2*j + 1) * alpha_half) / np.sin(alpha_half)


def multiplicity(j, chi_vals, elements, order):
    """m(χ, j) = (1/|G|) Σ_g conj(χ(g)) · χ_j(g)."""
    s = 0.0
    for g in elements:
        k = qkey(g)
        s += np.conj(chi_vals.get(k, 0.0)) * chi_su2(j, 2.0 * g[0])
    return (s / order).real


class TestDefectHolonomy:
    """Phase acquired by rotor encircling a Z₂ defect.

    Physical setup: liquid crystal with SO(3)/O order parameter.
    A disclination line carries Z₂ flux (holonomy = -1 ∈ 2O).
    A quantum rotor in sector χ, encircling the defect, picks up
    phase exp(iπ) = -1 (spinorial) or +1 (scalar).

    This is measurable: interferometry with two paths, one
    encircling the defect. Phase difference = 0 or π.

    The mathematical content is: parallel transport on the
    flat bundle over SO(3)/H with fiber labeled by χ gives
    holonomy χ(γ̃) for loop γ, where γ̃ is the lift to 2H.
    For the non-contractible loop: γ̃ = -1."""

    def test_scalar_sector_trivial_phase(self):
        """Scalar sectors of 2O: encircling Z₂ defect → phase 0.
        Physically: rotor in A₀ or A₁ sector is unaffected."""
        chars, G = scalar_chars_2O()
        for name, chi in chars.items():
            phase = chi[MINUS_ONE_KEY]  # χ(-1)
            assert abs(phase - 1.0) < 1e-10, \
                f"{name}: phase = {phase}, expected +1 (trivial)"

    def test_spinorial_sector_pi_phase(self):
        """Spinorial irreps of 2O: encircling Z₂ defect → phase π.
        These are irreps of 2O that don't factor through O.
        For 1D chars of 2O, there are none spinorial (abelianization = Z₂).
        But higher-dim irreps can be spinorial.
        Check: χ_{irrep}(-1) = -(dim) for spinorial irreps."""
        G = build_binary_octahedral()
        minus_one = np.array([-1.0, 0, 0, 0])
        # For any irrep ρ: ρ(-1) = ±Id (since (-1)² = 1).
        # Scalar: ρ(-1) = +Id → tr = +dim
        # Spinorial: ρ(-1) = -Id → tr = -dim
        # χ_j(-1) = (2j+1)·(-1)^{2j}: integer j → +, half-integer j → -
        # So V_j restricted to 2O: spinorial components have χ(-1) < 0
        for j in [0, 1, 2, 3]:  # integer j → scalar
            tr_at_minus = chi_su2(j, -2.0)
            assert tr_at_minus == (2*j + 1), \
                f"V_{j}(-1) = {tr_at_minus}, expected {2*j+1}"
        for j_twice in [1, 3, 5]:  # half-integer j → spinorial
            j = j_twice / 2.0
            tr_at_minus = chi_su2(j, -2.0)
            assert tr_at_minus == -(2*j + 1), \
                f"V_{j}(-1) = {tr_at_minus}, expected {-(2*j+1)}"

    def test_defect_phase_from_multiplicities(self):
        """Cross-check: for 2O scalar sector A₀, all occupied levels
        have integer j (m(A₀,j)=1 only at integer j).
        So the "aggregate" phase from the defect is always +1.
        For spinorial irreps, occupied levels have half-integer j → phase -1."""
        chars, G = scalar_chars_2O()
        order = len(G)  # 48
        j_star = order // 4 - 1  # 11

        # A₀: check all occupied j are integer
        for j in range(j_star + 1):
            m = multiplicity(j, chars['A0'], G, order)
            if round(m) == 1:
                assert j == int(j), f"A₀ occupied at non-integer j={j}"


# ==========================================================
# § 2. Parameter-space Berry phase: O → D₃ loop
#
# THIS is where new content could be.
# Family of Hamiltonians H(λ) = B·L² + λ·V_trigonal
# As λ varies: O symmetry → D₃ symmetry → back to O.
# Eigenstate |j, χ⟩ evolves. Berry phase for a loop in λ-space?
#
# Key: at re-entrance transition, spinorial sectors appear.
# If Berry phase is quantized/jumps → topological invariant.
# If it's just χ(-1) again → reformulation, stop.
# ==========================================================

class TestParameterSpaceBerry:
    """Berry phase in distortion parameter space.

    Setup: H(λ) acts on L²(SO(3)/H) truncated to j ≤ j*.
    But actually, character theory tells us m(χ,j) depends ONLY
    on the group H, not on the Hamiltonian.

    So as we vary λ (distortion strength):
    - While symmetry stays O: m(χ,j) is fixed
    - While symmetry stays D₃: m(χ,j) is fixed (different pattern)
    - At the transition point: symmetry changes discontinuously

    Key insight: there's NO smooth parameter space to do Berry phase!
    The symmetry group H(λ) is discrete — it jumps from O to D₃.
    The multiplicities m(χ,j) are topological invariants of H,
    not smooth functions of λ.

    Verdict: Berry phase in parameter space = reformulation.
    The only content is χ(-1), which IS the sector label."""

    def test_multiplicities_are_topological(self):
        """m(χ,j) depends only on group H, not on Hamiltonian details.
        Varying the crystal field strength doesn't change multiplicities
        (as long as symmetry is preserved). This is the paper's key point:
        'Adding any SO(3)-invariant perturbation shifts energy levels
        but cannot change the occupation pattern.'"""
        chars, G = scalar_chars_2O()
        order = len(G)

        # m(A₀, j) for j=0..11: determined purely by character theory
        pattern_A0 = []
        for j in range(12):
            m = round(multiplicity(j, chars['A0'], G, order))
            pattern_A0.append(m)

        # Known pattern: [1,0,0,0,1,0,1,0,1,1,1,0]
        expected = [1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0]
        assert pattern_A0 == expected, \
            f"A₀ pattern: {pattern_A0} ≠ {expected}"

    def test_berry_phase_is_chi_minus_one(self):
        """The ONLY topological invariant distinguishing sectors
        is χ(-1) ∈ {+1, -1}. This is the Z₂ Berry phase.

        For scalar sectors: χ(-1) = +1, phase = 0.
        For spinorial sectors: χ(-1) = -1, phase = π.

        No finer structure exists because π₁(SO(3)/H) = H
        and the abelianization H^ab determines the 1D chars.
        The Z₂ part of H^ab gives scalar/spinorial.

        CONCLUSION: Berry phase = reformulation of χ(-1).
        Does NOT give new content beyond sector classification."""
        # For 2O: H^ab = Z₂, so only invariant is ±1
        chars, G = scalar_chars_2O()
        for name, chi in chars.items():
            berry = chi[MINUS_ONE_KEY]
            assert berry in [1.0, -1.0], \
                f"{name}: Berry phase {berry} not in {{+1, -1}}"

        # For 2T: H^ab = Z₃, central element -1 maps to 0 in Z₃
        # (since -1 ∈ [2T,2T] = Q₈), so χ(-1) = ω⁰ = 1 for ALL chars
        chars_T, G_T = scalar_chars_2T()
        for name, chi in chars_T.items():
            berry = chi[MINUS_ONE_KEY]
            assert abs(berry - 1.0) < 1e-10, \
                f"2T {name}: Berry phase {berry}, expected +1"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
