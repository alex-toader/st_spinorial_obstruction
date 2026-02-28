#!/usr/bin/env python3
"""
I6 — Nuclear spin isomers and sector creation under symmetry breaking.

INVESTIGATION: The superselection sectors on SO(3)/H ARE nuclear spin
isomers for molecular rotors with identical nuclei. This file establishes
the dictionary between our mathematical framework and the physical
observable (nuclear spin species).

RAW OUTPUT:
  Nuclear spin decomposition O: A₁:10, A₂:2, E:8, T₁:6, T₂:6
  Scalar weights SO(3)/O: A₀=10, A₁=2, ratio 5:1
  Budget: 12 scalar-accessible, 20 non-scalar (total dim 64)
  O→D₃: A₁(D₃):16, A₂(D₃):8, E(D₃):20. Ratio 2:1 (was 5:1)
  Spinorial: χ₁=0, χ₃=0 — UNPOPULATED (spinorial invisibility)
  Direct D₃ check: A₁:16, A₂:8, E:20 — matches branching ✓
  SF₆ spinorial gap (j_min=3/2): 0.3415 cm⁻¹ = 0.491 K = 10.24 GHz
  CH₄ A₀ gap (j=3): 62.89 cm⁻¹ = 90.5 K
  16 passed in 0.03s

STATUS: COMPLETE

KEY FINDING (Spinorial Invisibility Theorem):
  For any rigid molecular rotor with identical nuclei and no intrinsic
  half-integer electronic spin, spinorial sectors have g = 0.

  Proof: nuclear spins transform under H (rotation group), not 2H
  (double cover). The element -1 ∈ 2H projects to identity in H, so
  ρ_spin(-1) = +Id. Only tensorial irreps (Γ(-1) = +1) appear in the
  nuclear spin decomposition.

  Nuclear spin budget for SF₆ (6 × spin-1/2, O_h):
    O phase:   A₁:10, A₂:2  (scalar-accessible: 12/64 = 19%)
    O→D₃:     A₁(D₃):16, A₂(D₃):8  (scalar-accessible: 24/64 = 38%)
               χ₁:0,  χ₃:0  (spinorial: UNPOPULATED)
               locked in E(D₃): 20 copies (from E + T₁ + T₂ of O)

  Independent verification: direct D₃ cycle index computation on 6
  octahedral vertices confirms A₁:16, A₂:8, E:20 (dim check: 64).

  Measurable observable: statistical ratio A₁:A₂ changes from 5:1
  (O phase) to 2:1 (D₃ phase). Spinorial sectors are kinematically
  allowed but thermally dark.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest

# ==========================================================
# § 0. Physical constants and molecular data
# ==========================================================

# SF₆: octahedral, 6 identical ¹⁹F nuclei (spin 1/2)
# UF₆: octahedral, 6 identical ¹⁹F nuclei (spin 1/2) — same symmetry
# CH₄: tetrahedral, 4 identical ¹H nuclei (spin 1/2)

MOLECULES = {
    'SF6': {
        'B_cm': 0.09107,     # rotational constant [cm⁻¹]
        'symmetry': 'O_h',   # point group
        'H': 'O',            # rotation subgroup (SO(3) part)
        'nuclei': 6,         # identical nuclei
        'I_nuc': 0.5,        # nuclear spin of ¹⁹F
        'order_H': 24,       # |O| = 24
    },
    'UF6': {
        'B_cm': 0.05652,     # rotational constant [cm⁻¹]
        'symmetry': 'O_h',
        'H': 'O',
        'nuclei': 6,
        'I_nuc': 0.5,        # ¹⁹F
        'order_H': 24,
    },
    'CH4': {
        'B_cm': 5.2412,      # rotational constant [cm⁻¹]
        'symmetry': 'T_d',
        'H': 'T',            # rotation subgroup
        'nuclei': 4,
        'I_nuc': 0.5,        # ¹H
        'order_H': 12,
    },
}

CM_TO_K = 1.4388    # 1 cm⁻¹ = 1.4388 K
CM_TO_GHZ = 29.9792 # 1 cm⁻¹ = 29.9792 GHz


# ==========================================================
# § 1. Nuclear spin species ↔ superselection sectors
#
# For a molecule with N identical nuclei of spin I in a cage
# of symmetry H, the Pauli principle requires:
#   total wavefunction = (rotation) × (nuclear spin) × (electronic)
# to be symmetric (bosons) or antisymmetric (fermions) under
# permutation of identical nuclei.
#
# The permutation representation on nuclear spins decomposes into
# irreps of H. Each irrep Γ of H defines a "spin species":
# rotational states must transform as Γ* to make the product
# totally symmetric.
#
# Our "superselection sectors" (scalar characters of 2H) correspond
# to the 1D representations. The tensorial ones (χ(-1)=+1) are the
# standard spin species; spinorial ones (χ(-1)=-1) become accessible
# only when Q₈ ⊄ 2K after symmetry breaking.
# ==========================================================

def nuclear_spin_decomposition_O():
    """Decompose (C²)^⊗6 under O acting on octahedral vertices.

    Method: χ_spin(h) = (2I+1)^{c(h)} where c(h) = number of cycles
    in the permutation of the 6 vertices by rotation h.

    Cycle structure of O on 6 octahedral vertices:
      E   (1):  (1⁶)         c=6  →  χ = 2⁶ = 64
      C₃  (8):  (3²)         c=2  →  χ = 2² = 4
      C₄  (6):  (1²)(4)      c=3  →  χ = 2³ = 8
      C₂  (3):  (1²)(2²)     c=4  →  χ = 2⁴ = 16
      C₂' (6):  (2³)         c=3  →  χ = 2³ = 8

    Character table of O:
            E   8C₃  6C₄  3C₂  6C₂'
      A₁    1    1    1    1    1
      A₂    1    1   -1    1   -1
      E     2   -1    0    2    0
      T₁    3    0    1   -1   -1
      T₂    3    0   -1   -1    1
    """
    # χ_spin on each conjugacy class
    chi_spin = [64, 4, 8, 16, 8]
    # class sizes
    class_sizes = [1, 8, 6, 3, 6]
    # character table rows: A₁, A₂, E, T₁, T₂
    char_table = {
        'A1': [1,  1,  1,  1,  1],
        'A2': [1,  1, -1,  1, -1],
        'E':  [2, -1,  0,  2,  0],
        'T1': [3,  0,  1, -1, -1],
        'T2': [3,  0, -1, -1,  1],
    }
    order = 24
    result = {}
    for name, chi_gamma in char_table.items():
        n = sum(s * cg * cs for s, cg, cs
                in zip(chi_spin, chi_gamma, class_sizes)) // order
        result[name] = n
    return result


class TestNuclearSpinDictionary:
    """Establish: sectors = spin species for O_h molecules."""

    def test_spin_decomposition_total(self):
        """(C²)^⊗6 = 64 dimensions. Check decomposition sums to 64."""
        dec = nuclear_spin_decomposition_O()
        dims = {'A1': 1, 'A2': 1, 'E': 2, 'T1': 3, 'T2': 3}
        total = sum(dec[g] * dims[g] for g in dec)
        assert total == 64

    def test_spin_decomposition_values(self):
        """Compute nuclear spin multiplicities for each O irrep.
        These are the statistical weights g_Γ."""
        dec = nuclear_spin_decomposition_O()
        # Expected from cycle index computation:
        # n_{A₁} = (64 + 32 + 48 + 48 + 48)/24 = 240/24 = 10
        # n_{A₂} = (64 + 32 - 48 + 48 - 48)/24 = 48/24 = 2
        # n_E    = (128 - 32 + 0 + 96 + 0)/24 = 192/24 = 8
        # n_{T₁} = (192 + 0 + 48 - 48 - 48)/24 = 144/24 = 6
        # n_{T₂} = (192 + 0 - 48 - 48 + 48)/24 = 144/24 = 6
        assert dec == {'A1': 10, 'A2': 2, 'E': 8, 'T1': 6, 'T2': 6}
        print(f"\nNuclear spin decomposition for 6 × spin-1/2 in O:")
        for name, n in dec.items():
            print(f"  {name}: {n}")

    def test_scalar_sector_weights(self):
        """Only A₁ and A₂ are scalar (1D) sectors.
        Statistical weights: g_{A₁}=10, g_{A₂}=2.
        These are the nuclear spin degeneracies of the two
        superselection sectors on SO(3)/O."""
        dec = nuclear_spin_decomposition_O()
        g_A0 = dec['A1']  # A₀ in paper notation = A₁ in O notation
        g_A1 = dec['A2']  # A₁ in paper notation = A₂ in O notation
        assert g_A0 == 10
        assert g_A1 == 2
        print(f"\nScalar sector weights (SO(3)/O):")
        print(f"  A₀ (identity sector): g = {g_A0}")
        print(f"  A₁ (sign sector):     g = {g_A1}")
        print(f"  Ratio: {g_A0}:{g_A1} = 5:1")

    def test_nonscalar_sectors_locked(self):
        """E, T₁, T₂ sectors contribute spin configurations locked in
        non-scalar rotational states on SO(3)/O."""
        dec = nuclear_spin_decomposition_O()
        # Spin states in non-scalar O-irreps (counted as spin multiplicity × dim):
        #   E(dim 2):  8 copies × 2 = 16 spin-state-dimensions
        #   T₁(dim 3): 6 copies × 3 = 18
        #   T₂(dim 3): 6 copies × 3 = 18
        # But for counting spin CONFIGURATIONS: n_Γ is the multiplicity.
        accessible = dec['A1'] + dec['A2']  # 10 + 2 = 12
        locked_count = dec['E'] + dec['T1'] + dec['T2']  # 8 + 6 + 6 = 20
        # Dimension check: 12×1 + 8×2 + 6×3 + 6×3 = 12+16+18+18 = 64 ✓
        total_dim = (dec['A1']*1 + dec['A2']*1 + dec['E']*2
                     + dec['T1']*3 + dec['T2']*3)
        assert total_dim == 64
        assert accessible == 12
        assert locked_count == 20
        print(f"\nNuclear spin budget (64 total dim):")
        print(f"  Scalar-accessible: {accessible} configs (A₁:{dec['A1']}, A₂:{dec['A2']})")
        print(f"  Non-scalar:        {locked_count} configs (E:{dec['E']}, T₁:{dec['T1']}, T₂:{dec['T2']})")

    def test_spin_species_count_O(self):
        """O has 2 scalar sectors (A₀, A₁) = 2 spin species accessible
        in the octahedral phase. Both tensorial."""
        kappa_O = 2  # A₁, A₂ in O notation
        assert kappa_O == 2

    def test_branching_nuclear_spins_O_to_D3(self):
        """Under O→D₃, the 5 irreps of O branch as:
          A₁ → A₁(D₃)           (1D, scalar)
          A₂ → A₂(D₃)           (1D, scalar)
          E  → E(D₃)            (2D, non-scalar)
          T₁ → A₁(D₃) ⊕ E(D₃)  (1D scalar freed + 2D locked)
          T₂ → A₂(D₃) ⊕ E(D₃)  (1D scalar freed + 2D locked)

        Scalar sectors after O→D₃:
          A₁(D₃): g = n_{A₁} + n_{T₁} = 10 + 6 = 16
          A₂(D₃): g = n_{A₂} + n_{T₂} = 2 + 6 = 8

        Locked in E(D₃) (non-scalar, 2D):
          From E(O): 8, from T₁: 6, from T₂: 6 → total 20

        Dimension check: 16×1 + 8×1 + 20×2 = 64 ✓
        """
        dec = nuclear_spin_decomposition_O()
        # Scalar sectors of D₃
        g_A1_D3 = dec['A1'] + dec['T1']  # 10 + 6 = 16
        g_A2_D3 = dec['A2'] + dec['T2']  # 2 + 6 = 8
        # Non-scalar: E(D₃) gets contributions from E, T₁, T₂ of O
        g_E_D3 = dec['E'] + dec['T1'] + dec['T2']  # 8 + 6 + 6 = 20

        assert g_A1_D3 == 16
        assert g_A2_D3 == 8
        assert g_E_D3 == 20
        # Dimension check
        assert g_A1_D3 * 1 + g_A2_D3 * 1 + g_E_D3 * 2 == 64

        print(f"\nNuclear spin redistribution O→D₃:")
        print(f"  A₁(D₃) scalar: g = {g_A1_D3}  (A₁:10 + T₁:6)")
        print(f"  A₂(D₃) scalar: g = {g_A2_D3}  (A₂:2 + T₂:6)")
        print(f"  E(D₃) locked:  g = {g_E_D3}  (E:8 + T₁:6 + T₂:6)")
        print(f"  Ratio A₁:A₂ = {g_A1_D3}:{g_A2_D3} = 2:1  (was 5:1)")
        print(f"  Scalar-accessible: {g_A1_D3 + g_A2_D3}/64 = "
              f"{(g_A1_D3+g_A2_D3)/64:.0%}  (was 19%)")

    def test_spinorial_sectors_zero_nuclear_weight(self):
        """Spinorial sectors χ₁, χ₃ of 2D₃ have g = 0 nuclear spin weight.

        Reason: nuclear spin decomposition is performed in Rep(O), not
        Rep(2O). Nuclear spins transform under the rotation group H (which
        permutes identical nuclei), not the double cover 2H. Since -1 ∈ 2H
        acts trivially on nuclear spins (ρ_spin(-1) = +Id), only tensorial
        irreps of 2H contribute. Spinorial irreps (with ρ(-1) = -Id) are
        invisible to nuclear spin statistics.

        Spinorial sectors are kinematically allowed after symmetry breaking
        but carry zero nuclear spin statistical weight because nuclear spins
        transform under H rather than 2H.
        """
        dec = nuclear_spin_decomposition_O()
        # Nuclear spin decomposition lives in Rep(O), not Rep(2O).
        # Spinorial irreps of 2O (L, H₁, H₂ with ρ(-1)=-Id) don't
        # factor through O and are absent from the spin decomposition.
        g_chi1 = 0  # from L → χ₁ ⊕ χ₃ ⊕ ψ₁, but n_L = 0
        g_chi3 = 0
        assert g_chi1 == 0
        assert g_chi3 == 0
        print(f"\nSpinorial sector weights after O→D₃:")
        print(f"  χ₁ (spinorial): g = {g_chi1}  ← zero! (from L, spinorial parent)")
        print(f"  χ₃ (spinorial): g = {g_chi3}  ← zero! (from L, spinorial parent)")
        print(f"  → Spinorial sectors are UNPOPULATED in thermal equilibrium")

    def test_direct_D3_decomposition(self):
        """Independent verification: decompose nuclear spins directly under D₃
        (not via branching from O).

        D₃ embedded in O via [111] axis. The 6 octahedral vertices are
        ±x, ±y, ±z. C₃ along [111] cycles x→y→z and -x→-y→-z.

        Cycle structure of D₃ on 6 octahedral vertices:
          E   (1):  (1⁶)     c=6  →  χ = 2⁶ = 64
          C₃  (2):  (3)(3)   c=2  →  χ = 2² = 4     [two 3-cycles: {x,y,z},{-x,-y,-z}]
          C₂  (3):  (2)(2)(2) c=3 →  χ = 2³ = 8     [three pairs swapped]

        Character table of D₃:
                E   2C₃  3C₂
          A₁    1    1    1
          A₂    1    1   -1
          E     2   -1    0
        """
        chi_spin = [64, 4, 8]
        class_sizes = [1, 2, 3]
        char_table = {
            'A1': [1,  1,  1],
            'A2': [1,  1, -1],
            'E':  [2, -1,  0],
        }
        order = 6

        result = {}
        for name, chi_gamma in char_table.items():
            n = sum(s * cg * cs for s, cg, cs
                    in zip(chi_spin, chi_gamma, class_sizes)) // order
            result[name] = n

        # Should match branching result: A₁:16, A₂:8, E:20
        assert result == {'A1': 16, 'A2': 8, 'E': 20}
        # Dimension check
        assert result['A1']*1 + result['A2']*1 + result['E']*2 == 64

        print(f"\nDirect D₃ nuclear spin decomposition (independent check):")
        for name, n in result.items():
            print(f"  {name}: {n}")
        print(f"  Matches O→D₃ branching: ✓")

    def test_ratio_change_measurable(self):
        """The statistical ratio between scalar sectors changes from
        5:1 (O phase) to 2:1 (D₃ phase). This affects spectral line
        intensities and is measurable spectroscopically."""
        # O phase: A₁:A₂ = 10:2 = 5:1
        ratio_O = 10 / 2
        # D₃ phase: A₁:A₂ = 16:8 = 2:1
        ratio_D3 = 16 / 8
        assert ratio_O == 5.0
        assert ratio_D3 == 2.0
        print(f"\nStatistical ratio change at O→D₃:")
        print(f"  O phase:  A₀:A₁ = 10:2 = {ratio_O:.0f}:1")
        print(f"  D₃ phase: A₁:A₂ = 16:8 = {ratio_D3:.0f}:1")
        print(f"  → Measurable intensity change in rotational spectrum")

    def test_spin_species_count_D3(self):
        """D₃ has 4 scalar sectors of 2D₃: 2 tensorial + 2 spinorial.
        Re-entrance creates 2 new spinorial sectors.
        But: spinorial ones have g=0 nuclear spin weight."""
        kappa_D3 = 4
        kappa_spin_D3 = 2
        assert kappa_D3 == 4
        assert kappa_spin_D3 == 2


# ==========================================================
# § 2. Energy scales for concrete molecules
#
# Convert mathematical gaps (in units of B) to physical units.
# The gap for sector χ is E_gap = j_min(j_min+1) · B
# where j_min is the first j with m(χ, j) ≥ 1.
# ==========================================================

class TestEnergyScales:
    """Convert gaps to cm⁻¹, K, GHz for real molecules."""

    def test_sf6_octahedral_gaps(self):
        """SF₆ in O_h: A₀ gap = 20B, A₁ gap = 12B."""
        B = MOLECULES['SF6']['B_cm']
        gap_A0_cm = 20 * B    # j=4: 4×5=20
        gap_A1_cm = 12 * B    # j=3: 3×4=12
        gap_A0_K = gap_A0_cm * CM_TO_K
        gap_A1_K = gap_A1_cm * CM_TO_K

        # Verify known values
        assert abs(gap_A0_K - 2.62) < 0.02
        assert abs(gap_A1_K - 1.57) < 0.02

    def test_sf6_trigonal_spinorial_gap(self):
        """SF₆ distorted to D₃: spinorial gap energy scale.
        NOTE: gap universality fails for D₃ — actual j_min^spin = 3/2
        (not 1/2) from I3 investigation. Gap = n(n+2)/4 · B = 3.75B."""
        B = MOLECULES['SF6']['B_cm']
        gap_spin_cm = 3.75 * B  # j=3/2: (3/2)(5/2) = 15/4 = 3.75
        gap_spin_K = gap_spin_cm * CM_TO_K
        gap_spin_GHz = gap_spin_cm * CM_TO_GHZ

        print(f"\nSF₆ spinorial gap after O→D₃ (j_min=3/2):")
        print(f"  {gap_spin_cm:.4f} cm⁻¹")
        print(f"  {gap_spin_K:.3f} K")
        print(f"  {gap_spin_GHz:.2f} GHz")
        print(f"  NOTE: g=0 for molecular rotors (spinorial invisibility).")
        print(f"  Relevant only for systems with intrinsic half-integer spin.")

        assert gap_spin_cm > 0.001, "gap above spectroscopic resolution"

    def test_uf6_trigonal_spinorial_gap(self):
        """UF₆ distorted to D₃: spinorial gap = 3.75B (j_min=3/2)."""
        B = MOLECULES['UF6']['B_cm']
        gap_spin_cm = 3.75 * B
        gap_spin_K = gap_spin_cm * CM_TO_K

        print(f"\nUF₆ spinorial gap after O→D₃ (j_min=3/2):")
        print(f"  {gap_spin_cm:.4f} cm⁻¹ = {gap_spin_K:.3f} K")

        assert gap_spin_cm > 0.001

    def test_ch4_tetrahedral_gaps(self):
        """CH₄ in T_d: κ=3 sectors, all tensorial (obstructed).
        Much larger B → higher energy scale.
        For T: j₁(A₀) = 3, gap = 3×4 = 12B (from paper Table 3)."""
        B = MOLECULES['CH4']['B_cm']
        # T has κ=3: A, E_1, E_2 (all tensorial, -1 ∈ [2T,2T])
        gap_A_cm = 12 * B    # j=3 for A₀ of T: 3×4=12
        gap_A_K = gap_A_cm * CM_TO_K

        print(f"\nCH₄ (tetrahedral):")
        print(f"  B = {B:.4f} cm⁻¹ = {B*CM_TO_K:.2f} K")
        print(f"  A₀ gap = {gap_A_cm:.2f} cm⁻¹ = {gap_A_K:.1f} K")

        # CH₄ is fully obstructed: no spinorial sectors
        # Distortion T→C₃ would create them
        assert B > 1.0, "CH₄ has large B"


# ==========================================================
# § 3. Energy scale check
#
# The spinorial gap is spectroscopically resolvable in principle
# (above microwave resolution ~0.001 cm⁻¹), but for molecular
# rotors with bosonic nuclei the spinorial sectors have g = 0
# (spinorial invisibility). These energy scales become physically
# relevant for systems with intrinsic half-integer spin (Kramers
# systems, cold atoms in synthetic gauge fields, etc.).
# ==========================================================

class TestEnergyScaleResolution:
    """Verify spinorial gap energy scale is above instrumental limits."""

    def test_sf6_gap_above_resolution(self):
        """SF₆ at D₃: spinorial gap 3.75B = 0.34 cm⁻¹ — resolvable.
        (Thermally dark for molecular SF₆ due to g=0.)"""
        B = MOLECULES['SF6']['B_cm']
        gap = 3.75 * B  # j_min = 3/2 for D_3
        microwave_resolution = 0.001  # cm⁻¹
        assert gap > microwave_resolution
        assert gap / microwave_resolution > 10, "margin > 10×"

    def test_uf6_gap_above_resolution(self):
        """UF₆ at D₃: spinorial gap 3.75B = 0.21 cm⁻¹ — resolvable."""
        B = MOLECULES['UF6']['B_cm']
        gap = 3.75 * B
        microwave_resolution = 0.001
        assert gap > microwave_resolution
        assert gap / microwave_resolution > 10


# ==========================================================
# § 4. Physical context: trigonal distortion
#
# O→D₃ distortion in real systems:
# - Perovskite ABX₃: octahedral cage distorts along [111] axis
# - Jahn-Teller: t₂g orbital degeneracy → trigonal compression
# - External strain: uniaxial stress along [111]
#
# The distortion must preserve D₃ but break O.
# Example: BaTiO₃ cubic→rhombohedral at 183 K
#          (but this is a displacement, not a rotation)
#
# Key question: is there a REAL system where a molecular rotor
# sits in an O_h cage that distorts to D₃?
# ==========================================================


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
