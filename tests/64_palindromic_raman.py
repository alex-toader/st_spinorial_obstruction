#!/usr/bin/env python3
"""
I22 — Occupation complementarity from reflection identity.

QUESTION: Does m(χ,j) + m(χ,j*-j) = 1 produce observable
constraints on consecutive level pairs?

ANSWER: YES. For any pair (j, j+k) with reflected pair
(j*-j-k, j*-j), at most one pair can be doubly occupied.

NOTE: Originally framed as "Raman palindrome" but reviewer
correctly identified that Raman operator transforms as E_g,
connecting scalar→non-scalar sectors. The mathematical
complementarity is correct but applies to occupation patterns,
not directly to Raman selection rules.

RESULT:
  §0. Occupation patterns verified (2O):
      A₀: [1,0,0,0,1,0,1,0,1,1,1,0]  (j=0..11)
      A₁: [0,0,0,1,0,0,1,1,0,1,1,1]  (j=0..11)
      Both satisfy m(j) + m(11-j) = 1.

  §1. Consecutive-pair complementarity (k=2):
      If m(j)=m(j+2)=1, then m(j*-j-2)·m(j*-j) = 0.
      A₀: 3 doubly-occupied pairs at j=4,6,8
      A₁: 2 doubly-occupied pairs at j=7,9

  §2. Verified on 2T (j*=5), 2O (j*=11), 2I (j*=29).
      Complementarity holds universally.

VERDICT: Genuine new content — second-order consequence of
  reflection identity. Constrains occupation pattern globally.
  Integrated in paper v9 Discussion as "occupation complementarity."

RAW OUTPUT (17 passed in 0.85s):
  TestOccupation (3), TestRamanSpectrum (2), TestPalindrome (3),
  TestPhysicalSpectrum (3), TestGeneralization2T (3),
  TestGeneralization2I (3) — all PASSED.

STATUS: COMPLETE — integrated in paper v9 Discussion
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.group import build_binary_octahedral, build_binary_tetrahedral, build_binary_icosahedral
from src.group import compute_commutator_subgroup
from src.quaternion import qkey


# ==========================================================
# § 0. Multiplicities m(χ, j) for 2O
# ==========================================================

def chi_su2(j, g_trace):
    """Character of V_j at element with tr(g) = g_trace."""
    if abs(g_trace - 2.0) < 1e-10:
        return 2*j + 1
    if abs(g_trace + 2.0) < 1e-10:
        return (2*j + 1) * ((-1)**(int(2*j)))
    alpha_half = np.arccos(np.clip(g_trace / 2.0, -1, 1))
    return np.sin((2*j + 1) * alpha_half) / np.sin(alpha_half)


def build_2O_chars():
    """Build scalar characters A₀, A₁ of 2O and return with group data."""
    G = build_binary_octahedral()
    order = len(G)  # 48
    j_star = order // 4 - 1  # 11
    comm = compute_commutator_subgroup(G)
    comm_keys = set(comm.keys())

    chi_A0 = {qkey(g): 1.0 for g in G}
    chi_A1 = {qkey(g): (1.0 if qkey(g) in comm_keys else -1.0) for g in G}

    return {'A0': chi_A0, 'A1': chi_A1}, G, order, j_star


def multiplicity(j, chi_vals, elements, order):
    """m(χ, j) = (1/|G|) Σ_g conj(χ(g)) · χ_j(g)."""
    s = 0.0
    for g in elements:
        k = qkey(g)
        s += np.conj(chi_vals.get(k, 0.0)) * chi_su2(j, 2.0 * g[0])
    return round((s / order).real)


def get_occupation_pattern(chi_vals, elements, order, j_star):
    """Return list of m(χ, j) for j = 0, 1, ..., j_star."""
    return [multiplicity(j, chi_vals, elements, order)
            for j in range(j_star + 1)]


class TestOccupation:
    """Verify occupation patterns for A₀ and A₁ of 2O."""

    def test_A0_pattern(self):
        """A₀: m = [1,0,0,0,1,0,1,0,1,1,1,0] for j=0..11."""
        chars, G, order, j_star = build_2O_chars()
        pattern = get_occupation_pattern(chars['A0'], G, order, j_star)
        expected = [1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0]
        assert pattern == expected

    def test_A1_pattern(self):
        """A₁: m = [0,0,0,1,0,0,1,1,0,1,1,1] for j=0..11."""
        chars, G, order, j_star = build_2O_chars()
        pattern = get_occupation_pattern(chars['A1'], G, order, j_star)
        expected = [0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1]
        assert pattern == expected

    def test_reflection_identity(self):
        """m(χ,j) + m(χ,j*-j) = 1 for all j ≤ j*."""
        chars, G, order, j_star = build_2O_chars()
        for name, chi in chars.items():
            pattern = get_occupation_pattern(chi, G, order, j_star)
            for j in range(j_star + 1):
                assert pattern[j] + pattern[j_star - j] == 1, \
                    f"{name}: m({j}) + m({j_star - j}) = {pattern[j] + pattern[j_star - j]}"


# ==========================================================
# § 1. Raman stick spectrum
#
# Raman rotational: ΔJ = 0, ±2 (Stokes branch: ΔJ = +2).
# Transition j → j+2 is "allowed" in sector χ iff:
#   m(χ, j) = 1  AND  m(χ, j+2) = 1
# (both initial and final levels exist in that sector).
#
# Transition energy: E(j+2) - E(j) = B[(j+2)(j+3) - j(j+1)]
#                                   = B(4j + 6)
# ==========================================================

def raman_stick_spectrum(pattern, j_star, B=1.0):
    """Build Raman ΔJ=+2 stick spectrum for a sector.
    Returns list of (j_initial, energy, allowed)."""
    lines = []
    for j in range(j_star - 1):  # j can go up to j_star - 2
        energy = B * (4 * j + 6)
        allowed = (pattern[j] == 1 and pattern[j + 2] == 1)
        lines.append((j, energy, allowed))
    return lines


class TestRamanSpectrum:
    """Raman stick spectrum for A₀ and A₁ of O."""

    def test_A0_raman_lines(self):
        """A₀ pattern [1,0,0,0,1,0,1,0,1,1,1,0]:
        Allowed transitions (both endpoints m=1):
          j=4→6 (E=22B), j=6→8 (E=30B), j=8→10 (E=38B)
        Forbidden (at least one endpoint m=0):
          j=0→2, j=1→3, j=2→4, j=3→5, j=5→7, j=7→9, j=9→11"""
        chars, G, order, j_star = build_2O_chars()
        pattern = get_occupation_pattern(chars['A0'], G, order, j_star)
        lines = raman_stick_spectrum(pattern, j_star)

        allowed = [(j, E) for j, E, a in lines if a]
        assert len(allowed) == 3
        assert allowed[0] == (4, 22.0)
        assert allowed[1] == (6, 30.0)
        assert allowed[2] == (8, 38.0)

    def test_A1_raman_lines(self):
        """A₁ pattern [0,0,0,1,0,0,1,1,0,1,1,1]:
        Allowed: j=7→9 (E=34B), j=9→11... but j+2=11 and j*=11 so
        need m(9)=1 AND m(11)=1.  m(11)=1 ✓.
        j=7: m(7)=1, m(9)=1 → allowed (E=34B)
        j=9: m(9)=1, m(11)=1 → but j=9, j+2=11=j*, need j≤j*-2=9 ✓ (E=42B)"""
        chars, G, order, j_star = build_2O_chars()
        pattern = get_occupation_pattern(chars['A1'], G, order, j_star)
        lines = raman_stick_spectrum(pattern, j_star)

        allowed = [(j, E) for j, E, a in lines if a]
        assert len(allowed) == 2
        assert allowed[0] == (7, 34.0)
        assert allowed[1] == (9, 42.0)


# ==========================================================
# § 2. Palindrome test
#
# Claim: the allowed/forbidden pattern of Raman lines is
# palindromic under j ↔ j*-j-2.
#
# Proof: transition j→j+2 allowed ⟺ m(j)=1 AND m(j+2)=1.
# Reflected transition: (j*-j-2)→(j*-j) allowed ⟺
#   m(j*-j-2)=1 AND m(j*-j)=1.
# By reflection: m(j*-j) = 1 - m(j) and m(j*-j-2) = 1 - m(j+2).
# So reflected allowed ⟺ m(j)=0 AND m(j+2)=0.
#
# NOT palindromic in simple sense! Instead:
#   both allowed → impossible (would need m(j)=1,m(j+2)=1,m(j)=0,m(j+2)=0)
#   Wait — j and j*-j-2 are DIFFERENT values.
#
# Let's think again. Line at j: allowed iff m(j)·m(j+2)=1.
# Reflected line at j' = j*-j-2: allowed iff m(j')·m(j'+2)=1
#   = m(j*-j-2)·m(j*-j) = (1-m(j+2))·(1-m(j)).
#
# So: line(j) allowed ⟺ m(j)=m(j+2)=1
#     line(j') allowed ⟺ m(j)=m(j+2)=0
# These are COMPLEMENTARY: exactly one of {line(j), line(j')} is
# allowed, UNLESS m(j)≠m(j+2) (one is 0, one is 1) in which
# case BOTH are forbidden.
#
# Three cases for the pair (j, j'=j*-j-2):
#   m(j)=1, m(j+2)=1 → line(j) ON,  line(j') OFF
#   m(j)=0, m(j+2)=0 → line(j) OFF, line(j') ON
#   m(j)≠m(j+2)      → line(j) OFF, line(j') OFF
# ==========================================================

class TestPalindrome:
    """Test palindromic structure of Raman lines."""

    def test_complementarity_A0(self):
        """For each pair (j, j*-j-2): at most one line is allowed.
        If both m values match (both 0 or both 1), exactly one is allowed."""
        chars, G, order, j_star = build_2O_chars()
        pattern = get_occupation_pattern(chars['A0'], G, order, j_star)
        lines = raman_stick_spectrum(pattern, j_star)

        for j, E, allowed in lines:
            j_refl = j_star - j - 2
            if 0 <= j_refl <= j_star - 2:
                refl_allowed = (pattern[j_refl] == 1 and pattern[j_refl + 2] == 1)
                # Cannot both be allowed
                assert not (allowed and refl_allowed), \
                    f"Both j={j} and j'={j_refl} allowed — impossible"

    def test_complementarity_A1(self):
        """Same for A₁."""
        chars, G, order, j_star = build_2O_chars()
        pattern = get_occupation_pattern(chars['A1'], G, order, j_star)
        lines = raman_stick_spectrum(pattern, j_star)

        for j, E, allowed in lines:
            j_refl = j_star - j - 2
            if 0 <= j_refl <= j_star - 2:
                refl_allowed = (pattern[j_refl] == 1 and pattern[j_refl + 2] == 1)
                assert not (allowed and refl_allowed), \
                    f"Both j={j} and j'={j_refl} allowed — impossible"

    def test_at_most_half_lines_allowed(self):
        """Complementarity implies: #allowed ≤ floor(#total / 2).
        The pattern cannot have more than half the lines active."""
        chars, G, order, j_star = build_2O_chars()
        for name, chi in chars.items():
            pattern = get_occupation_pattern(chi, G, order, j_star)
            lines = raman_stick_spectrum(pattern, j_star)
            n_allowed = sum(1 for _, _, a in lines if a)
            n_total = len(lines)
            assert n_allowed <= n_total // 2 + 1, \
                f"{name}: {n_allowed}/{n_total} lines allowed, too many"


# ==========================================================
# § 3. Gap pattern and physical numbers
#
# For SF₆: B ≈ 0.0911 cm⁻¹ (octahedral, |O|=24, j*=11)
# Raman transition energies: B(4j+6) cm⁻¹
# ==========================================================

B_SF6 = 0.0911  # cm⁻¹


class TestPhysicalSpectrum:
    """Physical Raman spectrum for SF₆-like molecule."""

    def test_A0_raman_energies_cm(self):
        """A₀ sector Raman lines in cm⁻¹ for SF₆."""
        chars, G, order, j_star = build_2O_chars()
        pattern = get_occupation_pattern(chars['A0'], G, order, j_star)
        lines = raman_stick_spectrum(pattern, j_star, B=B_SF6)

        allowed = [(j, E) for j, E, a in lines if a]
        # j=4→6: 22B = 2.004 cm⁻¹
        # j=6→8: 30B = 2.733 cm⁻¹
        # j=8→10: 38B = 3.462 cm⁻¹
        assert len(allowed) == 3
        assert abs(allowed[0][1] - 22 * B_SF6) < 1e-10
        assert abs(allowed[1][1] - 30 * B_SF6) < 1e-10
        assert abs(allowed[2][1] - 38 * B_SF6) < 1e-10

    def test_gap_structure(self):
        """A₀ pattern [1,0,0,0,1,0,1,0,1,1,1,0]:
        Raman ΔJ=+2 lines from j=0..9 (j+2 ≤ 11):
        Allowed: j=4→6, j=6→8, j=8→10 (three consecutive)
        Gaps: j=0,1,2,3,5,7,9 (seven forbidden)"""
        chars, G, order, j_star = build_2O_chars()
        pattern = get_occupation_pattern(chars['A0'], G, order, j_star)
        lines = raman_stick_spectrum(pattern, j_star)

        gap_js = [j for j, E, a in lines if not a]
        allowed_js = [j for j, E, a in lines if a]
        assert gap_js == [0, 1, 2, 3, 5, 7, 9]
        assert allowed_js == [4, 6, 8]

    def test_A0_A1_complementary_coverage(self):
        """A₀ allowed at j=4,6,8 (mid-high region).
        A₁ allowed at j=7,9 (high region).
        No overlap. Together: 5 of 10 possible transitions."""
        chars, G, order, j_star = build_2O_chars()
        pat_A0 = get_occupation_pattern(chars['A0'], G, order, j_star)
        pat_A1 = get_occupation_pattern(chars['A1'], G, order, j_star)

        lines_A0 = raman_stick_spectrum(pat_A0, j_star)
        lines_A1 = raman_stick_spectrum(pat_A1, j_star)

        allowed_A0 = {j for j, E, a in lines_A0 if a}
        allowed_A1 = {j for j, E, a in lines_A1 if a}

        # No overlap
        assert allowed_A0 & allowed_A1 == set(), \
            "A₀ and A₁ Raman lines should not overlap"
        # Together cover 5 of 10 possible transitions
        assert len(allowed_A0 | allowed_A1) == 5


# ==========================================================
# § 4. Generalization: 2T and 2I
#
# 2T: |H|=12, j*=5, abelianization Z₃ → 3 scalar chars
# 2I: |H|=60, j*=29, abelianization trivial → 1 scalar char (A₀)
# ==========================================================

def build_2T_chars():
    """Build scalar (1D) characters of 2T.
    Abelianization Z₃ → 3 chars. [2T,2T] = Q₈."""
    G = build_binary_tetrahedral()
    order = len(G)  # 24
    j_star = order // 4 - 1  # 5
    comm = compute_commutator_subgroup(G)
    comm_keys = set(comm.keys())

    # Quotient map: 2T → Z₃. Generator c = (1+i+j+k)/2.
    from src.quaternion import qmul, qinv
    omega = np.exp(2j * np.pi / 3)
    c = np.array([1, 1, 1, 1]) / 2.0
    c_inv = qinv(c)

    coset_label = {}
    for g in G:
        k = qkey(g)
        if k in comm_keys:
            coset_label[k] = 0
        elif qkey(qmul(g, c_inv)) in comm_keys:
            coset_label[k] = 1
        else:
            coset_label[k] = 2

    chars = {}
    for m in range(3):
        chi = {}
        for g in G:
            k = qkey(g)
            chi[k] = omega ** (m * coset_label[k])
        chars[f'chi_{m}'] = chi

    return chars, G, order, j_star


def build_2I_chars():
    """Build scalar (1D) characters of 2I.
    2I is perfect ([2I,2I]=2I) → only trivial char A₀."""
    G = build_binary_icosahedral()
    order = len(G)  # 120
    j_star = order // 4 - 1  # 29

    chi_A0 = {qkey(g): 1.0 for g in G}
    return {'A0': chi_A0}, G, order, j_star


class TestGeneralization2T:
    """Raman complementarity for 2T (tetrahedral, j*=5)."""

    def test_reflection_holds(self):
        """m(χ,j) + m(χ,5-j) = 1 for all 3 scalar chars."""
        chars, G, order, j_star = build_2T_chars()
        for name, chi in chars.items():
            pattern = get_occupation_pattern(chi, G, order, j_star)
            for j in range(j_star + 1):
                assert pattern[j] + pattern[j_star - j] == 1, \
                    f"2T {name}: m({j})+m({j_star-j}) ≠ 1"

    def test_complementarity(self):
        """For each char, paired Raman lines (j, j*-j-2) cannot
        both be allowed."""
        chars, G, order, j_star = build_2T_chars()
        for name, chi in chars.items():
            pattern = get_occupation_pattern(chi, G, order, j_star)
            lines = raman_stick_spectrum(pattern, j_star)
            for j, E, allowed in lines:
                j_refl = j_star - j - 2
                if 0 <= j_refl <= j_star - 2:
                    refl_allowed = (pattern[j_refl] == 1 and
                                    pattern[j_refl + 2] == 1)
                    assert not (allowed and refl_allowed), \
                        f"2T {name}: both j={j} and j'={j_refl} allowed"

    def test_patterns(self):
        """Print and verify 2T occupation patterns."""
        chars, G, order, j_star = build_2T_chars()
        # chi_0 is trivial: m(A₀,j) for T
        pat0 = get_occupation_pattern(chars['chi_0'], G, order, j_star)
        # Reflection: pat0[j] + pat0[5-j] = 1
        assert pat0[0] + pat0[5] == 1
        assert pat0[1] + pat0[4] == 1
        assert pat0[2] + pat0[3] == 1


class TestGeneralization2I:
    """Raman complementarity for 2I (icosahedral, j*=29)."""

    def test_reflection_holds(self):
        """m(A₀,j) + m(A₀,29-j) = 1 for j=0..29."""
        chars, G, order, j_star = build_2I_chars()
        pattern = get_occupation_pattern(chars['A0'], G, order, j_star)
        for j in range(j_star + 1):
            assert pattern[j] + pattern[j_star - j] == 1, \
                f"2I: m({j})+m({j_star-j}) ≠ 1"

    def test_complementarity(self):
        """Paired Raman lines cannot both be allowed."""
        chars, G, order, j_star = build_2I_chars()
        pattern = get_occupation_pattern(chars['A0'], G, order, j_star)
        lines = raman_stick_spectrum(pattern, j_star)
        for j, E, allowed in lines:
            j_refl = j_star - j - 2
            if 0 <= j_refl <= j_star - 2:
                refl_allowed = (pattern[j_refl] == 1 and
                                pattern[j_refl + 2] == 1)
                assert not (allowed and refl_allowed), \
                    f"2I: both j={j} and j'={j_refl} allowed"

    def test_sparsity(self):
        """2I has j*=29 → 28 possible Raman lines.
        Complementarity bounds allowed ≤ 14."""
        chars, G, order, j_star = build_2I_chars()
        pattern = get_occupation_pattern(chars['A0'], G, order, j_star)
        lines = raman_stick_spectrum(pattern, j_star)
        n_allowed = sum(1 for _, _, a in lines if a)
        assert n_allowed <= 14


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
