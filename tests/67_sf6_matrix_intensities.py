#!/usr/bin/env python3
"""
I29 — SF₆ in rare gas matrix: line-by-line intensity predictions.

QUESTION: How does the rotational spectrum of SF₆ change under
O → D₃ distortion (trigonal site in Ar matrix), and how does the
reflection identity constrain the intensity pattern?

SF₆: ¹⁹F (I=1/2), 6 equivalent nuclei, 2⁶=64 spin states.
B = 0.09107 cm⁻¹.

PLAN:
  §0. Nuclear spin statistical weights under O and D₃.
      Verify: g(A₀)=10, g(A₁)=2 under O; g(χ₀)=16, g(χ₂)=8 under D₃.
  §1. Line intensities I(j) ∝ g_χ · (2j+1) · m(χ,j) · exp(-E_j/kT).
      Compute stick spectrum for each sector, T = 5K, 10K, 20K.
  §2. Reflection constrains intensity ratios under O (obstructed);
      constraint BREAKS under D₃ (non-obstructed, D_n odd).
  §3. Observable predictions: intensity ratio changes O → D₃,
      new lines appearing, line spacing pattern.

STATUS: COMPLETE
RESULT: 21/21 pass.
  §0: O and D₃ multiplicities computed ab initio from build_binary_octahedral
      and build_dicyclic. Nuclear spin weights g(A₀)=10, g(A₁)=2, g(χ₀)=16,
      g(χ₂)=8 derived from cycle counting on octahedral atoms (not hardcoded).
      O recurrence m(j+12)=m(j)+1 verified against group computation.
  §1: Stick spectra for O (A₀,A₁) and D₃ (χ₀,χ₂) at T=5,10,20 K.
  §2: O satisfies anti-correlation (paired m sum = 1); D₃ BREAKS it
      (m(χ₀,0)+m(χ₀,2) = 2 ≠ 1). χ₂ fails trivially (both m=0).
  §3: j=2 line appears in D₃ χ₀ but absent in O A₀. Line density
      per period: D₃ 2/3 > O 1/2.

VERDICT: Concrete measurable predictions for SF₆ / Ar matrix isolation
spectroscopy. Reflection failure under D₃ → qualitatively different
intensity pattern, testable by comparing gas-phase and matrix spectra.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.group import (build_binary_octahedral, compute_conjugacy_classes,
                       compute_commutator_subgroup)
from src.quaternion import qkey, qmul, quat_to_so3


# ---- Physical constants ----
B_CM = 0.09107       # cm⁻¹
CM_TO_K = 1.4388     # 1 cm⁻¹ = 1.4388 K
B_K = B_CM * CM_TO_K  # rotational constant in Kelvin


def chi_su2(j, trace):
    """SU(2) spin-j character at given trace = 2cos(α)."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    return np.sin(d * ha) / np.sin(ha)


# ---- Group construction ----

def build_dicyclic(N):
    """Build Dic_N = 2D_N of order 4N. Standard construction."""
    elements = []
    for k in range(2 * N):
        angle = k * np.pi / N
        elements.append(np.array([np.cos(angle), 0, 0, np.sin(angle)]))
    b = np.array([0, 0, 1, 0], dtype=float)
    for k in range(2 * N):
        angle = k * np.pi / N
        ak = np.array([np.cos(angle), 0, 0, np.sin(angle)])
        elements.append(qmul(b, ak))
    return np.array(elements)


def build_2D3_in_2O():
    """Build 2D₃ ⊂ 2O by closure from generators.
    C₃ along [111], C₂ perpendicular. Order 12."""
    a = np.array([1, 1, 1, 1]) / 2.0       # C₃ along [111]
    x = np.array([0, 1, -1, 0]) / np.sqrt(2)  # C₂
    elements = [np.array([1, 0, 0, 0], dtype=float)]
    seen = {qkey(elements[0])}
    queue = list(elements)
    idx = 0
    while idx < len(queue):
        g = queue[idx]; idx += 1
        for gen in [a, x]:
            for h in [qmul(g, gen), qmul(gen, g)]:
                k = qkey(h)
                if k not in seen:
                    seen.add(k)
                    queue.append(h)
                    elements.append(h)
    return np.array(elements)


# SF₆ atom positions: vertices of octahedron
ATOM_POS = np.array([
    [1, 0, 0], [-1, 0, 0],
    [0, 1, 0], [0, -1, 0],
    [0, 0, 1], [0, 0, -1],
], dtype=float)


def nuclear_spin_character(q):
    """χ_spin(g) = 2^{cycles(σ(g))} for 6 spin-1/2 nuclei on octahedron."""
    R = quat_to_so3(q)
    rotated = ATOM_POS @ R.T
    n = len(ATOM_POS)
    perm = []
    for i in range(n):
        for j_idx in range(n):
            if np.allclose(rotated[i], ATOM_POS[j_idx], atol=1e-6):
                perm.append(j_idx)
                break
    visited = [False] * n
    n_cycles = 0
    for i in range(n):
        if not visited[i]:
            n_cycles += 1
            k = i
            while not visited[k]:
                visited[k] = True
                k = perm[k]
    return 2 ** n_cycles


def compute_nuclear_spin_weights(characters, classes, sizes, order):
    """Derive g(χ) for each 1D character from χ_spin = 2^{cycles}."""
    chi_spin = np.array([nuclear_spin_character(c['rep']) for c in classes],
                        dtype=float)
    weights = {}
    for ci, char in enumerate(characters):
        chi_k = char['chi_on_classes']
        g = round((np.sum(chi_k.conj() * chi_spin * sizes) / order).real)
        weights[ci] = g
    return weights


def compute_1d_characters(elements):
    """Compute all 1D characters via abelianization.
    Returns (characters, classes, sizes, |G|, obstructed)."""
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    G = int(sum(sizes))

    comm = compute_commutator_subgroup(elements)
    comm_keys = set(comm.keys())
    mk = qkey(np.array([-1, 0, 0, 0], dtype=float))
    obstructed = mk in comm_keys
    abel_size = G // len(comm)

    cosets = []
    elem_coset = {}
    for e in elements:
        k = qkey(e)
        if k in elem_coset:
            continue
        coset = set()
        for ck, cv in comm.items():
            coset.add(qkey(qmul(e, cv)))
        found = False
        for ci, existing in enumerate(cosets):
            if existing & coset:
                for kk in coset:
                    elem_coset[kk] = ci
                found = True
                break
        if not found:
            ci = len(cosets)
            cosets.append(coset)
            for kk in coset:
                elem_coset[kk] = ci

    gen_elem = None
    if abel_size > 1:
        for e in elements:
            k = qkey(e)
            if k in comm_keys:
                continue
            power = e.copy()
            for d in range(1, abel_size + 1):
                pk = qkey(power)
                if pk in comm_keys:
                    if d == abel_size:
                        gen_elem = e.copy()
                    break
                power = qmul(power, e)
            if gen_elem is not None:
                break

    coset_label = {}
    if abel_size == 1:
        coset_label[elem_coset[qkey(np.array([1, 0, 0, 0]))]] = 0
    else:
        power = np.array([1, 0, 0, 0], dtype=float)
        for d in range(abel_size):
            ci = elem_coset[qkey(power)]
            coset_label[ci] = d
            power = qmul(power, gen_elem)

    class_labels = np.array([coset_label[elem_coset[list(c['keys'])[0]]]
                             for c in classes])
    omega = np.exp(2j * np.pi / abel_size) if abel_size > 1 else 1.0

    minus_one_label = 0
    for i, c in enumerate(classes):
        if abs(c['trace'] + 2) < 1e-10:
            minus_one_label = class_labels[i]
            break

    characters = []
    for k in range(abel_size):
        chi_k = omega ** (k * class_labels)
        chi_m1 = (omega ** (k * minus_one_label)).real
        is_spin = chi_m1 < -0.5
        characters.append({
            'name': f'chi_{k}',
            'chi_on_classes': chi_k,
            'chi_minus_one': chi_m1,
            'is_spinorial': is_spin,
        })

    return characters, classes, sizes, G, obstructed


def compute_multiplicities(characters, classes, sizes, G, j_max=30):
    """Compute m(chi, j) for all characters and j up to j_max."""
    mult = {}
    for ci, char in enumerate(characters):
        chi_k = char['chi_on_classes']
        is_spin = char['is_spinorial']
        for two_j in range(0, int(2 * j_max) + 1):
            j = two_j / 2.0
            if is_spin and two_j % 2 == 0:
                continue
            if not is_spin and two_j % 2 == 1:
                continue
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            m = round((np.sum(chi_k.conj() * chi_Vj * sizes) / G).real)
            if m > 0:
                mult[(ci, j)] = m
    return mult


# ==========================================================
# § 0. Nuclear spin statistical weights
# ==========================================================

# Nuclear spin weights (from file 60, verified independently)
# O phase: g(A₀)=10, g(A₁)=2
# D₃ phase: g(χ₀)=16, g(χ₂)=8, g(χ₁)=g(χ₃)=0 (spinorial dark)
G_O_WEIGHTS = {'A0': 10, 'A1': 2}
G_D3_WEIGHTS = {'chi0': 16, 'chi2': 8, 'chi1': 0, 'chi3': 0}

# Multiplicity patterns (from tests 51, 64, verified)
# A₀ of O: m(A₀,j) for j=0..11 (j*=11)
M_A0_O = [1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0]
# A₁ of O: m(A₁,j) for j=0..11
M_A1_O = [0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1]


def m_sector_O(sector, j):
    """m(sector, j) for octahedral group. j integer, j* = 11."""
    base = M_A0_O if sector == 'A0' else M_A1_O
    if j < 0:
        return 0
    if j < len(base):
        return base[j]
    return m_sector_O(sector, j - 12) + 1


# Precompute O data
_O_ELEMS = np.array(build_binary_octahedral())
_O_CHARS, _O_CLASSES, _O_SIZES, _O_ORDER, _O_OBST = \
    compute_1d_characters(_O_ELEMS)
_O_MULT = compute_multiplicities(_O_CHARS, _O_CLASSES, _O_SIZES,
                                 _O_ORDER, j_max=50)
_O_TRIVIAL_IDX = None
for _i, _c in enumerate(_O_CHARS):
    if np.allclose(_c['chi_on_classes'], 1.0, atol=1e-8):
        _O_TRIVIAL_IDX = _i
        break
_O_SIGN_IDX = [i for i in range(len(_O_CHARS)) if i != _O_TRIVIAL_IDX][0]

# Nuclear spin weights for O (ab initio)
_O_SPIN_WEIGHTS = compute_nuclear_spin_weights(
    _O_CHARS, _O_CLASSES, _O_SIZES, _O_ORDER)

# Precompute D₃ data — use 2D₃⊂2O for spin weights (correct atom permutations)
_D3_PHYS_ELEMS = build_2D3_in_2O()
_O_KEYS = set(qkey(e) for e in _O_ELEMS)
assert all(qkey(e) in _O_KEYS for e in _D3_PHYS_ELEMS), \
    "2D₃ elements not all in 2O — embedding broken"
_D3_PHYS_CHARS, _D3_PHYS_CLASSES, _D3_PHYS_SIZES, _D3_PHYS_ORDER, _ = \
    compute_1d_characters(_D3_PHYS_ELEMS)
_D3_PHYS_TRIVIAL = None
for _i, _c in enumerate(_D3_PHYS_CHARS):
    if np.allclose(_c['chi_on_classes'], 1.0, atol=1e-8):
        _D3_PHYS_TRIVIAL = _i
        break
_D3_SPIN_WEIGHTS = compute_nuclear_spin_weights(
    _D3_PHYS_CHARS, _D3_PHYS_CLASSES, _D3_PHYS_SIZES, _D3_PHYS_ORDER)

# Precompute D₃ data — standard construction for multiplicities
_D3_ELEMS = build_dicyclic(3)
_D3_CHARS, _D3_CLASSES, _D3_SIZES, _D3_ORDER, _D3_OBST = \
    compute_1d_characters(_D3_ELEMS)
_D3_MULT = compute_multiplicities(_D3_CHARS, _D3_CLASSES, _D3_SIZES,
                                  _D3_ORDER, j_max=50)

# Identify tensorial character indices
_D3_TENS_INDICES = [i for i, c in enumerate(_D3_CHARS) if not c['is_spinorial']]
_D3_SPIN_INDICES = [i for i, c in enumerate(_D3_CHARS) if c['is_spinorial']]

# Identify trivial character: all class values ≈ 1
_D3_TRIVIAL_IDX = None
for _i, _c in enumerate(_D3_CHARS):
    if np.allclose(_c['chi_on_classes'], 1.0, atol=1e-8):
        _D3_TRIVIAL_IDX = _i
        break
_D3_SIGN_IDX = [i for i in _D3_TENS_INDICES if i != _D3_TRIVIAL_IDX][0]


def m_sector_D3(chi_idx, j):
    """m(χ_idx, j) for D₃ from precomputed table."""
    return _D3_MULT.get((chi_idx, float(j)), 0)


class TestNuclearWeights:
    """§0: Verify nuclear spin weights and multiplicities."""

    def test_O_group_structure(self):
        """2O: order 48, obstructed."""
        assert _O_ORDER == 48
        assert _O_OBST

    def test_O_multiplicities_ab_initio(self):
        """Build 2O, compute m(A₀,j) and m(A₁,j), verify hardcoded tables."""
        for j in range(12):
            m_triv = _O_MULT.get((_O_TRIVIAL_IDX, float(j)), 0)
            m_sign = _O_MULT.get((_O_SIGN_IDX, float(j)), 0)
            assert m_triv == M_A0_O[j], \
                f"A₀: m({j}) ab initio={m_triv}, table={M_A0_O[j]}"
            assert m_sign == M_A1_O[j], \
                f"A₁: m({j}) ab initio={m_sign}, table={M_A1_O[j]}"

    def test_O_recurrence_ab_initio(self):
        """Verify m(j+12) = m(j)+1 for O against ab initio computation."""
        for j in range(12, 25):
            for ci, label in [(_O_TRIVIAL_IDX, 'A0'), (_O_SIGN_IDX, 'A1')]:
                m_j = _O_MULT.get((ci, float(j)), 0)
                m_prev = _O_MULT.get((ci, float(j - 12)), 0)
                assert m_j == m_prev + 1, \
                    f"{label}: m({j})={m_j} ≠ m({j-12})+1={m_prev+1}"

    def test_O_reflection_A0(self):
        """A₀ reflection from ab initio: m(j)+m(11-j)=1 for j=0..11."""
        for j in range(12):
            m_j = _O_MULT.get((_O_TRIVIAL_IDX, float(j)), 0)
            m_jp = _O_MULT.get((_O_TRIVIAL_IDX, float(11 - j)), 0)
            assert m_j + m_jp == 1

    def test_O_reflection_A1(self):
        """A₁ reflection from ab initio: m(j)+m(11-j)=1 for j=0..11."""
        for j in range(12):
            m_j = _O_MULT.get((_O_SIGN_IDX, float(j)), 0)
            m_jp = _O_MULT.get((_O_SIGN_IDX, float(11 - j)), 0)
            assert m_j + m_jp == 1

    def test_D3_group_structure(self):
        """2D₃ = Dic₃: order 12, non-obstructed."""
        assert _D3_ORDER == 12
        assert not _D3_OBST  # D₃ is D_n with n=3 odd → non-obstructed

    def test_D3_four_characters(self):
        """D₃ has Abel = Z₄, so 4 scalar characters: 2 tensorial, 2 spinorial."""
        assert len(_D3_CHARS) == 4
        assert len(_D3_TENS_INDICES) == 2
        assert len(_D3_SPIN_INDICES) == 2

    def test_D3_tensorial_multiplicities(self):
        """D₃ tensorial m(χ₀, j) for j=0..5 matches analytic calculation.
        j*=2 for D₃ (|H|=6). Since non-obstructed, reflection FAILS."""
        ci = _D3_TRIVIAL_IDX  # χ₀ (trivial)
        expected = [1, 0, 1, 1, 2, 1]  # computed analytically
        for j in range(6):
            assert m_sector_D3(ci, j) == expected[j], \
                f"m(χ₀, {j}) = {m_sector_D3(ci, j)}, expected {expected[j]}"

    def test_O_spin_weights_ab_initio(self):
        """Derive g(A₀)=10, g(A₁)=2 from cycle counting on octahedral atoms."""
        assert _O_SPIN_WEIGHTS[_O_TRIVIAL_IDX] == 10
        assert _O_SPIN_WEIGHTS[_O_SIGN_IDX] == 2

    def test_D3_spin_weights_ab_initio(self):
        """Derive g(χ₀)=16, g(χ₂)=8 from cycle counting, 2D₃⊂2O embedding."""
        assert _D3_PHYS_ORDER == 12
        # Identify tensorial characters of the physical D₃
        tens = [i for i, c in enumerate(_D3_PHYS_CHARS)
                if not c['is_spinorial']]
        spin = [i for i, c in enumerate(_D3_PHYS_CHARS)
                if c['is_spinorial']]
        # Tensorial weights: 16 and 8
        tens_weights = sorted([_D3_SPIN_WEIGHTS[i] for i in tens], reverse=True)
        assert tens_weights == [16, 8]
        # Spinorial weights: both 0
        for i in spin:
            assert _D3_SPIN_WEIGHTS[i] == 0

    def test_spin_states_conservation(self):
        """Total spin states = 2⁶ = 64 for both groups.
        Σ_χ g(χ)·dim(χ) = 64. All 1D chars here, so dim=1."""
        # O: 2 chars, both 1D. But O also has higher-dim irreps carrying spin.
        # For 1D chars only: g(A₀)+g(A₁) = 12 (not 64, rest in 2D/3D irreps).
        assert _O_SPIN_WEIGHTS[_O_TRIVIAL_IDX] + \
               _O_SPIN_WEIGHTS[_O_SIGN_IDX] == 12
        # D₃: 4 chars (1D). g(χ₀)+g(χ₁)+g(χ₂)+g(χ₃) = 16+0+8+0 = 24.
        # Rest in 2D irrep E.
        assert sum(_D3_SPIN_WEIGHTS.values()) == 24


# ==========================================================
# § 1. Line intensities under O and D₃
# ==========================================================

def line_intensity(g_chi, m_chi_j, j, T):
    """I(χ, j, T) = g_χ · (2j+1) · m(χ,j) · exp(-B·j(j+1)/T)."""
    if m_chi_j == 0:
        return 0.0
    Ej = B_K * j * (j + 1)  # energy in Kelvin
    return g_chi * (2 * j + 1) * m_chi_j * np.exp(-Ej / T)


def stick_spectrum_O(sector, j_max, T):
    """Compute stick spectrum for O sector at temperature T."""
    g = G_O_WEIGHTS[sector]
    intensities = []
    for j in range(j_max + 1):
        m = m_sector_O(sector, j)
        intensities.append(line_intensity(g, m, j, T))
    return np.array(intensities)


def stick_spectrum_D3(chi_idx, g_chi, j_max, T):
    """Compute stick spectrum for D₃ sector at temperature T."""
    intensities = []
    for j in range(j_max + 1):
        m = m_sector_D3(chi_idx, j)
        intensities.append(line_intensity(g_chi, m, j, T))
    return np.array(intensities)


class TestLineIntensities:
    """§1: Line intensities I(j) ∝ g·(2j+1)·m·exp(-Ej/kT)."""

    J_MAX = 30

    def test_O_A0_spectrum_nonzero(self):
        """A₀ sector has lines at j = 0, 4, 6, 8, 9, 10, 12, ..."""
        spec = stick_spectrum_O('A0', self.J_MAX, 10.0)
        # j=0 present (m=1), j=1,2,3 absent (m=0), j=4 present (m=1)
        assert spec[0] > 0
        assert spec[1] == 0
        assert spec[4] > 0

    def test_O_A1_spectrum_nonzero(self):
        """A₁ sector has lines at j = 3, 6, 7, 9, 10, 11, ..."""
        spec = stick_spectrum_O('A1', self.J_MAX, 10.0)
        assert spec[0] == 0
        assert spec[3] > 0
        assert spec[6] > 0

    def test_Boltzmann_suppression(self):
        """Higher T populates more j; lower T concentrates on low j."""
        spec_5 = stick_spectrum_O('A0', self.J_MAX, 5.0)
        spec_20 = stick_spectrum_O('A0', self.J_MAX, 20.0)
        # At T=5K, j=0 should be relatively stronger vs high-j
        ratio_5 = spec_5[0] / (spec_5[0] + spec_5[12]) if spec_5[12] > 0 else 1.0
        ratio_20 = spec_20[0] / (spec_20[0] + spec_20[12])
        assert ratio_5 > ratio_20

    def test_total_weight_scales(self):
        """Total spectral weight ∝ g_χ at each T."""
        T = 10.0
        W_A0 = sum(stick_spectrum_O('A0', self.J_MAX, T))
        W_A1 = sum(stick_spectrum_O('A1', self.J_MAX, T))
        # Not exactly 5:1 because m patterns differ, but g_A0/g_A1 = 5
        # Total weight = g · Σ (2j+1)·m·exp(-Ej/T), so ratio ≈ 5 × (sum_A0/sum_A1)
        assert W_A0 > W_A1  # at minimum, A₀ dominates


# ==========================================================
# § 2. Reflection constrains intensity ratios
# ==========================================================

class TestReflectionConstraint:
    """§2: Under O (obstructed), reflection gives anti-correlated pairs.
    Under D₃ (non-obstructed), this constraint breaks."""

    def test_D3_reflection_FAILS(self):
        """D₃ χ₀: m(0)+m(2) = 2 ≠ 1. Reflection broken for tensorial."""
        ci = _D3_TRIVIAL_IDX
        j_star = 2  # |D₃|/2 - 1 = 2
        m0 = m_sector_D3(ci, 0)
        m2 = m_sector_D3(ci, j_star)
        assert m0 + m2 == 2, f"Expected reflection failure: m(0)+m(2)={m0+m2}"
        # Both lines PRESENT — no anti-correlation
        assert m0 == 1 and m2 == 1

    def test_D3_chi2_reflection_FAILS(self):
        """D₃ χ₂ (sign character): m(0)=m(2)=0, so m+m'=0 ≠ 1.
        Reflection fails trivially — neither level is occupied in χ₂."""
        ci = _D3_SIGN_IDX
        m0 = m_sector_D3(ci, 0)
        m2 = m_sector_D3(ci, 2)
        assert m0 == 0 and m2 == 0, \
            f"Expected both empty for χ₂: m(0)={m0}, m(2)={m2}"
        assert m0 + m2 != 1  # reflection says sum should be 1; it's 0

    def test_O_paired_intensity_ratio(self):
        """Under O, paired lines have infinite intensity ratio (one is zero)."""
        T = 10.0
        spec = stick_spectrum_O('A0', 11, T)
        for j in range(6):  # j and 11-j
            jp = 11 - j
            I_j = spec[j]
            I_jp = spec[jp]
            # Exactly one should be zero
            assert (I_j == 0) != (I_jp == 0), \
                f"j={j},{jp}: I={I_j:.4f},{I_jp:.4f} — both zero or both nonzero"

    def test_D3_paired_both_present(self):
        """Under D₃, j=0 and j=2 BOTH have intensity → ratio is finite."""
        ci = _D3_TRIVIAL_IDX
        T = 10.0
        spec = stick_spectrum_D3(ci, G_D3_WEIGHTS['chi0'], 5, T)
        # Both j=0 and j=2 lines present
        assert spec[0] > 0 and spec[2] > 0
        ratio = spec[0] / spec[2]
        assert 0 < ratio < np.inf


# ==========================================================
# § 3. Observable predictions
# ==========================================================

class TestObservablePredictions:
    """§3: Concrete predictions for SF₆ gas → Ar matrix."""

    J_MAX = 30

    def test_new_lines_in_D3(self):
        """Lines present in D₃ χ₀ that are absent in O A₀.
        j=2 is the key example: m(A₀,2)=0 but m(χ₀,2)=1."""
        m_O_j2 = M_A0_O[2]  # 0
        ci = _D3_TRIVIAL_IDX
        m_D3_j2 = m_sector_D3(ci, 2)  # 1
        assert m_O_j2 == 0, "j=2 should be dark in O A₀"
        assert m_D3_j2 == 1, "j=2 should appear in D₃ χ₀"

    def test_line_density_per_period(self):
        """D₃ has higher line density per period than O.
        O: period 12, 6 lines in A₀ → density 1/2.
        D₃: period 3, count active lines in one period → density 2/3."""
        # O: 6 active j in period [0..11]
        active_O = sum(1 for j in range(12) if M_A0_O[j] > 0)
        density_O = active_O / 12.0
        # D₃: count active j in one period [0..2]
        ci = _D3_TRIVIAL_IDX
        active_D3 = sum(1 for j in range(3) if m_sector_D3(ci, j) > 0)
        density_D3 = active_D3 / 3.0
        assert density_D3 > density_O, \
            f"D₃ density {density_D3:.2f} should exceed O density {density_O:.2f}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
