#!/usr/bin/env python3
"""
D1+D2 — Spectral zeta and heat kernel on S²/H.

QUESTION: Does the reflection identity produce analytic consequences
beyond the combinatorial level?

RAW OUTPUT:
  §0: Heat trace cross-check (7 groups: 3 polyhedral + 4 dihedral)
      — max |equivariant - mult| < 1e-10.
  §1: Seeley-DeWitt a₋₁ via Richardson (t=0.002,0.004):
      2T: a₋₁ = 1/12 ± O(1e-6),  2O: 1/24,  2I: 1/60
      D₃..D₆: all within 1e-4 of 1/(2n).
  §2: First-period polynomial ε(k):
      Polyhedral: binary ∈ {0,1}, reflection ε(k)+ε(P-1-k)=1, zero_frac=1/2.
      D_odd (n=3,5,7): binary, reflection fails, n_zeros=(n-1)/2 exact.
      D_even (n=4,6,8): binary, reflection holds, zero_frac = 1/2.
  §3: Spectral zeta Weyl ratio (j_max=2000):
      All groups: ratio(s) strictly ↗ toward 1/|H| as s → 1⁺.
      ratio < 1/|H| at all s (approaching from below).
      At s=1.2: 2T 66%, 2O 60%, 2I 52% of 1/|H| (slow convergence
      near abscissa — logarithmic in j_max).

ANSWER: YES (consolidation). Heat/zeta machinery confirms known structure:
  - a₋₁ = 1/|H| (orbifold volume) from Seeley-DeWitt.
  - First-period polynomial ε restates reflection ⟺ (B).
  - Spectral zeta ratio ζ_{A₀}/ζ_{S²} → 1/|H| (Weyl).
  NO new spectral invariant distinguishes (B) beyond the
  combinatorial level (first-period polynomial is a restatement).
  D2 direction is a negative result for novelty upgrade.

PLAN:
  §0. Heat trace cross-check: equivariant vs multiplicity sum [7 tests]
  §1. Leading Seeley-DeWitt coefficient a₋₁ = 1/|H| [7 tests]
  §2. First-period polynomial ε(k): binary, reflection iff (B) [9 tests]
  §3. Spectral zeta Weyl ratio [3 tests]

STATUS: COMPLETE — 26/26 tests pass
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
                       build_binary_icosahedral, compute_conjugacy_classes)
from src.quaternion import qmul


# ---- shared helpers (from file 66, with guards) ----

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
    x2 = qmul(j_quat, j_quat)
    a_n = np.array([np.cos(np.pi), 0, 0, np.sin(np.pi)])
    assert np.allclose(x2, a_n, atol=1e-12), "Dicyclic relation x²=aⁿ failed"
    return G


def chi_su2(j, trace):
    """SU(2) spin-j character at given trace = 2cos(α)."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def prepare_group(G, H_order):
    """Prepare group data from element array G."""
    classes = compute_conjugacy_classes(G)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    order = len(G)
    j_star = H_order // 2 - 1
    # Condition (B): all K = |H|/d even for non-identity classes.
    # Check ALL classes (no short-circuit) so d-validity guards run everywhere.
    cond_B = True
    for c in classes:
        tr = c['trace']
        if abs(tr - 2) < 1e-10 or abs(tr + 2) < 1e-10:
            continue
        alpha = np.arccos(np.clip(tr / 2, -1, 1))
        for d in range(2, 10000):
            if abs(np.sin(d * alpha)) < 1e-8:
                break
        else:
            raise ValueError(f"SO(3) order not found for trace={tr}")
        if H_order % d != 0:
            raise ValueError(f"d={d} does not divide |H|={H_order}")
        K = H_order // d
        if K % 2 == 1:
            cond_B = False
    return {
        'classes': classes, 'sizes': sizes,
        'order': order, 'H_order': H_order,
        'j_star': j_star, 'cond_B': cond_B,
    }


def multiplicity(j, grp):
    """m(A₀, j) — multiplicity of trivial character in V_j."""
    chi_Vj = np.array([chi_su2(j, c['trace']) for c in grp['classes']])
    m = np.sum(chi_Vj * grp['sizes']) / grp['order']
    val = m.real
    assert abs(val - round(val)) < 1e-6, \
        f"j={j}: non-integer multiplicity {val}"
    return int(round(val))


# ---- heat trace functions ----

def equivariant_heat_trace(trace_su2, t, j_max):
    """K_g(t) = Σ_j χ_j(g) · e^{-t·j(j+1)}, truncated at j_max."""
    # j_max ~ 5/√t is sufficient (exponential decay)
    s = 0.0
    for j in range(j_max + 1):
        s += chi_su2(j, trace_su2) * np.exp(-t * j * (j + 1))
    return s


def sector_heat_trace_equivariant(grp, t, j_max):
    """K_{A₀}(t) via equivariant formula: (1/|2H|) Σ_g K_g(t)."""
    s = 0.0
    for i, c in enumerate(grp['classes']):
        K_g = equivariant_heat_trace(c['trace'], t, j_max)
        s += K_g * grp['sizes'][i]
    return (s / grp['order']).real


def sector_heat_trace_multiplicities(grp, t, j_max):
    """K_{A₀}(t) via multiplicity sum: Σ_j m(A₀,j) e^{-tj(j+1)}."""
    s = 0.0
    for j in range(j_max + 1):
        m = multiplicity(j, grp)
        s += m * np.exp(-t * j * (j + 1))
    return s


# ==========================================================
# § 0. Heat trace cross-check
#
# Two independent computations of K_{A₀}(t):
#   Method A: (1/|2H|) Σ_g χ_j(g) e^{-tj(j+1)}  (equivariant)
#   Method B: Σ_j m(A₀,j) e^{-tj(j+1)}           (multiplicities)
# Must agree to machine precision.
# ==========================================================

class TestHeatTraceCrossCheck:
    """Equivariant heat trace must equal multiplicity-sum heat trace."""

    @pytest.mark.parametrize("name,builder,H_order", [
        ("2T", build_binary_tetrahedral, 12),
        ("2O", build_binary_octahedral, 24),
        ("2I", build_binary_icosahedral, 60),
    ])
    def test_polyhedral(self, name, builder, H_order):
        """Two independent computations of K_{A₀}(t) must agree."""
        grp = prepare_group(builder(), H_order)
        j_max = 200  # overkill for t ≥ 0.1 (need ~ 5/√t ≈ 16)
        for t in [0.1, 0.5, 1.0]:
            K_eq = sector_heat_trace_equivariant(grp, t, j_max)
            K_mult = sector_heat_trace_multiplicities(grp, t, j_max)
            assert abs(K_eq - K_mult) < 1e-10, \
                f"{name} t={t}: equivariant={K_eq}, mult={K_mult}"

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_dihedral(self, n):
        """Dihedral cross-check: equivariant vs multiplicity heat trace."""
        G = build_binary_dihedral(n)
        grp = prepare_group(G, 2 * n)
        j_max = 200
        for t in [0.1, 0.5, 1.0]:
            K_eq = sector_heat_trace_equivariant(grp, t, j_max)
            K_mult = sector_heat_trace_multiplicities(grp, t, j_max)
            assert abs(K_eq - K_mult) < 1e-10, \
                f"D_{n} t={t}: equivariant={K_eq}, mult={K_mult}"


# ==========================================================
# § 1. Leading Seeley-DeWitt coefficient
#
# K_{A₀}(t) ~ a₋₁/t + a₀ + O(t) as t→0.
# Volume term: a₋₁ = Area(S²)/(4π|H|) = 1/|H|.
# Extract a₋₁ = lim_{t→0} t·K(t) and verify.
# ==========================================================

class TestSeeleyDeWittLeading:
    """a₋₁ = 1/|H| for all groups (volume term)."""

    @pytest.mark.parametrize("name,builder,H_order", [
        ("2T", build_binary_tetrahedral, 12),
        ("2O", build_binary_octahedral, 24),
        ("2I", build_binary_icosahedral, 60),
    ])
    def test_polyhedral(self, name, builder, H_order):
        """t·K_{A₀}(t) → 1/|H| as t → 0 (volume term)."""
        grp = prepare_group(builder(), H_order)
        expected = 1.0 / H_order
        # Richardson extrapolation: f(t) = t·K(t) = a₋₁ + a₀·t + O(t²)
        # From two points: a₋₁ = (t₂·f(t₁) - t₁·f(t₂)) / (t₂ - t₁)
        # Residual O(a₁·t₁·t₂). With t₁=0.002, t₂=0.004: ~8e-6.
        t1, t2 = 0.002, 0.004
        j_max = 500  # need ~ 5/√0.002 ≈ 112, use 500 for safety
        K1 = sector_heat_trace_equivariant(grp, t1, j_max)
        K2 = sector_heat_trace_equivariant(grp, t2, j_max)
        f1, f2 = t1 * K1, t2 * K2
        a_minus1 = (t2 * f1 - t1 * f2) / (t2 - t1)
        assert abs(a_minus1 - expected) < 1e-4, \
            f"{name}: a₋₁={a_minus1}, expected 1/|H|={expected}"

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_dihedral(self, n):
        """t·K_{A₀}(t) → 1/|H| for dihedral groups."""
        G = build_binary_dihedral(n)
        H_order = 2 * n
        grp = prepare_group(G, H_order)
        expected = 1.0 / H_order
        t1, t2 = 0.002, 0.004
        j_max = 500
        K1 = sector_heat_trace_equivariant(grp, t1, j_max)
        K2 = sector_heat_trace_equivariant(grp, t2, j_max)
        f1, f2 = t1 * K1, t2 * K2
        a_minus1 = (t2 * f1 - t1 * f2) / (t2 - t1)
        assert abs(a_minus1 - expected) < 1e-4, \
            f"D_{n}: a₋₁={a_minus1}, expected 1/|H|={expected}"


# ==========================================================
# § 2. First-period polynomial ε(k) = m(A₀, k), k=0..P-1
#
# Known properties:
#   (a) ε(k) ∈ {0,1} for all groups (binary range)
#   (b) ε(k) + ε(P-1-k) = 1  ⟺  condition (B) holds
#   (c) Under (B): exactly P/2 zeros and P/2 ones
#   (d) Without (B): ε still binary but reflection fails
#
# These are RESTATEMENTS of known results (Theorem 3.3),
# verified here as consistency checks on the heat/zeta machinery.
# ==========================================================

class TestFirstPeriodPolynomial:
    """Properties of ε(k) = m(A₀, k) in the first period [0, P-1]."""

    @pytest.mark.parametrize("name,builder,H_order", [
        ("2T", build_binary_tetrahedral, 12),
        ("2O", build_binary_octahedral, 24),
        ("2I", build_binary_icosahedral, 60),
    ])
    def test_reflecting_groups(self, name, builder, H_order):
        """Polyhedral: ε binary, reflection holds, zero_frac = 1/2."""
        grp = prepare_group(builder(), H_order)
        P = grp['j_star'] + 1
        eps = [multiplicity(k, grp) for k in range(P)]
        # Binary
        assert all(e in (0, 1) for e in eps), \
            f"{name}: non-binary ε = {eps}"
        # Reflection
        for k in range(P):
            assert eps[k] + eps[P - 1 - k] == 1, \
                f"{name}: ε({k})+ε({P-1-k}) = {eps[k]+eps[P-1-k]} ≠ 1"
        # Zero fraction exactly 1/2
        n_zeros = sum(1 for e in eps if e == 0)
        assert n_zeros * 2 == P, \
            f"{name}: {n_zeros} zeros in {P} entries, expected {P}//2"

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_non_reflecting_dihedral(self, n):
        """D_n odd: ε binary but reflection FAILS, zero_frac = (n-1)/(2n)."""
        G = build_binary_dihedral(n)
        grp = prepare_group(G, 2 * n)
        assert not grp['cond_B'], f"D_{n} should not satisfy (B)"
        P = grp['j_star'] + 1
        eps = [multiplicity(k, grp) for k in range(P)]
        # Still binary
        assert all(e in (0, 1) for e in eps), \
            f"D_{n}: non-binary ε = {eps}"
        # Reflection FAILS for at least one k
        refl_fails = any(eps[k] + eps[P - 1 - k] != 1 for k in range(P))
        assert refl_fails, f"D_{n}: reflection holds but (B) fails"
        # Exact zero count: n_zeros = (n-1)/2 for D_n odd
        n_zeros = sum(1 for e in eps if e == 0)
        assert n_zeros == (n - 1) // 2, \
            f"D_{n}: {n_zeros} zeros, expected (n-1)/2 = {(n-1)//2}"

    @pytest.mark.parametrize("n", [4, 6, 8])
    def test_reflecting_dihedral(self, n):
        """D_n even: ε binary, reflection holds, zero_frac = 1/2."""
        G = build_binary_dihedral(n)
        grp = prepare_group(G, 2 * n)
        assert grp['cond_B'], f"D_{n} should satisfy (B)"
        P = grp['j_star'] + 1
        eps = [multiplicity(k, grp) for k in range(P)]
        assert all(e in (0, 1) for e in eps), \
            f"D_{n}: non-binary ε = {eps}"
        for k in range(P):
            assert eps[k] + eps[P - 1 - k] == 1, \
                f"D_{n}: ε({k})+ε({P-1-k}) = {eps[k]+eps[P-1-k]} ≠ 1"
        n_zeros = sum(1 for e in eps if e == 0)
        assert n_zeros * 2 == P, \
            f"D_{n}: {n_zeros} zeros in {P} entries, expected {P}//2"


# ==========================================================
# § 3. Spectral zeta Weyl ratio
#
# ζ_{A₀}(s) = Σ_{j≥1} m(A₀,j) / (j(j+1))^s
# ζ_{S²}(s) = Σ_{j≥1} (2j+1) / (j(j+1))^s
#
# Weyl law for orbifold S²/H:
#   ζ_{A₀}(s) / ζ_{S²}(s) → 1/|H|  as s → 1⁺
#
# because m(A₀,j) ~ (2j+1)/|H| for large j.
# Near s=1 (abscissa), ζ_{S²} diverges, drowning the
# finite correction from low-j sparsity.
#
# Test design: at s=1.01, ζ_{S²} converges like ln(j_max),
# so any finite truncation gives poor absolute accuracy.
# Instead test MONOTONE approach: |ratio − 1/|H|| shrinks
# as s decreases, plus direction and non-triviality checks.
# ==========================================================

class TestSpectralZetaWeyl:
    """Spectral zeta approaches Weyl prediction for s → 1⁺."""

    @pytest.mark.parametrize("name,builder,H_order", [
        ("2T", build_binary_tetrahedral, 12),
        ("2O", build_binary_octahedral, 24),
        ("2I", build_binary_icosahedral, 60),
    ])
    def test_weyl_ratio(self, name, builder, H_order):
        """ζ_{A₀}(s)/ζ_{S²}(s) → 1/|H| monotonically as s → 1⁺.

        At finite j_max, the partial-sum ratio converges slowly to 1/|H|
        near the abscissa (s=1). Instead of an unrealistic tight tolerance,
        we verify three structural Weyl properties:
          (a) errors strictly decrease as s → 1 (monotone convergence),
          (b) ratio < 1/|H| (approaching from below — low-j deficit),
          (c) ratio is non-trivially large (not stuck at 0)."""
        grp = prepare_group(builder(), H_order)
        j_max = 2000
        expected = 1.0 / H_order

        # Precompute multiplicities once
        m_vals = [multiplicity(j, grp) for j in range(j_max + 1)]

        def zeta_ratio(s):
            z_grp = sum(m_vals[j] / (j * (j + 1)) ** s
                        for j in range(1, j_max + 1) if m_vals[j] > 0)
            z_S2 = sum((2*j+1) / (j * (j + 1)) ** s
                       for j in range(1, j_max + 1))
            return z_grp / z_S2

        s_vals = [2.0, 1.5, 1.2]
        ratios = [zeta_ratio(s) for s in s_vals]
        errors = [abs(r - expected) for r in ratios]

        # (a) Strictly monotone approach to 1/|H|
        for i in range(len(errors) - 1):
            assert errors[i + 1] < errors[i], \
                f"{name}: err(s={s_vals[i+1]})={errors[i+1]:.6f} >= " \
                f"err(s={s_vals[i]})={errors[i]:.6f}"
        # (b) Approaching from below (low-j sparsity deficit).
        # Note: this is empirically true for polyhedral groups because
        # m(A₀,j) = 0 at low j creates a deficit. Not a general theorem —
        # a group with clustered high-multiplicity levels could overshoot.
        for i, s in enumerate(s_vals):
            assert ratios[i] < expected, \
                f"{name} s={s}: ratio={ratios[i]:.6f} >= 1/|H|={expected:.6f}"
        # (c) At s=1.2, ratio > 40% of 1/|H| (non-trivially large)
        assert ratios[-1] > 0.4 * expected, \
            f"{name}: ratio(1.2)={ratios[-1]:.6f} < 40% of 1/|H|"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
