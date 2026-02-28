#!/usr/bin/env python3
"""
I30 — Generalized sector-selective duality.

QUESTION: For which (group, character) pairs does the reflection
m(χ,j) + m(χ,j*-j) = 1 hold, even when condition (B) fails globally?

ANSWER: Duality holds for (H, χ) iff χ is immune to every
(B)-violating conjugacy class. A class β with SO(3) order d
is violating iff K=|H|/d is odd. χ is immune to β iff
f(α_β, j)=0 for all j in χ's support.

For D_n (n odd): only d=2 violates. Spinorial chars (half-int j)
are immune (sin((2j+1)π/2)=0 when 2j+1 even). Tensorial (int j)
are not immune. Result: duality holds for spinorial only.

For polyhedral (T,O,I): no violating classes → all chars immune
→ duality holds globally. Recovers Theorem 3.2.

PLAN:
  §0. Identify violating classes [12 tests]
  §1. Selective duality scan: all 4 chars of D_n [15 tests]
  §2. Immunity criterion: d=2 invisible at half-int j [6 tests]
  §3. Generalized criterion: iff verification [15 tests]

STATUS: COMPLETE — 48/48 tests pass
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from src.quaternion import qkey, qmul
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
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


def build_binary_dihedral(n):
    """2D_n = Dic_n in SU(2). Order 4n.
    Generators: a = e^{iπ/n}, x = j.  Defining relations: a^{2n}=1, x²=aⁿ."""
    j_quat = np.array([0.0, 0.0, 1.0, 0.0])
    result = []
    for k in range(2 * n):
        theta = k * np.pi / n
        q = np.array([np.cos(theta), 0, 0, np.sin(theta)])
        result.append(q)
        result.append(qmul(j_quat, q))
    G = np.array(result)
    # Verify dicyclic relation: x² = aⁿ
    x = j_quat
    x2 = qmul(x, x)
    a_n_angle = n * np.pi / n  # = π
    a_n = np.array([np.cos(a_n_angle), 0, 0, np.sin(a_n_angle)])
    assert np.allclose(x2, a_n, atol=1e-12), \
        f"Dicyclic relation x²=aⁿ failed: x²={x2}, aⁿ={a_n}"
    return G


def so3_order(trace_su2):
    """SO(3) order d of element with given SU(2) trace.
    trace = 2cos(α), α = half-angle. SO(3) angle = 2α.
    d = smallest k with 2kα ∈ 2πZ, i.e., kα ∈ πZ."""
    if abs(trace_su2 - 2) < 1e-10:
        return 1  # identity
    if abs(trace_su2 + 2) < 1e-10:
        return 1  # -1 in SU(2), projects to identity in SO(3)
    alpha = np.arccos(np.clip(trace_su2 / 2, -1, 1))
    for d in range(2, 10000):
        if abs(np.sin(d * alpha)) < 1e-8:
            return d
    raise ValueError(f"SO(3) order not found for trace={trace_su2} (α={alpha})")


def violating_classes(classes, H_order):
    """Return list of classes where K = |H|/d is odd (violating (B))."""
    viol = []
    for c in classes:
        d = so3_order(c['trace'])
        if d == 1:
            continue  # identity / central
        assert H_order % d == 0, \
            f"d={d} does not divide |H|={H_order} (Lagrange violated — numerical error?)"
        K = H_order // d
        if K % 2 == 1:
            viol.append((c, d, K))
    return viol


# ==========================================================
# § 0. Violating classes
# ==========================================================

class TestViolatingClasses:
    """Identify which groups have (B)-violating classes."""

    @pytest.mark.parametrize("name,builder,H_order", [
        ("T", build_binary_tetrahedral, 12),
        ("O", build_binary_octahedral, 24),
        ("I", build_binary_icosahedral, 60),
    ])
    def test_polyhedral_no_violating(self, name, builder, H_order):
        """T, O, I satisfy (B) globally — no violating classes."""
        G = builder()
        classes = compute_conjugacy_classes(G)
        viol = violating_classes(classes, H_order)
        assert len(viol) == 0, f"{name} has violating classes: {viol}"

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11])
    def test_Dn_odd_has_violating(self, n):
        """D_n (n odd): half-turns (d=2) are the ONLY violating class.
        K = |H|/2 = n, which is odd."""
        G = build_binary_dihedral(n)
        H_order = 2 * n
        classes = compute_conjugacy_classes(G)
        viol = violating_classes(classes, H_order)
        # All violating classes should have d=2
        assert len(viol) > 0, f"D_{n} should have violating classes"
        for c, d, K in viol:
            assert d == 2, f"D_{n}: unexpected violating class with d={d}"
            assert K == n, f"D_{n}: K should be {n}, got {K}"

    @pytest.mark.parametrize("n", [2, 4, 6, 8])
    def test_Dn_even_no_violating(self, n):
        """D_n (n even): satisfies (B), no violating classes."""
        G = build_binary_dihedral(n)
        H_order = 2 * n
        classes = compute_conjugacy_classes(G)
        viol = violating_classes(classes, H_order)
        assert len(viol) == 0, f"D_{n} (even) has violating classes"


# ==========================================================
# § 1. All characters of D_n (n odd): selective duality scan
# ==========================================================

def build_Dn_characters(n, elements, classes):
    """Build all 4 scalar characters of 2D_n (n odd).
    Abel(2D_n) ≅ Z₄ for n odd.
    χ₀: a→+1, x→+1   (tensorial)
    χ₁: a→+1, x→-1   (tensorial)
    χ₂: a→-1, x→+i   (spinorial)
    χ₃: a→-1, x→-i   (spinorial)
    """
    a_powers = {}
    for k in range(2 * n):
        angle = k * np.pi / n
        a_powers[k] = np.array([np.cos(angle), 0, 0, np.sin(angle)])
    x_elem = np.array([0, 0, 1, 0], dtype=float)
    xa_powers = {}
    for k in range(2 * n):
        xa_powers[k] = qmul(x_elem, a_powers[k])

    char_defs = {
        'chi0': (1.0, 1.0 + 0j),
        'chi1': (1.0, -1.0 + 0j),
        'chi2': (-1.0, 1j),
        'chi3': (-1.0, -1j),
    }
    chars = {}
    for cname, (chi_a, chi_x) in char_defs.items():
        elem_vals = {}
        for k in range(2 * n):
            elem_vals[qkey(a_powers[k])] = chi_a ** k
        for k in range(2 * n):
            elem_vals[qkey(xa_powers[k])] = chi_x * (chi_a ** k)
        assert len(elem_vals) == 4 * n, \
            f"{cname}: qkey collision — {len(elem_vals)} keys for {4*n} elements"
        chi_on_classes = []
        for c in classes:
            rep_key = list(c['keys'])[0]
            chi_on_classes.append(elem_vals[rep_key])
        chars[cname] = np.array(chi_on_classes)
    return chars


def prepare_Dn(n):
    """Build complete group data for D_n (n odd)."""
    elements = build_binary_dihedral(n)
    G_order = len(elements)
    H_order = 2 * n
    j_star = H_order // 2 - 1  # n - 1
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    chars = build_Dn_characters(n, elements, classes)
    return {
        'n': n, 'G': G_order, 'H': H_order,
        'j_star': j_star, 'classes': classes,
        'sizes': sizes, 'chars': chars,
    }


def m_chi(j, grp, chi_name):
    """Multiplicity of character chi_name in spin-j rep."""
    chi_Vj = np.array([chi_su2(j, c['trace']) for c in grp['classes']])
    m = np.sum(grp['chars'][chi_name].conj() * chi_Vj * grp['sizes']) / grp['G']
    val = m.real
    assert abs(val - round(val)) < 1e-6, \
        f"j={j}, {chi_name}: non-integer multiplicity {val}"
    return int(round(val))


def duality_holds(grp, chi_name):
    """Check if m(χ,j) + m(χ,j*-j) = 1 for all j in χ's support.
    Only checks each (j, j*-j) pair once."""
    j_star = grp['j_star']
    is_spinorial = chi_name in ('chi2', 'chi3')
    if is_spinorial:
        # Half-integers up to midpoint: j = 1/2, 3/2, ..., ≤ j*/2
        j_vals = [k / 2 for k in range(1, j_star + 1, 2)]
    else:
        # Integers up to midpoint: j = 0, 1, ..., ⌊j*/2⌋
        j_vals = list(range(0, j_star // 2 + 1))
    for j in j_vals:
        s = m_chi(j, grp, chi_name) + m_chi(j_star - j, grp, chi_name)
        if s != 1:
            return False, j, s
    return True, None, None


class TestSelectiveDualityScan:
    """For D_n (n odd), which characters have selective duality?"""

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11])
    def test_character_orthonormality(self, n):
        """Verify Σ_g |χ(g)|² = |G| for each character (norm check)."""
        grp = prepare_Dn(n)
        for chi_name, chi_vals in grp['chars'].items():
            norm_sq = np.sum(np.abs(chi_vals)**2 * grp['sizes'])
            assert abs(norm_sq - grp['G']) < 1e-10, \
                f"D_{n} {chi_name}: ||χ||² = {norm_sq}, expected {grp['G']}"

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11])
    def test_spinorial_duality_holds(self, n):
        """Spinorial χ₂, χ₃: duality holds (known from Prop 4.7)."""
        grp = prepare_Dn(n)
        for chi in ('chi2', 'chi3'):
            ok, j_fail, s = duality_holds(grp, chi)
            assert ok, f"D_{n} {chi}: fails at j={j_fail}, sum={s}"

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11])
    def test_tensorial_duality_fails(self, n):
        """Tensorial χ₀, χ₁: duality FAILS (sums are 0 or 2, not 1).
        In particular, j=0 gives sum=2 since m(χ₀,0)=1 and m(χ₀,j*)=1."""
        grp = prepare_Dn(n)
        j_star = grp['j_star']
        for chi in ('chi0', 'chi1'):
            sums = []
            for j in range(j_star + 1):
                s = m_chi(j, grp, chi) + m_chi(j_star - j, grp, chi)
                sums.append(s)
            # Sums should be 0 or 2 (not 1)
            assert all(s in (0, 2) for s in sums), \
                f"D_{n} {chi}: unexpected sums {sums}"
            # χ₀ is trivial → m(χ₀,0)=1 and m(χ₀,j*)=1 → sum=2 at j=0.
            # χ₁ has χ₁(x)=-1 → m(χ₁,0)=0 and m(χ₁,j*)=0 → sum=0 at j=0.
            expected_j0 = 2 if chi == 'chi0' else 0
            assert sums[0] == expected_j0, \
                f"D_{n} {chi}: sum at j=0 should be {expected_j0}, got {sums[0]}"


# ==========================================================
# § 2. Immunity criterion
#
# A character χ is IMMUNE to violating class β (SO(3) order d)
# iff f(π k'/d, j) = 0 for ALL j in the support of χ.
#
# For D_n (n odd), the only violating class has d=2.
# f(π/2, j) = sin((2j+1)π/2) / sin(π/2) = sin((2j+1)π/2).
# At half-integer j: 2j+1 is even → sin(even·π/2) = 0. IMMUNE.
# At integer j: 2j+1 is odd → sin(odd·π/2) = ±1 ≠ 0. NOT IMMUNE.
#
# So: spinorial (half-int j) → immune → duality holds.
#     tensorial (int j) → not immune → duality fails.
# ==========================================================

def mode_contribution_at_class(d, j):
    """f(π k'/d, j) for the class with SO(3) order d.
    For d=2, k'=1, half-angle = π/2."""
    alpha = np.pi / d  # simplest representative: k'=1
    return np.sin((2*j + 1) * alpha) / np.sin(alpha)


class TestImmunityCriterion:
    """Verify the immunity mechanism for d=2 violating class."""

    def test_d2_invisible_at_half_integer(self):
        """f(π/2, j) = 0 at all half-integer j."""
        for j2 in range(1, 40, 2):  # j = 1/2, 3/2, ..., 39/2
            j = j2 / 2
            f = mode_contribution_at_class(2, j)
            assert abs(f) < 1e-10, f"f(π/2, {j}) = {f} ≠ 0"

    def test_d2_visible_at_integer(self):
        """f(π/2, j) ≠ 0 at all integer j."""
        for j in range(0, 20):
            f = mode_contribution_at_class(2, j)
            assert abs(f) > 0.9, f"f(π/2, {j}) = {f} ≈ 0"

    @pytest.mark.parametrize("n", [3, 5, 7, 9])
    def test_non_violating_classes_antisymmetric(self, n):
        """Non-violating classes (d≠2) satisfy mode antisymmetry:
        f(α, j) + f(α, j*-j) = 0 for all j."""
        grp = prepare_Dn(n)
        j_star = grp['j_star']
        for c in grp['classes']:
            d = so3_order(c['trace'])
            if d <= 1:
                continue
            K = grp['H'] // d
            if K % 2 == 1:
                continue  # skip violating classes
            alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
            for j in range(j_star + 1):
                f_j = np.sin((2*j+1) * alpha) / np.sin(alpha)
                f_r = np.sin((2*(j_star-j)+1) * alpha) / np.sin(alpha)
                assert abs(f_j + f_r) < 1e-10, \
                    f"D_{n} d={d}: f({j}) + f({j_star-j}) = {f_j+f_r}"


# ==========================================================
# § 3. Generalized criterion and verification
#
# THEOREM (generalized selective duality):
# Let H ⊂ SO(3) with |H| even. A scalar character χ of 2H
# satisfies the spectral reflection m(χ,j) + m(χ,j*-j) = 1
# for all j in its support iff χ is immune to every
# (B)-violating conjugacy class of H.
#
# A class β with SO(3) order d is (B)-violating iff K=|H|/d is odd.
# χ is immune to β iff f(α_β, j) = 0 for all j with χ(g)·χ(g)* ≠ 0
# in that j-sector. For D_n (n odd), this reduces to:
#   spinorial (half-int j) → immune to d=2 → duality holds
#   tensorial (int j) → not immune to d=2 → duality fails
#
# For polyhedral groups (T, O, I): no violating classes exist,
# so ALL characters are trivially immune → duality holds globally.
# This recovers Theorem 3.2 as a special case.
# ==========================================================

class TestGeneralizedCriterion:
    """Verify the generalized selective duality criterion."""

    @pytest.mark.parametrize("n", [3, 5, 7, 9, 11, 13, 15])
    def test_criterion_iff(self, n):
        """For D_n (n odd): duality ⟺ immunity to d=2 class.
        Spinorial → immune → duality. Tensorial → not immune → no duality."""
        grp = prepare_Dn(n)
        j_star = grp['j_star']

        # Spinorial: immune to d=2, duality holds
        for chi in ('chi2', 'chi3'):
            ok, j_fail, s = duality_holds(grp, chi)
            assert ok, f"D_{n} {chi}: duality should hold but fails at j={j_fail}"

        # Tensorial: not immune to d=2, duality fails
        for chi in ('chi0', 'chi1'):
            ok, j_fail, s = duality_holds(grp, chi)
            assert not ok, f"D_{n} {chi}: duality should fail but holds"

    @pytest.mark.parametrize("n", [2, 4, 6, 8, 10])
    def test_Dn_even_all_duality(self, n):
        """D_n (n even): (B) holds, ALL 4 characters have duality.
        Abel(2D_n) ≅ Z₂×Z₂ for n even, all tensorial."""
        elements = build_binary_dihedral(n)
        G_order = len(elements)
        H_order = 2 * n
        j_star = H_order // 2 - 1
        classes = compute_conjugacy_classes(elements)
        classes.sort(key=lambda c: -c['trace'])
        sizes = np.array([c['size'] for c in classes], dtype=float)

        # Build all 4 characters of Z₂×Z₂ = Dic_n/[Dic_n,Dic_n]
        # For n even: a²ⁿ=1, x²=aⁿ, xax⁻¹=a⁻¹
        # [Dic_n, Dic_n] = ⟨a²⟩, Abel = {1,a,x,ax} with a²=1, x²=1
        # χ(a) ∈ {±1}, χ(x) ∈ {±1}
        a_powers = {}
        for k in range(2 * n):
            angle = k * np.pi / n
            a_powers[k] = np.array([np.cos(angle), 0, 0, np.sin(angle)])
        x_elem = np.array([0, 0, 1, 0], dtype=float)

        char_defs = [
            (1.0, 1.0),   # χ₀: a→+1, x→+1
            (1.0, -1.0),  # χ₁: a→+1, x→-1
            (-1.0, 1.0),  # χ₂: a→-1, x→+1
            (-1.0, -1.0), # χ₃: a→-1, x→-1
        ]

        for ci, (chi_a, chi_x) in enumerate(char_defs):
            # Build character values on elements
            elem_vals = {}
            for k in range(2 * n):
                elem_vals[qkey(a_powers[k])] = chi_a ** k
            for k in range(2 * n):
                xak = qmul(x_elem, a_powers[k])
                elem_vals[qkey(xak)] = chi_x * (chi_a ** k)
            assert len(elem_vals) == G_order, \
                f"D_{n} χ_{ci}: qkey collision — {len(elem_vals)} keys for {G_order} elements"

            # Character on classes
            chi_on_classes = np.array([
                elem_vals[list(c['keys'])[0]] for c in classes
            ])

            # Norm check: Σ_g |χ(g)|² = |G|
            norm_sq = np.sum(np.abs(chi_on_classes)**2 * sizes)
            assert abs(norm_sq - G_order) < 1e-10, \
                f"D_{n} χ_{ci}: ||χ||² = {norm_sq}, expected {G_order}"

            for j in range(j_star + 1):
                chi_Vj = np.array([chi_su2(j, c['trace'])
                                   for c in classes])
                chi_Vr = np.array([chi_su2(j_star - j, c['trace'])
                                   for c in classes])
                m_j = int(round((np.sum(chi_on_classes.conj() * chi_Vj
                                        * sizes) / G_order).real))
                m_r = int(round((np.sum(chi_on_classes.conj() * chi_Vr
                                        * sizes) / G_order).real))
                assert m_j + m_r == 1, \
                    f"D_{n} χ_{ci}: m({j})+m({j_star-j})={m_j+m_r}"

    @pytest.mark.parametrize("name,builder,H_order,n_chars", [
        ("T", build_binary_tetrahedral, 12, 3),   # Abel ≅ Z₃
        ("O", build_binary_octahedral, 24, 2),     # Abel ≅ Z₂
        ("I", build_binary_icosahedral, 60, 1),    # perfect group
    ])
    def test_polyhedral_all_characters(self, name, builder, H_order, n_chars):
        """All scalar characters of T, O, I satisfy duality.
        Builds chars from abelianization G/[G,G]."""
        G = builder()
        j_star = H_order // 2 - 1
        classes = compute_conjugacy_classes(G)
        classes.sort(key=lambda c: -c['trace'])
        sizes = np.array([c['size'] for c in classes], dtype=float)
        order = len(G)

        # Build all scalar chars via abelianization
        comm = compute_commutator_subgroup(G)
        comm_keys = set(comm.keys())

        # Assign coset labels: coset 0 = [G,G]
        from src.quaternion import qinv
        coset_reps = [G[0]]  # identity
        coset_of = {}
        coset_of[qkey(G[0])] = 0
        for g in G:
            gk = qkey(g)
            if gk in coset_of:
                continue
            found = False
            for ci, rep in enumerate(coset_reps):
                if qkey(qmul(g, qinv(rep))) in comm_keys:
                    coset_of[gk] = ci
                    found = True
                    break
            if not found:
                coset_of[gk] = len(coset_reps)
                coset_reps.append(g)
        n_ab = len(coset_reps)
        assert len(coset_of) == len(G), \
            f"{name}: only {len(coset_of)}/{len(G)} elements assigned cosets"
        assert n_ab == n_chars, \
            f"{name}: |G^ab|={n_ab}, expected {n_chars}"

        # Coset label per class (from representative)
        class_cosets = [coset_of[qkey(c['rep'])] for c in classes]

        # Characters: χ_k(class) = ω^{k·coset}, ω = e^{2πi/n_ab}.
        # This assumes G^ab is cyclic — true for 2T (Z₃), 2O (Z₂), 2I (trivial).
        # Would need multi-generator construction for non-cyclic (e.g. Z₂×Z₂).
        omega = np.exp(2j * np.pi / n_ab) if n_ab > 1 else 1.0

        for k in range(n_ab):
            chi = np.array([omega ** (k * cl) for cl in class_cosets])
            for j in range(j_star + 1):
                chi_j = np.array([chi_su2(j, c['trace']) for c in classes])
                chi_r = np.array([chi_su2(j_star - j, c['trace'])
                                  for c in classes])
                m_j = np.sum(chi.conj() * chi_j * sizes) / order
                m_r = np.sum(chi.conj() * chi_r * sizes) / order
                assert abs(m_j.real - round(m_j.real)) < 1e-6, \
                    f"{name} A_{k}: non-integer m({j})={m_j.real}"
                s = int(round(m_j.real)) + int(round(m_r.real))
                assert s == 1, \
                    f"{name} A_{k}: m({j})+m({j_star-j})={s}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
