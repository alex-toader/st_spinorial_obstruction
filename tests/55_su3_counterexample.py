#!/usr/bin/env python3
"""
Script 55: SU(3) negative result — Δ(27) has no complementary pairing.

Δ(27) = Heisenberg group over F₃, order 27, embedded in SU(3).
Elements: (p, q, r) ∈ Z₃³, multiplication (p₁+p₂, q₁+q₂+r₁p₂, r₁+r₂) mod 3.

Structure:
  - 11 conjugacy classes, sizes [1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3]
  - [G,G] = Z₃ (center), G^ab = Z₃ × Z₃ → 9 one-dimensional characters
  - 2 three-dimensional irreps (faithful)
  - Irrep dimensions: 9×1 + 2×3 = 9 + 18 = 27 ✓

SU(3) irreps V_{(p,q)} have dim = (p+1)(q+1)(p+q+2)/2 (degree 2 polynomial).
Characters computed via Jacobi-Trudi identity:
  χ_{(p,q)}(M) = h_{p+q} · h_q − h_{p+q+1} · h_{q−1}
where h_k are complete homogeneous symmetric polynomials of eigenvalues of M,
computed via recurrence h_k = e₁·h_{k-1} − e₂·h_{k-2} + e₃·h_{k-3}.

Central elements (0,q,0) → M = ω^q·I: handled analytically via N-ality,
  χ_{(p,q)}(ω^a·I) = ω^{a(p−q)} · dim V_{(p,q)}.

Key results:
  1. No complementary pairing m(χ,(p,q)) + m(χ,(P-p,Q-q)) = const exists
     at any threshold (P,Q) for any 1D character.
  2. All 8 non-trivial characters have identical multiplicities.
  3. N-ality selection: m ≠ 0 only when p−q ≡ 0 mod 3.
  4. Structural reason: dim V_{(p,q)} is degree 2, so no affine involution
     gives constant dimension sum.

References: main_v6.tex §5 (Discussion), tracker v5 item 0f.

Run: /usr/bin/python3 -m pytest tests/55_su3_counterexample.py -v
Date: Feb 2026

RAW OUTPUT (318 passed, 0.75s):
  TestGroupStructure::test_order PASSED
  TestGroupStructure::test_closure PASSED
  TestGroupStructure::test_inverses PASSED
  TestGroupStructure::test_associativity_sample PASSED
  TestGroupStructure::test_conjugacy_classes PASSED
  TestGroupStructure::test_center_is_Z3 PASSED
  TestGroupStructure::test_commutator_subgroup PASSED
  TestGroupStructure::test_nine_characters_orthonormal PASSED
  TestSU3Embedding::test_matrices_in_SU3 PASSED
  TestSU3Embedding::test_homomorphism PASSED
  TestSU3Embedding::test_faithful PASSED
  TestJacobiTrudiCharacter::test_dimension_at_identity PASSED
  TestJacobiTrudiCharacter::test_central_nality PASSED
  TestJacobiTrudiCharacter::test_multiplicities_are_nonneg_integers PASSED
  TestJacobiTrudiCharacter::test_dimension_decomposition PASSED
  TestJacobiTrudiCharacter::test_class_function PASSED
  TestHardcodedMultiplicities::test_trivial_table PASSED
  TestHardcodedMultiplicities::test_nontrivial_table PASSED
  TestHardcodedMultiplicities::test_trivial_pq_symmetry PASSED
  TestHardcodedMultiplicities::test_nontrivial_pq_symmetry PASSED
  TestAllNontrivialIdentical::test_equals_chi10[char_ab0..6] 7 PASSED
  TestNalitySelection::test_trivial_nality PASSED
  TestNalitySelection::test_nontrivial_nality PASSED
  TestNoPairing2D::test_trivial_no_constant_sum[threshold0..48] 49 PASSED
  TestNoPairing2D::test_nontrivial_no_constant_sum[threshold0..48] 49 PASSED
  TestNoPairing1D::test_trivial_along_p0[2..8] 7 PASSED
  TestNoPairing1D::test_trivial_along_0q[2..8] 7 PASSED
  TestNoPairing1D::test_trivial_along_diagonal[2..8] 7 PASSED
  TestNoPairing1D::test_nontrivial_along_p0[2..8] 7 PASSED
  TestNoPairing1D::test_nontrivial_along_0q[2..8] 7 PASSED
  TestNoPairing1D::test_nontrivial_along_diagonal[2..8] 7 PASSED
  TestNoAffineInvolution::test_trivial_no_affine_pairing[center0..48] 49 PASSED
  TestNoAffineInvolution::test_nontrivial_no_affine_pairing[center0..48] 49 PASSED
  TestDimensionSumNotConstant::test_dim_sum_varies[threshold0..48] 49 PASSED
  TestComparisonWithSU2::test_su2_dim_sum_constant PASSED
  TestComparisonWithSU2::test_su3_dim_sum_not_constant PASSED
"""
import numpy as np
import pytest


# ════════════════════════════════════════════════════════════
# Δ(27) group operations
# ════════════════════════════════════════════════════════════

ELEMENTS = [(p, q, r) for p in range(3) for q in range(3) for r in range(3)]
OMEGA = np.exp(2j * np.pi / 3)


def grp_mult(g1, g2):
    """Heisenberg multiplication: (p₁+p₂, q₁+q₂+r₁p₂, r₁+r₂) mod 3."""
    p1, q1, r1 = g1
    p2, q2, r2 = g2
    return ((p1 + p2) % 3, (q1 + q2 + r1 * p2) % 3, (r1 + r2) % 3)


def grp_inv(g):
    """Inverse in Δ(27)."""
    p, q, r = g
    return ((-p) % 3, (-(q - r * p)) % 3, (-r) % 3)


def conjugate(g, h):
    """h g h⁻¹."""
    return grp_mult(grp_mult(h, g), grp_inv(h))


def conjugacy_classes():
    """Compute conjugacy classes of Δ(27)."""
    assigned = set()
    classes = []
    for g in ELEMENTS:
        if g in assigned:
            continue
        cls = set()
        for h in ELEMENTS:
            cls.add(conjugate(g, h))
        classes.append(sorted(cls))
        assigned |= cls
    return classes


# ════════════════════════════════════════════════════════════
# 1D characters and SU(3) matrix representation
# ════════════════════════════════════════════════════════════

def chi_1d(a, b, g):
    """1D character χ_{a,b}(p,q,r) = ω^{ap+br}."""
    p, q, r = g
    return OMEGA ** ((a * p + b * r) % 3)


def matrix_rep(g):
    """Faithful 3D rep in SU(3): M = (ω^q · P^r · D^p)^T.

    Transpose converts the natural anti-homomorphism into a proper
    homomorphism: ρ(g₁g₂) = ρ(g₁)ρ(g₂). Traces are unchanged.
    """
    p, q, r = g
    D = np.diag([1.0 + 0j, OMEGA ** p, OMEGA ** (2 * p)])
    P = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        P[(i + r) % 3, i] = 1.0
    return (OMEGA ** q * (P @ D)).T


# ════════════════════════════════════════════════════════════
# SU(3) character via Jacobi-Trudi identity
# ════════════════════════════════════════════════════════════

def su3_character(pw, qw, g):
    """
    Character of V_{(pw, qw)} at group element g.

    Central elements (0,q,0): analytic formula via N-ality.
    Non-central: Jacobi-Trudi identity (no Vandermonde denominator).

    Jacobi-Trudi for partition (pw+qw, qw, 0):
      s = h_{pw+qw} · h_{qw} − h_{pw+qw+1} · h_{qw−1}
    """
    p, q, r = g
    dim = (pw + 1) * (qw + 1) * (pw + qw + 2) // 2

    # Central elements: M = ω^q · I
    if p == 0 and r == 0:
        return dim * OMEGA ** (q * (pw - qw))

    # Non-central: eigenvalues are distinct
    M = matrix_rep(g)
    eigs = np.linalg.eigvals(M)
    z1, z2, z3 = eigs

    e1 = z1 + z2 + z3
    e2 = z1 * z2 + z1 * z3 + z2 * z3
    e3 = z1 * z2 * z3  # = det(M) = 1

    max_k = pw + qw + 2
    h = [0.0 + 0j] * (max_k + 1)
    h[0] = 1.0
    for k in range(1, max_k + 1):
        h[k] = e1 * h[k - 1]
        if k >= 2:
            h[k] -= e2 * h[k - 2]
        if k >= 3:
            h[k] += e3 * h[k - 3]

    h_qm1 = h[qw - 1] if qw >= 1 else 0.0
    return h[pw + qw] * h[qw] - h[pw + qw + 1] * h_qm1


def multiplicity(a, b, pw, qw):
    """m(χ_{a,b}, V_{(pw,qw)}) = (1/27) Σ_g χ̄_{a,b}(g) · χ_{(pw,qw)}(g)."""
    total = 0.0 + 0j
    for g in ELEMENTS:
        total += np.conj(chi_1d(a, b, g)) * su3_character(pw, qw, g)
    m = total / 27.0
    assert abs(m.imag) < 1e-8, f"multiplicity has imaginary part {m.imag:.2e}"
    return int(round(m.real))


# ════════════════════════════════════════════════════════════
# Precomputed multiplicity tables (verified by computation)
# ════════════════════════════════════════════════════════════

# m(A₀, (p,q)) — trivial character, grid [0,8]²
MULT_TRIVIAL = {
    (0,0):1,  (0,1):0, (0,2):0, (0,3):2,  (0,4):0,  (0,5):0,  (0,6):4,  (0,7):0,  (0,8):0,
    (1,0):0,  (1,1):0, (1,2):0, (1,3):0,  (1,4):3,  (1,5):0,  (1,6):0,  (1,7):8,  (1,8):0,
    (2,0):0,  (2,1):0, (2,2):3, (2,3):0,  (2,4):0,  (2,5):9,  (2,6):0,  (2,7):0,  (2,8):18,
    (3,0):2,  (3,1):0, (3,2):0, (3,3):8,  (3,4):0,  (3,5):0,  (3,6):18, (3,7):0,  (3,8):0,
    (4,0):0,  (4,1):3, (4,2):0, (4,3):0,  (4,4):13, (4,5):0,  (4,6):0,  (4,7):28, (4,8):0,
    (5,0):0,  (5,1):0, (5,2):9, (5,3):0,  (5,4):0,  (5,5):24, (5,6):0,  (5,7):0,  (5,8):45,
    (6,0):4,  (6,1):0, (6,2):0, (6,3):18, (6,4):0,  (6,5):0,  (6,6):39, (6,7):0,  (6,8):0,
    (7,0):0,  (7,1):8, (7,2):0, (7,3):0,  (7,4):28, (7,5):0,  (7,6):0,  (7,7):56, (7,8):0,
    (8,0):0,  (8,1):0, (8,2):18,(8,3):0,  (8,4):0,  (8,5):45, (8,6):0,  (8,7):0,  (8,8):81,
}

# m(χ_{1,0}, (p,q)) — non-trivial character (identical for all 8 non-trivial)
MULT_NONTRIVIAL = {
    (0,0):0,  (0,1):0, (0,2):0, (0,3):1,  (0,4):0,  (0,5):0,  (0,6):3,  (0,7):0,  (0,8):0,
    (1,0):0,  (1,1):1, (1,2):0, (1,3):0,  (1,4):4,  (1,5):0,  (1,6):0,  (1,7):9,  (1,8):0,
    (2,0):0,  (2,1):0, (2,2):3, (2,3):0,  (2,4):0,  (2,5):9,  (2,6):0,  (2,7):0,  (2,8):18,
    (3,0):1,  (3,1):0, (3,2):0, (3,3):7,  (3,4):0,  (3,5):0,  (3,6):17, (3,7):0,  (3,8):0,
    (4,0):0,  (4,1):4, (4,2):0, (4,3):0,  (4,4):14, (4,5):0,  (4,6):0,  (4,7):29, (4,8):0,
    (5,0):0,  (5,1):0, (5,2):9, (5,3):0,  (5,4):0,  (5,5):24, (5,6):0,  (5,7):0,  (5,8):45,
    (6,0):3,  (6,1):0, (6,2):0, (6,3):17, (6,4):0,  (6,5):0,  (6,6):38, (6,7):0,  (6,8):0,
    (7,0):0,  (7,1):9, (7,2):0, (7,3):0,  (7,4):29, (7,5):0,  (7,6):0,  (7,7):57, (7,8):0,
    (8,0):0,  (8,1):0, (8,2):18,(8,3):0,  (8,4):0,  (8,5):45, (8,6):0,  (8,7):0,  (8,8):81,
}

P_MAX = 8


# ════════════════════════════════════════════════════════════
# Test classes
# ════════════════════════════════════════════════════════════

class TestGroupStructure:
    """Verify Δ(27) group properties."""

    def test_order(self):
        assert len(ELEMENTS) == 27

    def test_closure(self):
        for g1 in ELEMENTS:
            for g2 in ELEMENTS:
                assert grp_mult(g1, g2) in ELEMENTS

    def test_inverses(self):
        for g in ELEMENTS:
            assert grp_mult(g, grp_inv(g)) == (0, 0, 0)

    def test_associativity_sample(self):
        """Check associativity on a sample (full check is 27³ = 19683)."""
        sample = [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,1,1), (2,1,0)]
        for a in sample:
            for b in sample:
                for c in sample:
                    assert grp_mult(grp_mult(a, b), c) == grp_mult(a, grp_mult(b, c))

    def test_conjugacy_classes(self):
        classes = conjugacy_classes()
        assert len(classes) == 11
        sizes = sorted(len(c) for c in classes)
        assert sizes == [1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3]

    def test_center_is_Z3(self):
        """Center = {(0,q,0) : q ∈ Z₃}, isomorphic to Z₃."""
        center = []
        for g in ELEMENTS:
            if all(grp_mult(g, h) == grp_mult(h, g) for h in ELEMENTS):
                center.append(g)
        assert len(center) == 3
        assert all(g[0] == 0 and g[2] == 0 for g in center)

    def test_commutator_subgroup(self):
        """[G,G] = center = Z₃, so |G^ab| = 9."""
        commutators = set()
        for g in ELEMENTS:
            for h in ELEMENTS:
                c = grp_mult(grp_mult(g, h), grp_inv(grp_mult(h, g)))
                commutators.add(c)
        # Close under multiplication
        comm = set(commutators)
        changed = True
        while changed:
            changed = False
            new = set()
            for a in comm:
                for b in comm:
                    ab = grp_mult(a, b)
                    if ab not in comm:
                        new.add(ab)
            if new:
                comm |= new
                changed = True
        assert len(comm) == 3  # [G,G] = Z₃
        assert 27 // len(comm) == 9  # |G^ab| = 9 → 9 one-dim characters

    def test_nine_characters_orthonormal(self):
        """9 characters from G^ab ≅ Z₃ × Z₃, verify orthonormality."""
        classes = conjugacy_classes()
        sizes = [len(c) for c in classes]
        reps = [c[0] for c in classes]
        for a1, b1 in [(a, b) for a in range(3) for b in range(3)]:
            for a2, b2 in [(a, b) for a in range(3) for b in range(3)]:
                inner = sum(
                    np.conj(chi_1d(a1, b1, reps[i])) * chi_1d(a2, b2, reps[i]) * sizes[i]
                    for i in range(len(classes))
                ) / 27.0
                expected = 1.0 if (a1, b1) == (a2, b2) else 0.0
                assert abs(inner - expected) < 1e-10, (
                    f"⟨χ_{a1}{b1}, χ_{a2}{b2}⟩ = {inner:.6f}, expected {expected}")


class TestSU3Embedding:
    """Verify the faithful 3D representation in SU(3)."""

    def test_matrices_in_SU3(self):
        """All matrices have det = 1 and are unitary."""
        for g in ELEMENTS:
            M = matrix_rep(g)
            assert abs(np.linalg.det(M) - 1.0) < 1e-10, f"det ≠ 1 at {g}"
            assert np.allclose(M @ M.conj().T, np.eye(3), atol=1e-10), f"not unitary at {g}"

    def test_homomorphism(self):
        """ρ(g₁·g₂) = ρ(g₁)·ρ(g₂) for all pairs."""
        for g1 in ELEMENTS:
            M1 = matrix_rep(g1)
            for g2 in ELEMENTS:
                M2 = matrix_rep(g2)
                Mprod = matrix_rep(grp_mult(g1, g2))
                assert np.allclose(M1 @ M2, Mprod, atol=1e-10), (
                    f"ρ({g1})·ρ({g2}) ≠ ρ({grp_mult(g1, g2)})")

    def test_faithful(self):
        """Only identity maps to I."""
        for g in ELEMENTS:
            M = matrix_rep(g)
            if np.allclose(M, np.eye(3), atol=1e-10):
                assert g == (0, 0, 0), f"Non-identity {g} maps to I"


class TestJacobiTrudiCharacter:
    """Verify SU(3) character computation."""

    def test_dimension_at_identity(self):
        """χ_{(p,q)}(e) = dim V_{(p,q)}."""
        e = (0, 0, 0)
        for pw in range(P_MAX + 1):
            for qw in range(P_MAX + 1):
                chi = su3_character(pw, qw, e)
                dim = (pw + 1) * (qw + 1) * (pw + qw + 2) // 2
                assert abs(chi - dim) < 1e-8, (
                    f"χ({pw},{qw})(e) = {chi.real:.1f}, expected {dim}")

    def test_central_nality(self):
        """χ_{(p,q)}(ω^a·I) = ω^{a(p-q)} · dim."""
        for pw in range(P_MAX + 1):
            for qw in range(P_MAX + 1):
                dim = (pw + 1) * (qw + 1) * (pw + qw + 2) // 2
                for a in range(3):
                    g = (0, a, 0)
                    chi = su3_character(pw, qw, g)
                    expected = dim * OMEGA ** (a * (pw - qw))
                    assert abs(chi - expected) < 1e-8, (
                        f"χ({pw},{qw})(center a={a}): {chi:.4f} ≠ {expected:.4f}")

    def test_multiplicities_are_nonneg_integers(self):
        """All multiplicities are non-negative integers (all 9 characters)."""
        for pw in range(P_MAX + 1):
            for qw in range(P_MAX + 1):
                for a in range(3):
                    for b in range(3):
                        m = multiplicity(a, b, pw, qw)
                        assert m >= 0, f"Negative m({a},{b}, ({pw},{qw})) = {m}"

    def test_dimension_decomposition(self):
        """Sum of 1D mults × 1 + 3D mults × 3 = dim V."""
        for pw in range(6):
            for qw in range(6):
                dim = (pw + 1) * (qw + 1) * (pw + qw + 2) // 2
                sum_1d = sum(
                    multiplicity(a, b, pw, qw)
                    for a in range(3) for b in range(3))
                remainder = dim - sum_1d
                assert remainder >= 0, f"({pw},{qw}): 1D part {sum_1d} > dim {dim}"
                assert remainder % 3 == 0, (
                    f"({pw},{qw}): 3D remainder {remainder} not divisible by 3")

    def test_class_function(self):
        """Character is constant on each conjugacy class."""
        classes = conjugacy_classes()
        for pw in range(6):
            for qw in range(6):
                for cls in classes:
                    vals = [su3_character(pw, qw, g) for g in cls]
                    for v in vals[1:]:
                        assert abs(v - vals[0]) < 1e-8, (
                            f"χ({pw},{qw}) not constant on class of {cls[0]}: "
                            f"{vals[0]:.4f} vs {v:.4f}")


class TestHardcodedMultiplicities:
    """Verify precomputed tables against live computation."""

    def test_trivial_table(self):
        for (pw, qw), expected in MULT_TRIVIAL.items():
            m = multiplicity(0, 0, pw, qw)
            assert m == expected, (
                f"m(A₀, ({pw},{qw})) = {m}, expected {expected}")

    def test_nontrivial_table(self):
        for (pw, qw), expected in MULT_NONTRIVIAL.items():
            m = multiplicity(1, 0, pw, qw)
            assert m == expected, (
                f"m(χ_{{1,0}}, ({pw},{qw})) = {m}, expected {expected}")

    def test_trivial_pq_symmetry(self):
        """m(A₀, (p,q)) = m(A₀, (q,p)) — from V_{(p,q)}* ≅ V_{(q,p)}."""
        for pw in range(P_MAX + 1):
            for qw in range(pw + 1, P_MAX + 1):
                assert MULT_TRIVIAL[(pw, qw)] == MULT_TRIVIAL[(qw, pw)], (
                    f"m(A₀, ({pw},{qw})) = {MULT_TRIVIAL[(pw,qw)]} ≠ "
                    f"m(A₀, ({qw},{pw})) = {MULT_TRIVIAL[(qw,pw)]}")

    def test_nontrivial_pq_symmetry(self):
        """m(χ, (p,q)) = m(χ, (q,p)) — from duality + all non-trivial chars identical."""
        for pw in range(P_MAX + 1):
            for qw in range(pw + 1, P_MAX + 1):
                assert MULT_NONTRIVIAL[(pw, qw)] == MULT_NONTRIVIAL[(qw, pw)], (
                    f"m(χ_{{1,0}}, ({pw},{qw})) = {MULT_NONTRIVIAL[(pw,qw)]} ≠ "
                    f"m(χ_{{1,0}}, ({qw},{pw})) = {MULT_NONTRIVIAL[(qw,pw)]}")


class TestAllNontrivialIdentical:
    """All 8 non-trivial characters have the same multiplicities."""

    @pytest.fixture(params=[(a, b) for a in range(3) for b in range(3)
                            if (a, b) != (0, 0) and (a, b) != (1, 0)])
    def char_ab(self, request):
        return request.param

    def test_equals_chi10(self, char_ab):
        a, b = char_ab
        for pw in range(P_MAX + 1):
            for qw in range(P_MAX + 1):
                m_ab = multiplicity(a, b, pw, qw)
                m_10 = MULT_NONTRIVIAL[(pw, qw)]
                assert m_ab == m_10, (
                    f"m(χ_{{{a},{b}}}, ({pw},{qw})) = {m_ab} ≠ "
                    f"m(χ_{{1,0}}, ({pw},{qw})) = {m_10}")


class TestNalitySelection:
    """Multiplicities vanish outside the N-ality 0 sector."""

    def test_trivial_nality(self):
        """m(A₀, (p,q)) = 0 when (p-q) mod 3 ≠ 0."""
        for pw in range(P_MAX + 1):
            for qw in range(P_MAX + 1):
                if (pw - qw) % 3 != 0:
                    assert MULT_TRIVIAL[(pw, qw)] == 0, (
                        f"m(A₀, ({pw},{qw})) = {MULT_TRIVIAL[(pw,qw)]} ≠ 0 "
                        f"but (p-q) mod 3 = {(pw-qw)%3}")

    def test_nontrivial_nality(self):
        """m(χ_{1,0}, (p,q)) = 0 when (p-q) mod 3 ≠ 0."""
        for pw in range(P_MAX + 1):
            for qw in range(P_MAX + 1):
                if (pw - qw) % 3 != 0:
                    assert MULT_NONTRIVIAL[(pw, qw)] == 0, (
                        f"m(χ_{{1,0}}, ({pw},{qw})) = {MULT_NONTRIVIAL[(pw,qw)]} ≠ 0 "
                        f"but (p-q) mod 3 = {(pw-qw)%3}")


class TestNoPairing2D:
    """No complementary pairing on the 2D weight lattice."""

    @pytest.fixture(params=[(P, Q) for P in range(2, 9) for Q in range(2, 9)])
    def threshold(self, request):
        return request.param

    def test_trivial_no_constant_sum(self, threshold):
        """m(A₀,(p,q)) + m(A₀,(P-p,Q-q)) is NOT constant."""
        P, Q = threshold
        sums = set()
        for pw in range(P + 1):
            for qw in range(Q + 1):
                s = MULT_TRIVIAL[(pw, qw)] + MULT_TRIVIAL[(P - pw, Q - qw)]
                sums.add(s)
        assert len(sums) > 1, (
            f"Trivial char: UNEXPECTED constant sum {sums.pop()} at (P,Q)=({P},{Q})")

    def test_nontrivial_no_constant_sum(self, threshold):
        """m(χ,(p,q)) + m(χ,(P-p,Q-q)) is NOT constant."""
        P, Q = threshold
        sums = set()
        for pw in range(P + 1):
            for qw in range(Q + 1):
                s = MULT_NONTRIVIAL[(pw, qw)] + MULT_NONTRIVIAL[(P - pw, Q - qw)]
                sums.add(s)
        assert len(sums) > 1, (
            f"Non-trivial char: UNEXPECTED constant sum {sums.pop()} "
            f"at (P,Q)=({P},{Q})")


class TestNoPairing1D:
    """No complementary pairing along 1D slices."""

    @pytest.fixture(params=range(2, 9))
    def P(self, request):
        return request.param

    def test_trivial_along_p0(self, P):
        """m(A₀,(p,0)) + m(A₀,(P-p,0)) not constant."""
        sums = [MULT_TRIVIAL[(p, 0)] + MULT_TRIVIAL[(P - p, 0)]
                for p in range(P + 1)]
        assert len(set(sums)) > 1, (
            f"Trivial along (p,0): UNEXPECTED constant {sums[0]} at P={P}")

    def test_trivial_along_0q(self, P):
        """m(A₀,(0,q)) + m(A₀,(0,P-q)) not constant."""
        sums = [MULT_TRIVIAL[(0, q)] + MULT_TRIVIAL[(0, P - q)]
                for q in range(P + 1)]
        assert len(set(sums)) > 1, (
            f"Trivial along (0,q): UNEXPECTED constant {sums[0]} at Q={P}")

    def test_trivial_along_diagonal(self, P):
        """m(A₀,(p,p)) + m(A₀,(P-p,P-p)) not constant.

        Diagonal (p,p) = self-conjugate representations — the closest
        SU(3) analogue to the SU(2) weight parameter j.
        """
        sums = [MULT_TRIVIAL[(p, p)] + MULT_TRIVIAL[(P - p, P - p)]
                for p in range(P + 1)]
        assert len(set(sums)) > 1, (
            f"Trivial along diagonal: UNEXPECTED constant {sums[0]} at P={P}")

    def test_nontrivial_along_p0(self, P):
        """m(χ,(p,0)) + m(χ,(P-p,0)) not constant."""
        sums = [MULT_NONTRIVIAL[(p, 0)] + MULT_NONTRIVIAL[(P - p, 0)]
                for p in range(P + 1)]
        if any(s > 0 for s in sums):
            assert len(set(sums)) > 1, (
                f"Non-trivial along (p,0): UNEXPECTED constant {sums[0]} at P={P}")

    def test_nontrivial_along_0q(self, P):
        """m(χ,(0,q)) + m(χ,(0,P-q)) not constant."""
        sums = [MULT_NONTRIVIAL[(0, q)] + MULT_NONTRIVIAL[(0, P - q)]
                for q in range(P + 1)]
        if any(s > 0 for s in sums):
            assert len(set(sums)) > 1, (
                f"Non-trivial along (0,q): UNEXPECTED constant {sums[0]} at Q={P}")

    def test_nontrivial_along_diagonal(self, P):
        """m(χ,(p,p)) + m(χ,(P-p,P-p)) not constant."""
        sums = [MULT_NONTRIVIAL[(p, p)] + MULT_NONTRIVIAL[(P - p, P - p)]
                for p in range(P + 1)]
        if any(s > 0 for s in sums):
            assert len(set(sums)) > 1, (
                f"Non-trivial along diagonal: UNEXPECTED constant {sums[0]} at P={P}")


class TestDimensionSumNotConstant:
    """dim V_{(p,q)} + dim V_{(P-p,Q-q)} is never constant (degree 2)."""

    @pytest.fixture(params=[(P, Q) for P in range(2, 9) for Q in range(2, 9)])
    def threshold(self, request):
        return request.param

    def test_dim_sum_varies(self, threshold):
        P, Q = threshold
        sums = set()
        for pw in range(P + 1):
            for qw in range(Q + 1):
                d1 = (pw + 1) * (qw + 1) * (pw + qw + 2) // 2
                d2 = (P - pw + 1) * (Q - qw + 1) * (P - pw + Q - qw + 2) // 2
                sums.add(d1 + d2)
        assert len(sums) > 1, (
            f"UNEXPECTED constant dim sum at (P,Q)=({P},{Q})")


class TestNoAffineInvolution:
    """No affine involution (p,q) → (a-p, b-q) gives constant sum on the support.

    Stronger than TestNoPairing2D: tests only points where m > 0,
    so even a pairing that works on a sparse subset would be caught.
    """

    @pytest.fixture(params=[(a, b) for a in range(2, 9) for b in range(2, 9)])
    def center(self, request):
        return request.param

    def test_trivial_no_affine_pairing(self, center):
        a, b = center
        seen = set()
        sums = []
        for pw in range(P_MAX + 1):
            for qw in range(P_MAX + 1):
                if MULT_TRIVIAL[(pw, qw)] == 0:
                    continue
                rp, rq = a - pw, b - qw
                if 0 <= rp <= P_MAX and 0 <= rq <= P_MAX:
                    pair = frozenset([(pw, qw), (rp, rq)])
                    if pair not in seen:
                        seen.add(pair)
                        sums.append(MULT_TRIVIAL[(pw, qw)] + MULT_TRIVIAL[(rp, rq)])
        if len(sums) >= 2:
            assert len(set(sums)) > 1, (
                f"Trivial: UNEXPECTED constant sum {sums[0]} on support "
                f"under (p,q)→({a}-p,{b}-q)")

    def test_nontrivial_no_affine_pairing(self, center):
        a, b = center
        seen = set()
        sums = []
        for pw in range(P_MAX + 1):
            for qw in range(P_MAX + 1):
                if MULT_NONTRIVIAL[(pw, qw)] == 0:
                    continue
                rp, rq = a - pw, b - qw
                if 0 <= rp <= P_MAX and 0 <= rq <= P_MAX:
                    pair = frozenset([(pw, qw), (rp, rq)])
                    if pair not in seen:
                        seen.add(pair)
                        sums.append(MULT_NONTRIVIAL[(pw, qw)] + MULT_NONTRIVIAL[(rp, rq)])
        if len(sums) >= 2:
            assert len(set(sums)) > 1, (
                f"Non-trivial: UNEXPECTED constant sum {sums[0]} on support "
                f"under (p,q)→({a}-p,{b}-q)")


class TestComparisonWithSU2:
    """SU(2) duality works because dim V_j = 2j+1 is degree 1."""

    def test_su2_dim_sum_constant(self):
        """In SU(2): dim V_j + dim V_{j*-j} = 2j*+2 = |H|, always constant."""
        for j_star in [2, 5, 11, 29]:
            H = 2 * (j_star + 1)
            for j in range(j_star + 1):
                d1 = 2 * j + 1
                d2 = 2 * (j_star - j) + 1
                assert d1 + d2 == H, (
                    f"j*={j_star}, j={j}: {d1}+{d2}={d1+d2} ≠ {H}")

    def test_su3_dim_sum_not_constant(self):
        """In SU(3): dim V_{(p,q)} + dim V_{(P-p,Q-q)} is NOT constant."""
        P, Q = 6, 6
        sums = set()
        for pw in range(P + 1):
            for qw in range(Q + 1):
                d1 = (pw + 1) * (qw + 1) * (pw + qw + 2) // 2
                d2 = (P - pw + 1) * (Q - qw + 1) * (P - pw + Q - qw + 2) // 2
                sums.add(d1 + d2)
        assert len(sums) > 1
        assert max(sums) > 2 * min(sums), (
            f"Dim sum range [{min(sums)}, {max(sums)}] — wide variation confirms no pairing")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
