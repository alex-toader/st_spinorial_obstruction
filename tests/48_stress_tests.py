"""
Script 48: Comprehensive stress tests for the spectral duality paper.

Items covered: 48 (Q8 minimal), 49 (brute force), 53 (cl(-1) to n=200),
55 (duality defect), 56 (vector bundle binary), 60 (Hasse monotonicity).

Run: /usr/bin/python3 -m pytest tests/48_stress_tests.py -v -s
Date: Feb 2026

RAW OUTPUT (9 passed, 15.37s) — v3 after reviewer audit:

  v3 fixes: (a) find_idx guards (assert ≥ 0 everywhere), (b) Hasse uses BOTH
  direct Q8-subset check AND commutator test (no circularity), (c) VectorBundle
  computes per-irrep multiplicities m(H₁,j), m(H₂,j), m(L,j) individually.

  TestQ8Minimal: all 4 chars × 2 levels, duality perfect ✓

  TestBruteForce (ALL scalar chars):
    Q8..2D_12: 4 chars each, duality ✓, m≥1 above ✓
    2T: 3, 2O: 2, 2I: 1 chars. All pass. ✓

  TestCommutatorLength:
    2D_n structural: [a^{n/2}, b] = -1 for all even n in [2, 200] ✓
    Polyhedral: 2T, 2O, 2I all have pure imaginary orthogonal witness ✓

  TestDualityDefect: δ(j) = (-1)^j for n=3..51 (odd). EXACT. ✓

  TestVectorBundleBinary (PER-IRREP multiplicities on 2O):
    Irreps extracted: H₁=V_{1/2}, L=V_{3/2}, H₂ from V_{5/2}. All irreducible. ✓
    Duality m(ρ,j)+m(ρ,j*-j) = dim(ρ) verified for all half-integer j ≤ j*. ✓
    Binary question ANSWERED: intermediate values EXIST.
      H₁(dim 2): m=1 at j=0.5,3.5,4.5,5.5,6.5,7.5,10.5 (not binary)
      H₂(dim 2): m=1 at j=2.5,3.5,5.5,7.5,8.5 (not binary)
      L(dim 4):  m=1,2,3 at various j (not binary)
    Conclusion: binary {0,1} property is SPECIFIC TO SCALAR (rank-1) sectors.
    Higher-rank bundles have m(ρ,j)+m(ρ,j*-j)=dim(ρ) but m can be intermediate.

  TestHasseMonotone (TWO independent Q8 checks):
    Direct subset check: {±1,±i,±j,±k} ⊂ G tested element-by-element
    Commutator check: ∃[a,b]=-1 tested independently
    Both agree on all 9 subgroups. ✓
    2O, 2T, 2D4, Q8: OBSTRUCTED. 2D3, Z8, Z6, Z4, Z2: free. ✓
    Monotonicity: once free, stays free ✓

KEY RESULTS:
  1. cl(-1) = 1 structurally: witness is always [a^{n/2}, b] = [i, j] = -1
  2. Duality holds for ALL scalar chars — κ chars × j* levels
  3. δ(j) = (-1)^j for non-obstructed groups (proved analytically)
  4. Binary m∈{0,1} is rank-1 specific; bundles have intermediate multiplicities
  5. Q8 subset and commutator tests agree (Theorem 2.3 independently validated)
"""

import numpy as np
import pytest

# ============================================================
# Core utilities (from script 41, minimal)
# ============================================================

def qmul(a, b):
    return np.array([
        a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
        a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0],
    ])

def qinv(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

def commutator(a, b):
    return qmul(qmul(a, b), qmul(qinv(a), qinv(b)))

def is_minus_one(q, tol=1e-10):
    return abs(q[0] + 1) < tol and np.linalg.norm(q[1:]) < tol

def generate_group(generators, max_size=500):
    elements = [np.array([1, 0, 0, 0], dtype=float)]
    def is_in(q, lst):
        return any(np.linalg.norm(q - e) < 1e-10 for e in lst)
    for g in generators:
        for x in [g, qinv(g)]:
            if not is_in(x, elements):
                elements.append(x)
    changed = True
    while changed:
        changed = False
        new = []
        for a in elements:
            for b in elements:
                p = qmul(a, b)
                if not is_in(p, elements) and not is_in(p, new):
                    new.append(p)
                    if len(elements) + len(new) > max_size:
                        raise ValueError(f"Exceeds {max_size}")
        if new:
            elements.extend(new)
            changed = True
    return elements

def chi_j(g, j):
    """SU(2) character at spin j."""
    a0 = g[0]
    # Half-angle: cos(alpha) = a0, but handle ±1 specially
    if abs(abs(a0) - 1) < 1e-12:
        # g = ±1
        return ((-1)**(2*j)) * (2*j + 1) if a0 < 0 else (2*j + 1)
    alpha = np.arccos(np.clip(a0, -1, 1))
    s = np.sin(alpha)
    if abs(s) < 1e-14:
        return (2*j + 1)
    return np.sin((2*j + 1) * alpha) / s

def multiplicity(G, chi_vals, j):
    """m(chi, j) = (1/|G|) sum_g conj(chi(g)) * chi_j(g)."""
    n = len(G)
    total = 0.0
    for i, g in enumerate(G):
        total += np.conj(chi_vals[i]) * chi_j(g, j)
    return total / n

def trivial_char(G):
    return np.ones(len(G))

# ============================================================
# Group builders
# ============================================================

def make_2Dn(n):
    a = np.array([np.cos(np.pi/n), np.sin(np.pi/n), 0, 0])
    b = np.array([0, 0, 1, 0])
    return generate_group([a, b])

def make_2T():
    return generate_group([np.array([0, 1, 0, 0]), np.array([0.5, 0.5, 0.5, 0.5])])

def make_2O():
    return generate_group([np.array([0.5, 0.5, 0.5, 0.5]),
                           np.array([1/np.sqrt(2), 1/np.sqrt(2), 0, 0])])

def make_2I():
    phi = (1 + np.sqrt(5)) / 2
    return generate_group([np.array([0.5, 0.5, 0.5, 0.5]),
                           np.array([phi/2, 1/(2*phi), 0.5, 0])])


# ============================================================
# 48. Q8 minimal — complete verification
# ============================================================

class TestQ8Minimal:
    """j* = 1. Two levels: j=0, j=1. Everything must work."""

    def test_q8_complete(self):
        G = make_2Dn(2)
        assert len(G) == 8
        jstar = len(G) // 4 - 1
        assert jstar == 1

        A0 = trivial_char(G)

        # m(A0, 0) = 1, m(A0, 1) = 0
        m0 = round(multiplicity(G, A0, 0).real)
        m1 = round(multiplicity(G, A0, 1).real)
        print(f"Q8: m(A0, 0) = {m0}, m(A0, 1) = {m1}")
        assert m0 == 1
        assert m1 == 0
        assert m0 + m1 == 1, "Duality fails at Q8!"

        # Non-trivial scalar chars: Q8^ab = V4, so 4 chars total
        # Build all 1D chars by finding abelianization
        # Q8/{±1} = V4. Chars of V4 lifted to Q8.
        # Identify which element is which by half-angle
        pure_imag = [g for g in G if abs(g[0]) < 1e-10]
        assert len(pure_imag) == 6  # ±i, ±j, ±k

        # Group by axis
        axes = []
        for g in pure_imag:
            if g[1] > 0.5 or g[2] > 0.5 or g[3] > 0.5:
                axes.append(g)
        # Should be i, j, k (positive versions)

        # Build chi_i: +1 on ±1, ±i; -1 on ±j, ±k
        def make_chi(pos_axis_idx):
            """char that is +1 on elements along axis pos_axis_idx, -1 on others."""
            vals = np.ones(len(G))
            for idx, g in enumerate(G):
                if abs(g[0]) < 1e-10:  # pure imaginary
                    # Check which axis: g[1], g[2], g[3]
                    ax = np.argmax(np.abs(g[1:])) + 1
                    if ax != pos_axis_idx:
                        vals[idx] = -1.0
            return vals

        nontrivial_count = 0
        for ax in [1, 2, 3]:  # i, j, k axes
            chi = make_chi(ax)
            m_0 = round(multiplicity(G, chi, 0).real)
            m_1 = round(multiplicity(G, chi, 1).real)
            print(f"  chi_{ax}: m(0) = {m_0}, m(1) = {m_1}, sum = {m_0 + m_1}")
            assert m_0 == 0, f"chi_{ax} at j=0 should be 0"
            assert m_1 == 1, f"chi_{ax} at j=1 should be 1"
            assert m_0 + m_1 == 1, "Duality fails!"
            nontrivial_count += 1

        assert nontrivial_count == 3
        print("Q8 COMPLETE: all 4 chars × 2 levels, duality perfect ✓")


# ============================================================
# 49. Brute force: m(A0, j) for j = 0..2j*, all obstructed groups
# ============================================================

def scalar_characters(G):
    """Build ALL 1D characters of G by finding the abelianization.

    Method: compute all commutators, find normal closure = [G,G].
    Then G^ab = G/[G,G]. Characters of G^ab lift to 1D chars of G.
    We use the concrete approach: find [G,G] elements, then for each
    coset representative, assign a root of unity."""
    n = len(G)

    def find_idx(q):
        for i, e in enumerate(G):
            if np.linalg.norm(q - e) < 1e-10:
                return i
        return -1

    # Compute [G,G] as set of indices
    comm_elements = set()
    id_idx = find_idx(np.array([1, 0, 0, 0]))
    assert id_idx >= 0, "Identity not found in group!"
    comm_elements.add(id_idx)
    for a in G:
        for b in G:
            c = commutator(a, b)
            idx = find_idx(c)
            assert idx >= 0, f"Commutator not found in group! a={a}, b={b}, [a,b]={c}"
            comm_elements.add(idx)

    # Close under multiplication
    changed = True
    while changed:
        changed = False
        new = set()
        for i in comm_elements:
            for j_idx in comm_elements:
                p = qmul(G[i], G[j_idx])
                idx = find_idx(p)
                assert idx >= 0, f"Product not in group! G[{i}]·G[{j_idx}]={p}"
                if idx not in comm_elements and idx not in new:
                    new.add(idx)
        if new:
            comm_elements |= new
            changed = True

    comm_order = len(comm_elements)
    kappa = n // comm_order  # |G^ab|

    # Build coset labels: for each g, find which coset of [G,G] it's in
    coset_label = [-1] * n
    coset_reps = []
    for i in range(n):
        if coset_label[i] >= 0:
            continue
        label = len(coset_reps)
        coset_reps.append(i)
        # Mark all elements in this coset
        for c_idx in comm_elements:
            p = qmul(G[i], G[c_idx])
            j_idx = find_idx(p)
            assert j_idx >= 0, f"Coset product not in group! G[{i}]·G[{c_idx}]={p}"
            coset_label[j_idx] = label

    assert len(coset_reps) == kappa
    assert all(c >= 0 for c in coset_label)

    # Characters: each char sends coset k to ω^(m·k) for m = 0..κ-1
    # where ω = exp(2πi/κ). But G^ab may not be cyclic.
    # Safe approach: enumerate all homomorphisms G^ab → U(1).
    # For small κ (≤ 8 in our cases), brute force over κ-th roots of unity.

    # Find multiplication table of cosets
    coset_mult = [[0]*kappa for _ in range(kappa)]
    for a_label in range(kappa):
        for b_label in range(kappa):
            p = qmul(G[coset_reps[a_label]], G[coset_reps[b_label]])
            idx = find_idx(p)
            assert idx >= 0, f"Coset mult product not in group!"
            coset_mult[a_label][b_label] = coset_label[idx]

    # Find all homomorphisms: assign value to each coset, check multiplication
    chars = []
    # Try all κ-th roots of unity for each coset
    roots = [np.exp(2j * np.pi * k / kappa) for k in range(kappa)]

    def find_chars(partial, next_label):
        if next_label == kappa:
            # Verify homomorphism
            for a in range(kappa):
                for b in range(kappa):
                    ab = coset_mult[a][b]
                    if abs(partial[a] * partial[b] - partial[ab]) > 1e-10:
                        return
            chars.append(list(partial))
            return
        for r in roots:
            partial[next_label] = r
            find_chars(partial, next_label + 1)

    partial = [0] * kappa
    partial[0] = 1.0 + 0j  # identity coset → 1
    find_chars(partial, 1)

    assert len(chars) == kappa, f"Expected {kappa} chars, found {len(chars)}"

    # Convert to per-element character values
    result = []
    for char_coset in chars:
        vals = np.array([char_coset[coset_label[i]] for i in range(n)])
        result.append(vals)

    return result


class TestBruteForce:
    """Verify duality for ALL scalar characters, not just A₀."""

    def _test_group(self, name, G):
        n = len(G)
        jstar = n // 4 - 1
        chars = scalar_characters(G)
        kappa = len(chars)

        for ci, chi in enumerate(chars):
            # Sub-threshold: duality
            for j in range(jstar + 1):
                m_j = round(multiplicity(G, chi, j).real)
                m_dual = round(multiplicity(G, chi, jstar - j).real)
                assert m_j + m_dual == 1, \
                    f"{name} χ_{ci}: duality fails at j={j}, m={m_j}, m'={m_dual}"
                assert m_j in (0, 1), f"{name} χ_{ci}: non-binary m={m_j} at j={j}"

            # Above threshold first period: m >= 1
            P = n // 4
            for j in range(jstar + 1, jstar + P + 1):
                m_j = round(multiplicity(G, chi, j).real)
                assert m_j >= 1, \
                    f"{name} χ_{ci}: m={m_j} < 1 at j={j} (above threshold!)"

        print(f"{name}: duality ✓ (j≤{jstar}), m≥1 above ✓, all {kappa} chars ✓")

    def test_Q8(self):
        self._test_group("Q8", make_2Dn(2))

    def test_2D_even(self):
        for n in [4, 6, 8, 10, 12]:
            self._test_group(f"2D_{n}", make_2Dn(n))

    def test_polyhedral(self):
        self._test_group("2T", make_2T())
        self._test_group("2O", make_2O())
        self._test_group("2I", make_2I())


# ============================================================
# 53. Commutator length cl(-1) for 2D_n, n even, up to n=200
# ============================================================

class TestCommutatorLength:
    """cl(-1) = 1 for all obstructed groups."""

    def test_cl_minus_one_2Dn(self):
        """For 2D_n even: [a^{n/2}, b] = -1 (structural).
        We verify numerically for n up to 200."""
        fails = []
        for n in range(2, 202, 2):
            angle = np.pi / n
            a = np.array([np.cos(angle), np.sin(angle), 0, 0])
            b = np.array([0, 0, 1, 0])

            # Structural witness: a^{n/2} has half-angle (n/2)(π/n) = π/2.
            # So a^{n/2} is pure imaginary. [a^{n/2}, b] = -1.
            ap = np.array([1, 0, 0, 0], dtype=float)
            for _ in range(n // 2):
                ap = qmul(ap, a)
            c = commutator(ap, b)
            if not is_minus_one(c):
                fails.append(n)
            if n <= 20 or n % 50 == 0:
                print(f"2D_{n}: [a^{n//2}, b] = -1: {is_minus_one(c)} ✓")

        assert len(fails) == 0, f"cl(-1) > 1 for n = {fails}"
        print(f"cl(-1) = 1 for ALL even n in [2, 200] (structural: [a^n, b] = -1) ✓")

    def test_cl_minus_one_polyhedral(self):
        """2T, 2O, 2I: search all pairs for [a,b] = -1."""
        for name, G in [("2T", make_2T()), ("2O", make_2O()), ("2I", make_2I())]:
            found = False
            witness = None
            for a in G:
                for b in G:
                    if is_minus_one(commutator(a, b)):
                        found = True
                        witness = (a, b)
                        break
                if found:
                    break
            assert found, f"{name}: cl(-1) > 1!"
            # Verify witness is pure imaginary orthogonal pair
            a, b = witness
            a_pure = abs(a[0]) < 1e-10
            b_pure = abs(b[0]) < 1e-10
            ortho = abs(np.dot(a[1:], b[1:])) < 1e-10 if a_pure and b_pure else False
            print(f"{name}: [a,b]=-1 ✓ (pure_imag={a_pure and b_pure}, orthogonal={ortho})")


# ============================================================
# 55. Duality defect for non-obstructed groups
# ============================================================

class TestDualityDefect:
    """When duality fails (odd n), what is δ = m+m'-1?"""

    def test_defect_2Dn_odd(self):
        """Duality defect for non-obstructed 2D_n (odd n).

        ANALYTIC PROOF that δ(j) = (-1)^j:

        For 2D_n (odd n), |G| = 4n, j* = n-1, H = D_n, |H| = 2n.
        Non-central elements have SO(3) orders d | n (rotation type) and d = 2 (reflection type).

        For rotation elements (half-angle α = kπ/n, k = 1..n-1):
          K = |H|/d = 2n/d.  Since n is odd and d | n, d is odd, so K = 2n/d is even.
          Mode: f(α, j*-j) = (-1)^{K+1} f(α, j) = -f(α, j).  Antisymmetric. ✓

        For reflection elements (half-angle α = π/2, d = 2):
          K = |H|/2 = n.  Since n is odd, K is odd.
          Mode: f(π/2, j*-j) = (-1)^{n+1} f(π/2, j) = +f(π/2, j).  SYMMETRIC.

        So one mode (d=2) is symmetric, all others antisymmetric.
        The defect comes entirely from the d=2 mode:
          S(j) + S(j*-j) = 2 · A_{π/2} · f(π/2, j)
        where A_{π/2} = number of elements at half-angle π/2.

        Now f(π/2, j) = sin((2j+1)π/2)/sin(π/2) = sin((2j+1)π/2) = (-1)^j.

        And the central contribution: 2(2j+1)/|G| + 2(2(j*-j)+1)/|G| = 2(2j*+2)/|G| = 1.

        So: m(j) + m(j*-j) = 1 + (2A_{π/2}/|G|) · 2 · (-1)^j.

        For 2D_n: A_{π/2} = n (there are n elements of order 4 in SU(2)).
        |G| = 4n, so 2A_{π/2}/|G| · 2 = 2·n/(4n) · 2 = 1.

        Therefore: m(j) + m(j*-j) = 1 + (-1)^j = {2 if j even, 0 if j odd}.
        Equivalently: δ(j) = (-1)^j.  QED.
        """
        print("\nDuality defect δ(A₀, j) = m(j) + m(j*-j) - 1 for odd n:")
        print(f"{'n':>4} | {'j*':>3} | δ values | max|δ| | pattern")
        print("-" * 65)

        for n in [3, 5, 7, 9, 11, 13, 15, 21, 31, 51]:
            G = make_2Dn(n)
            jstar = len(G) // 4 - 1
            A0 = trivial_char(G)

            defects = []
            for j in range(jstar + 1):
                m_j = round(multiplicity(G, A0, j).real)
                m_dual = round(multiplicity(G, A0, jstar - j).real)
                delta = m_j + m_dual - 1
                defects.append(delta)

            max_abs = max(abs(d) for d in defects)
            # Compact representation
            d_str = ''.join(f"{d:+d}" if d != 0 else " 0" for d in defects[:15])
            if len(defects) > 15:
                d_str += "..."

            print(f"{n:>4} | {jstar:>3} | {d_str} | {max_abs:>5} | ", end="")

            # Check for periodicity
            nonzero = [j for j, d in enumerate(defects) if d != 0]
            if len(nonzero) == 0:
                print("all zero (duality holds!)")
            else:
                print(f"nonzero at j={nonzero[:8]}")

            # Exact pattern: δ(j) = (-1)^j for odd n
            expected = [(-1)**j for j in range(jstar + 1)]
            assert defects == expected, \
                f"n={n}: defect ≠ (-1)^j! got {defects}"


# ============================================================
# 56. Vector bundle: binary sub-structure?
# ============================================================

class TestVectorBundleBinary:
    """Per-irrep multiplicities for spinorial irreps of 2O.

    2O has 3 spinorial irreps: H₁(dim 2), H₂(dim 2), L(dim 4).
    Question: is m(ρ, j) ∈ {0, dim(ρ)} always? Or can intermediate values occur?

    Method: extract irrep characters from small-j representations.
      χ_{H₁} = χ_{1/2} (fundamental rep of SU(2) restricted to 2O)
      χ_L = χ_{3/2} (spin-3/2 restricted to 2O is irreducible, dim 4)
      χ_{H₂} = extracted from first j where it appears, by subtraction.
    Then compute m(ρ, j) = ⟨χ_ρ, χ_j⟩_G for each ρ individually.
    """

    def test_2O_per_irrep_duality(self):
        """Compute m(H₁, j), m(H₂, j), m(L, j) individually for all half-integer j ≤ j*."""
        G = make_2O()
        n = len(G)
        assert n == 48
        jstar = 11

        # Step 1: Extract spinorial irrep characters from small j.
        # H₁ = V_{1/2}|_{2O} (dim 2, the fundamental/defining representation)
        chi_H1 = np.array([chi_j(g, 0.5) for g in G])

        # L = V_{3/2}|_{2O} (dim 4, spin-3/2 is irreducible on 2O)
        chi_L = np.array([chi_j(g, 1.5) for g in G])

        # Verify these are irreducible: ⟨χ, χ⟩ = 1
        norm_H1 = np.sum(np.abs(chi_H1)**2).real / n
        norm_L = np.sum(np.abs(chi_L)**2).real / n
        assert abs(norm_H1 - 1) < 1e-10, f"H₁ not irreducible: ⟨χ,χ⟩ = {norm_H1}"
        assert abs(norm_L - 1) < 1e-10, f"L not irreducible: ⟨χ,χ⟩ = {norm_L}"

        # Verify H₁ ≠ L: ⟨H₁, L⟩ = 0
        cross_HL = np.sum(np.conj(chi_H1) * chi_L).real / n
        assert abs(cross_HL) < 1e-10, f"H₁ and L not orthogonal: {cross_HL}"

        # Step 2: Find H₂ from V_{5/2}|_{2O} (dim 6 = 2 + 4 or 2 + 2 + 2).
        chi_5half = np.array([chi_j(g, 2.5) for g in G])
        m_H1_5half = round((np.sum(np.conj(chi_H1) * chi_5half) / n).real)
        m_L_5half = round((np.sum(np.conj(chi_L) * chi_5half) / n).real)
        # Remaining must be H₂
        m_H2_5half_dim = 6 - 2 * m_H1_5half - 4 * m_L_5half
        assert m_H2_5half_dim > 0, f"H₂ not found at j=5/2: remaining dim = {m_H2_5half_dim}"
        m_H2_5half = m_H2_5half_dim // 2
        assert m_H2_5half * 2 == m_H2_5half_dim

        # Extract H₂ character
        chi_H2 = (chi_5half - m_H1_5half * chi_H1 - m_L_5half * chi_L) / m_H2_5half

        # Verify H₂ is irreducible
        norm_H2 = np.sum(np.abs(chi_H2)**2).real / n
        assert abs(norm_H2 - 1) < 1e-10, f"H₂ not irreducible: ⟨χ,χ⟩ = {norm_H2}"

        # Verify orthogonality
        cross_H1H2 = np.sum(np.conj(chi_H1) * chi_H2).real / n
        cross_H2L = np.sum(np.conj(chi_H2) * chi_L).real / n
        assert abs(cross_H1H2) < 1e-10, f"H₁ ⊥ H₂ fails: {cross_H1H2}"
        assert abs(cross_H2L) < 1e-10, f"H₂ ⊥ L fails: {cross_H2L}"

        print(f"\n2O spinorial irreps extracted:")
        print(f"  H₁ = V_{{1/2}}|_{{2O}} (dim 2, ⟨χ,χ⟩={norm_H1:.6f})")
        print(f"  H₂ from V_{{5/2}} (dim 2, ⟨χ,χ⟩={norm_H2:.6f})")
        print(f"  L  = V_{{3/2}}|_{{2O}} (dim 4, ⟨χ,χ⟩={norm_L:.6f})")

        # Step 3: Compute per-irrep multiplicities for all half-integer j ≤ j*
        irreps = [("H₁", chi_H1, 2), ("H₂", chi_H2, 2), ("L", chi_L, 4)]

        print(f"\n{'j':>5} | {'H₁':>3} {'H₂':>3} {'L':>3} | {'dim':>3} | duality pair")
        print("-" * 55)

        has_intermediate = False
        for j_half in range(1, 2*jstar + 2, 2):  # half-integer j from 1/2 to j*+1/2
            j = j_half / 2.0
            jd = jstar - j
            if jd < 0:
                break

            chi_j_vals = np.array([chi_j(g, j) for g in G])
            chi_jd_vals = np.array([chi_j(g, jd) for g in G])

            mults = {}
            mults_dual = {}
            for rho_name, chi_rho, dim_rho in irreps:
                m_j = round((np.sum(np.conj(chi_rho) * chi_j_vals) / n).real)
                m_jd = round((np.sum(np.conj(chi_rho) * chi_jd_vals) / n).real)
                mults[rho_name] = m_j
                mults_dual[rho_name] = m_jd

                # Duality check: m + m' = dim(ρ)
                assert m_j + m_jd == dim_rho, \
                    f"j={j}, {rho_name}: m({j}) + m({jd}) = {m_j} + {m_jd} ≠ {dim_rho}"

                # Check binary question: is m ∈ {0, dim(ρ)}?
                if m_j not in (0, dim_rho):
                    has_intermediate = True

            dim_check = 2*mults["H₁"] + 2*mults["H₂"] + 4*mults["L"]
            assert dim_check == int(2*j + 1), \
                f"j={j}: dim mismatch {dim_check} ≠ {int(2*j+1)}"

            pair_str = f"({mults['H₁']}+{mults_dual['H₁']}) ({mults['H₂']}+{mults_dual['H₂']}) ({mults['L']}+{mults_dual['L']})"
            marker = ""
            for rho_name, _, dim_rho in irreps:
                if mults[rho_name] not in (0, dim_rho):
                    marker = " ← INTERMEDIATE"
                    break
            print(f"{j:>5.1f} | {mults['H₁']:>3} {mults['H₂']:>3} {mults['L']:>3} | {int(2*j+1):>3} | {pair_str}{marker}")

        if has_intermediate:
            print("\nRESULT: intermediate multiplicities EXIST — NOT fully binary at higher rank.")
        else:
            print("\nRESULT: m(ρ,j) ∈ {0, dim(ρ)} always — fully binary even for bundles!")


# ============================================================
# 60. Hasse monotonicity — obstruction can't return
# ============================================================

class TestHasseMonotone:
    """Compute Q8 containment dynamically for all subgroups of 2O."""

    def _contains_Q8_commutator(self, G):
        """Check if ∃ a,b ∈ G with [a,b] = -1 (Theorem 2.3 criterion)."""
        for a in G:
            for b in G:
                if is_minus_one(commutator(a, b)):
                    return True
        return False

    def _contains_Q8_subset(self, G):
        """Direct subset check: generate Q8 = {±1, ±i, ±j, ±k}, verify all 8 ∈ G."""
        q8 = [
            np.array([1, 0, 0, 0]),
            np.array([-1, 0, 0, 0]),
            np.array([0, 1, 0, 0]),
            np.array([0, -1, 0, 0]),
            np.array([0, 0, 1, 0]),
            np.array([0, 0, -1, 0]),
            np.array([0, 0, 0, 1]),
            np.array([0, 0, 0, -1]),
        ]
        def is_in(q, lst):
            return any(np.linalg.norm(q - e) < 1e-10 for e in lst)
        return all(is_in(q, G) for q in q8)

    def _make_cyclic(self, n):
        """Z_{2n} ⊂ SU(2), order 2n."""
        angle = np.pi / n
        a = np.array([np.cos(angle), np.sin(angle), 0, 0])
        return generate_group([a])

    def test_O_subgroups_computed(self):
        """Dynamically generate each binary lift 2H' and check Q8 containment."""
        groups = [
            ("2O",  make_2O,         48),
            ("2T",  make_2T,         24),
            ("2D4", lambda: make_2Dn(4), 16),
            ("2D3", lambda: make_2Dn(3), 12),
            ("Q8",  lambda: make_2Dn(2),  8),
            ("Z8",  lambda: self._make_cyclic(4), 8),
            ("Z6",  lambda: self._make_cyclic(3), 6),
            ("Z4",  lambda: self._make_cyclic(2), 4),
            ("Z2",  lambda: self._make_cyclic(1), 2),
        ]

        print("\n2O subgroups — Q8 containment (two independent checks):")
        results = {}
        for name, builder, expected_order in groups:
            G = builder()
            assert len(G) == expected_order, f"{name}: |G|={len(G)} ≠ {expected_order}"
            has_q8_comm = self._contains_Q8_commutator(G)
            has_q8_sub = self._contains_Q8_subset(G)
            assert has_q8_comm == has_q8_sub, \
                f"{name}: commutator test ({has_q8_comm}) ≠ subset test ({has_q8_sub})!"
            results[name] = has_q8_sub  # use direct subset as primary
            status = "OBSTRUCTED" if has_q8_sub else "free"
            print(f"  {name:>4} (|G|={expected_order:>3}): subset={has_q8_sub}, comm={has_q8_comm}  [{status}]")

        # Verify expected values
        assert results["2O"] == True
        assert results["2T"] == True
        assert results["2D4"] == True
        assert results["Q8"] == True
        assert results["2D3"] == False
        assert results["Z8"] == False
        assert results["Z6"] == False
        assert results["Z4"] == False
        assert results["Z2"] == False

        # Monotonicity: subgroups of free groups are free
        # 2D3 is free. Its subgroups (Z6, Z3, Z2) must also be free. ✓
        # Z8 is free. Its subgroups (Z4, Z2) must also be free. ✓
        print("\nMonotonicity: all subgroups of free groups are free ✓")
        print("Proof: Q8 ⊄ H' ∧ H'' ⊂ H' ⟹ Q8 ⊄ H'' (subgroup transitivity)")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
