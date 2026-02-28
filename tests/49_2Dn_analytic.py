"""
Script 49: Analytic proof of divisibility condition (B) for the binary dihedral family.

Goal: prove that for 2D_n (order 4n), the divisibility condition holds iff n is even.

The binary dihedral group 2D_n (also called dicyclic Dic_n or generalized quaternion
for n = 2^k) has order 4n and sits inside SU(2). Its SO(3) image is the dihedral
group D_n of order 2n.

Conjugacy class structure of D_n (n ≥ 3):
  - {e}: identity, order 1
  - Rotations: a^k for k = 1, ..., n-1 grouped by order
    * For each k: a^k and a^{n-k} are conjugate (since bab⁻¹ = a⁻¹ in D_n)
    * Class {a^k, a^{n-k}} has SO(3) order d = n/gcd(k, n)
    * Number of rotation classes: ⌊(n-1)/2⌋
  - Reflections: elements of the form a^k b
    * n even: two reflection classes of size n/2 each
      {b, a²b, a⁴b, ...} and {ab, a³b, ...}
    * n odd: one reflection class of size n
      {b, ab, a²b, ..., a^{n-1}b}
    * Reflection elements have SO(3) order 2

Divisibility condition (B): |D_n|/ord(h) is even for all h in D_n, h != e.
  |D_n| = 2n, so need: 2n/d is even for all non-identity h of order d.

For reflections: d = 2, so 2n/2 = n. Need n even. ✓/✗

For rotations: element a^k has order d = n/gcd(k,n).
  K = 2n/d = 2n·gcd(k,n)/n = 2·gcd(k,n). ALWAYS EVEN! ✓

So the ONLY obstruction to divisibility comes from the reflection classes:
  K_refl = 2n/2 = n.
  n even → K_refl even → duality holds.
  n odd → K_refl odd → duality breaks.

This is the complete analytic proof: (B) holds for 2D_n ⟺ n even.

Run: /usr/bin/python3 -m pytest tests/49_2Dn_analytic.py -v -s
Date: Feb 2026

RAW OUTPUT (12 passed, 0.13s) — v3 after second reviewer pass:

  v3 fixes: (a) rotation index k added to dihedral_classes tuple (6-element),
  eliminates fragile string parsing of desc field,
  (b) polyhedral element count verification (assert sum = |H|-1),
  (c) n=2 edge case documented in docstring.

  v2 fixes: (a) rotation test reduced (K=2·gcd is trivially even),
  (b) polyhedral test rewritten with clear element counts (no string hacks),
  (c) tolerance guard on _multiplicity_via_classes catches non-integer m,
  (d) half-angle bug FOUND AND FIXED: was using α=π/d instead of α=kπ/n.
      The old code gave m=0.354 for n=8,j=1 — round() silently produced 0.

  TestDivisibilityAnalytic:
    test_rotation_K_formula — K = 2·gcd(k,n) formula verified on 10 values of n. ✓
    test_reflection_class_parity_determines_duality — K_refl = n, n=2..100. ✓
    test_full_divisibility_condition — (B) holds ⟺ n even, n=2..100. ✓
    test_class_count_consistency — class sizes sum to |D_n|-1 for all n. ✓
    test_print_structure — D_2..D_20 with actual half-angles. ✓

  TestModeAntisymmetry (now uses actual half-angles and k'K sign):
    test_antisymmetry_obstructed — k'K always even for n even → all modes antisymmetric. ✓
    test_antisymmetry_non_obstructed — reflections: k'K=n odd → mode symmetric. ✓

  TestK_Formula:
    test_K_rotation_formula — K_rot = 2·gcd(k,n). Proof statement. ✓
    test_polyhedral_K_values — T,O,I: all K values even (with element counts). ✓

  TestDualityNumerical (class-based, independent of script 48):
    test_duality_via_class_formula — m(j)+m(j*-j)=1, n=2..50. Exact integers. ✓
    test_defect_via_class_formula — δ(j)=(-1)^j, n=3..11. ✓

  TestSummary — Complete proof statement. ✓

KEY RESULT:
  For 2D_n: K_rotation = 2·gcd(k,n) (ALWAYS even), K_reflection = n.
  Therefore (B) ⟺ n even. Eliminates infinite family from case-by-case.
  Only 2T, 2O, 2I remain as 3 finite checks (all K values even).

  Refined antisymmetry: the sign at half-angle α = k'π/d is (-1)^{k'K+1}.
  For rotations: k'K = k'·2·gcd(k,n) is always even (factor of 2). ✓
  For reflections: k'=1, K=n. Even iff n even. Same conclusion.
"""

import numpy as np
import pytest


# ============================================================
# Analytic computation of K = |H|/d for all conjugacy classes
# ============================================================

def dihedral_classes(n):
    """Return all non-identity conjugacy classes of D_n with their SO(3) orders.

    Returns list of (desc, size, d, K, alpha, k) where:
      desc: class description string
      size: number of elements in the class
      d: SO(3) order of elements in the class
      K: |D_n|/d = 2n/d
      alpha: actual SU(2) half-angle (kπ/n for rotations, π/2 for reflections)
      k: rotation index (k for a^k class, None for reflections)

    Note: D_2 = V₄ is abelian, so "rotation" vs "reflection" is conventional.
    The generic D_n analysis assumes n ≥ 3 for the class structure, but the
    formulas remain correct for n = 2 (Q₈ case).
    """
    classes = []

    # Rotation classes: a^k and a^{n-k} are conjugate for k=1..n-1
    seen = set()
    for k in range(1, n):
        if k in seen:
            continue
        partner = n - k
        d = n // np.gcd(k, n)  # SO(3) order of a^k
        K = 2 * n // d
        alpha = k * np.pi / n  # actual half-angle of a^k
        if k == partner:
            size = 1  # self-conjugate (only when n even, k = n/2)
        else:
            size = 2
            seen.add(partner)
        classes.append((f"a^{k}", size, d, K, alpha, k))

    # Reflection classes: all reflections have half-angle π/2 (order 4 in SU(2))
    if n % 2 == 0:
        classes.append(("b-type even", n // 2, 2, n, np.pi / 2, None))
        classes.append(("b-type odd", n // 2, 2, n, np.pi / 2, None))
    else:
        classes.append(("b-type all", n, 2, n, np.pi / 2, None))

    return classes


class TestDivisibilityAnalytic:
    """Prove (B) for 2D_n analytically via conjugacy class structure."""

    def test_rotation_K_formula(self):
        """K = 2·gcd(k,n) for rotation class a^k.

        Proof: d = ord(a^k) = n/gcd(k,n). K = 2n/d = 2·gcd(k,n).
        This is 2 times an integer, hence always even. QED.

        We verify the formula K = 2·gcd(k,n) on a few values (the evenness
        is trivial and doesn't need 10000-iteration verification).
        """
        for n in [2, 3, 5, 6, 10, 12, 15, 20, 60, 100]:
            for k in range(1, n):
                d = n // np.gcd(k, n)
                K = 2 * n // d
                assert K == 2 * np.gcd(k, n), f"n={n}, k={k}: K formula mismatch"

    def test_reflection_class_parity_determines_duality(self):
        """K_refl = n. Even iff n even."""
        for n in range(2, 101):
            K_refl = n
            if n % 2 == 0:
                assert K_refl % 2 == 0, f"n={n}: K_refl should be even"
            else:
                assert K_refl % 2 == 1, f"n={n}: K_refl should be odd"

    def test_full_divisibility_condition(self):
        """(B) holds ⟺ ALL non-identity classes have K even ⟺ n even."""
        for n in range(2, 101):
            classes = dihedral_classes(n)
            all_K_even = all(K % 2 == 0 for (_, _, _, K, _, _) in classes)

            if n % 2 == 0:
                assert all_K_even, f"n={n} (even): (B) should hold but some K is odd!"
            else:
                assert not all_K_even, f"n={n} (odd): (B) should fail but all K even!"

            # Check that it's ONLY the reflections that can fail
            rotation_K_even = all(K % 2 == 0 for (_, _, _, K, _, k_val) in classes
                                  if k_val is not None)
            assert rotation_K_even, f"n={n}: some rotation class has odd K!"

    def test_class_count_consistency(self):
        """Verify total element count from classes equals |D_n| - 1."""
        for n in range(2, 101):
            classes = dihedral_classes(n)
            total_non_id = sum(size for (_, size, _, _, _, _) in classes)
            assert total_non_id == 2 * n - 1, \
                f"n={n}: class sizes sum to {total_non_id}, expected {2*n - 1}"

    def test_print_structure(self):
        """Print the class structure for small n to verify."""
        for n in [2, 3, 4, 5, 6, 7, 8, 12, 15, 20]:
            classes = dihedral_classes(n)
            all_even = all(K % 2 == 0 for (_, _, _, K, _, _) in classes)
            print(f"\nD_{n} (|D_n|={2*n}):")
            for desc, size, d, K, alpha, k_val in classes:
                parity = "even" if K % 2 == 0 else "ODD"
                print(f"  {desc:>12}: size={size}, d={d}, K=2n/d={K} ({parity}), α={alpha:.4f}")
            print(f"  (B) holds: {all_even} ({'n even' if n%2==0 else 'n odd'})")


class TestModeAntisymmetry:
    """Verify mode antisymmetry f(α, j*-j) = (-1)^{K+1} f(α, j) numerically."""

    def test_antisymmetry_obstructed(self):
        """For n even: all modes antisymmetric, duality holds.

        The antisymmetry sign at half-angle α = k'π/d is (-1)^{k'K+1}.
        For rotations: K = 2·gcd(k,n), so k'K = k'·2·gcd(k,n) is always even.
        For reflections: k' = 1, K = n. k'K = n, even when n even.
        So all modes are antisymmetric when n is even.
        """
        for n in [2, 4, 6, 8, 10, 20, 50, 100]:
            jstar = n - 1  # |G|/4 - 1 = 4n/4 - 1 = n - 1
            classes = dihedral_classes(n)
            for desc, size, d, K, alpha, k_val in classes:
                # k' = reduced numerator: k/gcd(k,n) for rotations, 1 for reflections
                kprime = 1 if k_val is None else k_val // np.gcd(k_val, n)
                sign_exp = kprime * K  # (-1)^{sign_exp+1} is the antisymmetry sign

                for j in range(jstar + 1):
                    sa = np.sin(alpha)
                    if abs(sa) < 1e-14:
                        f_j = 2*j + 1
                        f_paired = 2*(jstar-j) + 1
                    else:
                        f_j = np.sin((2*j + 1) * alpha) / sa
                        f_paired = np.sin((2*(jstar-j) + 1) * alpha) / sa

                    # f(α, j*-j) = (-1)^{k'K+1} f(α, j)
                    expected = (-1)**(sign_exp + 1) * f_j
                    assert abs(f_paired - expected) < 1e-10, \
                        f"n={n}, {desc}, j={j}: f_paired={f_paired:.6f}, expected={expected:.6f}"

                    # For n even: sign_exp is always even, so (-1)^{sign_exp+1} = -1
                    assert sign_exp % 2 == 0, \
                        f"n={n}, {desc}: k'K={sign_exp} should be even!"
                    assert abs(f_j + f_paired) < 1e-10, \
                        f"n={n}, {desc}, j={j}: modes should cancel!"

    def test_antisymmetry_non_obstructed(self):
        """For n odd: reflection mode is symmetric (k'K = n odd → sign = +1)."""
        for n in [3, 5, 7, 9, 11, 21, 51]:
            jstar = n - 1
            classes = dihedral_classes(n)
            for desc, size, d, K, alpha, k_val in classes:
                kprime = 1 if k_val is None else k_val // np.gcd(k_val, n)
                sign_exp = kprime * K

                for j in range(jstar + 1):
                    sa = np.sin(alpha)
                    if abs(sa) < 1e-14:
                        continue
                    f_j = np.sin((2*j + 1) * alpha) / sa
                    f_paired = np.sin((2*(jstar-j) + 1) * alpha) / sa

                    expected = (-1)**(sign_exp + 1) * f_j
                    assert abs(f_paired - expected) < 1e-10, \
                        f"n={n}, {desc}, j={j}: antisymmetry formula fails"

                    if k_val is None:  # reflection class
                        # k'=1, K=n (odd): sign_exp=n odd, (-1)^{n+1} = +1: SYMMETRIC
                        assert sign_exp % 2 == 1
                        assert abs(f_paired - f_j) < 1e-10, \
                            f"n={n}, {desc}, j={j}: reflection mode should be SYMMETRIC"


class TestK_Formula:
    """The key formula: K = 2·gcd(k, n) for rotations, K = n for reflections."""

    def test_K_rotation_formula(self):
        """K = |D_n|/d = 2n/d. For a^k: d = n/gcd(k,n), so K = 2·gcd(k,n)."""
        print("\nKey analytic result:")
        print("  For rotation a^k in D_n: K = 2·gcd(k, n). ALWAYS EVEN.")
        print("  Proof: d = ord(a^k) = n/gcd(k,n). K = 2n/d = 2n·gcd(k,n)/n = 2·gcd(k,n).")
        print()
        print("  For reflection a^k·b in D_n: d = 2. K = 2n/2 = n.")
        print("  K even ⟺ n even.")
        print()
        print("  Therefore: (B) holds for 2D_n ⟺ n even. QED.")
        print()
        print("  This eliminates the ENTIRE infinite dicyclic family from")
        print("  case-by-case verification. Only 3 polyhedral groups (2T, 2O, 2I)")
        print("  remain as finite checks.")

    def test_polyhedral_K_values(self):
        """Cross-check: K = |H|/d for all polyhedral SO(3) element orders.

        For each polyhedral group H, we only need to check that K is even
        for each distinct SO(3) order d that appears among non-identity elements.
        The collapsed coefficients from the paper (element counts in the binary
        group 2G, not in H) are shown for reference but not used in assertions.
        """
        # Format: (group_name, |H|, list of (d, count_in_H, description))
        polyhedral_data = [
            ("T", 12, [
                (3, 8, "vertex rotations"),
                (2, 3, "edge rotations"),
            ]),
            ("O", 24, [
                (4, 6, "face rotations 90/270"),
                (3, 8, "vertex rotations 120/240"),
                (2, 9, "face 180 (3) + edge 180 (6)"),
            ]),
            ("I", 60, [
                (5, 24, "vertex rotations 72/144/216/288"),
                (3, 20, "face rotations 120/240"),
                (2, 15, "edge rotations 180"),
            ]),
        ]

        for name, H_order, orders in polyhedral_data:
            print(f"\n{name} (|H|={H_order}):")
            total_non_id = 0
            for d, count, desc in orders:
                K = H_order // d
                print(f"  d={d}: {count} elements, K={H_order}/{d}={K} (even={K%2==0})  [{desc}]")
                assert K % 2 == 0, f"{name}: K={K} odd at d={d}!"
                total_non_id += count
            assert total_non_id == H_order - 1, \
                f"{name}: element counts sum to {total_non_id}, expected {H_order - 1}"


class TestDualityNumerical:
    """Numerical verification that the analytic K analysis gives correct duality."""

    def _multiplicity_via_classes(self, n, j):
        """Compute m(A₀, j) for 2D_n using class structure.

        m(A₀, j) = (1/|G|) [2(2j+1) + S(j)]
        where S(j) = sum over non-central classes of A_α · f(α, j).
        For 2D_n: |G| = 4n.
        Central contribution: g=1 gives χ_j(1)=2j+1, g=-1 gives χ_j(-1)=2j+1 (j integer).
        So central = 2(2j+1).
        Non-central: each element g with half-angle α contributes f(α, j).
        But we need elements of the BINARY group 2D_n, not just D_n.
        Each non-identity h ∈ D_n lifts to a pair {g, -g} in 2D_n.
        For integer j: χ_j(-g) = χ_j(g). So each pair contributes 2·f(α, j).
        Total non-central contribution = 2 * Σ_{h in D_n, h!=e} f(α_h, j).
        """
        G_order = 4 * n
        jstar = n - 1

        # Central contribution
        central = 2 * (2 * j + 1)

        # Non-central: sum over D_n elements
        S = 0.0
        classes = dihedral_classes(n)
        for desc, size, d, K, alpha, k_val in classes:
            sa = np.sin(alpha)
            if abs(sa) < 1e-14:
                f_val = 2 * j + 1
            else:
                f_val = np.sin((2 * j + 1) * alpha) / sa
            # Each class has 'size' elements in D_n, each contributes f_val
            # In the binary group, each lifts to a pair contributing 2·f_val
            S += 2 * size * f_val

        m = (central + S) / G_order
        m_rounded = round(m)
        assert abs(m - m_rounded) < 1e-8, \
            f"n={n}, j={j}: m={m} not close to integer (rounding would mask error)"
        return m_rounded

    def test_duality_via_class_formula(self):
        """Verify m(A₀, j) + m(A₀, j*-j) = 1 using class-based computation."""
        for n in [2, 4, 6, 8, 10, 12, 20, 50]:
            jstar = n - 1
            for j in range(jstar + 1):
                m_j = self._multiplicity_via_classes(n, j)
                m_paired = self._multiplicity_via_classes(n, jstar - j)
                assert m_j + m_paired == 1, \
                    f"2D_{n}: m({j}) + m({jstar-j}) = {m_j} + {m_paired} ≠ 1"

        print("Duality via class formula: verified for n = 2,4,6,8,10,12,20,50")

    def test_defect_via_class_formula(self):
        """For n odd: compute defect δ(j) = m(j) + m(j*-j) - 1."""
        for n in [3, 5, 7, 9, 11]:
            jstar = n - 1
            print(f"\n2D_{n} (n odd, non-obstructed):")
            for j in range(jstar + 1):
                m_j = self._multiplicity_via_classes(n, j)
                m_paired = self._multiplicity_via_classes(n, jstar - j)
                delta = m_j + m_paired - 1
                expected = (-1)**j
                print(f"  j={j}: m={m_j}, m'={m_paired}, δ={delta} (expected (-1)^j={expected})")
                assert delta == expected, f"n={n}, j={j}: δ={delta} ≠ (-1)^j"


class TestSummary:
    """Summary of the analytic proof."""

    def test_summary(self):
        """Print the complete analytic argument."""
        print("""
╔══════════════════════════════════════════════════════════════════╗
║  ANALYTIC PROOF: (B) for 2D_n ⟺ n even                        ║
╠══════════════════════════════════════════════════════════════════╣
║                                                                  ║
║  D_n has two types of non-identity elements:                     ║
║                                                                  ║
║  1. ROTATIONS a^k (k = 1, ..., n-1):                            ║
║     SO(3) order: d = n/gcd(k, n)                                ║
║     K = |D_n|/d = 2n/d = 2·gcd(k, n)                           ║
║     K is ALWAYS EVEN. ✓                                          ║
║                                                                  ║
║  2. REFLECTIONS a^k·b:                                           ║
║     SO(3) order: d = 2                                           ║
║     K = |D_n|/d = 2n/2 = n                                      ║
║     K even ⟺ n even.                                            ║
║                                                                  ║
║  Result: (B) holds ⟺ n even. QED.                              ║
║                                                                  ║
║  Combined with (A) ⟺ (B) via classification:                   ║
║  -1 ∈ [G,G] ⟺ n even (for 2D_n).                              ║
║                                                                  ║
║  Paper impact: eliminates infinite dicyclic family from          ║
║  case-by-case check. Only 2T, 2O, 2I need finite verification. ║
║                                                                  ║
║  The key identity: K_rot = 2·gcd(k, n) is structural.           ║
║  It follows from Lagrange's theorem applied to ⟨a^k⟩ ⊂ ⟨a⟩.    ║
╚══════════════════════════════════════════════════════════════════╝
""")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
