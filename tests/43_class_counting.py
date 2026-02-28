"""
Script 43: Conjugacy class counting for G → H = G/{±1}.

Key question: When does 2·#classes(H) > #classes(G)?
(equivalently: n_+ > n_-, equivalently: -1 is a single commutator)

For the 2:1 cover G → H:
- Each conjugacy class C of G projects to a conjugacy class of H.
- Two classes C and -C (where -C = {-g : g ∈ C}) project to the same class in H.
- Cases:
  (a) C = -C (self-conjugate under ±1): C projects to one class in H
  (b) C ≠ -C: the pair {C, -C} projects to one class in H

Let k = #classes(G), s = # self-conjugate classes, p = # paired classes.
Then k = s + 2p (each non-self-conjugate class has a partner).
And #classes(H) = s + p (each self-conjugate class → 1, each pair → 1).

So 2·#classes(H) - #classes(G) = 2(s+p) - (s+2p) = s.

Therefore: n_+ - n_- = s = # self-conjugate classes of G under g ↦ -g.

A class C is self-conjugate iff -g ∈ C for every g ∈ C.
For the central classes {1} and {-1}: these are ALWAYS self-conjugate (each is its own class, and -1 ∈ {-1}).
Wait, {1} is not self-conjugate: -{1} = {-1} ≠ {1}. So {1} and {-1} form a PAIR.

Hmm, unless 1 and -1 are in the same conjugacy class, which they're not (conjugacy preserves trace in SU(2), and tr(1) ≠ tr(-1)).

So {1} and {-1} always form a pair (p → 1), NOT self-conjugate.

Self-conjugate classes are those where C = -C = {g, -g, ...}.
In SU(2), g and -g have the same axis but opposite angles (θ vs 2π-θ).
They're in the same conjugacy class iff they have the same trace, but tr(g) = 2cos(θ/2) and tr(-g) = -2cos(θ/2). These are equal iff cos(θ/2) = 0, i.e., θ/2 = π/2, i.e., θ = π.

So the ONLY self-conjugate classes are those at half-angle α = π/2 (rotation angle π).

Elements of order 4 in SU(2) (i.e., half-angle π/2, rotation angle π in SO(3), order 2 in SO(3)): these have g and -g in the same conjugacy class of G.

Let me verify this for each ADE group.
"""

import numpy as np
import pytest

# ============================================================
# Count self-conjugate classes for each ADE group
# ============================================================

class TestClassCounting:
    """Count self-conjugate classes (C = -C) for each finite SU(2) subgroup."""

    def test_cyclic(self):
        """Z_{2n}: classes are {exp(ikπ/n)} for k = 0, ..., 2n-1.
        -g = exp(i(k+n)π/n). So C_k = -C_k iff k ≡ k+n (mod 2n) iff n ≡ 0.
        That's impossible for n > 0. So s = 0 for all cyclic groups.
        """
        for n in [1, 2, 3, 4, 5, 10]:
            s = 0  # no self-conjugate classes in cyclic groups
            print(f"Z_{2*n}: s = {s}, n_+ - n_- = {s}, "
                  f"N(-1) = {2*n * s} → {'single comm' if s > 0 else 'not single comm'}")
            # For cyclic: -1 ∉ [G,G] (abelian), consistent with s = 0

    def test_binary_dihedral(self):
        """2D_n (order 4n). Conjugacy classes:
        - {1}, {-1}: 2 central classes, form a pair, NOT self-conjugate
        - For k = 1, ..., n-1: {a^k, a^{-k}} and {-a^k, -a^{-k}}
          These are paired: C_k and -C_k. Self-conjugate iff a^k = -a^{-k} etc.
          a^k has half-angle kπ/n. -a^k has half-angle π - kπ/n.
          Same class iff kπ/n = π - kπ/n, i.e., k = n/2.
          So C_{n/2} is self-conjugate iff n is even.
        - {b, -b, ab, -ab, ...}: classes of "b-type" elements
          For n ≥ 3: two classes of b-type if n even, one if n odd.
          Each b-type element has half-angle π/2 (order 4 in SU(2)).
          -b is in the same b-type class: self-conjugate!

        For n even: self-conjugate classes =
          - C_{n/2} (the "equatorial" rotation class): 1
          - two b-type classes: 2
          Total s = 3.
        For n odd: self-conjugate classes =
          - no equatorial class (n/2 not integer)
          - one b-type class: but it may or may not be self-conjugate...
          Actually for n odd: the b-type elements have half-angle π/2.
          There's one class of b-type: {b, ab, a²b, ..., a^{n-1}b} ∪ {-b, -ab, ...}
          Wait, I need to be more careful.

        Let me just count s = 2·#classes(H) - #classes(G) and verify.
        """
        for n in range(2, 20):
            # Known from script 42:
            k_G = n + 3  # conj classes of 2D_n (for n ≥ 2)
            if n == 2:
                k_G = 5  # Q₈ special case

            if n % 2 == 0:
                k_H = n // 2 + 3
            else:
                k_H = (n + 3) // 2

            s = 2 * k_H - k_G

            obstructed = (n % 2 == 0)
            print(f"2D_{n:>2}: k_G={k_G:>2}, k_H={k_H:>2}, "
                  f"s = 2×{k_H}-{k_G} = {s}, "
                  f"obstructed={'Yes' if obstructed else 'No'}")

            if obstructed:
                assert s == 3, f"Expected s=3 for even n={n}"
            else:
                assert s == 0, f"Expected s=0 for odd n={n}"

    def test_polyhedral(self):
        """Polyhedral groups: count s directly."""
        data = [
            ('2T', 7, 4, True),   # k_G=7, k_H=4
            ('2O', 8, 5, True),   # k_G=8, k_H=5
            ('2I', 9, 5, True),   # k_G=9, k_H=5
        ]

        for name, k_G, k_H, obstructed in data:
            s = 2 * k_H - k_G
            print(f"{name}: k_G={k_G}, k_H={k_H}, s={s}")
            if obstructed:
                assert s > 0, f"{name}: s should be > 0!"

    def test_structural_theorem(self):
        """THE KEY RESULT:

        n_+ - n_- = s = # self-conjugate conjugacy classes of G (under g ↦ -g).

        A class C is self-conjugate iff it contains both g and -g.
        In SU(2), this happens iff the half-angle is π/2
        (elements of order 4 in SU(2), i.e., order 2 in SO(3)).

        For obstructed groups (4 | |G|, -1 ∈ [G,G]):
        - There exist elements of half-angle π/2 (involutions in H)
        - These form self-conjugate classes
        - Therefore s > 0, hence N(-1) > 0, hence -1 is a single commutator

        The structural argument:
        1. -1 ∈ [G,G] implies obstruction (by definition)
        2. Frobenius formula: N(-1) = |G| · s where s = # self-conjugate classes
        3. s = # conjugacy classes at half-angle π/2
        4. If G has elements of SU(2)-order 4 (half-angle π/2),
           these form self-conjugate classes, giving s > 0.
        5. Question: does -1 ∈ [G,G] guarantee existence of order-4 elements?

        For SU(2) subgroups: -1 has order 2. If -1 ∈ [G,G] = [a,b]·...
        This requires elements that don't commute, which requires elements
        of order ≥ 4 in SU(2) (order ≥ 2 in SO(3)).
        """
        print("Structural summary:")
        print("  n_+ - n_- = s = # self-conjugate classes")
        print("  Self-conjugate ⟺ half-angle = π/2 (SU(2)-order 4)")
        print()

        all_groups = [
            # (name, k_G, k_H, has_order4_elements)
            ('Z_2',  2, 1, False),   # only ±1
            ('Z_4',  4, 2, True),    # has i
            ('Z_6',  6, 3, False),   # no pure imaginary elements
            ('2D_2 (Q₈)', 5, 4, True),
            ('2D_3', 6, 3, True),    # has elements at π/2
            ('2D_4', 7, 5, True),
            ('2D_5', 8, 4, True),
            ('2D_6', 9, 6, True),
            ('2T', 7, 4, True),
            ('2O', 8, 5, True),
            ('2I', 9, 5, True),
        ]

        for name, k_G, k_H, has_ord4 in all_groups:
            s = 2 * k_H - k_G
            print(f"  {name:>12}: s = {s}, has order-4 elements: {has_ord4}")

        # Key observation: 2D_3 and 2D_5 HAVE order-4 elements (b ∈ 2D_n has order 4)
        # but s = 0 for odd n. Why?
        # Because for odd n, the b-type class is NOT self-conjugate.
        # Actually wait: -b = a^n · b (since a^n = -1). So -b is in the b-orbit.
        # If n is odd, a^n b is NOT conjugate to b in 2D_n.
        # Hmm, this needs more careful analysis.

        print()
        print("Note: 2D_3 has order-4 elements but s=0.")
        print("This is because the b-type class is not self-conjugate for odd n.")
        print("Specifically: -b = a^n·b, and a^n·b ∉ same conjugacy class as b when n odd.")
        print()
        print("KEY: s > 0 ⟺ n_+ > n_- ⟺ N(-1) > 0 ⟺ -1 is single commutator")
        print("And from ADE data: s > 0 ⟺ -1 ∈ [G,G] (obstruction)")
        print()
        print("Can we prove s > 0 ⟹ -1 ∈ [G,G] structurally? YES: Frobenius formula.")
        print("Can we prove -1 ∈ [G,G] ⟹ s > 0 structurally? THIS IS THE QUESTION.")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
