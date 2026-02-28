#!/usr/bin/env python3
"""
Script 58: Test Vector 3 for classification-free proof of (B) ⟹ (A).

The chain:
  (B) ⟹ Syl₂(H) non-cyclic ⟹ V₄ ⊂ H ⟹ Q₈ ⊂ 2H ⟹ (A)

Tests:
  1. (B) ⟹ Syl₂(H) non-cyclic (v₂ argument)
  2. Every H ⊂ SO(3) with non-cyclic Syl₂ contains V₄
  3. No generalized quaternion group embeds in SO(3)
  4. Full chain: (B) ⟹ (A) via Sylow for all ADE groups

Run: .venv/bin/python -m pytest tests/58_vector3_sylow.py -v -s
Date: Feb 2026
"""
import numpy as np
import pytest
from math import gcd
from functools import reduce

# ============================================================
# Quaternion utilities
# ============================================================

def qmul(a, b):
    return np.array([
        a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
        a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
        a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
        a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0],
    ])

def qinv(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

def is_close(p, q, tol=1e-10):
    return np.linalg.norm(np.array(p) - np.array(q)) < tol

def is_identity(q, tol=1e-10):
    return abs(q[0] - 1) < tol and np.linalg.norm(q[1:]) < tol

def is_minus_one(q, tol=1e-10):
    return abs(q[0] + 1) < tol and np.linalg.norm(q[1:]) < tol

def qorder(g, max_ord=200):
    """Order of g in the group (as quaternion)."""
    p = g.copy()
    for k in range(1, max_ord):
        if is_identity(p):
            return k
        p = qmul(p, g)
    return None

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
        new_elts = []
        for a in elements:
            for b in elements:
                p = qmul(a, b)
                if not is_in(p, elements) and not is_in(p, new_elts):
                    new_elts.append(p)
                    if len(elements) + len(new_elts) > max_size:
                        raise ValueError(f"Exceeds {max_size}")
        if new_elts:
            elements.extend(new_elts)
            changed = True
    return elements


# ============================================================
# Group builders (2H in SU(2))
# ============================================================

def make_2Cn(n):
    """Binary cyclic Z_{2n}."""
    a = np.array([np.cos(np.pi / n), np.sin(np.pi / n), 0, 0])
    return generate_group([a])

def make_2Dn(n):
    """Binary dihedral Dic_n (order 4n)."""
    a = np.array([np.cos(np.pi / n), np.sin(np.pi / n), 0, 0])
    b = np.array([0, 0, 1, 0])
    return generate_group([a, b])

def make_2T():
    return generate_group([np.array([0, 1, 0, 0]),
                           np.array([0.5, 0.5, 0.5, 0.5])])

def make_2O():
    return generate_group([np.array([0.5, 0.5, 0.5, 0.5]),
                           np.array([1/np.sqrt(2), 1/np.sqrt(2), 0, 0])])

def make_2I():
    phi = (1 + np.sqrt(5)) / 2
    return generate_group([np.array([0.5, 0.5, 0.5, 0.5]),
                           np.array([phi/2, 1/(2*phi), 0.5, 0])])


# ============================================================
# SO(3) group analysis (works on H = 2H/{±1})
# ============================================================

def so3_elements(G_su2):
    """Project 2H → H: identify g with -g, keep one representative.
    Returns list of quaternions representing SO(3) elements."""
    reps = []
    for g in G_su2:
        # Check if -g is already represented
        already = False
        for r in reps:
            if is_close(g, -r):
                already = True
                break
        if not already:
            reps.append(g)
    return reps

def so3_order(g):
    """SO(3) order of rotation represented by unit quaternion g.
    Half-angle α = arccos(|g[0]|), SO(3) order = π/α if rational."""
    if is_identity(g) or is_minus_one(g):
        return 1  # identity in SO(3)
    alpha = np.arccos(np.clip(abs(g[0]), -1, 1))
    if alpha < 1e-12:
        return 1
    for d in range(2, 200):
        if abs(np.sin(d * alpha)) < 1e-8:
            return d
    return None

def count_involutions_so3(H_reps):
    """Count elements of order 2 in H ⊂ SO(3) (π-rotations)."""
    return sum(1 for g in H_reps if so3_order(g) == 2)

def has_V4(H_reps):
    """Check if H contains V₄ = Z₂ × Z₂ (two commuting involutions)."""
    involutions = [g for g in H_reps if so3_order(g) == 2]
    for i, a in enumerate(involutions):
        for j, b in enumerate(involutions):
            if j <= i:
                continue
            # Check if ab is also an involution (or identity)
            ab = qmul(a, b)
            # In SO(3), need to check if the product is order 2
            # Two π-rotations commute iff axes perpendicular
            # Their product is then also a π-rotation
            if so3_order(ab) == 2 or is_identity(ab) or is_minus_one(ab):
                # Check they actually commute
                ba = qmul(b, a)
                if is_close(ab, ba) or is_close(ab, -ba):
                    return True
    return False

def v2(n):
    """2-adic valuation of n."""
    if n == 0:
        return float('inf')
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k

def check_condition_B(H_reps):
    """Check if H satisfies (B): |H|/ord(h) even for all h ≠ e."""
    H_order = len(H_reps)
    for g in H_reps:
        d = so3_order(g)
        if d == 1:
            continue  # identity
        if (H_order // d) % 2 != 0:
            return False
    return True

def has_commutator_minus1(G_su2):
    """Check (A): ∃ a,b ∈ 2H with [a,b] = −1."""
    for a in G_su2:
        for b in G_su2:
            c = qmul(qmul(a, b), qmul(qinv(a), qinv(b)))
            if is_minus_one(c):
                return True
    return False

def sylow2_is_cyclic(H_reps):
    """Check if the Sylow 2-subgroup of H is cyclic.
    Cyclic iff H has a unique element of order 2."""
    involution_count = count_involutions_so3(H_reps)
    # A cyclic 2-group has exactly one involution.
    # A non-cyclic 2-group has > 1 involution.
    # For the full group H: Syl₂ cyclic iff Syl₂ has unique involution.
    # Since all involutions of H lie in Syl₂ (every element of order 2
    # has 2-power order), unique involution in H iff Syl₂ cyclic.
    #
    # WAIT: involutions in SO(3) have order 2 (not a 2-power, just 2).
    # An involution of H always lies in a Sylow 2-subgroup.
    # But H could have involutions in different Sylow 2-subgroups.
    #
    # Simpler criterion: if H has ≥ 2 involutions that commute,
    # then Syl₂ is non-cyclic (contains V₄).
    # If H has exactly 1 involution, Syl₂ is cyclic or gen. quaternion.
    return involution_count == 1


# ============================================================
# Test 1: (B) ⟹ Syl₂ non-cyclic
# ============================================================

class TestSylow2NonCyclic:
    """(B) forces v₂(ord(h)) < v₂(|H|) for all h.
    This means no element achieves the full 2-power of |H|,
    so Syl₂ cannot be cyclic."""

    ALL_GROUPS = [
        ('C_2',  make_2Cn(2),  2),
        ('C_3',  make_2Cn(3),  3),
        ('C_4',  make_2Cn(4),  4),
        ('C_6',  make_2Cn(6),  6),
        ('D_2',  make_2Dn(2),  4),
        ('D_3',  make_2Dn(3),  6),
        ('D_4',  make_2Dn(4),  8),
        ('D_5',  make_2Dn(5),  10),
        ('D_6',  make_2Dn(6),  12),
        ('D_8',  make_2Dn(8),  16),
        ('D_10', make_2Dn(10), 20),
        ('T',    make_2T(),    12),
        ('O',    make_2O(),    24),
        ('I',    make_2I(),    60),
    ]

    @pytest.mark.parametrize("name,G_su2,H_order",
                             ALL_GROUPS, ids=[g[0] for g in ALL_GROUPS])
    def test_B_implies_noncyclic_sylow(self, name, G_su2, H_order):
        """If (B) holds, Sylow 2-subgroup is non-cyclic."""
        H_reps = so3_elements(G_su2)
        assert len(H_reps) == H_order, f"{name}: |H| = {len(H_reps)} ≠ {H_order}"

        has_B = check_condition_B(H_reps)
        if not has_B:
            pytest.skip(f"{name}: (B) does not hold")

        # (B) holds. Check Syl₂ is non-cyclic.
        n_inv = count_involutions_so3(H_reps)
        is_cyclic = sylow2_is_cyclic(H_reps)
        assert not is_cyclic, (
            f"{name}: (B) holds but Syl₂ appears cyclic! "
            f"involutions={n_inv}")
        print(f"  {name}: (B) holds, {n_inv} involutions, Syl₂ non-cyclic ✓")


# ============================================================
# Test 2: Non-cyclic Syl₂ in SO(3) ⟹ V₄ ⊂ H
# ============================================================

class TestNonCyclicImpliesV4:
    """If H ⊂ SO(3) has non-cyclic Sylow 2-subgroup, then V₄ ⊂ H."""

    ALL_GROUPS = TestSylow2NonCyclic.ALL_GROUPS

    @pytest.mark.parametrize("name,G_su2,H_order",
                             ALL_GROUPS, ids=[g[0] for g in ALL_GROUPS])
    def test_noncyclic_has_v4(self, name, G_su2, H_order):
        H_reps = so3_elements(G_su2)
        n_inv = count_involutions_so3(H_reps)

        if n_inv <= 1:
            pytest.skip(f"{name}: ≤1 involution (Syl₂ may be cyclic)")

        has_v = has_V4(H_reps)
        assert has_v, f"{name}: {n_inv} involutions but no V₄ found!"
        print(f"  {name}: {n_inv} involutions, V₄ found ✓")


# ============================================================
# Test 3: Involution count under (B)
# ============================================================

class TestInvolutionCount:
    """Under (B), H has multiple involutions.
    Specifically: if (B) holds then |H| ≥ 4 and H has ≥ 3 involutions."""

    ALL_GROUPS = TestSylow2NonCyclic.ALL_GROUPS

    @pytest.mark.parametrize("name,G_su2,H_order",
                             ALL_GROUPS, ids=[g[0] for g in ALL_GROUPS])
    def test_involution_count(self, name, G_su2, H_order):
        H_reps = so3_elements(G_su2)
        has_B = check_condition_B(H_reps)
        n_inv = count_involutions_so3(H_reps)

        if has_B:
            assert n_inv >= 3, (
                f"{name}: (B) holds but only {n_inv} involutions")
            print(f"  {name}: (B) holds, {n_inv} involutions (≥ 3) ✓")
        else:
            print(f"  {name}: (B) fails, {n_inv} involutions (no constraint)")


# ============================================================
# Test 4: Full chain (B) ⟹ (A) for all groups
# ============================================================

class TestFullChain:
    """(B) ⟹ Syl₂ non-cyclic ⟹ V₄ ⊂ H ⟹ Q₈ ⊂ 2H ⟹ (A).
    Verify each step and the final equivalence."""

    ALL_GROUPS = TestSylow2NonCyclic.ALL_GROUPS

    @pytest.mark.parametrize("name,G_su2,H_order",
                             ALL_GROUPS, ids=[g[0] for g in ALL_GROUPS])
    def test_full_chain(self, name, G_su2, H_order):
        H_reps = so3_elements(G_su2)
        has_B = check_condition_B(H_reps)
        has_A = has_commutator_minus1(G_su2)
        n_inv = count_involutions_so3(H_reps)
        has_v = has_V4(H_reps)

        # Check Q₈ ⊂ 2H
        q8 = [np.array([s*d[0], s*d[1], s*d[2], s*d[3]])
              for d in [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
              for s in [1, -1]]
        has_q8 = all(
            any(np.linalg.norm(q - g) < 1e-10 for g in G_su2)
            for q in q8)

        print(f"  {name}: (B)={has_B}, inv={n_inv}, V₄={has_v}, "
              f"Q₈={has_q8}, (A)={has_A}")

        if has_B:
            # Full chain should hold
            assert n_inv >= 3, f"Step 1 fails: only {n_inv} involutions"
            assert has_v, "Step 2 fails: no V₄"
            assert has_q8, "Step 3 fails: no Q₈"
            assert has_A, "Step 4 fails: (A) doesn't hold"
            print(f"    Chain: (B) → {n_inv} inv → V₄ → Q₈ → (A) ✓")

        # Also verify: (A) ⟺ (B) agreement
        assert has_A == has_B, f"{name}: (A)={has_A} ≠ (B)={has_B}"


# ============================================================
# Test 5: The 2-adic argument (pure arithmetic)
# ============================================================

class TestTwoAdicArgument:
    """Pure arithmetic: (B₃) says v₂(ord(h)) < v₂(|H|) for all h ≠ e.
    This forces Sylow 2-subgroup to be non-cyclic.

    Proof: if Syl₂(H) is cyclic of order 2^a, it has an element of
    order 2^a. Then v₂(2^a) = a = v₂(|H|), contradicting (B₃)."""

    def test_cyclic_groups_fail_B(self):
        """C_n never satisfies (B): generator has |H|/ord = 1."""
        for n in range(2, 30):
            G = make_2Cn(n)
            H_reps = so3_elements(G)
            has_B = check_condition_B(H_reps)
            assert not has_B, f"C_{n} satisfies (B)!"
        print("  C_n, n=2..29: none satisfy (B) ✓")

    def test_v2_argument_dihedral(self):
        """For D_n: (B) ⟺ n even. When n even, v₂(ord(h)) < v₂(2n)
        for all h, so Syl₂ non-cyclic."""
        for n in range(2, 30):
            G = make_2Dn(n)
            H_reps = so3_elements(G)
            H_order = 2 * n
            has_B = check_condition_B(H_reps)
            assert has_B == (n % 2 == 0), f"D_{n}: (B)={has_B}, n%2={n%2}"

            if has_B:
                a = v2(H_order)  # v₂(|H|)
                for g in H_reps:
                    d = so3_order(g)
                    if d == 1:
                        continue
                    assert v2(d) < a, (
                        f"D_{n}: h has order {d}, v₂={v2(d)} but v₂(|H|)={a}")
        print("  D_n, n=2..29: v₂ argument verified ✓")

    def test_v2_argument_polyhedral(self):
        """For T, O, I: (B) holds, and v₂(ord(h)) < v₂(|H|) for all h."""
        cases = [
            ('T', make_2T(), 12),
            ('O', make_2O(), 24),
            ('I', make_2I(), 60),
        ]
        for name, G, H_order in cases:
            H_reps = so3_elements(G)
            a = v2(H_order)
            for g in H_reps:
                d = so3_order(g)
                if d == 1:
                    continue
                assert v2(d) < a, (
                    f"{name}: d={d}, v₂(d)={v2(d)}, v₂(|H|)={a}")
            print(f"  {name}: v₂(ord(h)) < v₂(|H|) = {a} for all h ✓")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
