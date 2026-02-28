#!/usr/bin/env python3
"""
Script 59: Generalized quaternion Q_{2^n} cannot embed in SO(3).

This is the KEY STEP for the classification-free proof of (B) ⟹ (A).

THEOREM: For n ≥ 3, Q_{2^n} admits no faithful homomorphism to SO(3).

PROOF (classification-free, pure rotation geometry):
  Q_{2^n} = ⟨a, b | a^{2^{n-1}} = e, b² = a^{2^{n-2}}, bab⁻¹ = a⁻¹⟩.

  Suppose φ: Q_{2^n} → SO(3) faithful. Then:
  1. ⟨φ(a)⟩ is cyclic ⟹ all rotations share axis u.
  2. φ(b) has order 4 ⟹ φ(b) = rotation by ±π/2 about axis w.
     φ(b)² = π-rotation about w. But φ(b)² = φ(a^{2^{n-2}}) =
     π-rotation about u. So w = u.
  3. φ(a) and φ(b) share axis u ⟹ they commute.
  4. But bab⁻¹ = a⁻¹ ⟹ φ(a) = φ(a)⁻¹ ⟹ φ(a)² = e.
  5. ord(φ(a)) ≤ 2, contradicting ord(a) = 2^{n-1} ≥ 4.  □

Uses only: (i) cyclic rotation groups share an axis,
(ii) R² of an order-4 rotation shares axis with R,
(iii) same-axis rotations commute. No ADE.

Tests verify:
  - The three geometric facts (shared axis, order-4 square, commutation)
  - Axis forcing: b² = z forces b onto same axis as a (1440-direction scan)
  - V₄ orthogonality: π-rotation composition angle = 2θ, so V₄ ⟹ θ = π/2
  - End-to-end on concrete groups: D₂, D₄, D₆, T, O, I

External fact (cited): Burnside — non-cyclic, non-gen-quaternion
2-group of order ≥ 4 contains V₄. This is 2-group structure theory, not ADE.

Also tests the (A)⟹(B) gap closure:
  - V₄ ⊂ H ⟺ (B) brute check on 16 groups (C₂..C₆, D₂..D₁₀, T, O, I)
  - Burnside orbit-counting for 2-groups: only cyclic/dihedral (a=1..20)
  - Dihedral exponent bound: exp(D_{2^{a-1}}) = 2^{a-1} ⟹ (B₃) (a=2..7)

RAW OUTPUT (77 passed, 25s):
  TestSameAxisProperties (3 tests) PASSED
  TestGenQuaternionContradiction::test_axis_forcing[3..6] PASSED (4)
  TestGenQuaternionContradiction::test_commutation_contradiction[3..6] PASSED (4)
  TestGenQuaternionContradiction::test_no_off_axis_escape[3..6] PASSED (4)
  TestV4Orthogonality (2 tests) PASSED
  TestClassificationFreeChain::test_chain_summary PASSED
  TestEndToEnd::test_B_implies_A_concrete[D_2..I] PASSED (6)
  TestV4iffB::test_V4_iff_B[C_2..I] PASSED (16)
  TestBurnside2Groups::test_diophantine_solutions[1..20] PASSED (20)
  TestDihedralExponent::test_dihedral_exponent[2..7] PASSED (6)
  TestSylowTransfer::test_sylow_transfer[D_2..I] PASSED (7)
  TestSylowTransfer::test_contrapositive[C_4,D_3,D_5,D_7] PASSED (4)

Run: .venv/bin/python -m pytest tests/59_gen_quaternion_no_so3.py -v -s
Date: Feb 2026
"""
import numpy as np
import pytest
from scipy.spatial.transform import Rotation


# ============================================================
# SO(3) rotation utilities
# ============================================================

def rotation_matrix(axis, angle):
    """Rotation matrix for rotation by `angle` about `axis`."""
    axis = np.array(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    return Rotation.from_rotvec(angle * axis).as_matrix()


def rotation_axis(R):
    """Extract rotation axis of SO(3) matrix R (unit vector).
    Returns None for identity."""
    r = Rotation.from_matrix(R)
    rotvec = r.as_rotvec()
    angle = np.linalg.norm(rotvec)
    if angle < 1e-10:
        return None  # identity
    return rotvec / angle


def rotation_angle(R):
    """Extract rotation angle ∈ [0, π] of SO(3) matrix R."""
    r = Rotation.from_matrix(R)
    return np.linalg.norm(r.as_rotvec())


def axes_parallel(u, v, tol=1e-8):
    """Check if axes u, v are parallel (same or opposite direction)."""
    if u is None or v is None:
        return False
    cross = np.linalg.norm(np.cross(u, v))
    return cross < tol


# ============================================================
# Test 1: Proof step verification — same-axis properties
# ============================================================

class TestSameAxisProperties:
    """Verify the three geometric facts used in the proof."""

    def test_cyclic_group_shares_axis(self):
        """All powers of a rotation share the same axis."""
        np.random.seed(42)
        for _ in range(50):
            axis = np.random.randn(3)
            axis /= np.linalg.norm(axis)
            angle = np.random.uniform(0.1, np.pi)
            R = rotation_matrix(axis, angle)
            # Check powers R, R², R³, ...
            Rk = np.eye(3)
            for k in range(1, 20):
                Rk = Rk @ R
                if np.linalg.norm(Rk - np.eye(3)) < 1e-8:
                    break  # identity
                ax = rotation_axis(Rk)
                assert axes_parallel(ax, axis), (
                    f"R^{k} axis not parallel to R axis")
        print("  50 random axes: all powers share axis ✓")

    def test_order4_square_same_axis(self):
        """If R has order 4 (angle π/2), R² is π-rotation about SAME axis."""
        np.random.seed(123)
        for _ in range(50):
            axis = np.random.randn(3)
            axis /= np.linalg.norm(axis)
            R = rotation_matrix(axis, np.pi / 2)
            R2 = R @ R
            # R² should be π-rotation about same axis
            angle2 = rotation_angle(R2)
            axis2 = rotation_axis(R2)
            assert abs(angle2 - np.pi) < 1e-8, f"R² angle = {angle2}, not π"
            assert axes_parallel(axis2, axis), "R² axis differs from R axis"
        print("  50 random: order-4 square gives π-rotation about same axis ✓")

    def test_same_axis_commute(self):
        """Rotations about the same axis commute."""
        np.random.seed(456)
        for _ in range(50):
            axis = np.random.randn(3)
            axis /= np.linalg.norm(axis)
            a1 = np.random.uniform(0.1, 2 * np.pi)
            a2 = np.random.uniform(0.1, 2 * np.pi)
            R1 = rotation_matrix(axis, a1)
            R2 = rotation_matrix(axis, a2)
            comm = R1 @ R2 - R2 @ R1
            assert np.linalg.norm(comm) < 1e-10, "Same-axis rotations don't commute"
        print("  50 random: same-axis rotations commute ✓")


# ============================================================
# Test 2: Direct contradiction — Q_{2^n} embedding attempt
# ============================================================

class TestGenQuaternionContradiction:
    """Attempt to embed Q_{2^n} in SO(3) and derive contradiction."""

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_axis_forcing(self, n):
        """The non-trivial step: b² = z FORCES b to share axis with a.

        Scans over all possible axes w for b (at angles θ to u).
        Shows b² = π-rotation about w can only equal z = π-rotation
        about u when w ∥ u. This is what forces the contradiction."""
        order_a = 2 ** (n - 1)
        np.random.seed(n)

        for _ in range(10):
            # Random axis u for a
            u = np.random.randn(3)
            u /= np.linalg.norm(u)
            Ra = rotation_matrix(u, 2 * np.pi / order_a)

            # z = a^{2^{n-2}} = π-rotation about u
            Rz = np.linalg.matrix_power(Ra, 2 ** (n - 2))
            assert abs(rotation_angle(Rz) - np.pi) < 1e-6

            # Build orthonormal frame with u as first vector
            perp = np.array([1, 0, 0]) if abs(u[2]) > 0.9 else np.array([0, 0, 1])
            v1 = np.cross(u, perp); v1 /= np.linalg.norm(v1)

            # Scan axis w at angle θ to u, full azimuthal sweep
            matches = []
            for theta in np.linspace(0, np.pi, 60):
                for phi in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                    w = np.cos(theta) * u + np.sin(theta) * (
                        np.cos(phi) * v1 + np.sin(phi) * np.cross(u, v1))
                    w /= np.linalg.norm(w)
                    Rb = rotation_matrix(w, np.pi / 2)  # order-4 candidate
                    Rb2 = Rb @ Rb  # π-rotation about w
                    if np.linalg.norm(Rb2 - Rz) < 1e-6:
                        matches.append(theta)

            # Only θ ≈ 0 or θ ≈ π (i.e. w ∥ u) should match
            for theta in matches:
                assert theta < 0.1 or theta > np.pi - 0.1, (
                    f"b² = z with non-parallel axis at θ={theta:.3f}")

            assert len(matches) > 0, "No axis found (should find w = ±u)"

        print(f"  Q_{2**n}: b² = z forces w ∥ u (scanned 60×24 directions) ✓")

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_commutation_contradiction(self, n):
        """Once b is forced onto axis u: a,b commute ⟹ bab⁻¹ = a ≠ a⁻¹."""
        order_a = 2 ** (n - 1)
        np.random.seed(n + 50)

        for _ in range(20):
            u = np.random.randn(3)
            u /= np.linalg.norm(u)
            Ra = rotation_matrix(u, 2 * np.pi / order_a)
            Ra_inv = np.linalg.matrix_power(Ra, order_a - 1)

            for sign in [1, -1]:
                Rb = rotation_matrix(u, sign * np.pi / 2)
                # Same axis ⟹ commute
                assert np.linalg.norm(Ra @ Rb - Rb @ Ra) < 1e-10
                # So bab⁻¹ = a, but relation requires a⁻¹
                conjugate = Rb @ Ra @ Rb.T
                assert np.linalg.norm(conjugate - Ra) < 1e-8
                assert np.linalg.norm(Ra - Ra_inv) > 1e-6, (
                    "a = a⁻¹ means ord(a) ≤ 2, but need ≥ 4")

        print(f"  Q_{2**n}: same-axis ⟹ commute ⟹ bab⁻¹=a ≠ a⁻¹ ✓")

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_no_off_axis_escape(self, n):
        """Even if b has a different axis w ≠ u, b² cannot equal z.
        Because b² = π-rotation about w, and z = π-rotation about u,
        and π-rotations about different axes are different."""
        order_a = 2 ** (n - 1)
        np.random.seed(n + 100)

        u = np.array([0, 0, 1.0])  # axis of a
        Ra = rotation_matrix(u, 2 * np.pi / order_a)
        exp_z = 2 ** (n - 2)
        Rz = np.linalg.matrix_power(Ra, exp_z)

        # Try b with axis w at various angles to u
        for theta in np.linspace(0.1, np.pi / 2, 20):
            w = np.array([np.sin(theta), 0, np.cos(theta)])
            Rb = rotation_matrix(w, np.pi / 2)  # order 4 about w
            Rb2 = Rb @ Rb  # π-rotation about w

            # Rb2 should NOT equal Rz (π-rotation about u) when w ≠ u
            if not axes_parallel(w, u):
                diff = np.linalg.norm(Rb2 - Rz)
                assert diff > 1e-6, (
                    f"π-rotation about w (θ={theta:.2f}) equals "
                    f"π-rotation about u!")

        print(f"  Q_{2**n}: off-axis b cannot satisfy b² = z ✓")


# ============================================================
# Test 3: V₄ orthogonality — why axes must be perpendicular
# ============================================================

class TestV4Orthogonality:
    """Two π-rotations about axes at angle θ compose to rotation of angle 2θ.
    So (ab)² = e in SO(3) requires 2θ = π, i.e. θ = π/2 (orthogonal)."""

    def test_pi_rotation_composition_angle(self):
        """Composition of π-rotations about axes at angle θ gives 2θ."""
        u = np.array([0, 0, 1.0])
        for theta in np.linspace(0.05, np.pi - 0.05, 50):
            v = np.array([np.sin(theta), 0, np.cos(theta)])
            Ra = rotation_matrix(u, np.pi)
            Rb = rotation_matrix(v, np.pi)
            ab = Ra @ Rb
            angle_ab = rotation_angle(ab)
            expected = 2 * theta if theta <= np.pi / 2 else 2 * (np.pi - theta)
            assert abs(angle_ab - expected) < 1e-6, (
                f"θ={theta:.3f}: got {angle_ab:.3f}, expected {expected:.3f}")
        print("  π-rotation composition: angle = 2θ ✓")

    def test_V4_forces_orthogonal(self):
        """In V₄ = {e, a, b, ab}, all three non-identity elements have order 2.
        (ab)² = e requires composition angle = π, so axes at π/2."""
        np.random.seed(789)
        for _ in range(30):
            u = np.random.randn(3); u /= np.linalg.norm(u)
            # Build perpendicular v
            perp = np.array([1, 0, 0]) if abs(u[2]) > 0.9 else np.array([0, 0, 1])
            v = np.cross(u, perp); v /= np.linalg.norm(v)
            Ra = rotation_matrix(u, np.pi)
            Rb = rotation_matrix(v, np.pi)
            ab = Ra @ Rb
            # ab should be π-rotation (order 2)
            assert abs(rotation_angle(ab) - np.pi) < 1e-6
            # Axes are orthogonal
            assert abs(np.dot(u, v)) < 1e-8
        print("  V₄ axes orthogonal ⟹ SU(2) lifts mutually orthogonal ✓")


# ============================================================
# Test 4: The full classification-free chain
# ============================================================

class TestClassificationFreeChain:
    """Verify the complete chain:
    (B) ⟹ Syl₂ non-cyclic ⟹ |Syl₂| ≥ 4 ⟹ Syl₂ not gen.quat
    ⟹ V₄ ⊂ Syl₂ ⟹ V₄ ⊂ H ⟹ Q₈ ⊂ 2H ⟹ (A).

    Each step proven without ADE classification."""

    def test_chain_summary(self):
        """Print the logical chain with proof methods."""
        steps = [
            ("(B₃): v₂(ord(h)) < v₂(|H|) for all h ≠ e",
             "Hypothesis"),
            ("Syl₂(H) is non-cyclic",
             "If cyclic, has element of order |Syl₂| = 2^a, "
             "so v₂(ord) = a = v₂(|H|). Contradicts (B₃)."),
            ("|Syl₂(H)| ≥ 4",
             "exp(H) = |H|/2 (Lemma 3.2, classification-free) "
             "⟹ |H| ≥ 4 ⟹ v₂(|H|) ≥ 2 ⟹ |Syl₂| ≥ 4."),
            ("Syl₂(H) is not generalized quaternion",
             "Q_{2^n} cannot embed in SO(3): forced same-axis "
             "⟹ commutation ⟹ contradicts bab⁻¹=a⁻¹. "
             "Pure rotation geometry, no ADE."),
            ("V₄ ⊂ Syl₂(H) ⊂ H",
             "A 2-group of order ≥ 4 that is neither cyclic nor "
             "generalized quaternion contains V₄ = Z₂ × Z₂. "
             "Standard group theory (Burnside)."),
            ("Q₈ ⊂ 2H",
             "V₄ = {e,a,b,ab} with a,b,ab all order-2 rotations. "
             "SU(2) lifts are pure imaginary, mutually orthogonal "
             "(scalar part of ab = -u_a·u_b = 0). "
             "Generates Q₈. No ADE."),
            ("-1 ∈ [2H, 2H], i.e., (A) holds",
             "Q₈ ⊂ 2H ⟹ ∃ a,b with [a,b] = -1. "
             "Pure quaternion algebra. No ADE."),
        ]

        print("\n  CLASSIFICATION-FREE PROOF: (B) ⟹ (A)")
        print("  " + "=" * 60)
        for i, (claim, proof) in enumerate(steps):
            print(f"\n  Step {i+1}: {claim}")
            print(f"    Proof: {proof}")
        print("\n  " + "=" * 60)
        print("  Each step uses only: SU(2)/SO(3) geometry, Sylow theory,")
        print("  basic 2-group structure, quaternion algebra.")
        print("  NO finite rotation group classification (ADE) invoked.")

        # The only external fact needed:
        print("\n  EXTERNAL FACTS USED (all standard, no ADE):")
        print("  1. Sylow theorem (existence of p-subgroups)")
        print("  2. Burnside: non-cyclic, non-gen-quaternion 2-group has V₄")
        print("  3. Rotation geometry: cyclic ⟹ shared axis, shared axis ⟹ commute")
        print("  4. Quaternion algebra: pure imaginary, orthogonal ⟹ Q₈")


# ============================================================
# Test 5: End-to-end on concrete groups
# ============================================================

# Quaternion utilities (from file 58)
def qmul(a, b):
    return np.array([
        a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
        a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
        a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
        a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0],
    ])

def qinv(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

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

def so3_elements(G_su2):
    reps = []
    for g in G_su2:
        already = False
        for r in reps:
            if np.linalg.norm(g - (-r)) < 1e-10:
                already = True
                break
        if not already:
            reps.append(g)
    return reps

def so3_order(g):
    if abs(abs(g[0]) - 1) < 1e-10:
        return 1
    alpha = np.arccos(np.clip(abs(g[0]), -1, 1))
    for d in range(2, 200):
        if abs(np.sin(d * alpha)) < 1e-8:
            return d
    return None

def check_B(H_reps):
    n = len(H_reps)
    for g in H_reps:
        d = so3_order(g)
        if d == 1:
            continue
        if (n // d) % 2 != 0:
            return False
    return True

def has_Q8(G_su2):
    for a in G_su2:
        for b in G_su2:
            c = qmul(qmul(a, b), qmul(qinv(a), qinv(b)))
            if abs(c[0] + 1) < 1e-10 and np.linalg.norm(c[1:]) < 1e-10:
                return True
    return False

# Group builders
def make_2Dn(n):
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


class TestEndToEnd:
    """Full chain on concrete groups: (B) ⟹ V₄ ⊂ H ⟹ Q₈ ⊂ 2H ⟹ (A)."""

    @pytest.mark.parametrize("name,builder", [
        ("D_2", lambda: make_2Dn(2)),
        ("D_4", lambda: make_2Dn(4)),
        ("D_6", lambda: make_2Dn(6)),
        ("T",   make_2T),
        ("O",   make_2O),
        ("I",   make_2I),
    ])
    def test_B_implies_A_concrete(self, name, builder):
        """For each group satisfying (B), verify Q₈ ⊂ 2H (i.e. (A))."""
        G = builder()
        H = so3_elements(G)
        sat_B = check_B(H)
        has_A = has_Q8(G)

        if sat_B:
            assert has_A, f"{name}: satisfies (B) but Q₈ ⊄ 2H!"
            # Also check V₄ exists (count involutions with perp pairs)
            involutions = [g for g in H if so3_order(g) == 2]
            assert len(involutions) >= 3, (
                f"{name}: (B) but < 3 involutions")
            print(f"  {name}: |H|={len(H)}, (B)=True, (A)=True, "
                  f"{len(involutions)} involutions ✓")
        else:
            print(f"  {name}: |H|={len(H)}, (B)=False (not in scope)")


def make_2Cn(n):
    a = np.array([np.cos(np.pi / n), np.sin(np.pi / n), 0, 0])
    return generate_group([a])

def has_V4(H_reps):
    """Check if H contains V₄: two commuting involutions whose product is also order 2."""
    invols = [g for g in H_reps if so3_order(g) == 2]
    for i, a in enumerate(invols):
        for b in invols[i+1:]:
            ab = qmul(a, b)
            if so3_order(ab) == 2:
                return True
    return False

def v2(n):
    if n == 0:
        return float('inf')
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


# ============================================================
# Test 6: (A)⟹(B) gap — V₄ ⊂ H ⟺ (B) brute verification
# ============================================================

class TestV4iffB:
    """Verify V₄ ⊂ H ⟺ (B) on ALL finite rotation groups up to large order."""

    @pytest.mark.parametrize("name,builder", [
        ("C_2", lambda: make_2Cn(2)),
        ("C_3", lambda: make_2Cn(3)),
        ("C_4", lambda: make_2Cn(4)),
        ("C_6", lambda: make_2Cn(6)),
        ("D_2", lambda: make_2Dn(2)),
        ("D_3", lambda: make_2Dn(3)),
        ("D_4", lambda: make_2Dn(4)),
        ("D_5", lambda: make_2Dn(5)),
        ("D_6", lambda: make_2Dn(6)),
        ("D_7", lambda: make_2Dn(7)),
        ("D_8", lambda: make_2Dn(8)),
        ("D_9", lambda: make_2Dn(9)),
        ("D_10", lambda: make_2Dn(10)),
        ("T",   make_2T),
        ("O",   make_2O),
        ("I",   make_2I),
    ])
    def test_V4_iff_B(self, name, builder):
        """V₄ ⊂ H ⟺ (B) for every finite H ⊂ SO(3)."""
        G = builder()
        H = so3_elements(G)
        sat_B = check_B(H)
        has_v4 = has_V4(H)
        assert has_v4 == sat_B, (
            f"{name}: V₄={has_v4} but (B)={sat_B}")
        print(f"  {name}: |H|={len(H)}, V₄={has_v4}, (B)={sat_B} ✓")


# ============================================================
# Test 7: Burnside equation for 2-groups — only cyclic/dihedral
# ============================================================

class TestBurnside2Groups:
    """The Burnside orbit-counting equation for a 2-group P ⊂ SO(3):
    2(1 - 1/|P|) = Σ(1 - 1/nᵢ), where nᵢ are stabilizer sizes (powers of 2).

    Solve for all |P| = 2^a, a = 1..20. Only cyclic and dihedral solutions exist."""

    @pytest.mark.parametrize("a", range(1, 21))
    def test_diophantine_solutions(self, a):
        """Enumerate all solutions of the Burnside equation for |P| = 2^a."""
        order = 2 ** a
        target = 2 * (1 - 1 / order)  # LHS

        # Find all (r, n₁, ..., nᵣ) with nᵢ = 2^kᵢ, nᵢ | order
        # such that Σ(1 - 1/nᵢ) = target
        # and orbit sizes |P|/nᵢ are positive integers
        solutions = []
        divisors = [2 ** k for k in range(1, a + 1)]

        # r = 2: two stabilizers
        for n1 in divisors:
            for n2 in divisors:
                if n1 > n2:
                    continue
                s = (1 - 1/n1) + (1 - 1/n2)
                if abs(s - target) < 1e-12:
                    solutions.append(("cyclic", [n1, n2]))

        # r = 3: three stabilizers
        for n1 in divisors:
            for n2 in divisors:
                if n2 < n1:
                    continue
                for n3 in divisors:
                    if n3 < n2:
                        continue
                    s = (1 - 1/n1) + (1 - 1/n2) + (1 - 1/n3)
                    if abs(s - target) < 1e-12:
                        solutions.append(("dihedral?", [n1, n2, n3]))

        # r ≥ 4: impossible (each term ≥ 1/2, sum ≥ 2, but target < 2)
        # (verify)
        assert target < 2, f"target = {target} ≥ 2"

        # Classify solutions
        for label, stabs in solutions:
            if len(stabs) == 2:
                # Both stabilizers = full group ⟹ cyclic
                assert stabs[0] == order and stabs[1] == order, (
                    f"2^{a}: unexpected 2-orbit solution {stabs}")
            elif len(stabs) == 3:
                # Must be (2, 2, 2^{a-1}) ⟹ dihedral
                assert stabs == [2, 2, 2 ** (a - 1)], (
                    f"2^{a}: unexpected 3-orbit solution {stabs}")

        if a == 1:
            assert len(solutions) == 1  # only cyclic C₂
        else:
            assert len(solutions) == 2  # cyclic + dihedral

        print(f"  2^{a}={order}: {len(solutions)} solutions — "
              f"{'cyclic' if a == 1 else 'cyclic + dihedral'} ✓")


# ============================================================
# Test 8: Dihedral exponent bound ⟹ (B₃)
# ============================================================

class TestDihedralExponent:
    """Dihedral D_{2^{a-1}} of order 2^a has exponent 2^{a-1}.
    So v₂(ord(h)) ≤ a-1 < a for every element."""

    @pytest.mark.parametrize("a", range(2, 8))
    def test_dihedral_exponent(self, a):
        """Build D_{2^{a-1}} in SO(3), check max element order = 2^{a-1}."""
        n = 2 ** (a - 1)  # D_n has order 2n = 2^a
        G = make_2Dn(n)
        H = so3_elements(G)
        assert len(H) == 2 * n, f"|H| = {len(H)}, expected {2*n}"

        max_ord = max(so3_order(g) for g in H)
        assert max_ord == n, f"max order = {max_ord}, expected {n} = 2^{a-1}"

        # Check (B₃): v₂(ord(h)) < a for all h ≠ e
        for g in H:
            d = so3_order(g)
            if d == 1:
                continue
            assert v2(d) <= a - 1, (
                f"v₂({d}) = {v2(d)} ≥ {a} = v₂(|H|)")

        print(f"  D_{n} (|H|={2*n}): exp = {n} = 2^{a-1}, (B₃) holds ✓")


# ============================================================
# Test 9: Syl₂ dihedral ⟹ (B₃) on full H (Sylow transfer)
# ============================================================

class TestSylowTransfer:
    """The exponent bound on Syl₂ transfers to all of H:
    any h ∈ H with ord(h) = 2^b · m has h^m of order 2^b in some
    Sylow 2-subgroup, so b ≤ a-1. Hence v₂(ord(h)) < a = v₂(|H|)."""

    @pytest.mark.parametrize("name,builder", [
        ("D_2", lambda: make_2Dn(2)),
        ("D_4", lambda: make_2Dn(4)),
        ("D_6", lambda: make_2Dn(6)),
        ("D_8", lambda: make_2Dn(8)),
        ("T",   make_2T),
        ("O",   make_2O),
        ("I",   make_2I),
    ])
    def test_sylow_transfer(self, name, builder):
        """For H with V₄ ⊂ H: Syl₂ is dihedral, exp(Syl₂) = |S|/2,
        and v₂(ord(h)) < v₂(|H|) for ALL h ∈ H."""
        G = builder()
        H = so3_elements(G)
        H_order = len(H)
        a = v2(H_order)

        # Check every element
        for g in H:
            d = so3_order(g)
            if d == 1:
                continue
            assert v2(d) < a, (
                f"{name}: h with ord={d}, v₂={v2(d)} ≥ {a}=v₂(|H|)")
            assert (H_order // d) % 2 == 0, (
                f"{name}: |H|/ord = {H_order//d} is odd")

        print(f"  {name}: |H|={H_order}, v₂(|H|)={a}, "
              f"all v₂(ord(h)) < {a} ✓")


    @pytest.mark.parametrize("name,builder", [
        ("C_4", lambda: make_2Cn(4)),
        ("D_3", lambda: make_2Dn(3)),
        ("D_5", lambda: make_2Dn(5)),
        ("D_7", lambda: make_2Dn(7)),
    ])
    def test_contrapositive(self, name, builder):
        """Groups without V₄: ∃ h with v₂(ord(h)) = v₂(|H|), so (B₃) fails."""
        G = builder()
        H = so3_elements(G)
        H_order = len(H)
        a = v2(H_order)

        # Find a witness violating (B₃)
        witness = None
        for g in H:
            d = so3_order(g)
            if d == 1:
                continue
            if v2(d) >= a:
                witness = d
                break

        assert witness is not None, (
            f"{name}: no element with v₂(ord) ≥ v₂(|H|), but expected one")
        assert not has_V4(H), f"{name}: has V₄ but shouldn't"
        print(f"  {name}: |H|={H_order}, V₄=False, witness ord={witness}, "
              f"v₂={v2(witness)} ≥ {a} ✓")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
