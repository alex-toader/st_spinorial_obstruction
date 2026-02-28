"""
Test Suite: Binary Octahedral Group — No 1D Spinorial Character
================================================================

Last run: 44/44 PASSED, 0.39s (Feb 2026, Python 3.14.3, macOS ARM64)

Core claim verified here:
  -1 ∈ [2O, 2O]
  ⇒ every 1D character χ of 2O satisfies χ(-1) = +1
  ⇒ no scalar (1D) quantization on SO(3)/O can produce spinors
  ⇒ spinorial sector requires multi-component wavefunction (dim ≥ 2)

By "scalar quantization" we mean one-dimensional unitary representations
of π₁(SO(3)/O) = 2O acting as phase boundary conditions (Laidlaw-DeWitt
1971, Schulman 1968). SU(2) is not postulated as an internal symmetry;
it appears as the universal cover of SO(3) in this quantization scheme.

What these tests DO verify:
  - 2O construction (48 unit quaternions, no duplicates, group closure)
  - |[2O,2O]| = 24, |2O/[2O,2O]| = 2
  - -1 ∈ [2O,2O]  (THE KEY RESULT)
  - Character table cross-check: χ(-1) = +dim for non-spin, -dim for spin irreps
  - H₁ and H₂ are chirality-conjugate: H₁ × A₁ = H₂ (A₁ = sign rep of O ≅ S₄)
  - Eg(O) lifts to non-spinorial E(2O), not to H₁ or H₂
  - Full character table orthogonality (all 8 irreps, 8 classes, from GAP/Groupprops)
  - Bilinear analysis: ∧²(H₁) = A₀ (symplectic), Sym²(H₁) = T_vec (vector rep)
  - Explicit SU(2) matrices: fundamental rep of 2O ⊂ SU(2) = H₁
  - E vs H₁ boundary condition: E(-1) = +I₂, H₁(-1) = -I₂
  - Contrast: Z₂ is abelian, so -1 ∉ [Z₂,Z₂] (1D spinorial exists)
  - Same property for 2T: -1 ∈ [2T,2T] (not unique to 2O)

What these tests do NOT verify:
  - Which physical DOF carries the spinorial representation
  - Whether the shear polarization space is the carrier
  - Why nature selects the spinorial sector
  These are physical identifications beyond group theory.

References:
- Schulman (1968): Phys. Rev. 176, 1558
- Laidlaw & DeWitt (1971): Phys. Rev. D 3, 1375
- Landsman (2013): arXiv:1302.3637
- Groupprops: binary octahedral group character table
"""

import numpy as np
import pytest

from src.quaternion import qmul, qinv, qkey, quat_to_so3, quat_to_su2
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
                       verify_closure, compute_commutator_subgroup,
                       MINUS_ONE_KEY)


# ============================================================
# TESTS
# ============================================================

class TestBinaryOctahedralConstruction:
    """Verify 2O is correctly constructed as a group of 48 unit quaternions."""

    def test_order_48(self):
        """2O has exactly 48 elements."""
        elements = build_binary_octahedral()
        assert len(elements) == 48

    def test_all_unit_quaternions(self):
        """All elements are unit quaternions."""
        elements = build_binary_octahedral()
        norms = np.linalg.norm(elements, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-12)

    def test_contains_minus_one(self):
        """2O contains the central element -1."""
        elements = build_binary_octahedral()
        keys = {qkey(e) for e in elements}
        assert MINUS_ONE_KEY in keys

    def test_no_duplicates(self):
        """All 48 elements are distinct (canonical keys are unique)."""
        elements = build_binary_octahedral()
        keys = {qkey(e) for e in elements}
        assert len(keys) == 48, (
            f"Expected 48 distinct elements, got {len(keys)}"
        )

    def test_elements_well_separated(self):
        """Minimum pairwise distance >> rounding tolerance (ROUND_DIGITS = 12).

        All 48 elements lie on S³. This test proves that rounding at 10⁻¹²
        cannot conflate distinct elements, validating qkey()-based identity
        checks used throughout the test suite.
        """
        elements = build_binary_octahedral()
        min_dist = float('inf')
        for i in range(len(elements)):
            for j in range(i + 1, len(elements)):
                d = np.linalg.norm(elements[i] - elements[j])
                if d < min_dist:
                    min_dist = d
        # Actual minimum ≈ 0.765 (adjacent vertices on S³)
        assert min_dist > 0.5, (
            f"Min pairwise distance {min_dist:.6f} too small for safe rounding"
        )

    def test_element_orders(self):
        """Element orders in 2O match known distribution.

        In SU(2), the order of q is the smallest n with q^n = +1.
        Expected distribution for 2O: {1:1, 2:1, 3:8, 4:12, 6:8, 8:18}.
        This identifies conjugacy class structure independently of trace.
        """
        elements = build_binary_octahedral()
        identity = np.array([1, 0, 0, 0], dtype=float)
        id_key = qkey(identity)
        order_counts = {}
        for e in elements:
            q = e.copy()
            for n in range(1, 25):
                if qkey(q) == id_key:
                    order_counts[n] = order_counts.get(n, 0) + 1
                    break
                q = qmul(q, e)
        expected = {1: 1, 2: 1, 3: 8, 4: 18, 6: 8, 8: 12}
        assert order_counts == expected, (
            f"Order distribution {order_counts} != expected {expected}"
        )


class TestGroupClosure:
    """Verify constructed sets are groups under quaternion multiplication."""

    def test_2O_closure(self):
        """2O is closed under quaternion multiplication."""
        elements = build_binary_octahedral()
        assert verify_closure(elements, "2O")

    def test_2T_closure(self):
        """2T is closed under quaternion multiplication."""
        elements = build_binary_tetrahedral()
        assert verify_closure(elements, "2T")


class TestNoOneDimSpinorialCharacter:
    """THE KEY RESULT: no 1D spinorial character exists for 2O.

    Method: compute [2O,2O], show -1 ∈ [2O,2O].
    Since any 1D rep χ sends commutators to 1
    (χ(ghg⁻¹h⁻¹) = χ(g)χ(h)χ(g)⁻¹χ(h)⁻¹ = 1),
    and -1 is in [2O,2O], we get χ(-1) = 1 for all 1D reps.
    """

    def test_abelianization_order_2(self):
        """2O has exactly 2 one-dimensional representations.
        |2O/[2O,2O]| = 2, |[2O,2O]| = 24.
        """
        elements = build_binary_octahedral()
        comm = compute_commutator_subgroup(elements)
        comm_order = len(comm)
        ab_order = len(elements) // comm_order

        assert comm_order == 24, (
            f"Expected |[2O,2O]| = 24, got {comm_order}"
        )
        assert ab_order == 2, (
            f"Expected |2O/[2O,2O]| = 2, got {ab_order}"
        )

    def test_minus_one_in_commutator(self):
        """-1 ∈ [2O,2O] ⇒ all 1D reps have χ(-1) = +1.
        This is why no scalar quantization on SO(3)/O produces spinors.
        """
        elements = build_binary_octahedral()
        comm = compute_commutator_subgroup(elements)

        assert MINUS_ONE_KEY in comm, (
            "-1 NOT in [2O,2O] — this would mean a 1D spinorial "
            "character exists, contradicting the known result!"
        )

    def test_explicit_commutator_witness(self):
        """Direct witness: [i, j] = ij(-i)(-j) = -1.

        g = i, h = j are quaternion units in Q₈ ⊂ 2O.
        ghg⁻¹h⁻¹ = ij·(-i)·(-j) = k·(-i)·(-j) = j·(-j) = -1.
        This provides -1 ∈ [2O,2O] without closure computation.
        """
        g = np.array([0, 1, 0, 0], dtype=float)  # i
        h = np.array([0, 0, 1, 0], dtype=float)  # j
        comm = qmul(
            qmul(g, h),
            qmul(qinv(g), qinv(h))
        )
        np.testing.assert_allclose(
            comm, [-1, 0, 0, 0], atol=1e-12,
            err_msg="[i,j] should equal -1"
        )


class TestCommutatorSubgroupIs2T:
    """[2O, 2O] = 2T as a concrete subgroup (not just order 24).

    Verifies set equality: the commutator subgroup of 2O, computed
    from scratch via quaternion arithmetic, is exactly the set of
    24 elements of 2T constructed independently.
    """

    def test_commutator_equals_2T(self):
        """[2O, 2O] = 2T (set equality as quaternion subgroups of SU(2))."""
        elements_2O = build_binary_octahedral()
        elements_2T = build_binary_tetrahedral()

        comm_keys = set(compute_commutator_subgroup(elements_2O).keys())
        t_keys = {qkey(e) for e in elements_2T}

        assert comm_keys == t_keys, (
            f"[2O,2O] ≠ 2T. "
            f"In [2O,2O] but not 2T: {comm_keys - t_keys}. "
            f"In 2T but not [2O,2O]: {t_keys - comm_keys}."
        )


class TestPresentation:
    """2O identified at the presentation level, not just by enumeration.

    Standard presentation: ⟨r, s | r⁴ = s³ = (rs)² = -1⟩
    where -1 is the unique central element of order 2.

    Generators: r = (1+i)/√2 (order 8), s = (1+i+j+k)/2 (order 6).
    These generate all 48 elements of 2O.
    """

    def test_presentation_relations(self):
        """r⁴ = s³ = (rs)² = -1 for explicit generators of 2O."""
        sq2 = 1.0 / np.sqrt(2)
        r = np.array([sq2, sq2, 0, 0])       # (1+i)/√2
        s = np.array([0.5, 0.5, 0.5, 0.5])   # (1+i+j+k)/2

        def qpow(q, n):
            result = np.array([1, 0, 0, 0], dtype=float)
            for _ in range(n):
                result = qmul(result, q)
            return result

        minus_one = np.array([-1, 0, 0, 0], dtype=float)

        r4 = qpow(r, 4)
        s3 = qpow(s, 3)
        rs = qmul(r, s)
        rs2 = qpow(rs, 2)

        np.testing.assert_allclose(r4, minus_one, atol=1e-12,
                                   err_msg="r⁴ ≠ -1")
        np.testing.assert_allclose(s3, minus_one, atol=1e-12,
                                   err_msg="s³ ≠ -1")
        np.testing.assert_allclose(rs2, minus_one, atol=1e-12,
                                   err_msg="(rs)² ≠ -1")

    def test_generators_produce_full_group(self):
        """r = (1+i)/√2 and s = (1+i+j+k)/2 generate all 48 elements of 2O."""
        sq2 = 1.0 / np.sqrt(2)
        r = np.array([sq2, sq2, 0, 0])
        s = np.array([0.5, 0.5, 0.5, 0.5])

        generated = {qkey(np.array([1, 0, 0, 0], dtype=float)):
                     np.array([1, 0, 0, 0], dtype=float)}

        def add(q):
            k = qkey(q)
            if k not in generated:
                generated[k] = q.copy()
                return True
            return False

        add(r); add(s)
        changed = True
        while changed:
            changed = False
            for a in list(generated.values()):
                for b in list(generated.values()):
                    if add(qmul(a, b)):
                        changed = True

        assert len(generated) == 48, (
            f"Generators produce {len(generated)} elements, expected 48"
        )

        # Verify it's the same set as the enumerated 2O
        enum_keys = {qkey(e) for e in build_binary_octahedral()}
        assert set(generated.keys()) == enum_keys, (
            "Generated group ≠ enumerated 2O"
        )


class TestCoveringMap:
    """2O is the preimage of O under the universal cover SU(2) → SO(3).

    This connects the abstract group theory to the topology of SO(3)/O:
    the covering map q → R(q) sends 48 quaternions to exactly 24 distinct
    rotations, with kernel {+1, -1} at every fiber.
    """

    @staticmethod
    def _so3_key(R):
        return tuple(round(float(x), 10) for x in R.flatten())

    def test_48_to_24_rotations(self):
        """48 quaternions in 2O map to exactly 24 distinct SO(3) rotations."""
        elements = build_binary_octahedral()
        so3_keys = {self._so3_key(quat_to_so3(q)) for q in elements}
        assert len(so3_keys) == 24, (
            f"Expected 24 distinct rotations, got {len(so3_keys)}"
        )

    def test_kernel_is_pm1(self):
        """Every rotation has exactly 2 preimages: +q and -q."""
        elements = build_binary_octahedral()
        from collections import Counter
        fiber_sizes = Counter()
        for q in elements:
            fiber_sizes[self._so3_key(quat_to_so3(q))] += 1
        sizes = set(fiber_sizes.values())
        assert sizes == {2}, (
            f"Expected all fibers of size 2, got sizes {sizes}"
        )

    def test_preimages_are_pm_q(self):
        """The two preimages of each rotation are exactly +q and -q."""
        elements = build_binary_octahedral()
        fibers = {}
        for q in elements:
            k = self._so3_key(quat_to_so3(q))
            fibers.setdefault(k, []).append(q)
        for k, pair in fibers.items():
            assert len(pair) == 2
            np.testing.assert_allclose(
                pair[0] + pair[1], [0, 0, 0, 0], atol=1e-12,
                err_msg="Preimage pair should be +q and -q"
            )

    def test_so3_matrices_orthogonal(self):
        """All SO(3) images are proper rotations: R^T R = I, det R = +1."""
        elements = build_binary_octahedral()
        for q in elements:
            R = quat_to_so3(q)
            np.testing.assert_allclose(R.T @ R, np.eye(3), atol=1e-12)
            np.testing.assert_allclose(np.linalg.det(R), 1.0, atol=1e-12)


class TestCharacterTableCrossCheck:
    """Structural checks on the published character table of 2O.

    Source: Groupprops / GAP SmallGroup(48,28).
    These tests verify structural properties of the hardcoded table
    (sum of squares, spin/non-spin split, no 1D spinorial irrep).
    They do not derive the table from first principles; they guard
    against accidental edits and confirm consistency with the
    commutator-subgroup result (-1 ∈ [2O,2O] ⇒ no 1D spinorial).
    """

    def test_character_at_minus_one(self):
        """Verify χ(-1) = ±dim from character table (Groupprops/GAP).

        For -1 ∈ center of 2O, Schur's lemma gives ρ(-1) = λ·I
        with λ² = 1, so λ = ±1, hence χ(-1) = ±dim.

        Non-spin irreps (factor through O):
          A₀(1): χ(-1) = +1, A₁(1): χ(-1) = +1
          E(2): χ(-1) = +2, T₁(3): χ(-1) = +3, T₂(3): χ(-1) = +3

        Spin irreps (genuine 2O, do NOT factor through O):
          H₁(2): χ(-1) = -2, H₂(2): χ(-1) = -2, L(4): χ(-1) = -4
        """
        # Character values at -1, from Groupprops/GAP (external source)
        # Format: (name, dim, chi_at_minus_one)
        character_table_at_minus_one = [
            ('A0', 1,  +1),
            ('A1', 1,  +1),
            ('E',  2,  +2),
            ('T1', 3,  +3),
            ('T2', 3,  +3),
            ('H1', 2,  -2),
            ('H2', 2,  -2),
            ('L',  4,  -4),
        ]

        # Verify: sum of squares = |2O|
        dims = [d for _, d, _ in character_table_at_minus_one]
        assert sum(d**2 for d in dims) == 48, "Sum of squares check"

        # Verify: χ(-1) = ±dim (Schur's lemma for central element)
        for name, dim, chi in character_table_at_minus_one:
            assert abs(chi) == dim, (
                f"{name}: |χ(-1)| = {abs(chi)} ≠ dim = {dim}"
            )

        # Extract spin vs non-spin
        non_spin = [(n, d) for n, d, c in character_table_at_minus_one if c > 0]
        spin = [(n, d) for n, d, c in character_table_at_minus_one if c < 0]

        # Non-spin: all 1D irreps are here (consistent with no 1D spinorial)
        one_dim_irreps = [n for n, d in non_spin if d == 1]
        assert len(one_dim_irreps) == 2, (
            f"Expected 2 one-dim non-spin irreps, got {one_dim_irreps}"
        )

        # Spin: no 1D irreps exist
        one_dim_spin = [n for n, d in spin if d == 1]
        assert len(one_dim_spin) == 0, (
            f"Found 1D spin irreps: {one_dim_spin} — should not exist!"
        )

        # Min spinorial dimension = 2
        spin_dims = sorted([d for _, d in spin])
        assert spin_dims == [2, 2, 4], (
            f"Expected spin dims [2,2,4], got {spin_dims}"
        )


class TestH1vsH2ChiralityConjugate:
    """H₁ and H₂ are chirality-conjugate: H₁ ⊗ A₁ = H₂.

    A₁(2O) is the lift of the sign representation of O ≅ S₄,
    which distinguishes even permutations (T ≅ A₄, sign = +1)
    from odd permutations (O\\T, sign = -1).

    Consequence: The choice between H₁ and H₂ is a chirality convention,
    analogous to (½,0) vs (0,½) in Lorentz group representations.
    No classical observable distinguishes them.

    Full character table source: GAP SmallGroup(48,28) / Groupprops.
    These tests are regression fixtures for the published table: they verify
    internal consistency (orthogonality) and derived identities (chirality
    conjugation, tensor products, Eg lift). They do not derive the table
    from the quaternion construction.
    """

    # Full character table of 2O (GAP/Groupprops)
    # Conjugacy classes: 1a, 4a, 3a, 4b, 2a, 8a, 6a, 8b
    # Class sizes:        1,  12,  8,   6,  1,  6,   8,  6
    # NOTE: T₁/T₂ labels follow GAP convention. Convention-independent:
    # T_vec  (χ(C₄)=+1, vector)      = T₂(GAP) = T₁(cryst)
    # T_pvec (χ(C₄)=-1, pseudo-vec)  = T₁(GAP) = T₂(cryst)
    SQ2 = np.sqrt(2)
    CLASS_SIZES = np.array([1, 12, 8, 6, 1, 6, 8, 6])
    CLASS_NAMES = ['1a', '4a', '3a', '4b', '2a', '8a', '6a', '8b']
    CHAR_TABLE = {
        'A0': np.array([1,  1,  1,  1,  1,  1,  1,  1], dtype=float),
        'A1': np.array([1, -1,  1,  1,  1, -1,  1, -1], dtype=float),
        'E':  np.array([2,  0, -1,  2,  2,  0, -1,  0], dtype=float),
        'H1': np.array([2,  0, -1,  0, -2, -SQ2, 1,  SQ2]),
        'H2': np.array([2,  0, -1,  0, -2,  SQ2, 1, -SQ2]),
        'T1': np.array([3,  1,  0, -1,  3, -1,  0, -1], dtype=float),
        'T2': np.array([3, -1,  0, -1,  3,  1,  0,  1], dtype=float),
        'L':  np.array([4,  0,  1,  0, -4,  0, -1,  0], dtype=float),
    }

    def _inner(self, chi_i, chi_j):
        """Weighted inner product: (1/|G|) Σ |c| χ_i(c) χ̄_j(c)."""
        return np.sum(self.CLASS_SIZES * chi_i * np.conj(chi_j)) / 48.0

    def _decompose(self, chi):
        """Decompose a character into irrep multiplicities."""
        return {name: round(self._inner(chi, irrep).real)
                for name, irrep in self.CHAR_TABLE.items()}

    def test_full_orthogonality(self):
        """All 8 irreps are orthonormal under weighted inner product."""
        names = list(self.CHAR_TABLE.keys())
        for i, ni in enumerate(names):
            for j, nj in enumerate(names):
                inner = self._inner(self.CHAR_TABLE[ni], self.CHAR_TABLE[nj])
                expected = 1.0 if i == j else 0.0
                assert abs(inner - expected) < 1e-10, (
                    f"<{ni},{nj}> = {inner}, expected {expected}"
                )

    def test_h1_times_a1_equals_h2(self):
        """H₁ ⊗ A₁ = H₂ (chirality conjugation).

        A₁(2O) = lift of sign representation of O ≅ S₄.
        This means H₁ and H₂ differ by a chirality sign,
        analogous to left-handed vs right-handed spinors.
        """
        h1_times_a1 = self.CHAR_TABLE['H1'] * self.CHAR_TABLE['A1']
        np.testing.assert_allclose(h1_times_a1, self.CHAR_TABLE['H2'], atol=1e-12)

    def test_eg_lifts_to_e_not_spinorial(self):
        """Eg(O) lifts to E(2O), which is non-spinorial (ρ(-1) = +I₂).

        The classical shear doublet transforms as Eg under O.
        Its lift to 2O is E, NOT H₁ or H₂.
        The spinorial property is a quantum sector choice, not classical.

        Projection: 2O → O maps classes as:
          1a,2a → 1;  3a,6a → C₃;  4a → C₂ₑ;  4b → C₂f;  8a,8b → C₄
        """
        # O character table: E = [2, -1, 0, 2, 0] at classes [1, C3, C2e, C2f, C4]
        # Lift to 2O by composing with projection:
        e_o = {'1': 2, 'C3': -1, 'C2e': 0, 'C2f': 2, 'C4': 0}
        projection = {
            '1a': '1', '2a': '1',
            '3a': 'C3', '6a': 'C3',
            '4a': 'C2e', '4b': 'C2f',
            '8a': 'C4', '8b': 'C4',
        }
        e_lifted = np.array([e_o[projection[c]] for c in self.CLASS_NAMES],
                            dtype=float)

        np.testing.assert_allclose(e_lifted, self.CHAR_TABLE['E'], atol=1e-12,
            err_msg="Lift of Eg(O) does not match E(2O)")

        # Verify E is non-spinorial: χ(-1) = +dim
        idx_2a = self.CLASS_NAMES.index('2a')
        assert self.CHAR_TABLE['E'][idx_2a] == +2, "E should have χ(-1) = +2"

    def test_tensor_product_h1_h1(self):
        """H₁ ⊗ H₁ = A₀ ⊕ T_vec.

        Two spinors from same sector form a scalar (A₀) + vector (T_vec = T₂ in GAP).
        The A₀ component = spinor singlet, analogous to spin-0 from two spin-½.
        """
        product = self.CHAR_TABLE['H1'] * self.CHAR_TABLE['H1']
        decomp = self._decompose(product)
        assert decomp == {'A0': 1, 'A1': 0, 'E': 0, 'H1': 0, 'H2': 0,
                          'T1': 0, 'T2': 1, 'L': 0}

    def test_tensor_product_h1_h2(self):
        """H₁ ⊗ H₂ = A₁ ⊕ T_pvec.

        Two spinors from opposite chirality sectors form a pseudo-scalar (A₁)
        + pseudo-vector (T_pvec = T₁ in GAP). No scalar singlet —
        chirality must match for pairing.
        """
        product = self.CHAR_TABLE['H1'] * self.CHAR_TABLE['H2']
        decomp = self._decompose(product)
        assert decomp == {'A0': 0, 'A1': 1, 'E': 0, 'H1': 0, 'H2': 0,
                          'T1': 1, 'T2': 0, 'L': 0}

    def test_L_is_H1_tensor_E(self):
        """L = H₁ ⊗ E: the 4D spinorial irrep derived from explicit constructions.

        H₁ character verified from SU(2) matrices (TestFundamentalRepIsH1).
        E character verified from Eg matrices (TestFundamentalRepIsH1).
        H₂ derived as H₁ ⊗ A₁ (test_h1_times_a1_equals_h2).

        Product H₁·E is shown to be irreducible and orthogonal to ALL other
        irreps — proving it is the 8th irrep without using GAP's L row.
        This closes the R2 gap: all spinorial dims (H₁=2, H₂=2, L=4)
        are derived from explicit matrix constructions.
        """
        L_derived = self.CHAR_TABLE['H1'] * self.CHAR_TABLE['E']
        H2_derived = self.CHAR_TABLE['H1'] * self.CHAR_TABLE['A1']

        # Irreducible: <L_derived, L_derived> = 1
        norm = self._inner(L_derived, L_derived)
        assert abs(norm - 1.0) < 1e-10, f"Not irreducible: <L_derived, L_derived> = {norm}"

        # Orthogonal to all 7 other irreps (using derived H₂, not hardcoded)
        others = {
            'A0': self.CHAR_TABLE['A0'],
            'A1': self.CHAR_TABLE['A1'],
            'E': self.CHAR_TABLE['E'],
            'H1': self.CHAR_TABLE['H1'],
            'H2_derived': H2_derived,
            'T1': self.CHAR_TABLE['T1'],
            'T2': self.CHAR_TABLE['T2'],
        }
        for name, chi in others.items():
            inner = self._inner(L_derived, chi)
            assert abs(inner) < 1e-10, (
                f"<L_derived, {name}> = {inner}, expected 0"
            )

    def test_h1_h2_same_at_o_classes(self):
        """H₁ and H₂ are identical when restricted to classes factoring through O.

        They differ ONLY at classes 8a, 8b (order-8 elements, lifts of C₄).
        No classical observable can distinguish them.
        """
        # Classes that project to distinct O elements (classical):
        # 1a, 4a, 3a, 4b — each maps to a unique O class
        classical_indices = [self.CLASS_NAMES.index(c)
                            for c in ['1a', '4a', '3a', '4b']]
        for idx in classical_indices:
            assert abs(self.CHAR_TABLE['H1'][idx] - self.CHAR_TABLE['H2'][idx]) < 1e-12, (
                f"H1 ≠ H2 at class {self.CLASS_NAMES[idx]} — should be equal"
            )

        # Classes where they differ: 8a, 8b (order-8 elements)
        for cn in ['8a', '8b']:
            idx = self.CLASS_NAMES.index(cn)
            assert abs(self.CHAR_TABLE['H1'][idx] - self.CHAR_TABLE['H2'][idx]) > 1.0, (
                f"H1 = H2 at class {cn} — should differ"
            )


class TestBilinearCarrierCompatibility:
    """Verify that H₁ sector adds SU(2) spinor structure to 2D fiber.

    Compares Sym² and ∧² decompositions of E (classical) vs H₁ (spinorial).
    Uses the power map (squaring) to compute symmetric/antisymmetric parts.
    Character table data: GAP SmallGroup(48,28) / Groupprops.

    Key results:
      ∧²(H₁) = A₀ → SU(2) symplectic structure (ε_{ab} invariant)
      Sym²(H₁) = T₂ → spin density is a 3D vector
      Classical observables (A₀) identical in both sectors
    """

    SQ2 = np.sqrt(2)
    CLASS_SIZES = np.array([1, 12, 8, 6, 1, 6, 8, 6])
    CLASS_NAMES = ['1a', '4a', '3a', '4b', '2a', '8a', '6a', '8b']
    # Power map: squaring each conjugacy class → target class
    POWER2 = {'1a': '1a', '4a': '2a', '3a': '3a', '4b': '2a',
              '2a': '1a', '8a': '4b', '6a': '3a', '8b': '4b'}
    CHARS = {
        'A0': np.array([1,  1,  1,  1,  1,  1,  1,  1], dtype=float),
        'A1': np.array([1, -1,  1,  1,  1, -1,  1, -1], dtype=float),
        'E':  np.array([2,  0, -1,  2,  2,  0, -1,  0], dtype=float),
        'T2': np.array([3, -1,  0, -1,  3,  1,  0,  1], dtype=float),
        'H1': np.array([2,  0, -1,  0, -2, -SQ2, 1,  SQ2]),
    }

    def _idx(self, cn):
        return self.CLASS_NAMES.index(cn)

    def _sym2(self, chi):
        """Sym²(V) character: χ_{S²}(g) = (χ(g)² + χ(g²)) / 2."""
        return np.array([(chi[self._idx(cn)]**2 +
                          chi[self._idx(self.POWER2[cn])]) / 2
                         for cn in self.CLASS_NAMES])

    def _alt2(self, chi):
        """∧²(V) character: χ_{∧²}(g) = (χ(g)² - χ(g²)) / 2."""
        return np.array([(chi[self._idx(cn)]**2 -
                          chi[self._idx(self.POWER2[cn])]) / 2
                         for cn in self.CLASS_NAMES])

    def _inner(self, chi_i, chi_j):
        return np.sum(self.CLASS_SIZES * chi_i * chi_j) / 48.0

    def test_alt2_h1_is_trivial(self):
        """∧²(H₁) = A₀: invariant antisymmetric form exists (SU(2) symplectic)."""
        alt2 = self._alt2(self.CHARS['H1'])
        np.testing.assert_allclose(alt2, self.CHARS['A0'], atol=1e-12)

    def test_alt2_e_is_pseudoscalar(self):
        """∧²(E) = A₁: antisymmetric form is a pseudoscalar, not invariant."""
        alt2 = self._alt2(self.CHARS['E'])
        np.testing.assert_allclose(alt2, self.CHARS['A1'], atol=1e-12)

    def test_sym2_h1_is_vector(self):
        """Sym²(H₁) = T_vec (= T₂ in GAP): spin density transforms as 3D vector.

        T_vec = SO(3) vector rep restricted to O (verified in TestFundamentalRepIsH1).
        See convention note in test_vector_rep_is_t2 for T_vec/T_pvec definitions.
        """
        sym2 = self._sym2(self.CHARS['H1'])
        np.testing.assert_allclose(sym2, self.CHARS['T2'], atol=1e-12)

    def test_sym2_e_is_scalar_plus_tensor(self):
        """Sym²(E) = A₀ ⊕ E: classical bilinears are scalar + 2D tensor."""
        sym2 = self._sym2(self.CHARS['E'])
        expected = self.CHARS['A0'] + self.CHARS['E']
        np.testing.assert_allclose(sym2, expected, atol=1e-12)

    def test_power_map_from_quaternions(self):
        """Derive the squaring map from quaternion arithmetic, verify POWER2.

        POWER2 (which class does g² land in?) is used to compute Sym² and ∧².
        If POWER2 were wrong, R4 and R5 would pass with incorrect results.
        This test derives the map from scratch via conjugacy class identification.
        """
        elements = build_binary_octahedral()

        # Build conjugacy classes by conjugation
        visited = set()
        classes = []
        for e in elements:
            k = qkey(e)
            if k in visited:
                continue
            cls = set()
            for g in elements:
                conj = qmul(qmul(g, e), qinv(g))
                cls.add(qkey(conj))
            for ck in cls:
                visited.add(ck)
            classes.append({
                'keys': cls, 'size': len(cls),
                'tr': 2 * e[0],  # SU(2) trace = 2w
                'rep': e.copy(),
            })

        assert len(classes) == 8, f"Expected 8 classes, got {len(classes)}"

        # Identify classes by (size, trace) — same logic as _build_class_traces
        from collections import defaultdict
        by_size = defaultdict(list)
        for ci in classes:
            by_size[ci['size']].append(ci)

        for ci in by_size[1]:
            ci['name'] = '1a' if ci['tr'] > 0 else '2a'
        by_size[12][0]['name'] = '4a'
        for ci in by_size[8]:
            ci['name'] = '3a' if ci['tr'] < 0 else '6a'
        for ci in by_size[6]:
            tr = ci['tr']
            if abs(tr) < 1e-10:
                ci['name'] = '4b'
            elif tr < 0:
                ci['name'] = '8a'
            else:
                ci['name'] = '8b'

        # Map each element key to its class name
        key_to_class = {}
        for ci in classes:
            for ck in ci['keys']:
                key_to_class[ck] = ci['name']

        # Compute power map: square each class representative
        computed = {}
        for ci in classes:
            sq = qmul(ci['rep'], ci['rep'])
            computed[ci['name']] = key_to_class[qkey(sq)]

        assert computed == self.POWER2, (
            f"Computed power map: {computed}\nHardcoded POWER2: {self.POWER2}"
        )


class TestFundamentalRepIsH1:
    """The defining representation of 2O ⊂ SU(2) is the irrep H₁.

    Method: map all 48 unit quaternions q → SU(2) via
      q = w + xi + yj + zk → [[w+iz, ix-y],[ix+y, w-iz]]
    Classify into conjugacy classes, compute tr(ρ(g)) per class,
    match to the H₁ character from GAP/Groupprops.

    Also verifies the classical representation E (Eg block of l=2)
    and the boundary condition: E(-1) = +I₂ vs H₁(-1) = -I₂.

    """

    SQ2 = np.sqrt(2)
    CLASS_NAMES = ['1a', '4a', '3a', '4b', '2a', '8a', '6a', '8b']
    CHAR_H1 = np.array([2, 0, -1, 0, -2, -np.sqrt(2), 1, np.sqrt(2)])
    CHAR_E = np.array([2, 0, -1, 2, 2, 0, -1, 0])

    @staticmethod
    def _so3_to_eg(R):
        """Extract the 2×2 Eg block from the l=2 rep of SO(3)."""
        sq2 = 1.0 / np.sqrt(2)
        sq6 = 1.0 / np.sqrt(6)
        B = [
            np.diag([1, -1, 0]) * sq2,
            np.diag([-1, -1, 2]) * sq6,
        ]
        D = np.zeros((2, 2))
        for i in range(2):
            Bp = R @ B[i] @ R.T
            for j in range(2):
                D[j, i] = np.trace(B[j] @ Bp)
        return D

    CHAR_T1 = np.array([3, 1, 0, -1, 3, -1, 0, -1], dtype=float)
    CHAR_T2 = np.array([3, -1, 0, -1, 3, 1, 0, 1], dtype=float)

    def _build_class_traces(self):
        """Build 2O, compute SU(2), Eg, and SO(3) traces per conjugacy class."""
        elements = build_binary_octahedral()
        su2_mats = [quat_to_su2(q) for q in elements]
        so3_mats = [quat_to_so3(q) for q in elements]
        eg_mats = [self._so3_to_eg(R) for R in so3_mats]

        # Classify by conjugation
        key_to_idx = {qkey(e): i for i, e in enumerate(elements)}
        visited = set()
        classes = []
        for i in range(len(elements)):
            ki = qkey(elements[i])
            if ki in visited:
                continue
            cls = set()
            for j in range(len(elements)):
                inv_j = elements[j].copy()
                inv_j[1:] *= -1
                conj = qmul(qmul(elements[j], elements[i]), inv_j)
                cls.add(qkey(conj))
            for k in cls:
                visited.add(k)
            rep_idx = key_to_idx[min(cls)]
            classes.append({
                'size': len(cls),
                'tr_fund': np.trace(su2_mats[rep_idx]).real,
                'tr_eg': np.trace(eg_mats[rep_idx]).real,
                'tr_so3': np.trace(so3_mats[rep_idx]).real,
                'su2_mat': su2_mats[rep_idx],
                'eg_mat': eg_mats[rep_idx],
            })

        # Identify classes by (size, SU(2) trace) → standard names.
        # In 2O, conjugacy classes are uniquely determined by this pair:
        # sizes {1,1,6,6,6,8,8,12} with distinct traces within each size group.
        from collections import defaultdict
        by_size = defaultdict(list)
        for ci in classes:
            by_size[ci['size']].append(ci)

        for ci in by_size[1]:
            ci['name'] = '1a' if ci['tr_fund'] > 0 else '2a'
        by_size[12][0]['name'] = '4a'
        for ci in by_size[8]:
            ci['name'] = '3a' if ci['tr_fund'] < 0 else '6a'
        for ci in by_size[6]:
            tr = ci['tr_fund']
            if abs(tr) < 1e-10:
                ci['name'] = '4b'
            elif tr < 0:
                ci['name'] = '8a'
            else:
                ci['name'] = '8b'

        return {ci['name']: ci for ci in classes}

    @pytest.fixture
    def class_data(self):
        return self._build_class_traces()

    def test_fund_character_matches_h1(self, class_data):
        """tr(ρ_fund(g)) = χ_{H₁}(g) at all 8 conjugacy classes."""
        for i, cn in enumerate(self.CLASS_NAMES):
            np.testing.assert_allclose(
                class_data[cn]['tr_fund'], self.CHAR_H1[i], atol=1e-10,
                err_msg=f"fund ≠ H₁ at class {cn}"
            )

    def test_eg_character_matches_e(self, class_data):
        """tr(ρ_Eg(g)) = χ_E(g) at all 8 conjugacy classes."""
        for i, cn in enumerate(self.CLASS_NAMES):
            np.testing.assert_allclose(
                class_data[cn]['tr_eg'], self.CHAR_E[i], atol=1e-10,
                err_msg=f"Eg ≠ E at class {cn}"
            )

    def test_su2_matrices_unitary(self):
        """All SU(2) matrices satisfy U†U = I, det U = +1."""
        elements = build_binary_octahedral()
        for q in elements:
            U = quat_to_su2(q)
            np.testing.assert_allclose(U.conj().T @ U, np.eye(2), atol=1e-12)
            det = U[0, 0] * U[1, 1] - U[0, 1] * U[1, 0]  # 2x2 det, exact
            np.testing.assert_allclose(det, 1.0, atol=1e-12)

    def test_h1_at_minus1_is_neg_identity(self, class_data):
        """H₁(-1) = -I₂ (spinorial boundary condition)."""
        mat = class_data['2a']['su2_mat']
        np.testing.assert_allclose(mat, -np.eye(2), atol=1e-12)

    def test_e_at_minus1_is_pos_identity(self, class_data):
        """E(-1) = +I₂ (factors through O, classical boundary condition)."""
        mat = class_data['2a']['eg_mat']
        np.testing.assert_allclose(mat, np.eye(2), atol=1e-12)

    def test_vector_rep_is_t2(self, class_data):
        """SO(3) vector representation restricted to O matches T₂ in our labeling.

        Verified by computing tr(R(g)) for the 3×3 SO(3) rotation matrix
        at each conjugacy class and matching to CHAR_T2.

        Convention note: In our identification of classes (via geometric
        analysis in check_class_geometry.py), the GAP irrep labeled T₂
        corresponds to the vector rep (x,y,z) — i.e. T₁ in standard
        crystallographic convention. The GAP T₁ corresponds to the other
        3D irrep (cryst T₂, pseudo-vector). This is a naming convention
        difference, not a physics difference.

        Convention-independent labels:
          T_vec  = 3D irrep with χ(C₄) = +1 = vector (x,y,z)
                 = T₂(GAP) = T₁(cryst) in our class mapping
          T_pvec = 3D irrep with χ(C₄) = -1 = pseudo-vector (yz,xz,xy)
                 = T₁(GAP) = T₂(cryst) in our class mapping

        Class labeling is derived from geometric identification:
          4a = C₂_edge (π about 110), 4b = C₂_face (π about 100),
          8a/8b = C₄ (±90° about 100), 3a/6a = C₃ (120° about 111).

        Consequence: Sym²(H₁) = T_vec, confirming that spin density
        in the H₁ sector transforms as a physical 3D vector.
        """
        # Geometric sanity checks (independent of character table):
        # C₄ classes (8a, 8b): 90° rotation → tr(SO3) = 1 + 2cos(90°) = 1
        for cn in ['8a', '8b']:
            np.testing.assert_allclose(
                class_data[cn]['tr_so3'], 1.0, atol=1e-10,
                err_msg=f"C₄ class {cn}: tr(SO3) should be 1")
        # C₂ classes (4a, 4b): 180° rotation → tr(SO3) = 1 + 2cos(180°) = -1
        for cn in ['4a', '4b']:
            np.testing.assert_allclose(
                class_data[cn]['tr_so3'], -1.0, atol=1e-10,
                err_msg=f"C₂ class {cn}: tr(SO3) should be -1")
        # C₃ classes (3a, 6a): 120° rotation → tr(SO3) = 1 + 2cos(120°) = 0
        for cn in ['3a', '6a']:
            np.testing.assert_allclose(
                class_data[cn]['tr_so3'], 0.0, atol=1e-10,
                err_msg=f"C₃ class {cn}: tr(SO3) should be 0")

        # Full character match against T₂(GAP) = T_vec
        for i, cn in enumerate(self.CLASS_NAMES):
            np.testing.assert_allclose(
                class_data[cn]['tr_so3'], self.CHAR_T2[i], atol=1e-10,
                err_msg=f"SO(3) trace ≠ T₂ at class {cn}"
            )


class TestExplicit1DCharacters:
    """Construct both 1D characters of 2O from quaternion arithmetic (no GAP table).

    2O/[2O,2O] ≅ Z₂, so exactly two 1D characters exist:
      χ₀ (trivial): all elements → +1
      χ₁ (sign on coset): 2T elements → +1, non-2T elements → -1

    Both satisfy χ(-1) = +1 (since -1 ∈ 2T = [2O,2O]).
    This is the constructive proof that no 1D spinorial character exists.
    """

    def test_sign_character_is_homomorphism(self):
        """χ₁(ab) = χ₁(a)·χ₁(b) for all a,b ∈ 2O."""
        elements = build_binary_octahedral()
        elements_2T = build_binary_tetrahedral()
        t_keys = {qkey(e) for e in elements_2T}

        def chi(q):
            return +1 if qkey(q) in t_keys else -1

        for i in range(len(elements)):
            for j in range(len(elements)):
                prod = qmul(elements[i], elements[j])
                assert chi(prod) == chi(elements[i]) * chi(elements[j]), (
                    f"χ₁ not multiplicative at elements {i}, {j}"
                )

    def test_both_characters_send_minus1_to_plus1(self):
        """χ₀(-1) = +1 and χ₁(-1) = +1: no 1D spinorial character exists."""
        elements_2T = build_binary_tetrahedral()
        t_keys = {qkey(e) for e in elements_2T}

        minus_one = np.array([-1, 0, 0, 0], dtype=float)

        # χ₀ (trivial) always returns +1
        chi_0 = +1

        # χ₁ (sign on coset): -1 ∈ 2T → +1
        chi_1 = +1 if qkey(minus_one) in t_keys else -1

        assert chi_0 == +1, "Trivial character should give +1 at -1"
        assert chi_1 == +1, "Sign character should give +1 at -1 (since -1 ∈ 2T)"

    def test_A1_matches_coset_sign_character(self):
        """A₁ from character table = coset sign character from quaternion construction.

        Glue test: links the topological result (2T = [2O,2O], coset structure)
        to the representation-theoretic fixture (A₁ row in GAP character table).
        For each conjugacy class: verifies uniform 2T membership (normality)
        and that the coset sign matches the hardcoded A₁ value.
        """
        elements = build_binary_octahedral()
        t_keys = {qkey(e) for e in build_binary_tetrahedral()}

        # Build conjugacy classes
        visited = set()
        classes = []
        for e in elements:
            k = qkey(e)
            if k in visited:
                continue
            cls = set()
            for g in elements:
                conj = qmul(qmul(g, e), qinv(g))
                cls.add(qkey(conj))
            for ck in cls:
                visited.add(ck)
            classes.append({
                'keys': cls, 'size': len(cls),
                'tr': 2 * e[0], 'rep': e.copy(),
            })

        # Identify classes by (size, trace)
        from collections import defaultdict
        by_size = defaultdict(list)
        for ci in classes:
            by_size[ci['size']].append(ci)
        for ci in by_size[1]:
            ci['name'] = '1a' if ci['tr'] > 0 else '2a'
        by_size[12][0]['name'] = '4a'
        for ci in by_size[8]:
            ci['name'] = '3a' if ci['tr'] < 0 else '6a'
        for ci in by_size[6]:
            if abs(ci['tr']) < 1e-10:
                ci['name'] = '4b'
            elif ci['tr'] < 0:
                ci['name'] = '8a'
            else:
                ci['name'] = '8b'

        # A₁ from GAP table (sign representation of O ≅ S₄)
        A1_TABLE = {
            '1a': 1, '4a': -1, '3a': 1, '4b': 1,
            '2a': 1, '8a': -1, '6a': 1, '8b': -1,
        }

        for ci in classes:
            rep_sign = +1 if qkey(ci['rep']) in t_keys else -1

            # Normality: all elements in class have same coset sign
            for ck in ci['keys']:
                elem_sign = +1 if ck in t_keys else -1
                assert elem_sign == rep_sign, (
                    f"Class {ci['name']}: non-uniform coset sign "
                    f"(normality of [2O,2O] violated)"
                )

            # Glue: coset sign matches A₁ from table
            assert rep_sign == A1_TABLE[ci['name']], (
                f"Class {ci['name']}: coset sign {rep_sign} != "
                f"A1 = {A1_TABLE[ci['name']]}"
            )

    def test_coset_structure(self):
        """2O splits into exactly 2 cosets of 2T, each of size 24."""
        elements = build_binary_octahedral()
        elements_2T = build_binary_tetrahedral()
        t_keys = {qkey(e) for e in elements_2T}

        in_2T = sum(1 for e in elements if qkey(e) in t_keys)
        not_in_2T = len(elements) - in_2T

        assert in_2T == 24, f"Expected 24 elements in 2T, got {in_2T}"
        assert not_in_2T == 24, f"Expected 24 elements outside 2T, got {not_in_2T}"


class TestComparisonZ2:
    """Contrast with Z₂ = π₁(SO(3)): the standard case.

    Z₂ is abelian ⇒ [Z₂,Z₂] = {1} ⇒ -1 ∉ [Z₂,Z₂].
    Therefore a 1D spinorial character EXISTS for Z₂.
    This is why scalar LDW works for SO(3) but not for SO(3)/O.
    """

    def test_z2_commutator_trivial(self):
        """[Z₂,Z₂] = {1} (Z₂ is abelian)."""
        z2 = np.array([[1, 0, 0, 0], [-1, 0, 0, 0]], dtype=float)
        comm = compute_commutator_subgroup(z2)

        assert len(comm) == 1, (
            f"[Z₂,Z₂] should be trivial, got order {len(comm)}"
        )
        assert MINUS_ONE_KEY not in comm, (
            "-1 should NOT be in [Z₂,Z₂]"
        )


class TestComparison2T:
    """Comparative check: 2T (binary tetrahedral) has the same property.

    -1 ∈ [2T,2T] ⇒ no 1D spinorial character for 2T either.
    This shows the property is NOT unique to octahedral symmetry.
    """

    def test_2T_abelianization(self):
        """|2T/[2T,2T]| = 3 (T ≅ A₄, [A₄,A₄] = V₄)."""
        elements = build_binary_tetrahedral()
        comm = compute_commutator_subgroup(elements)
        ab_order = len(elements) // len(comm)

        assert ab_order == 3, (
            f"Expected |2T/[2T,2T]| = 3, got {ab_order}"
        )

    def test_2T_minus_one_in_commutator(self):
        """-1 ∈ [2T,2T] ⇒ same no-1D-spinorial property as 2O."""
        elements = build_binary_tetrahedral()
        comm = compute_commutator_subgroup(elements)

        assert MINUS_ONE_KEY in comm, (
            "-1 should be in [2T,2T]"
        )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
