"""
Script 41: Commutator width of -1 in finite SU(2) subgroups.

Question: For each finite G ⊂ SU(2) with -1 ∈ [G,G],
is -1 always a SINGLE commutator [a,b], or does it require
a product of multiple commutators?

If width = 1 always: strong evidence that (3)⇒(1) in Theorem 2.4
might be provable structurally.

If width > 1 for some G: classification is genuinely needed.

Groups tested: Q₈, 2D_n (n even, n=2..20), 2T, 2O, 2I.
"""

import numpy as np
import pytest

# ============================================================
# SU(2) element utilities
# ============================================================

def su2_element(angle, axis):
    """Create SU(2) element from rotation angle and axis."""
    axis = np.array(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    half = angle / 2
    return np.array([np.cos(half),
                     np.sin(half) * axis[0],
                     np.sin(half) * axis[1],
                     np.sin(half) * axis[2]])

def qmul(a, b):
    """Quaternion multiplication."""
    return np.array([
        a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
        a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0],
    ])

def qinv(a):
    """Quaternion inverse (conjugate for unit quaternions)."""
    return np.array([a[0], -a[1], -a[2], -a[3]])

def commutator(a, b):
    """[a,b] = a b a⁻¹ b⁻¹"""
    return qmul(qmul(a, b), qmul(qinv(a), qinv(b)))

def is_minus_one(q, tol=1e-10):
    """Check if quaternion is -1 = (-1, 0, 0, 0)."""
    return abs(q[0] + 1) < tol and np.linalg.norm(q[1:]) < tol

def elements_equal(a, b, tol=1e-10):
    """Check if two quaternions represent same SU(2) element."""
    return np.linalg.norm(a - b) < tol or np.linalg.norm(a + b) < tol

# ============================================================
# Group generators
# ============================================================

def generate_group(generators, max_size=500):
    """Generate finite group from generators by closure.
    Uses saturation: keep multiplying until nothing new appears."""
    elements = [np.array([1, 0, 0, 0], dtype=float)]

    def is_in_list(q, lst):
        for e in lst:
            if np.linalg.norm(q - e) < 1e-10:
                return True
        return False

    # Add generators
    for g in generators:
        if not is_in_list(g, elements):
            elements.append(g)
        ginv = qinv(g)
        if not is_in_list(ginv, elements):
            elements.append(ginv)

    # Saturate: multiply all pairs until no new elements
    changed = True
    while changed:
        changed = False
        new_elements = []
        for a in elements:
            for b in elements:
                p = qmul(a, b)
                if not is_in_list(p, elements) and not is_in_list(p, new_elements):
                    new_elements.append(p)
                    if len(elements) + len(new_elements) > max_size:
                        raise ValueError(f"Group exceeds {max_size} elements")
        if new_elements:
            elements.extend(new_elements)
            changed = True

    return elements

def make_Q8():
    """Q₈ = {±1, ±i, ±j, ±k}"""
    return generate_group([
        np.array([0, 1, 0, 0]),  # i
        np.array([0, 0, 1, 0]),  # j
    ])

def make_2Dn(n):
    """Binary dihedral 2D_n, order 4n.
    Generators: a = exp(iπ/n), b = j (in quaternion terms).
    a has order 2n, b has order 4, b a b⁻¹ = a⁻¹.
    """
    angle = np.pi / n
    a = np.array([np.cos(angle), np.sin(angle), 0, 0])
    b = np.array([0, 0, 1, 0])
    return generate_group([a, b])

def make_2T():
    """Binary tetrahedral 2T, order 24."""
    # Generators: i, j, and (1+i+j+k)/2
    omega = np.array([0.5, 0.5, 0.5, 0.5])
    i = np.array([0, 1, 0, 0])
    return generate_group([i, omega])

def make_2O():
    """Binary octahedral 2O, order 48."""
    # Add (1+i)/√2 to 2T generators
    omega = np.array([0.5, 0.5, 0.5, 0.5])
    s = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0, 0])
    return generate_group([omega, s])

def make_2I():
    """Binary icosahedral 2I, order 120."""
    phi = (1 + np.sqrt(5)) / 2  # golden ratio
    # Generators: (1+i+j+k)/2 and (φ + φ⁻¹i + j)/2
    omega = np.array([0.5, 0.5, 0.5, 0.5])
    tau = np.array([phi/2, 1/(2*phi), 0.5, 0])
    return generate_group([omega, tau])


# ============================================================
# Commutator width computation
# ============================================================

def find_minus_one_as_commutator(elements):
    """Search for a, b ∈ G such that [a,b] = -1.
    Returns (a, b) if found, None otherwise."""
    for a in elements:
        for b in elements:
            c = commutator(a, b)
            if is_minus_one(c):
                return (a, b)
    return None

def minus_one_in_commutator_subgroup(elements):
    """Check if -1 is in [G,G] (as product of commutators)."""
    minus_one = np.array([-1, 0, 0, 0], dtype=float)
    # First check: is -1 in the group?
    found = False
    for e in elements:
        if np.linalg.norm(e - minus_one) < 1e-10:
            found = True
            break
    if not found:
        return False

    # Compute all single commutators
    comm_set = []
    for a in elements:
        for b in elements:
            c = commutator(a, b)
            is_new = True
            for existing in comm_set:
                if np.linalg.norm(c - existing) < 1e-10:
                    is_new = False
                    break
            if is_new:
                comm_set.append(c)

    # Check if -1 is a single commutator
    for c in comm_set:
        if is_minus_one(c):
            return True
    return False


# ============================================================
# Tests
# ============================================================

class TestCommutatorWidth:
    """Test commutator width of -1 for all obstructed SU(2) subgroups."""

    def test_Q8(self):
        """Q₈: smallest obstructed group. [i,j] = -1."""
        G = make_Q8()
        assert len(G) == 8
        result = find_minus_one_as_commutator(G)
        assert result is not None, "Q₈: -1 is NOT a single commutator!"
        a, b = result
        c = commutator(a, b)
        assert is_minus_one(c)
        print(f"Q₈ (|G|=8): -1 = [a,b] with a={a}, b={b}")

    def test_2D_even(self):
        """2D_n for even n: all should have -1 as single commutator."""
        for n in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]:
            G = make_2Dn(n)
            assert len(G) == 4 * n, f"2D_{n}: expected {4*n}, got {len(G)}"
            result = find_minus_one_as_commutator(G)
            assert result is not None, f"2D_{n}: -1 is NOT a single commutator!"
            a, b = result
            c = commutator(a, b)
            assert is_minus_one(c)
            print(f"2D_{n} (|G|={4*n}): -1 = [a,b] ✓")

    def test_2D_odd_no_obstruction(self):
        """2D_n for odd n: -1 should NOT be in [G,G]."""
        for n in [3, 5, 7, 9]:
            G = make_2Dn(n)
            assert len(G) == 4 * n
            result = find_minus_one_as_commutator(G)
            assert result is None, f"2D_{n}: -1 IS a commutator, but n is odd!"
            print(f"2D_{n} (|G|={4*n}): -1 ∉ single commutators ✓")

    def test_2T(self):
        """2T: binary tetrahedral, order 24."""
        G = make_2T()
        assert len(G) == 24
        result = find_minus_one_as_commutator(G)
        assert result is not None, "2T: -1 is NOT a single commutator!"
        a, b = result
        c = commutator(a, b)
        assert is_minus_one(c)
        # Count how many pairs give -1
        count = sum(1 for a in G for b in G if is_minus_one(commutator(a, b)))
        print(f"2T (|G|=24): -1 = [a,b] ✓ ({count} pairs out of {24*24})")

    def test_2O(self):
        """2O: binary octahedral, order 48."""
        G = make_2O()
        assert len(G) == 48
        result = find_minus_one_as_commutator(G)
        assert result is not None, "2O: -1 is NOT a single commutator!"
        a, b = result
        c = commutator(a, b)
        assert is_minus_one(c)
        count = sum(1 for a in G for b in G if is_minus_one(commutator(a, b)))
        print(f"2O (|G|=48): -1 = [a,b] ✓ ({count} pairs out of {48*48})")

    def test_2I(self):
        """2I: binary icosahedral, order 120."""
        G = make_2I()
        assert len(G) == 120
        result = find_minus_one_as_commutator(G)
        assert result is not None, "2I: -1 is NOT a single commutator!"
        a, b = result
        c = commutator(a, b)
        assert is_minus_one(c)
        count = sum(1 for a in G for b in G if is_minus_one(commutator(a, b)))
        print(f"2I (|G|=120): -1 = [a,b] ✓ ({count} pairs out of {120*120})")


class TestCommutatorSubgroupStructure:
    """Deeper analysis: what fraction of commutators equal -1?"""

    def test_commutator_density(self):
        """For each obstructed group, compute:
        - Total commutators [a,b]
        - How many equal -1
        - Fraction
        """
        groups = {
            'Q₈': make_Q8(),
            '2D₄': make_2Dn(4),
            '2D₆': make_2Dn(6),
            '2T': make_2T(),
            '2O': make_2O(),
            '2I': make_2I(),
        }

        print("\nCommutator density of -1:")
        print(f"{'Group':>6} | {'|G|':>4} | {'|G|²':>6} | "
              f"{'[a,b]=-1':>8} | {'fraction':>10}")
        print("-" * 50)

        for name, G in groups.items():
            n = len(G)
            count = sum(1 for a in G for b in G
                       if is_minus_one(commutator(a, b)))
            frac = count / (n * n)
            print(f"{name:>6} | {n:>4} | {n*n:>6} | "
                  f"{count:>8} | {frac:>10.4f}")
            assert count > 0, f"{name}: -1 not found as commutator!"

    def test_commutator_subgroup_equals_single_commutators(self):
        """Check: is the set of single commutators already a subgroup?
        If yes: [G,G] = {single commutators} and width is always 1.
        If no: there exist products of commutators that are not single commutators.
        """
        groups = {
            'Q₈': make_Q8(),
            '2D₄': make_2Dn(4),
            '2T': make_2T(),
            '2O': make_2O(),
            '2I': make_2I(),
        }

        print("\nSingle commutators vs commutator subgroup:")

        for name, G in groups.items():
            # Compute all single commutators
            single_comms = []
            for a in G:
                for b in G:
                    c = commutator(a, b)
                    is_new = True
                    for existing in single_comms:
                        if np.linalg.norm(c - existing) < 1e-10:
                            is_new = False
                            break
                    if is_new:
                        single_comms.append(c)

            # Compute commutator subgroup (closure of single commutators)
            comm_subgroup = generate_group(single_comms, max_size=len(G))

            n_single = len(single_comms)
            n_subgroup = len(comm_subgroup)

            is_subgroup = (n_single == n_subgroup)
            status = "= subgroup" if is_subgroup else "⊊ subgroup"

            print(f"{name:>6}: |single comms| = {n_single:>3}, "
                  f"|[G,G]| = {n_subgroup:>3}  {status}")

            # For the groups we care about, single commutators
            # should already form the commutator subgroup
            # (this is the key question!)


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
