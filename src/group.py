"""
Binary polyhedral group construction and algorithms.

Groups: 2T (order 24), 2O (order 48), 2I (order 120).
All realized as finite subgroups of SU(2) via unit quaternions.
"""
import numpy as np
from .quaternion import qmul, qinv, qkey


# ============================================================
# Group construction
# ============================================================

def build_binary_tetrahedral():
    """Construct all 24 elements of 2T as unit quaternions.
    2T = Q8 ∪ {(±1±i±j±k)/2}
    """
    elements = []
    for s in [1, -1]:
        elements.append(np.array([s, 0, 0, 0], dtype=float))
        elements.append(np.array([0, s, 0, 0], dtype=float))
        elements.append(np.array([0, 0, s, 0], dtype=float))
        elements.append(np.array([0, 0, 0, s], dtype=float))
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                for s4 in [1, -1]:
                    elements.append(np.array([s1, s2, s3, s4]) / 2.0)
    return np.array(elements)


def build_binary_octahedral():
    """Construct all 48 elements of 2O as unit quaternions.
    2O = Q8 ∪ {(±1±i±j±k)/2} ∪ {(±a±b)/√2, a≠b}
    """
    elements = []

    # Q8: {±1, ±i, ±j, ±k}
    for s in [1, -1]:
        elements.append(np.array([s, 0, 0, 0], dtype=float))
        elements.append(np.array([0, s, 0, 0], dtype=float))
        elements.append(np.array([0, 0, s, 0], dtype=float))
        elements.append(np.array([0, 0, 0, s], dtype=float))

    # 2T \ Q8: (±1±i±j±k)/2, all 16 sign combos
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                for s4 in [1, -1]:
                    elements.append(np.array([s1, s2, s3, s4]) / 2.0)

    # 2O \ 2T: (±a±b)/√2 where a,b distinct from {1,i,j,k}
    sq2 = 1.0 / np.sqrt(2)
    for i, j in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
        for si in [1, -1]:
            for sj in [1, -1]:
                q = np.zeros(4)
                q[i] = si * sq2
                q[j] = sj * sq2
                elements.append(q)

    return np.array(elements)


def build_binary_icosahedral():
    """Construct all 120 elements of 2I as unit quaternions.

    2I = 2T ∪ {even permutations of (0, ±1, ±φ⁻¹, ±φ)/2}
    where φ = (1+√5)/2 (golden ratio).

    Normalization check: (0² + 1² + φ⁻² + φ²)/4 = (0+1+3)/4 = 1.
    """
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    inv_phi = phi - 1.0  # = 1/φ = (√5-1)/2

    elements = list(build_binary_tetrahedral())  # 24 elements

    # Base values to be permuted
    vals = [0.0, 1.0, inv_phi, phi]

    # All 12 even permutations of {0,1,2,3}
    # Identity + 8 three-cycles + 3 double transpositions
    even_perms = [
        [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2],
        [1, 0, 3, 2], [1, 2, 0, 3], [1, 3, 2, 0],
        [2, 0, 1, 3], [2, 1, 3, 0], [2, 3, 0, 1],
        [3, 0, 2, 1], [3, 1, 0, 2], [3, 2, 1, 0],
    ]

    for perm in even_perms:
        base_q = np.array([vals[perm[i]] for i in range(4)])
        zero_pos = int(np.argmin(np.abs(base_q)))
        nonzero = [i for i in range(4) if i != zero_pos]
        # 2³ = 8 sign combinations for 3 nonzero entries
        for s0 in [1, -1]:
            for s1 in [1, -1]:
                for s2 in [1, -1]:
                    q = base_q.copy()
                    q[nonzero[0]] *= s0
                    q[nonzero[1]] *= s1
                    q[nonzero[2]] *= s2
                    elements.append(q / 2.0)

    return np.array(elements)


# ============================================================
# Group algorithms
# ============================================================

def verify_closure(elements, name="group"):
    """Verify closure under quaternion multiplication.
    Returns True or raises AssertionError.
    """
    key_set = {qkey(e) for e in elements}
    n = len(elements)
    for i in range(n):
        for j in range(n):
            prod = qmul(elements[i], elements[j])
            assert qkey(prod) in key_set, (
                f"{name}: product of elements {i} and {j} not in group"
            )
    return True


def compute_commutator_subgroup(elements):
    """Compute [G,G] by generating all commutators and closing.
    Returns dict {qkey: quaternion_array}.
    """
    n = len(elements)
    key_to_q = {}

    def add(q):
        k = qkey(q)
        if k not in key_to_q:
            key_to_q[k] = q.copy()
            return True
        return False

    add(np.array([1, 0, 0, 0], dtype=float))

    for i in range(n):
        for j in range(n):
            comm = qmul(
                qmul(elements[i], elements[j]),
                qmul(qinv(elements[i]), qinv(elements[j]))
            )
            add(comm)

    changed = True
    while changed:
        changed = False
        current_qs = list(key_to_q.values())
        for a in current_qs:
            for b in current_qs:
                prod = qmul(a, b)
                if add(prod):
                    changed = True

    return key_to_q


def compute_conjugacy_classes(elements):
    """Compute conjugacy classes by conjugation.
    Returns list of dicts: {'keys': set, 'size': int, 'rep': array, 'trace': float}.
    trace = 2*w = SU(2) trace.
    """
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
            'keys': cls,
            'size': len(cls),
            'rep': e.copy(),
            'trace': 2 * e[0],
        })
    return classes


def identify_classes(classes):
    """Assign standard names to conjugacy classes by (size, SU(2) trace).
    Works for 2O (8 classes). 2T (7 classes) partial — two class pairs share
    (size, trace) and need additional invariant to distinguish.
    Returns dict {class_name: class_dict}.
    """
    from collections import defaultdict
    by_size = defaultdict(list)
    for ci in classes:
        by_size[ci['size']].append(ci)

    n_classes = len(classes)

    if n_classes == 8:
        # 2O: sizes {1,1,6,6,6,8,8,12}
        for ci in by_size[1]:
            ci['name'] = '1a' if ci['trace'] > 0 else '2a'
        by_size[12][0]['name'] = '4a'
        for ci in by_size[8]:
            ci['name'] = '3a' if ci['trace'] < 0 else '6a'
        for ci in by_size[6]:
            tr = ci['trace']
            if abs(tr) < 1e-10:
                ci['name'] = '4b'
            elif tr < 0:
                ci['name'] = '8a'
            else:
                ci['name'] = '8b'
    elif n_classes == 7:
        # 2T: sizes {1,1,4,4,4,4,6}
        for ci in by_size[1]:
            ci['name'] = '1a' if ci['trace'] > 0 else '2a'
        by_size[6][0]['name'] = '4a'
        for ci in by_size[4]:
            tr = ci['trace']
            if tr > 0.5:
                ci['name'] = '6a'
            elif tr < -0.5:
                ci['name'] = '3a'
            elif tr > -0.5 and tr < 0.5:
                ci['name'] = '6b' if tr > 0 else '3b'
    else:
        raise ValueError(f"identify_classes: {n_classes} classes not handled yet")

    return {ci['name']: ci for ci in classes}


# ============================================================
# Constants
# ============================================================

MINUS_ONE = np.array([-1, 0, 0, 0], dtype=float)
MINUS_ONE_KEY = qkey(MINUS_ONE)
