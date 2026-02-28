#!/usr/bin/env python3
"""
Sector-resolved spectral comparison under symmetry reduction O -> D3.

PHYSICAL SETUP:
  Rigid rotor on SO(3)/H, quantization sectors = Hom(pi_1, U(1)).
  Sectors are superselected: a system is in ONE sector, fixed by
  boundary conditions, not by dynamics or temperature.

SECTOR CORRESPONDENCE under O -> D3:
  2O has 2 scalar sectors: A0(T), A1(T). Both tensorial.
  2D3 has 4 scalar sectors: chi_0(T), chi_1(S), chi_2(T), chi_3(S).

  Restriction of 2O characters to 2D3 subgroup:
    A0(2O)|_{2D3} = chi_0(T)  [trivial -> trivial]
    A1(2O)|_{2D3} = chi_2(T)  [sign -> sign]

  chi_1(S), chi_3(S) have no ancestor in O.
  These spinorial sectors EXIST on SO(3)/D3 but not on SO(3)/O.

  NOTE: Sector matching uses trace + commutator-subgroup membership to
  identify conjugacy classes. Pure trace lookup would fail here (2D3
  has two classes with trace 0). The commutator membership disambiguates.

SPECTRAL DENSIFICATION (same sector, different symmetry):
  A0(O):    gap = 20B, j = 0, 4, 6, 8, 9, 10, ...
  chi_0(D3): gap = 6B,  j = 0, 2, 3, 4, 5, 6, ...  (3.3x denser)

  A1(O):    gap = 12B, j = 3, 6, 7, 9, 10, 11, ...
  chi_2(D3): gap = 2B,  j = 1, 3, 4, 5, 6, 7, ...  (6x denser)

  Levels forbidden by O symmetry become allowed under D3 symmetry.

SPINORIAL SECTORS (forbidden on SO(3)/O, allowed on SO(3)/D3):
  chi_1(S), chi_3(S): ground state at j=3/2, E=3.75B
  These sectors are topologically forbidden on SO(3)/O.

ENERGY SCALE: SF6 rotational constant B = 0.09107 cm-1 = 0.1310 K.
  SF6 is used as an energy scale for octahedral rotors. SF6 itself
  does not undergo O -> D3; a material with trigonal distortion of an
  octahedral cage (e.g., BaTiO3-type perovskite) would be needed.

Asserts: 22
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.quaternion import qkey, qmul, qconj
from src.group import (build_binary_octahedral, build_binary_tetrahedral,
                       compute_conjugacy_classes, compute_commutator_subgroup)

# ==========================================================
# Physical constants
# ==========================================================
B_CM = 0.09107      # SF6 rotational constant [cm-1]
CM_K = 1.4388       # 1 cm-1 = 1.4388 K
B_K  = B_CM * CM_K  # B in Kelvin


# ==========================================================
# Mathematical infrastructure
# ==========================================================

def chi_su2(j, trace):
    """SU(2) spin-j character at group element with given trace."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def build_dicyclic(N):
    """Build Dic_N = 2D_N of order 4N. Analytic construction."""
    elements = []
    for k in range(2 * N):
        angle = k * np.pi / N
        elements.append(np.array([np.cos(angle), 0, 0, np.sin(angle)]))
    b = np.array([0, 0, 1, 0], dtype=float)
    for k in range(2 * N):
        angle = k * np.pi / N
        ak = np.array([np.cos(angle), 0, 0, np.sin(angle)])
        elements.append(qmul(b, ak))
    return np.array(elements)


def generate_group(generators, max_order=300):
    """Generate group from generators by closure."""
    IDENTITY = np.array([1, 0, 0, 0], dtype=float)
    elements = [IDENTITY.copy()]
    seen = {qkey(elements[0])}
    queue = list(elements)
    idx = 0
    while idx < len(queue):
        g = queue[idx]; idx += 1
        for gen in generators:
            for h in [qmul(g, gen), qmul(gen, g)]:
                k = qkey(h)
                if k not in seen:
                    seen.add(k)
                    queue.append(h)
                    elements.append(h)
        if len(elements) > max_order:
            raise RuntimeError(f"Group too large: {len(elements)}")
    return elements


def compute_1d_characters(elements):
    """Compute all 1D characters. Returns (characters, classes, sizes, |G|, obstructed)."""
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    G = int(sum(sizes))

    comm = compute_commutator_subgroup(elements)
    comm_keys = set(comm.keys())
    mk = qkey(np.array([-1, 0, 0, 0], dtype=float))
    obstructed = mk in comm_keys
    abel_size = G // len(comm)

    cosets = []
    elem_coset = {}
    for e in elements:
        k = qkey(e)
        if k in elem_coset:
            continue
        coset = set()
        for ck, cv in comm.items():
            coset.add(qkey(qmul(e, cv)))
        found = False
        for ci, existing in enumerate(cosets):
            if existing & coset:
                for kk in coset:
                    elem_coset[kk] = ci
                found = True
                break
        if not found:
            ci = len(cosets)
            cosets.append(coset)
            for kk in coset:
                elem_coset[kk] = ci

    gen_elem = None
    if abel_size > 1:
        for e in elements:
            k = qkey(e)
            if k in comm_keys:
                continue
            power = e.copy()
            for d in range(1, abel_size + 1):
                pk = qkey(power)
                if pk in comm_keys:
                    if d == abel_size:
                        gen_elem = e.copy()
                    break
                power = qmul(power, e)
            if gen_elem is not None:
                break

    coset_label = {}
    if abel_size == 1:
        coset_label[elem_coset[qkey(np.array([1, 0, 0, 0]))]] = 0
    else:
        power = np.array([1, 0, 0, 0], dtype=float)
        for d in range(abel_size):
            ci = elem_coset[qkey(power)]
            coset_label[ci] = d
            power = qmul(power, gen_elem)

    class_labels = np.array([coset_label[elem_coset[list(c['keys'])[0]]]
                             for c in classes])
    omega = np.exp(2j * np.pi / abel_size) if abel_size > 1 else 1.0

    minus_one_label = 0
    for i, c in enumerate(classes):
        if abs(c['trace'] + 2) < 1e-10:
            minus_one_label = class_labels[i]
            break

    characters = []
    for k in range(abel_size):
        chi_k = omega ** (k * class_labels)
        chi_m1 = (omega ** (k * minus_one_label)).real
        is_spin = chi_m1 < -0.5
        characters.append({
            'name': f'chi_{k}',
            'chi_on_classes': chi_k,
            'chi_minus_one': chi_m1,
            'is_spinorial': is_spin,
        })

    return characters, classes, sizes, G, obstructed


def compute_multiplicities(characters, classes, sizes, G, j_max=30):
    """Compute m(chi, j) for all characters and j up to j_max."""
    mult = {}
    for ci, char in enumerate(characters):
        chi_k = char['chi_on_classes']
        is_spin = char['is_spinorial']
        for two_j in range(0, int(2 * j_max) + 1):
            j = two_j / 2.0
            if is_spin and two_j % 2 == 0:
                continue
            if not is_spin and two_j % 2 == 1:
                continue
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            m = round((np.sum(chi_k.conj() * chi_Vj * sizes) / G).real)
            if m > 0:
                mult[(ci, j)] = m
    return mult


def Z_sector(mult, ci, T):
    """Partition function for sector ci at temperature T [Kelvin]."""
    Z = 0.0
    for (c, j), m in mult.items():
        if c == ci:
            Z += m * (2 * j + 1) * np.exp(-B_K * j * (j + 1) / T)
    return Z


def C_sector(mult, ci, T):
    """Specific heat C/k_B for sector ci at temperature T [Kelvin]."""
    Z = 0.0; E1 = 0.0; E2 = 0.0
    for (c, j), m in mult.items():
        if c == ci:
            Ej = B_K * j * (j + 1)
            w = m * (2 * j + 1) * np.exp(-Ej / T)
            Z += w; E1 += w * Ej; E2 += w * Ej**2
    if Z < 1e-30:
        return 0.0
    return (E2 / Z - (E1 / Z)**2) / T**2


# ==========================================================
# Main
# ==========================================================

def main():
    n_asserts = 0

    # Build groups
    elements_O = np.array(build_binary_octahedral())
    elements_D3 = build_dicyclic(3)
    elements_2T = build_binary_tetrahedral()
    keys_2T = set(qkey(e) for e in elements_2T)

    # Characters and multiplicities
    chars_O, cls_O, sz_O, G_O, obst_O = compute_1d_characters(elements_O)
    chars_D3, cls_D3, sz_D3, G_D3, obst_D3 = compute_1d_characters(elements_D3)

    J_MAX = 50
    mult_O = compute_multiplicities(chars_O, cls_O, sz_O, G_O, j_max=J_MAX)
    mult_D3 = compute_multiplicities(chars_D3, cls_D3, sz_D3, G_D3, j_max=J_MAX)

    assert obst_O and not obst_D3; n_asserts += 1
    assert len(chars_O) == 2 and len(chars_D3) == 4; n_asserts += 1

    # Build element->class index map for 2D3 (avoids trace-based lookup,
    # which fails when two classes share the same trace, as here: classes
    # 2 and 3 both have trace 0)
    elem_to_class_D3 = {}
    for ci_c, c in enumerate(cls_D3):
        for k in c['keys']:
            elem_to_class_D3[k] = ci_c
    assert len(elem_to_class_D3) == G_D3; n_asserts += 1

    # ===========================================================
    # PART 1: SECTOR CORRESPONDENCE O -> D3
    # ===========================================================
    print("=" * 70)
    print("PART 1: SECTOR CORRESPONDENCE O -> D3")
    print("=" * 70)

    # Build 2D3 as subgroup of 2O (using standard generators)
    SQ2 = np.sqrt(2)
    a_D3 = np.array([0.5, 0.5, 0.5, 0.5])
    b_D3 = np.array([0, 1/SQ2, -1/SQ2, 0])
    elems_2D3_in_2O = generate_group([a_D3, b_D3, qconj(a_D3)])
    assert len(elems_2D3_in_2O) == 12; n_asserts += 1

    # Verify all elements of 2D3 are in 2O
    keys_2O = set(qkey(e) for e in elements_O)
    for g in elems_2D3_in_2O:
        assert qkey(g) in keys_2O
    n_asserts += 1

    # A0(2O) = trivial character -> trivially restricts to chi_0(2D3)
    # A1(2O) = sign character: +1 on [2O,2O]=2T, -1 outside
    # Check: A1 restricted to 2D3 elements
    n_in_2T = sum(1 for g in elems_2D3_in_2O if qkey(g) in keys_2T)
    n_out_2T = len(elems_2D3_in_2O) - n_in_2T
    print(f"\n  2D3 elements in 2T (= [2O,2O]): {n_in_2T}/12")
    print(f"  2D3 elements NOT in 2T: {n_out_2T}/12")
    assert n_in_2T == 6 and n_out_2T == 6; n_asserts += 1

    # Identify which chi_k of 2D3 equals A1|_{2D3}
    # Build the coset structure of 2D3 via compute_1d_characters
    comm_D3 = compute_commutator_subgroup(elements_D3)
    comm_D3_keys = set(comm_D3.keys())
    abel_D3 = G_D3 // len(comm_D3)

    # We need the coset labels for elements of 2D3-as-subgroup-of-2O
    # Use the character values directly: evaluate each chi_k on 2D3 elements
    # and compare with A1 values
    # Map 2D3-as-subgroup-of-2O elements to 2D3-standalone element keys.
    # Both constructions produce the same abstract group, but element keys
    # differ. Match by trace + norm of imaginary part to build the map.
    # For 2D3 (|G|=12), we can match exactly via the isomorphism.
    #
    # Simpler approach: the inner product sum over elements uses A1 values
    # (from 2T membership) and chi_k values (from class index). We need
    # to map each 2D3-in-2O element to the correct 2D3 class index.
    # Since elements_D3 and elems_2D3_in_2O are different embeddings of
    # the same group, match via trace (= 2*q0) which is conjugation-invariant.
    # When traces collide (classes 2,3 both trace=0), further distinguish
    # by whether the element is in [2D3,2D3].
    comm_D3 = compute_commutator_subgroup(elements_D3)
    comm_D3_keys = set(comm_D3.keys())

    # Build 2D3-in-2O commutator subgroup to distinguish trace-degenerate classes
    comm_D3_in_2O = compute_commutator_subgroup(elems_2D3_in_2O)
    comm_D3_in_2O_keys = set(qkey(e) if isinstance(e, np.ndarray) else e
                              for e in comm_D3_in_2O.keys())

    print("\n  Matching A1(2O)|_{2D3} to chi_k(2D3):")
    a1_match = None
    for k in range(4):
        chi_k_vals = chars_D3[k]['chi_on_classes']
        # Compute (1/|G|) sum chi_k* A1 over 2D3 elements.
        # Map each element to its class in the standalone 2D3.
        inner = 0.0
        for g in elems_2D3_in_2O:
            a1_val = 1.0 if qkey(g) in keys_2T else -1.0
            tr = 2 * g[0]
            in_comm = qkey(g) in comm_D3_in_2O_keys
            # Find matching class: same trace AND same commutator membership
            matched = False
            for ci_c, c in enumerate(cls_D3):
                if abs(c['trace'] - tr) < 1e-6:
                    # Check commutator membership consistency
                    rep_key = list(c['keys'])[0]
                    c_in_comm = rep_key in comm_D3_keys
                    if c_in_comm == in_comm:
                        chi_k_val = chi_k_vals[ci_c].real
                        inner += a1_val * chi_k_val
                        matched = True
                        break
            assert matched, f"Element {qkey(g)} not matched to any class"
        inner /= G_D3
        is_spin = chars_D3[k]['is_spinorial']
        typ = 'S' if is_spin else 'T'
        print(f"    <A1, chi_{k}({typ})> = {inner:.4f}")
        if abs(inner - 1.0) < 0.01:
            a1_match = k

    assert a1_match is not None; n_asserts += 1
    print(f"\n  Result: A0(2O) -> chi_0(T),  A1(2O) -> chi_{a1_match}(T)")
    # Property-based check: a1_match must be tensorial and non-trivial
    assert not chars_D3[a1_match]['is_spinorial']; n_asserts += 1
    assert a1_match != 0; n_asserts += 1  # non-trivial (A1 is the sign char)

    print("\n  Spinorial sectors chi_1(S), chi_3(S): no ancestor in 2O")
    print("  These exist on SO(3)/D3 but are topologically forbidden on SO(3)/O")

    # ===========================================================
    # PART 2: SPECTRAL DENSIFICATION PER SECTOR
    # ===========================================================
    print("\n" + "=" * 70)
    print("PART 2: SPECTRAL DENSIFICATION PER SECTOR")
    print("=" * 70)

    B_str = f"B = {B_K:.4f} K"

    # A0(O) -> chi_0(D3)
    print(f"\n  A0(O) -> chi_0(D3)  [{B_str}]")
    print(f"  {'j':>5} {'E [K]':>8}  {'m(A0)':>6} {'m(chi0)':>8}")
    print("  " + "-" * 35)
    j_list_0 = set()
    for (ci, j), m in mult_O.items():
        if ci == 0:
            j_list_0.add(j)
    for (ci, j), m in mult_D3.items():
        if ci == 0:
            j_list_0.add(j)
    for j in sorted(j_list_0):
        if j > 12:
            break
        m_O = mult_O.get((0, j), 0)
        m_D3 = mult_D3.get((0, j), 0)
        if m_O > 0 or m_D3 > 0:
            E = B_K * j * (j + 1)
            print(f"  {j:5.1f} {E:8.3f}  {m_O:6d} {m_D3:8d}")

    # Compute gaps from data
    first_A0 = min(j for (ci, j), m in mult_O.items() if ci == 0 and j > 0)
    first_chi0 = min(j for (ci, j), m in mult_D3.items() if ci == 0 and j > 0)
    gap_A0_O = first_A0 * (first_A0 + 1)      # j(j+1) in units of B
    gap_chi0_D3 = first_chi0 * (first_chi0 + 1)
    assert gap_A0_O == 20; n_asserts += 1      # j=4: 4*5=20
    assert gap_chi0_D3 == 6; n_asserts += 1    # j=2: 2*3=6
    print(f"\n  Gap: {gap_A0_O:.0f}B = {gap_A0_O*B_K:.3f} K  ->  "
          f"{gap_chi0_D3:.0f}B = {gap_chi0_D3*B_K:.3f} K  "
          f"(ratio {gap_A0_O/gap_chi0_D3:.1f}x)")

    # A1(O) -> chi_2(D3)
    print(f"\n  A1(O) -> chi_2(D3)  [{B_str}]")
    print(f"  {'j':>5} {'E [K]':>8}  {'m(A1)':>6} {'m(chi2)':>8}")
    print("  " + "-" * 35)
    j_list_1 = set()
    for (ci, j), m in mult_O.items():
        if ci == 1:
            j_list_1.add(j)
    for (ci, j), m in mult_D3.items():
        if ci == a1_match:
            j_list_1.add(j)
    for j in sorted(j_list_1):
        if j > 12:
            break
        m_O = mult_O.get((1, j), 0)
        m_D3 = mult_D3.get((a1_match, j), 0)
        if m_O > 0 or m_D3 > 0:
            E = B_K * j * (j + 1)
            print(f"  {j:5.1f} {E:8.3f}  {m_O:6d} {m_D3:8d}")

    # Compute gaps from data
    first_A1 = min(j for (ci, j), m in mult_O.items() if ci == 1 and j > 0)
    first_chi2 = min(j for (ci, j), m in mult_D3.items()
                     if ci == a1_match and j > 0)
    gap_A1_O = first_A1 * (first_A1 + 1)
    gap_chi2_D3 = first_chi2 * (first_chi2 + 1)
    assert gap_A1_O == 12; n_asserts += 1      # j=3: 3*4=12
    assert gap_chi2_D3 == 2; n_asserts += 1    # j=1: 1*2=2
    print(f"\n  Gap: {gap_A1_O:.0f}B = {gap_A1_O*B_K:.3f} K  ->  "
          f"{gap_chi2_D3:.0f}B = {gap_chi2_D3*B_K:.3f} K  "
          f"(ratio {gap_A1_O/gap_chi2_D3:.1f}x)")

    # New spinorial sectors
    print(f"\n  NEW spinorial sectors chi_1(S), chi_3(S):")
    for ci in [1, 3]:
        first_j = min(j for (c, j), m in mult_D3.items() if c == ci)
        E_first = B_K * first_j * (first_j + 1)
        print(f"    chi_{ci}(S): ground state at j = {first_j}, "
              f"E = {first_j*(first_j+1):.2f}B = {E_first:.3f} K")
    first_spin = min(j for (c, j), m in mult_D3.items()
                     if chars_D3[c]['is_spinorial'])
    assert first_spin == 1.5; n_asserts += 1

    # ===========================================================
    # PART 3: PARTITION FUNCTION AND SPECIFIC HEAT PER SECTOR
    # ===========================================================
    print("\n" + "=" * 70)
    print("PART 3: PARTITION FUNCTION PER SECTOR (physical units)")
    print("=" * 70)

    print(f"\n  Z_chi(T) = sum_j m(chi,j) (2j+1) exp(-B j(j+1)/T)")
    print(f"  B = {B_K:.4f} K (SF6)")

    # A0(O) vs chi_0(D3) — the ground-state sector
    print(f"\n  Ground-state sector: A0(O) vs chi_0(D3)")
    print(f"  {'T [K]':>7} {'Z(A0,O)':>10} {'Z(chi0,D3)':>12} {'ratio':>7} "
          f"{'C(A0)':>8} {'C(chi0)':>9}")
    print("  " + "-" * 58)

    temps = [0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0]
    for T in temps:
        Z_O = Z_sector(mult_O, 0, T)
        Z_D3 = Z_sector(mult_D3, 0, T)
        C_O = C_sector(mult_O, 0, T)
        C_D3 = C_sector(mult_D3, 0, T)
        ratio = Z_D3 / Z_O if Z_O > 1e-30 else float('inf')
        print(f"  {T:7.2f} {Z_O:10.4f} {Z_D3:12.4f} {ratio:7.2f} "
              f"{C_O:8.4f} {C_D3:9.4f}")

    # High-T ratio: |2O|/|2D3| = 48/12 = 4 (for trivial sector)
    Z_O_hi = Z_sector(mult_O, 0, 10.0)
    Z_D3_hi = Z_sector(mult_D3, 0, 10.0)
    ratio_hi = Z_D3_hi / Z_O_hi
    print(f"\n  High-T ratio Z(chi_0)/Z(A0) -> {ratio_hi:.2f}  "
          f"(expected: |2O|/|2D3| = {G_O//G_D3})")
    assert abs(ratio_hi - G_O / G_D3) < 0.1; n_asserts += 1

    # A1(O) vs chi_2(D3)
    print(f"\n  Excited sector: A1(O) vs chi_2(D3)")
    print(f"  {'T [K]':>7} {'Z(A1,O)':>10} {'Z(chi2,D3)':>12} {'ratio':>7}")
    print("  " + "-" * 40)
    for T in temps:
        Z_O = Z_sector(mult_O, 1, T)
        Z_D3 = Z_sector(mult_D3, a1_match, T)
        ratio = Z_D3 / Z_O if Z_O > 1e-10 else float('inf')
        if Z_O > 1e-10:
            print(f"  {T:7.2f} {Z_O:10.4f} {Z_D3:12.4f} {ratio:7.2f}")
        else:
            print(f"  {T:7.2f} {Z_O:10.2e} {Z_D3:12.4f}     {'---':>4}")

    # High-T ratio for A1 sector: same Weyl law
    Z_A1_hi = Z_sector(mult_O, 1, 10.0)
    Z_chi2_hi = Z_sector(mult_D3, a1_match, 10.0)
    ratio_A1_hi = Z_chi2_hi / Z_A1_hi
    print(f"\n  High-T ratio Z(chi_2)/Z(A1) -> {ratio_A1_hi:.2f}  "
          f"(expected: {G_O//G_D3})")
    assert abs(ratio_A1_hi - G_O / G_D3) < 0.1; n_asserts += 1

    # Spinorial sector content
    print(f"\n  Spinorial sector chi_1(S) (new, no O ancestor):")
    print(f"  {'T [K]':>7} {'Z(chi1)':>10} {'C(chi1)':>10}")
    print("  " + "-" * 30)
    for T in temps:
        Z_s = Z_sector(mult_D3, 1, T)
        C_s = C_sector(mult_D3, 1, T)
        print(f"  {T:7.2f} {Z_s:10.4f} {C_s:10.4f}")

    # ===========================================================
    # PART 4: CONVERGENCE CHECK
    # ===========================================================
    print("\n" + "=" * 70)
    print("PART 4: CONVERGENCE CHECK (J_MAX = %d)" % J_MAX)
    print("=" * 70)

    # Compare Z at J_MAX vs J_MAX-10 to verify truncation is negligible
    mult_O_low = compute_multiplicities(chars_O, cls_O, sz_O, G_O, j_max=J_MAX-10)
    mult_D3_low = compute_multiplicities(chars_D3, cls_D3, sz_D3, G_D3, j_max=J_MAX-10)
    T_check = 10.0
    Z_O_full = Z_sector(mult_O, 0, T_check)
    Z_O_trunc = Z_sector(mult_O_low, 0, T_check)
    trunc_err = abs(Z_O_full - Z_O_trunc) / Z_O_full
    print(f"\n  Z(A0,O) at T={T_check}K: J_MAX={J_MAX}: {Z_O_full:.4f}, "
          f"J_MAX={J_MAX-10}: {Z_O_trunc:.4f}")
    print(f"  Truncation error: {trunc_err:.2e}")
    assert trunc_err < 0.01, f"Truncation error {trunc_err} too large"; n_asserts += 1

    Z_D3_full = Z_sector(mult_D3, 0, T_check)
    Z_D3_trunc = Z_sector(mult_D3_low, 0, T_check)
    trunc_err_D3 = abs(Z_D3_full - Z_D3_trunc) / Z_D3_full
    print(f"  Z(chi0,D3) at T={T_check}K: J_MAX={J_MAX}: {Z_D3_full:.4f}, "
          f"J_MAX={J_MAX-10}: {Z_D3_trunc:.4f}")
    print(f"  Truncation error: {trunc_err_D3:.2e}")
    assert trunc_err_D3 < 0.01; n_asserts += 1

    # ===========================================================
    # PART 5: SCHOTTKY PEAK VERIFICATION
    # ===========================================================
    print("\n" + "=" * 70)
    print("PART 5: SCHOTTKY PEAK VERIFICATION")
    print("=" * 70)

    # Scan C(T) on fine grid to find peak location
    T_fine = np.linspace(0.05, 5.0, 500)
    C_A0_fine = [C_sector(mult_O, 0, T) for T in T_fine]
    C_chi0_fine = [C_sector(mult_D3, 0, T) for T in T_fine]

    peak_A0_idx = np.argmax(C_A0_fine)
    peak_chi0_idx = np.argmax(C_chi0_fine)
    T_peak_A0 = T_fine[peak_A0_idx]
    T_peak_chi0 = T_fine[peak_chi0_idx]
    C_peak_A0 = C_A0_fine[peak_A0_idx]
    C_peak_chi0 = C_chi0_fine[peak_chi0_idx]

    print(f"\n  A0(O) Schottky peak: T = {T_peak_A0:.3f} K, C/kB = {C_peak_A0:.4f}")
    print(f"  chi_0(D3) Schottky peak: T = {T_peak_chi0:.3f} K, C/kB = {C_peak_chi0:.4f}")
    print(f"  Peak shift ratio: {T_peak_A0/T_peak_chi0:.2f}x")

    # Schottky peak for a multi-level system: T_peak ~ 0.3-0.5 * gap_in_K
    # (shifted from the two-level 0.42*gap by degeneracy and higher levels)
    gap_A0_K = gap_A0_O * B_K    # 20 * 0.131 = 2.62 K
    gap_chi0_K = gap_chi0_D3 * B_K  # 6 * 0.131 = 0.786 K
    assert 0.25 * gap_A0_K < T_peak_A0 < 0.5 * gap_A0_K; n_asserts += 1
    assert 0.25 * gap_chi0_K < T_peak_chi0 < 0.65 * gap_chi0_K; n_asserts += 1
    # Peak ordering tracks gap ordering
    assert T_peak_A0 > T_peak_chi0; n_asserts += 1
    # Peak ratio should track gap ratio (3.3x) within a factor
    peak_ratio = T_peak_A0 / T_peak_chi0
    gap_ratio = gap_A0_O / gap_chi0_D3
    print(f"  Gap ratio: {gap_ratio:.1f}x, peak ratio: {peak_ratio:.2f}x")
    assert 1.5 < peak_ratio < 5.0; n_asserts += 1

    # ===========================================================
    # PART 6: PHYSICAL SUMMARY
    # ===========================================================
    print("\n" + "=" * 70)
    print("PHYSICAL SUMMARY")
    print("=" * 70)

    # Compute j values forbidden by O but allowed by D3 (from data, not hardcoded)
    j_A0 = set(j for (ci, j), m in mult_O.items() if ci == 0 and j <= 12)
    j_chi0 = set(j for (ci, j), m in mult_D3.items() if ci == 0 and j <= 12)
    new_j_chi0 = sorted(j_chi0 - j_A0)
    j_A1 = set(j for (ci, j), m in mult_O.items() if ci == 1 and j <= 12)
    j_chi2 = set(j for (ci, j), m in mult_D3.items() if ci == a1_match and j <= 12)
    new_j_chi2 = sorted(j_chi2 - j_A1)

    new_j_chi0_str = ', '.join(str(int(j)) for j in new_j_chi0)
    new_j_chi2_str = ', '.join(str(int(j)) for j in new_j_chi2)

    print(f"""
Symmetry reduction O -> D3 (energy scale: SF6, B = {B_K:.4f} K).
SF6 is an octahedral rotor used for the energy scale. A material
with actual O -> D3 distortion (e.g., perovskite) would be needed
to observe these effects.

1. SECTOR CORRESPONDENCE (via character restriction)
   A0(2O) -> chi_0(2D3)  [trivial -> trivial]
   A1(2O) -> chi_2(2D3)  [sign -> sign]
   chi_1, chi_3: spinorial sectors with no O ancestor

   Sectors are superselected. A system in A0(O) remains in chi_0(D3).
   Spinorial sectors exist but require spinorial boundary conditions.

2. SPECTRAL DENSIFICATION
   Ground sector gap: {gap_A0_O:.0f}B -> {gap_chi0_D3:.0f}B = {gap_A0_O*B_K:.3f} -> {gap_chi0_D3*B_K:.3f} K  ({gap_A0_O/gap_chi0_D3:.1f}x)
   Excited sector gap: {gap_A1_O:.0f}B -> {gap_chi2_D3:.0f}B = {gap_A1_O*B_K:.3f} -> {gap_chi2_D3*B_K:.3f} K  ({gap_A1_O/gap_chi2_D3:.1f}x)

   Levels forbidden by O symmetry, allowed under D3:
     chi_0: j = {new_j_chi0_str} ({len(new_j_chi0)} levels below j=12)
     chi_2: j = {new_j_chi2_str} ({len(new_j_chi2)} levels below j=12)

3. SECTOR AVAILABILITY UNDER SYMMETRY REDUCTION
   On SO(3)/O: no spinorial scalar sector exists (topological obstruction).
   On SO(3)/D3: two spinorial scalar sectors exist, ground state at {first_spin*(first_spin+1):.2f}B.

   A system prepared with spinorial boundary conditions
   CAN propagate on SO(3)/D3 as a scalar field. On SO(3)/O it cannot.

4. THERMODYNAMIC SIGNATURES
   Schottky peak (ground sector): {T_peak_A0:.2f} K (O) -> {T_peak_chi0:.2f} K (D3)
   High-T ratio: Z(chi_0)/Z(A0) -> {G_O//G_D3} (Weyl law: |2O|/|2D3|)

   For the spinorial sector (if accessible):
   - Ground state at j = 3/2, E = {first_spin*(first_spin+1):.2f}B = {first_spin*(first_spin+1)*B_K:.3f} K
   - Half-integer j spectrum: j = 3/2, 5/2, 7/2, ...
""")

    print(f"All checks passed ({n_asserts} asserts).")


if __name__ == '__main__':
    main()
