#!/usr/bin/env python3
"""
Sector-resolved thermal predictions for SF₆-scale octahedral rotor.

Physical setup: rigid rotor on SO(3)/H with H = O (octahedral).
SF₆ rotational constant B = 0.09107 cm⁻¹ = 0.1310 K.
(SF₆ sets the energy scale; actual O→D₃ distortion requires a
trigonally distorted octahedral cage, e.g. perovskite.)

In sector χ, the partition function is
  Z_χ(T) = Σ_j m(χ,j) · (2j+1) · exp(-B·j(j+1)/T)
and the orientational specific heat is
  C_χ(T)/k_B = (1/T²) · [<E²>_χ - <E>_χ²].

Sector populations are fixed by preparation (superselection), not
by thermal equilibrium across sectors. Each C_χ(T) is the specific
heat of a system prepared in sector χ.

Key predictions (SF₆ energy scale):
  A₀(O):  gap = 20B = 2.62 K, Schottky peak at T ≈ 0.83 K
  χ₀(D₃): gap = 6B  = 0.79 K, Schottky peak at T ≈ 0.34 K
  χ₁(D₃): gap = 3.75B = 0.49 K (spinorial, forbidden on SO(3)/O)
           C(T) rises monotonically to 3/2 without Schottky overshoot

The spinorial sector becomes thermally active below ~0.5 K,
providing a thermal signature with no counterpart in the
undistorted octahedral phase.

Asserts: 17
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.quaternion import qkey, qmul, qconj
from src.group import (build_binary_octahedral,
                       compute_conjugacy_classes, compute_commutator_subgroup)

# ==========================================================
# Physical constants
# ==========================================================
B_CM = 0.09107      # SF₆ rotational constant [cm⁻¹]
CM_TO_K = 1.4388    # 1 cm⁻¹ = 1.4388 K
B_K = B_CM * CM_TO_K  # B in Kelvin

# ==========================================================
# Infrastructure (from file 36, minimal)
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
    """Build Dic_N = 2D_N of order 4N."""
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


def compute_multiplicities(characters, classes, sizes, G, j_max=50):
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

    elements_O = np.array(build_binary_octahedral())
    elements_D3 = build_dicyclic(3)

    chars_O, cls_O, sz_O, G_O, obst_O = compute_1d_characters(elements_O)
    chars_D3, cls_D3, sz_D3, G_D3, obst_D3 = compute_1d_characters(elements_D3)

    J_MAX = 50
    mult_O = compute_multiplicities(chars_O, cls_O, sz_O, G_O, j_max=J_MAX)
    mult_D3 = compute_multiplicities(chars_D3, cls_D3, sz_D3, G_D3, j_max=J_MAX)

    assert obst_O and not obst_D3; n_asserts += 1

    # ==========================================================
    # D2: SCHOTTKY ANOMALY PER SECTOR
    # ==========================================================
    print("=" * 70)
    print("D2: SECTOR-RESOLVED SCHOTTKY ANOMALY (SF6 energy scale)")
    print("=" * 70)

    print(f"\n  B = {B_CM} cm-1 = {B_K:.4f} K")

    # Energy scales
    gaps = {
        'A0(O)': 20,      # j=4: 4*5=20
        'A1(O)': 12,      # j=3: 3*4=12
        'chi0(D3)': 6,    # j=2: 2*3=6
        'chi2(D3)': 2,    # j=1: 1*2=2
        'chi1(D3)': 3.75, # j=3/2: 3/2*5/2=3.75
    }

    print("\n  Spectral gaps in physical units:")
    print(f"  {'Sector':>12} {'gap/B':>7} {'gap [cm-1]':>12} {'gap [K]':>8}")
    print("  " + "-" * 45)
    for name, g in gaps.items():
        print(f"  {name:>12} {g:7.2f} {g*B_CM:12.5f} {g*B_K:8.3f}")

    # Verify gaps from data
    j1_A0 = min(j for (ci, j), m in mult_O.items() if ci == 0 and j > 0)
    j1_chi0 = min(j for (ci, j), m in mult_D3.items() if ci == 0 and j > 0)
    j1_spin = min(j for (ci, j), m in mult_D3.items()
                  if chars_D3[ci]['is_spinorial'])
    assert j1_A0 * (j1_A0 + 1) == 20; n_asserts += 1
    assert j1_chi0 * (j1_chi0 + 1) == 6; n_asserts += 1
    assert abs(j1_spin * (j1_spin + 1) - 3.75) < 1e-10; n_asserts += 1

    # Schottky peaks via fine scan
    T_fine = np.linspace(0.05, 3.0, 1000)
    C_A0 = [C_sector(mult_O, 0, T) for T in T_fine]
    C_chi0 = [C_sector(mult_D3, 0, T) for T in T_fine]

    # Spinorial sector: chi_1 (index 1)
    C_chi1 = [C_sector(mult_D3, 1, T) for T in T_fine]

    peak_A0 = T_fine[np.argmax(C_A0)]
    peak_chi0 = T_fine[np.argmax(C_chi0)]
    Cmax_A0 = max(C_A0)
    Cmax_chi0 = max(C_chi0)

    print(f"\n  Schottky peak (ground-state sector):")
    print(f"    A0(O):     T_peak = {peak_A0:.2f} K,  C_max/kB = {Cmax_A0:.3f}")
    print(f"    chi0(D3):  T_peak = {peak_chi0:.2f} K,  C_max/kB = {Cmax_chi0:.3f}")
    print(f"    Ratio: {peak_A0/peak_chi0:.2f}x")

    # Schottky peak for A1(O)
    C_A1 = [C_sector(mult_O, 1, T) for T in T_fine]
    peak_A1 = T_fine[np.argmax(C_A1)]
    Cmax_A1 = max(C_A1)
    print(f"    A1(O):     T_peak = {peak_A1:.2f} K,  C_max/kB = {Cmax_A1:.3f}")

    # Verify Schottky peaks are in expected range
    assert 0.6 < peak_A0 < 1.2; n_asserts += 1   # ~0.83 K
    assert 0.2 < peak_chi0 < 0.6; n_asserts += 1  # ~0.34 K

    # Key prediction: thermal contrast at 0.5 K
    T_probe = 0.5
    C_A0_probe = C_sector(mult_O, 0, T_probe)
    C_chi0_probe = C_sector(mult_D3, 0, T_probe)
    C_chi1_probe = C_sector(mult_D3, 1, T_probe)
    Z_chi1_probe = Z_sector(mult_D3, 1, T_probe)

    print(f"\n  Thermal contrast at T = {T_probe} K:")
    print(f"    C(A0,O)/kB   = {C_A0_probe:.4f}  (still climbing)")
    print(f"    C(chi0,D3)/kB = {C_chi0_probe:.4f}  (past peak)")
    print(f"    C(chi1,D3)/kB = {C_chi1_probe:.4f}  (spinorial, new)")
    print(f"    Z(chi1,D3)   = {Z_chi1_probe:.4f}  (spinorial active)")

    # The spinorial onset temperature: where Z > 1.1 (10% thermal activation)
    T_onset_idx = next(i for i, T in enumerate(T_fine)
                       if Z_sector(mult_D3, 1, T) > 1.1)
    T_onset = T_fine[T_onset_idx]
    print(f"\n  Spinorial onset (Z > 1.1): T = {T_onset:.2f} K")
    print(f"    = {T_onset/B_K:.1f} B")
    assert T_onset < 1.0; n_asserts += 1  # well below 1 K

    # Verify chi_1 has no Schottky overshoot (monotonic approach to 3/2)
    assert max(C_chi1) < 1.505; n_asserts += 1

    # Summary table: C/kB at key temperatures
    print(f"\n  Specific heat comparison C/kB:")
    print(f"  {'T [K]':>7} {'A0(O)':>8} {'chi0(D3)':>10} {'chi1(D3,S)':>12}")
    print("  " + "-" * 42)
    for T in [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
        c_a0 = C_sector(mult_O, 0, T)
        c_chi0 = C_sector(mult_D3, 0, T)
        c_chi1 = C_sector(mult_D3, 1, T)
        print(f"  {T:7.1f} {c_a0:8.4f} {c_chi0:10.4f} {c_chi1:12.4f}")

    # ==========================================================
    # VERIFY D3 MULTIPLICITIES AGAINST PAPER TABLE
    # ==========================================================
    print("\n" + "=" * 70)
    print("D3 MULTIPLICITY VERIFICATION (paper §4.2)")
    print("=" * 70)

    # Paper table: chi_0(T), chi_1(S), chi_2(T), chi_3(S)
    expected_D3 = {
        (0, 0.0): 1, (2, 0.0): 0,  # chi_0: j=0 → 1, j=0 (chi_2) → 0
        (2, 1.0): 1,                # chi_2: j=1 → 1
        (1, 1.5): 1, (3, 1.5): 1,  # chi_1, chi_3: j=3/2 → 1
        (0, 2.0): 1, (2, 2.0): 0,  # chi_0: j=2 → 1, chi_2: j=2 → 0
        (1, 2.5): 1, (3, 2.5): 1,  # chi_1, chi_3: j=5/2 → 1
    }
    for (ci, j), expected in expected_D3.items():
        actual = mult_D3.get((ci, j), 0)
        name = chars_D3[ci]['name']
        ok = "✓" if actual == expected else "✗"
        print(f"  {ok} m({name}, {j}) = {actual} (expected {expected})")
        assert actual == expected, f"m({name},{j})={actual}, expected {expected}"
        n_asserts += 1

    print(f"All checks passed ({n_asserts} asserts).")


if __name__ == '__main__':
    main()
