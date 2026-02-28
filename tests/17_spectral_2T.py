#!/usr/bin/env python3
"""
Sector-resolved spectrum of the rigid rotor on SO(3)/T.

2T (binary tetrahedral, order 24) is obstructed: -1 ∈ [2T,2T] = Q₈.
Abel(2T) ≅ Z₃, so 3 scalar sectors — all tensorial.

Characters: χ_k(a) = ω^k where ω = e^{2πi/3}, k = 0,1,2.
Since χ_k(-1) = χ_k(a^? ) — need to check which power of generator gives -1.

Compare with 2O (Abel ≅ Z₂, 2 sectors) and 2D₃ (Abel ≅ Z₄, 4 sectors).

RAW OUTPUT:
===========
======================================================================
SECTOR-RESOLVED SPECTRUM: RIGID ROTOR ON SO(3)/T
======================================================================

|2T| = 24
|[2T,2T]| = 8
Abel(2T) ≅ Z_{3}
-1 ∈ [2T,2T]: True → obstructed

2T: 7 classes, |G| = 24

Character values at -1:
  χ_0(-1) = 1.0000+0.0000j → tensorial
  χ_1(-1) = 1.0000+0.0000j → tensorial
  χ_2(-1) = 1.0000+0.0000j → tensorial

Orthogonality: all ⟨χ_i, χ_j⟩ = δ_ij ✓

    j   dim   type     χ_0     χ_1     χ_2          E
-----------------------------------------------------
   0.0     1    (T)       1       0       0       0.00
   0.5     2    (S)       0       0       0       0.75
   1.0     3    (T)       0       0       0       2.00
   1.5     4    (S)       0       0       0       3.75
   2.0     5    (T)       0       1       1       6.00
   2.5     6    (S)       0       0       0       8.75
   3.0     7    (T)       1       0       0      12.00
   3.5     8    (S)       0       0       0      15.75
   4.0     9    (T)       1       1       1      20.00
   4.5    10    (S)       0       0       0      24.75
   5.0    11    (T)       0       1       1      30.00
   5.5    12    (S)       0       0       0      35.75
   6.0    13    (T)       2       1       1      42.00
   6.5    14    (S)       0       0       0      48.75

Sector χ_0 (trivial):   j = 0, 3, 4, 6     (5 states)
Sector χ_1 (ω):         j = 2, 4, 5, 6     (4 states)
Sector χ_2 (ω²):        j = 2, 4, 5, 6     (4 states)
Total scalar states (j ≤ 6): 13

Comparison: 2T(13 states) vs 2O(5 states) vs 2D₃(35 states)
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
from src.quaternion import qkey, qmul
from src.group import (build_binary_tetrahedral, compute_conjugacy_classes,
                       compute_commutator_subgroup)


def chi_su2(j, trace):
    """SU(2) spin-j character at given trace."""
    d = int(2 * j + 1)
    if abs(trace - 2) < 1e-10:
        return float(d)
    if abs(trace + 2) < 1e-10:
        return float(d * ((-1) ** int(round(2 * j))))
    ha = np.arccos(np.clip(trace / 2, -1, 1))
    if abs(np.sin(ha)) < 1e-14:
        return float(d)
    return np.sin(d * ha) / np.sin(ha)


def main():
    print("=" * 70)
    print("SECTOR-RESOLVED SPECTRUM: RIGID ROTOR ON SO(3)/T")
    print("=" * 70)

    elements = build_binary_tetrahedral()
    G_order = len(elements)

    # Commutator subgroup
    comm = compute_commutator_subgroup(elements)
    comm_keys = set(comm.keys())
    mk = qkey(np.array([-1, 0, 0, 0], dtype=float))
    print(f"\n|2T| = {G_order}")
    print(f"|[2T,2T]| = {len(comm)}")
    print(f"Abel(2T) ≅ Z_{{{G_order // len(comm)}}}")
    print(f"-1 ∈ [2T,2T]: {mk in comm_keys} → obstructed")

    abel_size = G_order // len(comm)  # = 3

    # Conjugacy classes
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    G = sum(sizes)
    sizes_arr = np.array(sizes, dtype=float)

    print(f"\n2T: {len(classes)} classes, |G| = {G}")
    for i, c in enumerate(classes):
        in_comm = list(c['keys'])[0] in comm_keys
        print(f"  C{i}: size={c['size']}, trace={c['trace']:.4f}, in [G,G]: {in_comm}")

    # Build quotient map G → Z₃
    # Find generator of Z₃ = G/[G,G]
    # Need element NOT in [G,G] with order 3 in quotient

    # Coset identification
    cosets = []
    elem_coset = {}
    for e in elements:
        k = qkey(e)
        if k in elem_coset:
            continue
        coset = set()
        for ck, cv in comm.items():
            prod = qmul(e, cv)
            coset.add(qkey(prod))
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

    print(f"\nCosets: {len(cosets)} (should be {abel_size})")

    # Find generator
    gen_elem = None
    for e in elements:
        k = qkey(e)
        if k in comm_keys:
            continue  # in [G,G], maps to 0
        # Check if e has order 3 in quotient (e³ ∈ [G,G])
        e2 = qmul(e, e)
        e3 = qmul(e2, e)
        if qkey(e3) in comm_keys:
            # Check e² ∉ [G,G]
            if qkey(e2) not in comm_keys:
                gen_elem = e.copy()
                break

    if gen_elem is None:
        print("ERROR: no generator found!")
        return

    print(f"Generator of Z₃ quotient: {gen_elem}, trace={2*gen_elem[0]:.4f}")

    # Assign Z₃ labels
    coset_label = {}
    power = np.array([1, 0, 0, 0], dtype=float)
    for d in range(abel_size):
        ci = elem_coset[qkey(power)]
        coset_label[ci] = d
        power = qmul(power, gen_elem)

    # Build 3 characters: χ_k(g) = ω^{k·label(g)}
    omega = np.exp(2j * np.pi / 3)

    # Character values per class
    class_labels = []
    for c in classes:
        rep_key = list(c['keys'])[0]
        ci = elem_coset[rep_key]
        class_labels.append(coset_label[ci])
    class_labels = np.array(class_labels)

    # Check character values at -1
    print("\nCharacter values at -1:")
    for i, c in enumerate(classes):
        if abs(c['trace'] + 2) < 1e-10:
            for k in range(abel_size):
                val = omega ** (k * class_labels[i])
                typ = "tensorial" if val.real > 0.5 else "SPINORIAL"
                print(f"  χ_{k}(-1) = {val:.4f} → {typ}")
            break

    # Verify orthogonality
    print("\nOrthogonality:")
    for k1 in range(abel_size):
        for k2 in range(k1, abel_size):
            chi1 = omega ** (k1 * class_labels)
            chi2 = omega ** (k2 * class_labels)
            inner = np.sum(chi1.conj() * chi2 * sizes_arr) / G
            if abs(inner) > 1e-10 or k1 == k2:
                expected = 1.0 if k1 == k2 else 0.0
                status = "✓" if abs(inner - expected) < 1e-10 else "✗"
                print(f"  ⟨χ_{k1}, χ_{k2}⟩ = {inner.real:.4f} {status}")

    # Full multiplicity table
    print("\n" + "=" * 70)
    print("FULL MULTIPLICITY TABLE")
    print("=" * 70)

    header = f"{'j':>5}  {'dim':>4}  {'type':>5}"
    for k in range(abel_size):
        header += f"  {'χ_'+str(k):>6}"
    header += f"  {'E':>9}"
    print(f"\n{header}")
    print("-" * len(header))

    for two_j in range(0, 14):
        j = two_j / 2.0
        chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
        typ = "(T)" if two_j % 2 == 0 else "(S)"
        E = j * (j + 1)

        line = f"  {j:4.1f}  {int(2*j+1):4d}  {typ:>5}"
        for k in range(abel_size):
            chi_k = omega ** (k * class_labels)
            m = np.sum(chi_k.conj() * chi_Vj * sizes_arr) / G
            line += f"  {round(m.real):6d}"
        line += f"  {E:9.2f}"
        print(line)

    # Physical spectrum per sector
    print("\n" + "=" * 70)
    print("PHYSICAL SPECTRUM PER SECTOR (all tensorial)")
    print("=" * 70)

    for k in range(abel_size):
        chi_k_vals = omega ** (k * class_labels)
        # All tensorial → integer j only
        label = "trivial" if k == 0 else f"ω^{k}"
        print(f"\nSector χ_{k} ({label}):")
        print(f"  {'j':>5}  {'E=j(j+1)':>9}  {'mult':>5}")
        print(f"  {'-'*25}")
        total = 0
        for j_int in range(7):
            j = float(j_int)
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            m = np.sum(chi_k_vals.conj() * chi_Vj * sizes_arr) / G
            m_r = round(m.real)
            if m_r > 0:
                E = j * (j + 1)
                print(f"  {j:5.1f}  {E:9.2f}  {m_r:5d}")
                total += m_r
        print(f"  Total states (j ≤ 6): {total}")

    # Comparison
    print("\n" + "=" * 70)
    print("COMPARISON: SO(3)/T vs SO(3)/O vs SO(3)/D₃")
    print("=" * 70)
    print("""
SO(3)/T (obstructed, Abel ≅ Z₃):
  3 scalar sectors, ALL tensorial (integer j only)
  More sectors than O, but still no spinorial scalars
  κ = 3 (most scalar sectors among obstructed groups)

SO(3)/O (obstructed, Abel ≅ Z₂):
  2 scalar sectors, both tensorial
  κ = 2

SO(3)/D₃ (non-obstructed, Abel ≅ Z₄):
  4 scalar sectors: 2 tensorial + 2 spinorial
  κ = 4, with spinorial scalars at j = 3/2, 5/2, ...

KEY: More scalar sectors ≠ spinorial access.
  T has MORE sectors than O (3 vs 2) but STILL no spinorial.
  D₃ has 4 sectors INCLUDING spinorial — obstruction is the criterion.""")


if __name__ == '__main__':
    main()
