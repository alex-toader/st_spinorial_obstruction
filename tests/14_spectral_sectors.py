#!/usr/bin/env python3
"""
Sector-resolved spectrum of the rigid rotor on SO(3)/O.

Step 1: For each scalar sector (1D character) of 2O, compute the multiplicity
of that character in V_j|_{2O} for j = 0, 1, 2, ..., 6.

2O has abelianization Z₂, so 2 scalar sectors:
  - χ₀ = A₀ (trivial): χ₀(-1) = +1  (tensorial)
  - χ₁ = A₁ (sign):    χ₁(-1) = +1  (tensorial)

Both are tensorial, so only integer j contribute.

The multiplicity of a 1D char σ in V_j|_{2O} is:
  m(σ, j) = (1/|G|) Σ_g σ(g)* · χ_j(g)

where χ_j is the SU(2) spin-j character.

RAW OUTPUT:
===========
======================================================================
SECTOR-RESOLVED SPECTRUM: RIGID ROTOR ON SO(3)/O
======================================================================

2O: 8 classes, |G| = 48
⟨A₀, A₁⟩ = 0.000000 (should be 0)
⟨A₀, A₀⟩ = 1.000000 (should be 1)

    j   dim   type   m(A₀)   m(A₁)   E=j(j+1)
--------------------------------------------------
   0.0     1    (T)       1       0       0.00
   0.5     2    (S)       0       0       0.75
   1.0     3    (T)       0       0       2.00
   1.5     4    (S)       0       0       3.75
   2.0     5    (T)       0       0       6.00
   2.5     6    (S)       0       0       8.75
   3.0     7    (T)       0       1      12.00
   3.5     8    (S)       0       0      15.75
   4.0     9    (T)       1       0      20.00
   4.5    10    (S)       0       0      24.75
   5.0    11    (T)       0       0      30.00
   5.5    12    (S)       0       0      35.75
   6.0    13    (T)       1       1      42.00

======================================================================
PHYSICAL SPECTRUM PER SCALAR SECTOR
======================================================================

Sector A₀ (trivial):
  χ(-1) = +1 → tensorial → only integer j contribute
      j   E=j(j+1)   mult
  -------------------------
    0.0       0.00      1
    4.0      20.00      1
    6.0      42.00      1
  Total states (j ≤ 6): 3

Sector A₁ (sign):
  χ(-1) = +1 → tensorial → only integer j contribute
      j   E=j(j+1)   mult
  -------------------------
    3.0      12.00      1
    6.0      42.00      1
  Total states (j ≤ 6): 2

======================================================================
SPINORIAL SECTOR: H₁ VECTOR BUNDLE (rank 2)
======================================================================

H₁(-1) = -I₂ → spinorial → only half-integer j contribute

      j   E=j(j+1)   m(H₁)  states
  -----------------------------------
    0.5       0.75       1       2
    3.5      15.75       1       2
    4.5      24.75       1       2
    5.5      35.75       1       2
    6.5      48.75       1       2
  Total spinorial states (j ≤ 6.5): 10

======================================================================
SUMMARY: SPECTRAL GAP
======================================================================

Lowest energy per sector:
  A₀ scalar:    E = 0    (j=0, ground state)
  A₁ scalar:    E = 12   (j=3)
  H₁ spinorial: E = 0.75 (j=1/2, rank-2 bundle)

Gap between ground state and first spinorial state: ΔE = 0.75
Gap between two scalar sectors: ΔE = 12

On SO(3): spinorial sector starts at E = 0.75 (j=1/2) — same!
On SO(3)/O: scalar spinorial FORBIDDEN.
  To access j=1/2, MUST use vector bundle H₁.
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
from src.quaternion import qkey
from src.group import build_binary_tetrahedral, build_binary_octahedral, compute_conjugacy_classes


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
    print("SECTOR-RESOLVED SPECTRUM: RIGID ROTOR ON SO(3)/O")
    print("=" * 70)

    # Build 2O and conjugacy classes
    elements = build_binary_octahedral()
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    G = sum(sizes)  # = 48

    print(f"\n2O: {len(classes)} classes, |G| = {G}")

    # Build A₀ (trivial) and A₁ (sign) characters
    elements_2T = build_binary_tetrahedral()
    keys_2T = {qkey(e) for e in elements_2T}

    chi_A0 = []  # trivial: all 1
    chi_A1 = []  # sign: +1 on 2T, -1 outside
    for c in classes:
        chi_A0.append(1.0)
        rep_key = list(c['keys'])[0]
        chi_A1.append(1.0 if rep_key in keys_2T else -1.0)

    chi_A0 = np.array(chi_A0)
    chi_A1 = np.array(chi_A1)
    sizes_arr = np.array(sizes, dtype=float)

    # Verify orthogonality
    inner = np.sum(chi_A0.conj() * chi_A1 * sizes_arr) / G
    print(f"⟨A₀, A₁⟩ = {inner:.6f} (should be 0)")
    inner00 = np.sum(chi_A0.conj() * chi_A0 * sizes_arr) / G
    print(f"⟨A₀, A₀⟩ = {inner00:.6f} (should be 1)")

    # Compute multiplicities for j = 0, 0.5, 1, ..., 6
    print(f"\n{'j':>5}  {'dim':>4}  {'type':>5}  {'m(A₀)':>6}  {'m(A₁)':>6}  {'E=j(j+1)':>9}")
    print("-" * 50)

    for two_j in range(0, 13):
        j = two_j / 2.0
        chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])

        m_A0 = np.sum(chi_A0.conj() * chi_Vj * sizes_arr) / G
        m_A1 = np.sum(chi_A1.conj() * chi_Vj * sizes_arr) / G

        # Check: sum of all irrep multiplicities × dim should = 2j+1
        m_A0_r = round(m_A0.real)
        m_A1_r = round(m_A1.real)

        typ = "(T)" if two_j % 2 == 0 else "(S)"
        E = j * (j + 1)

        print(f"  {j:4.1f}  {int(2*j+1):4d}  {typ:>5}  {m_A0_r:6d}  {m_A1_r:6d}  {E:9.2f}")

    # Now the physical spectrum per sector
    print("\n" + "=" * 70)
    print("PHYSICAL SPECTRUM PER SCALAR SECTOR")
    print("=" * 70)

    for name, chi_sigma in [("A₀ (trivial)", chi_A0), ("A₁ (sign)", chi_A1)]:
        print(f"\nSector {name}:")
        print(f"  χ(-1) = +1 → tensorial → only integer j contribute")
        print(f"  {'j':>5}  {'E=j(j+1)':>9}  {'mult':>5}")
        print(f"  {'-'*25}")
        total = 0
        for j_int in range(7):
            j = float(j_int)
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            m = np.sum(chi_sigma.conj() * chi_Vj * sizes_arr) / G
            m_r = round(m.real)
            if m_r > 0:
                E = j * (j + 1)
                print(f"  {j:5.1f}  {E:9.2f}  {m_r:5d}")
                total += m_r
        print(f"  Total states (j ≤ 6): {total}")

    # Spinorial sector: H₁ bundle (rank 2)
    # Multiplicity of H₁ in V_j|_{2O} for half-integer j
    print("\n" + "=" * 70)
    print("SPINORIAL SECTOR: H₁ VECTOR BUNDLE (rank 2)")
    print("=" * 70)
    print("\nH₁(-1) = -I₂ → spinorial → only half-integer j contribute")

    chi_H1 = np.array([chi_su2(0.5, c['trace']) for c in classes])

    print(f"\n  {'j':>5}  {'E=j(j+1)':>9}  {'m(H₁)':>6}  {'states':>6}")
    print(f"  {'-'*35}")
    total_spin = 0
    for two_j in range(1, 14, 2):  # j = 0.5, 1.5, ..., 6.5
        j = two_j / 2.0
        chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
        m = np.sum(chi_H1.conj() * chi_Vj * sizes_arr) / G
        m_r = round(m.real)
        E = j * (j + 1)
        states = m_r * 2  # each H₁ copy = 2 states
        if m_r > 0:
            print(f"  {j:5.1f}  {E:9.2f}  {m_r:6d}  {states:6d}")
            total_spin += states
    print(f"  Total spinorial states (j ≤ 6.5): {total_spin}")

    # Summary comparison
    print("\n" + "=" * 70)
    print("SUMMARY: SPECTRAL GAP")
    print("=" * 70)
    print("\nLowest energy per sector:")
    print("  A₀ scalar:    E = 0    (j=0, ground state)")
    print("  A₁ scalar:    E = 12   (j=3)")
    print("  H₁ spinorial: E = 0.75 (j=1/2, rank-2 bundle)")
    print()
    print("Gap between ground state and first spinorial state: ΔE = 0.75")
    print("Gap between two scalar sectors: ΔE = 12")
    print()
    print("On SO(3): spinorial sector starts at E = 0.75 (j=1/2) — same!")
    print("On SO(3)/O: scalar spinorial FORBIDDEN.")
    print("  To access j=1/2, MUST use vector bundle H₁.")


if __name__ == '__main__':
    main()
