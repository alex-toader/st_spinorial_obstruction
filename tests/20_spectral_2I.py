#!/usr/bin/env python3
"""
Sector-resolved spectrum of the rigid rotor on SO(3)/I.

2I (binary icosahedral, order 120) is the strongest obstruction:
  - 2I is PERFECT: [2I,2I] = 2I
  - Abel(2I) ≅ Z₁ → only 1 scalar sector (trivial)
  - κ = 1: the minimum possible

Only integer j contribute, and only those where the trivial character
appears in V_j|_{2I}. Expected: extreme sparsity.

RAW OUTPUT:
===========

======================================================================
SECTOR-RESOLVED SPECTRUM: RIGID ROTOR ON SO(3)/I
======================================================================

|2I| = 120
|[2I,2I]| = 120 (perfect group: [2I,2I] = 2I)
Abel(2I) ≅ Z_{1}
-1 ∈ [2I,2I]: True → OBSTRUCTED (total)

2I: 9 classes, |G| = 120
  C0: size=1, trace=2.000000
  C1: size=12, trace=1.618034
  C2: size=20, trace=1.000000
  C3: size=12, trace=0.618034
  C4: size=30, trace=0.000000
  C5: size=12, trace=-0.618034
  C6: size=20, trace=-1.000000
  C7: size=12, trace=-1.618034
  C8: size=1, trace=-2.000000

======================================================================
SCALAR SPECTRUM (only sector: trivial)
======================================================================

      j   dim   m(A)   E=j(j+1)
  ------------------------------
    0.0     1      1       0.00
    6.0    13      1      42.00
   10.0    21      1     110.00
   12.0    25      1     156.00
   15.0    31      1     240.00

  Total scalar states (j ≤ 15): 5
  First excited state: j=6, E=42
  Spectral gap: ΔE = 42

======================================================================
VERIFICATION: half-integer j multiplicities (should all be 0)
======================================================================
  All half-integer j: m(A) = 0 ✓

======================================================================
COMPARISON: SPECTRAL SPARSITY OF OBSTRUCTED GROUPS
======================================================================

  2T (SO(3)/T): |G|=24, κ=3
    Trivial sector j values: 0, 3, 4, 6, 7, 8, 9, 10, 11, 12
    Trivial states (j≤15): 22
    First excited: j=3, ΔE=12
    Density (trivial): 22/16 = 1.375 per integer j

  2O (SO(3)/O): |G|=48, κ=2
    Trivial sector j values: 0, 4, 6, 8, 9, 10, 12, 13, 14, 15
    Trivial states (j≤15): 11
    First excited: j=4, ΔE=20
    Density (trivial): 11/16 = 0.688 per integer j

  2I (SO(3)/I): |G|=120, κ=1
    Trivial sector j values: 0, 6, 10, 12, 15
    Trivial states (j≤15): 5
    First excited: j=6, ΔE=42
    Density (trivial): 5/16 = 0.312 per integer j

======================================================================
SUMMARY: ICOSAHEDRAL = MAXIMUM OBSTRUCTION
======================================================================

2I is the most extreme case:
  - Perfect group: [2I,2I] = 2I, Abel = Z₁
  - Only ONE scalar sector (trivial)
  - κ = 1: literally the minimum possible
  - Spectral density ~ 1/120 per (2j+1)

First excited state (trivial sector):
  SO(3):   j=1  (E=2)
  SO(3)/T: j=3  (E=12)
  SO(3)/O: j=4  (E=20)
  SO(3)/I: j=6  (E=42)    ← maximum gap

The bigger the group, the sparser the scalar spectrum.
Icosahedral symmetry produces the most "silent" quantum rotor.
"""
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))

import numpy as np
from src.group import build_binary_icosahedral, compute_conjugacy_classes
from src.group import compute_commutator_subgroup
from src.quaternion import qkey


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
    print("SECTOR-RESOLVED SPECTRUM: RIGID ROTOR ON SO(3)/I")
    print("=" * 70)

    elements = build_binary_icosahedral()
    G_order = len(elements)

    # Verify
    comm = compute_commutator_subgroup(elements)
    mk = qkey(np.array([-1, 0, 0, 0], dtype=float))
    print(f"\n|2I| = {G_order}")
    print(f"|[2I,2I]| = {len(comm)} (perfect group: [2I,2I] = 2I)")
    print(f"Abel(2I) ≅ Z_{{{G_order // len(comm)}}}")
    print(f"-1 ∈ [2I,2I]: {mk in set(comm.keys())} → OBSTRUCTED (total)")

    # Conjugacy classes
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    G = int(sum(sizes))

    print(f"\n2I: {len(classes)} classes, |G| = {G}")
    for i, c in enumerate(classes):
        print(f"  C{i}: size={c['size']}, trace={c['trace']:.6f}")

    # Only 1 scalar sector: trivial character (all 1s)
    chi_triv = np.ones(len(classes))

    # Multiplicity of trivial in V_j = number of times V_j|_{2I} contains A
    # This is the number of icosahedral-invariant harmonics at angular momentum j
    print("\n" + "=" * 70)
    print("SCALAR SPECTRUM (only sector: trivial)")
    print("=" * 70)
    print(f"\n  {'j':>5}  {'dim':>4}  {'m(A)':>5}  {'E=j(j+1)':>9}")
    print(f"  {'-'*30}")

    total_states = 0
    first_excited = None
    for two_j in range(0, 31):  # j up to 15
        j = two_j / 2.0
        # Only integer j (tensorial)
        if two_j % 2 != 0:
            continue
        chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
        m = round((np.sum(chi_triv * chi_Vj * sizes) / G).real)
        if m > 0:
            E = j * (j + 1)
            print(f"  {j:5.1f}  {int(2*j+1):4d}  {m:5d}  {E:9.2f}")
            total_states += m
            if j > 0 and first_excited is None:
                first_excited = j

    print(f"\n  Total scalar states (j ≤ 15): {total_states}")
    if first_excited:
        print(f"  First excited state: j={first_excited:.0f}, E={first_excited*(first_excited+1):.0f}")
        print(f"  Spectral gap: ΔE = {first_excited*(first_excited+1):.0f}")

    # Full table including half-integer (all zero)
    print("\n" + "=" * 70)
    print("VERIFICATION: half-integer j multiplicities (should all be 0)")
    print("=" * 70)
    any_nonzero = False
    for two_j in range(1, 31, 2):
        j = two_j / 2.0
        chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
        m = round((np.sum(chi_triv * chi_Vj * sizes) / G).real)
        if m != 0:
            print(f"  j={j}: m={m} ← ERROR!")
            any_nonzero = True
    if not any_nonzero:
        print("  All half-integer j: m(A) = 0 ✓")

    # Comparison with other obstructed groups
    print("\n" + "=" * 70)
    print("COMPARISON: SPECTRAL SPARSITY OF OBSTRUCTED GROUPS")
    print("=" * 70)

    # Count states up to j=15 for each
    from src.group import build_binary_tetrahedral, build_binary_octahedral

    for name, builder in [("2T (SO(3)/T)", build_binary_tetrahedral),
                          ("2O (SO(3)/O)", build_binary_octahedral),
                          ("2I (SO(3)/I)", build_binary_icosahedral)]:
        elems = builder()
        cls = compute_conjugacy_classes(elems)
        cls.sort(key=lambda c: -c['trace'])
        sz = np.array([c['size'] for c in cls], dtype=float)
        g = int(sum(sz))
        cm = compute_commutator_subgroup(elems)
        abel = g // len(cm)

        # Count trivial multiplicities at integer j
        total = 0
        first_j = None
        j_list = []
        for j_int in range(16):
            j = float(j_int)
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in cls])
            # Sum over ALL 1D chars (= abel many)
            # For trivial: just check trivial
            # Actually count ALL scalar states = sum over all 1D chars
            # For obstructed: all 1D chars are tensorial
            omega = np.exp(2j * np.pi / abel) if abel > 1 else 1.0

            # Quick: use Peter-Weyl. Total 1D content = sum of all 1D char multiplicities
            total_at_j = 0
            for k in range(abel):
                if abel == 1:
                    chi_k = np.ones(len(cls))
                else:
                    # Need proper character construction... too complex here
                    # Just count trivial for simplicity
                    chi_k = np.ones(len(cls))
                    if k > 0:
                        continue
                m = round((np.sum(chi_k * chi_Vj * sz) / g).real)
                total_at_j += m

            # For obstructed groups with known abel:
            # 2T: abel=3, 2O: abel=2, 2I: abel=1
            # Approximate: multiply trivial count by abel
            m_triv = round((np.sum(np.ones(len(cls)) * chi_Vj * sz) / g).real)
            if m_triv > 0:
                j_list.append(f"{j_int}")
                if j_int > 0 and first_j is None:
                    first_j = j_int
            total += m_triv

        print(f"\n  {name}: |G|={g}, κ={abel}")
        print(f"    Trivial sector j values: {', '.join(j_list[:10])}")
        print(f"    Trivial states (j≤15): {total}")
        print(f"    First excited: j={first_j}, ΔE={first_j*(first_j+1) if first_j else '?'}")
        print(f"    Density (trivial): {total}/16 = {total/16:.3f} per integer j")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY: ICOSAHEDRAL = MAXIMUM OBSTRUCTION")
    print("=" * 70)
    print("""
2I is the most extreme case:
  - Perfect group: [2I,2I] = 2I, Abel = Z₁
  - Only ONE scalar sector (trivial)
  - κ = 1: literally the minimum possible
  - Spectral density ~ 1/120 per (2j+1)

First excited state (trivial sector):
  SO(3):   j=1  (E=2)
  SO(3)/T: j=3  (E=12)
  SO(3)/O: j=4  (E=20)
  SO(3)/I: j=6  (E=42)    ← maximum gap

The bigger the group, the sparser the scalar spectrum.
Icosahedral symmetry produces the most "silent" quantum rotor.
""")


if __name__ == '__main__':
    main()
