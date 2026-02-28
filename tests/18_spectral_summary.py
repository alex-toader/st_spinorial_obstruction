#!/usr/bin/env python3
"""
Spectral comparison table: SO(3) vs SO(3)/T vs SO(3)/O vs SO(3)/D₃.

Compares:
  1. Number of scalar sectors and their type (tensorial/spinorial)
  2. Scalar state count at each energy level
  3. Asymptotic density (Weyl law check)
  4. Spectral gap between first tensorial and first spinorial states

This is the central result for the "Direction B: Spectral content" paper section.

RAW OUTPUT:
===========
======================================================================
SPECTRAL COMPARISON TABLE
======================================================================

      j          E     SO(3)   SO(3)/T   SO(3)/O  SO(3)/D₃
  --------------------------------------------------------
    0.0       0.00         1         1         1         1
    0.5       0.75         2         0         0         0
    1.0       2.00         3         0         0         1
    1.5       3.75         4         0         0         2
    2.0       6.00         5         2         0         1
    2.5       8.75         6         0         0         2
    3.0      12.00         7         1         1         3
    3.5      15.75         8         0         0         2
    4.0      20.00         9         3         1         3
    4.5      24.75        10         0         0         4
    5.0      30.00        11         2         0         3
    5.5      35.75        12         0         0         4
    6.0      42.00        13         4         2         5
    6.5      48.75        14         0         0         4

                      SO(3)   SO(3)/T   SO(3)/O  SO(3)/D₃
    Total (j≤6.5)       105        13         5        35
   Tensorial only        49        13         5        17
   Spinorial only        56         0         0        18

SPECTRAL GAPS:
  SO(3):   T→ΔE=2.00 (j=1), S→ΔE=0.75 (j=0.5)
  SO(3)/T: T→ΔE=6.00 (j=2), S→FORBIDDEN
  SO(3)/O: T→ΔE=12.0 (j=3), S→FORBIDDEN
  SO(3)/D₃: T→ΔE=2.00 (j=1), S→ΔE=3.75 (j=1.5)

KEY: Obstruction decimates spectrum. 2O has 5 vs D₃'s 35 scalar states.
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
from src.quaternion import qkey, qmul
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       compute_conjugacy_classes, compute_commutator_subgroup)


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


def build_dicyclic(N):
    """Build Dic_N (order 4N) as unit quaternions."""
    elements = []
    for k in range(2 * N):
        angle = k * np.pi / N
        elements.append(np.array([np.cos(angle), 0, 0, np.sin(angle)]))
    x = np.array([0, 0, 1, 0], dtype=float)
    for k in range(2 * N):
        angle = k * np.pi / N
        ak = np.array([np.cos(angle), 0, 0, np.sin(angle)])
        elements.append(qmul(x, ak))
    return np.array(elements)


def compute_scalar_spectrum(elements, max_j=10):
    """Compute scalar sector spectrum for a group.

    Returns dict with:
      - 'sectors': list of (name, type, spectrum)
      - 'total_per_j': dict j -> total scalar states
    """
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    G = int(sum(sizes))

    comm = compute_commutator_subgroup(elements)
    comm_keys = set(comm.keys())
    mk = qkey(np.array([-1, 0, 0, 0], dtype=float))
    obstructed = mk in comm_keys
    abel_size = G // len(comm)

    # Build coset structure
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

    # Find generator of abelianization
    gen_elem = None
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

    if gen_elem is None and abel_size > 1:
        # Fallback for weird cases
        return None

    # Assign labels
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

    # Find -1 class label
    minus_one_label = None
    for i, c in enumerate(classes):
        if abs(c['trace'] + 2) < 1e-10:
            minus_one_label = class_labels[i]
            break

    # Build sectors
    sectors = []
    n_tensorial = 0
    n_spinorial = 0
    for k in range(abel_size):
        chi_k = omega ** (k * class_labels) if abel_size > 1 else np.ones(len(classes))
        chi_minus = omega ** (k * minus_one_label) if minus_one_label is not None else 1.0
        is_spinorial = chi_minus.real < -0.5
        typ = 'S' if is_spinorial else 'T'
        if is_spinorial:
            n_spinorial += 1
        else:
            n_tensorial += 1

        spectrum = {}
        for two_j in range(0, 2 * max_j + 1):
            j = two_j / 2.0
            if is_spinorial and two_j % 2 == 0:
                continue
            if not is_spinorial and two_j % 2 == 1:
                continue
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            m = round((np.sum(chi_k.conj() * chi_Vj * sizes) / G).real)
            if m > 0:
                spectrum[j] = m

        sectors.append((f'χ_{k}', typ, spectrum))

    # Total scalar states per j
    total_per_j = {}
    for two_j in range(0, 2 * max_j + 1):
        j = two_j / 2.0
        total = 0
        for _, typ, spec in sectors:
            total += spec.get(j, 0)
        if total > 0:
            total_per_j[j] = total

    return {
        'G': G,
        'abel_size': abel_size,
        'obstructed': obstructed,
        'n_tensorial': n_tensorial,
        'n_spinorial': n_spinorial,
        'sectors': sectors,
        'total_per_j': total_per_j,
        'comm_size': len(comm),
    }


def main():
    print("=" * 70)
    print("SPECTRAL COMPARISON TABLE")
    print("=" * 70)

    # Groups to compare
    groups = [
        ("SO(3) [Z₂]", np.array([[1,0,0,0],[-1,0,0,0]], dtype=float)),
        ("SO(3)/T [2T]", build_binary_tetrahedral()),
        ("SO(3)/O [2O]", build_binary_octahedral()),
        ("SO(3)/D₃ [Dic₃]", build_dicyclic(3)),
    ]

    results = {}
    for name, elements in groups:
        result = compute_scalar_spectrum(elements, max_j=10)
        results[name] = result

        print(f"\n{'─'*50}")
        print(f"  {name}")
        print(f"{'─'*50}")
        print(f"  |G| = {result['G']}, |[G,G]| = {result['comm_size']}")
        print(f"  Abel ≅ Z_{result['abel_size']}")
        print(f"  Obstructed: {'YES' if result['obstructed'] else 'NO'}")
        print(f"  Sectors: {result['n_tensorial']}T + {result['n_spinorial']}S"
              f" = {result['abel_size']} total")

        for sname, typ, spec in result['sectors']:
            if spec:
                j_vals = sorted(spec.keys())
                E_min = j_vals[0] * (j_vals[0] + 1)
                print(f"    {sname}({typ}): first j={j_vals[0]:.1f} (E={E_min:.2f})"
                      f", states(j≤10): {sum(spec.values())}")

    # Comparison table
    print("\n" + "=" * 70)
    print("SCALAR STATES PER ENERGY LEVEL")
    print("=" * 70)

    # Collect all j values
    all_j = set()
    for r in results.values():
        all_j.update(r['total_per_j'].keys())
    all_j = sorted(all_j)

    # Header
    names_short = ["SO(3)", "SO(3)/T", "SO(3)/O", "SO(3)/D₃"]
    header = f"  {'j':>5}  {'E':>9}"
    for n in names_short:
        header += f"  {n:>8}"
    print(f"\n{header}")
    print("  " + "-" * (len(header) - 2))

    name_keys = list(results.keys())
    for j in all_j:
        if j > 6.5:
            break
        E = j * (j + 1)
        line = f"  {j:5.1f}  {E:9.2f}"
        for nk in name_keys:
            val = results[nk]['total_per_j'].get(j, 0)
            line += f"  {val:8d}"
        print(line)

    # Cumulative counts
    print(f"\n{'':>17}", end="")
    for n in names_short:
        print(f"  {n:>8}", end="")
    print()

    print(f"  {'Total (j≤6.5)':>15}", end="")
    for nk in name_keys:
        total = sum(v for j, v in results[nk]['total_per_j'].items() if j <= 6.5)
        print(f"  {total:8d}", end="")
    print()

    print(f"  {'Tensorial only':>15}", end="")
    for nk in name_keys:
        total_t = 0
        for _, typ, spec in results[nk]['sectors']:
            if typ == 'T':
                total_t += sum(v for j, v in spec.items() if j <= 6.5)
        print(f"  {total_t:8d}", end="")
    print()

    print(f"  {'Spinorial only':>15}", end="")
    for nk in name_keys:
        total_s = 0
        for _, typ, spec in results[nk]['sectors']:
            if typ == 'S':
                total_s += sum(v for j, v in spec.items() if j <= 6.5)
        print(f"  {total_s:8d}", end="")
    print()

    # Asymptotic density check
    print("\n" + "=" * 70)
    print("ASYMPTOTIC DENSITY CHECK (Weyl law)")
    print("=" * 70)
    print("\nAsymptotic: m(χ_k, j) → (2j+1)/|G| for each 1D char χ_k")
    print("  → total scalar at j ~ κ·(2j+1)/|G|")
    print()

    for nk, name in zip(name_keys, names_short):
        r = results[nk]
        print(f"\n  {name}: κ={r['abel_size']}, |G|={r['G']}")
        print(f"    Prediction: density ~ {r['abel_size']}/{r['G']}"
              f" = {r['abel_size']/r['G']:.4f} per (2j+1)")

        # Check at high j
        for j in [8, 9, 10]:
            total = r['total_per_j'].get(j, 0) + r['total_per_j'].get(j + 0.5, 0)
            expected = r['abel_size'] * (2 * j + 1) / r['G']
            if total > 0:
                ratio = total / (2 * j + 1)
                print(f"    j={j}: actual={total}, dim={2*j+1}"
                      f", ratio={ratio:.4f} (expect {r['abel_size']/r['G']:.4f})")

    # Spectral gap summary
    print("\n" + "=" * 70)
    print("SPECTRAL GAPS")
    print("=" * 70)

    for nk, name in zip(name_keys, names_short):
        r = results[nk]
        print(f"\n  {name}:")

        # First tensorial excited state
        first_T = None
        for _, typ, spec in r['sectors']:
            if typ == 'T':
                for j in sorted(spec.keys()):
                    if j > 0:
                        if first_T is None or j < first_T:
                            first_T = j
                        break

        # First spinorial state
        first_S = None
        for _, typ, spec in r['sectors']:
            if typ == 'S':
                for j in sorted(spec.keys()):
                    if first_S is None or j < first_S:
                        first_S = j
                    break

        E_ground = 0
        if first_T is not None:
            print(f"    Ground → first excited tensorial: ΔE = {first_T*(first_T+1):.2f}"
                  f" (j={first_T:.1f})")
        if first_S is not None:
            print(f"    Ground → first spinorial scalar:  ΔE = {first_S*(first_S+1):.2f}"
                  f" (j={first_S:.1f})")
        elif r['obstructed']:
            print(f"    Spinorial scalar: FORBIDDEN (obstructed)")

    # Final summary
    print("\n" + "=" * 70)
    print("PHYSICAL SUMMARY")
    print("=" * 70)
    print("""
The Q₈ obstruction has three spectral signatures:

1. SECTOR ELIMINATION: Obstructed groups have NO spinorial scalar sectors.
   All scalar states are at integer j only. Half-integer j excluded entirely.

2. SPECTRAL SPARSITY: Obstructed groups have dramatically fewer scalar states.
   At j ≤ 6.5: SO(3)/O has ~5 states, SO(3)/D₃ has ~35 states.
   The obstruction "decimates" the scalar spectrum.

3. SPINORIAL ENERGY SHIFT: For non-obstructed groups, the first spinorial
   scalar state appears at j_min = |[G,G]|/2, not at j = 1/2.
   The commutator subgroup delays spinorial onset:
     Cyclic (|[G,G]|=1): j_min = 1/2 (E = 0.75)
     Dic₃  (|[G,G]|=3): j_min = 3/2 (E = 3.75)
     Dic₅  (|[G,G]|=5): j_min = 5/2 (E = 8.75)

These are measurable predictions for quantum rotors with discrete symmetry.
""")


if __name__ == '__main__':
    main()
