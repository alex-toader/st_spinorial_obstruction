#!/usr/bin/env python3
"""
Weyl law: average spectral density converges to κ/|G|.

The number of scalar states at a given j fluctuates, but the AVERAGE
over a window of j values converges to a well-defined density.

Total scalar states up to angular momentum J:
  N_scalar(J) = Σ_{j=0}^{J} Σ_k m(χ_k, j)

For large J, N_scalar(J) → κ · (J+1)² / |G|

where (J+1)² ≈ Σ_{j=0}^{J} (2j+1) counts all states on SO(3).

The ratio N_scalar(J) / N_total(J) → κ/|G| is the "scalar fraction."

Note: for obstructed groups, only integer j contribute.
For non-obstructed, both integer and half-integer.
The sum counts both.

RAW OUTPUT:
===========

======================================================================
WEYL LAW: CUMULATIVE SCALAR STATE COUNT
======================================================================

N_scalar(J) / N_total(J) → κ/|G|  as J → ∞
where N_total(J) = Σ_{j=0,1/2,1,...}^{J} (2j+1) = (2J+1)(J+1)

───────────────────────────────────────────────────────
  SO(3): |G|=2, κ=2, target κ/|G|=1.000000
───────────────────────────────────────────────────────
      J  N_scalar   N_total       ratio      target    err%
   50.0      5151      5151    1.000000    1.000000    0.0%

───────────────────────────────────────────────────────
  SO(3)/D₃: |G|=12, κ=4, target κ/|G|=0.333333
───────────────────────────────────────────────────────
      J  N_scalar   N_total       ratio      target    err%
   50.0      1717      5151    0.333333    0.333333    0.0%

───────────────────────────────────────────────────────
  SO(3)/T: |G|=24, κ=3, target κ/|G|=0.125000
───────────────────────────────────────────────────────
      J  N_scalar   N_total       ratio      target    err%
   50.0       651      5151    0.126383    0.125000    1.1%

───────────────────────────────────────────────────────
  SO(3)/O: |G|=48, κ=2, target κ/|G|=0.041667
───────────────────────────────────────────────────────
      J  N_scalar   N_total       ratio      target    err%
   50.0       217      5151    0.042128    0.041667    1.1%

───────────────────────────────────────────────────────
  SO(3)/I: |G|=120, κ=1, target κ/|G|=0.008333
───────────────────────────────────────────────────────
      J  N_scalar   N_total       ratio      target    err%
   50.0        44      5151    0.008542    0.008333    2.5%

SUMMARY AT J=50:
         Space    |G|    κ   N_scl   N_tot     ratio     κ/|G|    err%
         SO(3)      2    2    5151    5151   1.00000   1.00000    0.0%
      SO(3)/D₃     12    4    1717    5151   0.33333   0.33333    0.0%
       SO(3)/T     24    3     651    5151   0.12638   0.12500    1.1%
       SO(3)/O     48    2     217    5151   0.04213   0.04167    1.1%
       SO(3)/I    120    1      44    5151   0.00854   0.00833    2.5%
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
from src.quaternion import qkey, qmul
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       build_binary_icosahedral,
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


def total_scalar_at_j(j, classes, sizes, G, abel_size):
    """Total scalar states at angular momentum j.

    Uses: Σ_k m(χ_k, j) = (|Abel|/|G|) Σ_{g∈[G,G]} χ_j(g)

    Equivalently: for each class C_i, if C_i ⊂ [G,G], its contribution is
    weighted by |Abel|, otherwise 0.

    Proof: Σ_k χ_k(g)* = |Abel| if g∈[G,G], else 0 (character orthogonality
    on quotient group G/[G,G]).
    """
    chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])

    # Total = (1/|G|) Σ_g (Σ_k χ_k(g)*) · χ_j(g)
    # = (|Abel|/|G|) Σ_{g∈[G,G]} χ_j(g)
    # = (|Abel|/|G|) Σ_{classes in [G,G]} size_i · χ_j(trace_i)

    # For this, we need to know which classes are in [G,G].
    # Use: class is in [G,G] iff representative is in [G,G].
    # We mark classes during setup.
    raise NotImplementedError("Use precomputed version")


def main():
    print("=" * 70)
    print("WEYL LAW: CUMULATIVE SCALAR STATE COUNT")
    print("=" * 70)
    print("\nN_scalar(J) / N_total(J) → κ/|G|  as J → ∞")
    print("where N_total(J) = Σ_{j=0,1/2,1,...}^{J} (2j+1) = (2J+1)(J+1)")

    groups = [
        ("SO(3)", np.array([[1,0,0,0],[-1,0,0,0]], dtype=float)),
        ("SO(3)/D₃", build_dicyclic(3)),
        ("SO(3)/T", build_binary_tetrahedral()),
        ("SO(3)/O", build_binary_octahedral()),
        ("SO(3)/I", build_binary_icosahedral()),
    ]

    J_max = 50

    for name, elements in groups:
        classes = compute_conjugacy_classes(elements)
        classes.sort(key=lambda c: -c['trace'])
        sizes = np.array([c['size'] for c in classes], dtype=float)
        G = int(sum(sizes))
        comm = compute_commutator_subgroup(elements)
        comm_keys = set(comm.keys())
        abel = G // len(comm)
        target = abel / G

        # Mark which classes are in [G,G]
        class_in_comm = []
        for c in classes:
            rep_key = list(c['keys'])[0]
            class_in_comm.append(rep_key in comm_keys)
        class_in_comm = np.array(class_in_comm)

        # Compute cumulative scalar count
        N_scalar = 0
        N_total = 0

        print(f"\n{'─'*55}")
        print(f"  {name}: |G|={G}, κ={abel}, target κ/|G|={target:.6f}")
        print(f"{'─'*55}")
        print(f"  {'J':>5}  {'N_scalar':>8}  {'N_total':>8}  {'ratio':>10}  {'target':>10}  {'err%':>6}")

        checkpoints = [5, 10, 15, 20, 25, 30, 40, 50]

        for two_j in range(0, 2 * J_max + 1):
            j = two_j / 2.0
            dim = int(2 * j + 1)
            N_total += dim

            # Scalar count at this j:
            # = (abel/G) * Σ_{classes in [G,G]} size * chi_j(trace)
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            scalar_j = round((abel / G * np.sum(
                sizes * chi_Vj * class_in_comm)).real)
            N_scalar += scalar_j

            j_val = two_j / 2.0
            if j_val in [float(c) for c in checkpoints]:
                ratio = N_scalar / N_total if N_total > 0 else 0
                err = abs(ratio - target) / target * 100
                print(f"  {j_val:5.1f}  {N_scalar:8d}  {N_total:8d}  "
                      f"{ratio:10.6f}  {target:10.6f}  {err:5.1f}%")

    # Comparison table at J=50
    print("\n" + "=" * 70)
    print("SUMMARY AT J=50")
    print("=" * 70)
    print(f"\n  {'Space':>12}  {'|G|':>5}  {'κ':>3}  {'N_scl':>6}  {'N_tot':>6}"
          f"  {'ratio':>8}  {'κ/|G|':>8}  {'err%':>6}")
    print(f"  {'-'*62}")

    for name, elements in groups:
        classes = compute_conjugacy_classes(elements)
        classes.sort(key=lambda c: -c['trace'])
        sizes = np.array([c['size'] for c in classes], dtype=float)
        G = int(sum(sizes))
        comm = compute_commutator_subgroup(elements)
        comm_keys = set(comm.keys())
        abel = G // len(comm)
        target = abel / G

        class_in_comm = np.array([list(c['keys'])[0] in comm_keys for c in classes])

        N_scalar = 0
        N_total = 0
        for two_j in range(0, 101):
            j = two_j / 2.0
            dim = int(2 * j + 1)
            N_total += dim
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            scalar_j = round((abel / G * np.sum(
                sizes * chi_Vj * class_in_comm)).real)
            N_scalar += scalar_j

        ratio = N_scalar / N_total
        err = abs(ratio - target) / target * 100
        print(f"  {name:>12}  {G:5d}  {abel:3d}  {N_scalar:6d}  {N_total:6d}"
              f"  {ratio:8.5f}  {target:8.5f}  {err:5.1f}%")

    print("""
INTERPRETATION:
  N_scalar(J) / N_total(J) converges to κ/|G|.
  This is the fraction of the full SO(3) Hilbert space that survives
  as scalar states on SO(3)/H.

  For SO(3)/I: only 1/120 ≈ 0.83% of states survive as scalars.
  For SO(3)/O: 2/48 ≈ 4.2%.
  For SO(3)/D₃: 4/12 ≈ 33%.

  The remaining states exist in higher-rank vector bundle sectors.
  The obstruction determines WHICH states are scalar and which require bundles.
""")


if __name__ == '__main__':
    main()
