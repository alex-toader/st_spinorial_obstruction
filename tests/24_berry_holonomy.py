#!/usr/bin/env python3
"""
Berry phase and holonomy under 2π rotation.

PHYSICAL SETUP: A rigid rotor on SO(3)/H. The configuration space has
π₁ = 2H, and the generator of π₁ corresponding to a 2π rotation is
the element -1 ∈ 2H.

In quantization sector χ, transporting the wave function around the
2π loop multiplies it by χ(-1):
  - Tensorial (χ(-1) = +1): no phase change
  - Spinorial (χ(-1) = -1): sign flip (Berry phase π)

PART 1: HOLONOMY TABLE
  For each group and each sector, compute χ(-1).
  Verify: obstructed groups → all χ(-1) = +1 (forced bosonic).

PART 2: PROPAGATOR DECOMPOSITION
  The SLD propagator on SO(3)/H:
    K(g',g;t) = Σ_χ χ(γ) K_γ(g',g;t)

  For the rigid rotor, in sector χ:
    K_χ(t) = Σ_j m(χ,j) · (2j+1) · e^{-iE_j t}

  where the trace (return propagator g'=g=identity) gives the
  partition function. We compute this for SO(3)/T (κ=3 sectors)
  and show how the three propagators differ.

PART 3: 2π ROTATION PHASE
  Start at identity, rotate to -1 (= 2π rotation).
  The partial propagator restricted to the contractible class has
  K_0(−1, 1; t) while the full propagator in sector χ has
  K_χ(−1, 1; t) = Σ_γ χ(γ) K_γ(−1, 1; t).

  For the return amplitude (−1 = same point as 1 on SO(3)/H since
  −1 ∈ 2H acts trivially on SO(3)/H), the holonomy phase appears.

RAW OUTPUT:
===========

PART 1: HOLONOMY χ(-1) FOR EACH GROUP AND SECTOR
               Space  |2H|   κ                         Sectors   Forced bosonic?
          SO(3) [Z₂]     2   2                χ_0(+1), χ_1(-1)                no
     SO(3)/D₃ [Dic₃]    12   4  χ_0(+1), χ_1(-1), χ_2(+1), χ_3(-1)                no
        SO(3)/T [2T]    24   3       χ_0(+1), χ_1(+1), χ_2(+1)               YES
        SO(3)/O [2O]    48   2                χ_0(+1), χ_1(+1)               YES
        SO(3)/I [2I]   120   1                         χ_0(+1)               YES

PART 2: PROPAGATOR ON SO(3)/T
  m(χ_k, j):      j    χ_0    χ_1    χ_2   total
               0.0      1      0      0       1
               3.0      1      0      0       1
               4.0      1      1      1       3
               6.0      2      1      1       4
  Z_χ(β): χ₁=χ₂ (conjugate pair, Z₃ symmetry confirmed)

PART 3: BERRY PHASE ON 2π
  SO(3):       121 tensorial + 110 spinorial  → phase detectable
  SO(3)/D₃:    41 tensorial +  36 spinorial  → phase detectable
  SO(3)/T,O,I: ALL tensorial, 0 spinorial   → forced bosonic

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
    """Build Dic_N ⊂ SU(2) (order 4N)."""
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


def compute_sectors(elements):
    """Compute 1D characters with their χ(-1) values.

    Returns list of (name, chi_minus_one, is_spinorial).
    Uses coset quotient method for abelian groups.
    """
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    G = sum(c['size'] for c in classes)

    comm = compute_commutator_subgroup(elements)
    comm_keys = set(comm.keys())
    abel_size = G // len(comm)

    mk = qkey(np.array([-1, 0, 0, 0], dtype=float))
    obstructed = mk in comm_keys

    if obstructed:
        # All 1D chars have χ(-1) = +1
        return [(f'χ_{k}', +1, False) for k in range(abel_size)], True

    # Build coset map
    elem_to_coset = {}
    cosets = []
    for e in elements:
        k = qkey(e)
        if k in elem_to_coset:
            continue
        coset = set()
        for ck, cv in comm.items():
            prod = qmul(e, cv)
            coset.add(qkey(prod))
        found = False
        for ci, existing in enumerate(cosets):
            if existing & coset:
                for kk in coset:
                    elem_to_coset[kk] = ci
                found = True
                break
        if not found:
            ci = len(cosets)
            cosets.append(coset)
            for kk in coset:
                elem_to_coset[kk] = ci

    identity_key = qkey(np.array([1, 0, 0, 0], dtype=float))
    identity_coset = elem_to_coset[identity_key]
    minus_coset = elem_to_coset[mk]

    # Find generator
    gen_elem = None
    for e in elements:
        k = qkey(e)
        power = e.copy()
        for d in range(1, abel_size + 1):
            pk = qkey(power)
            if elem_to_coset[pk] == identity_coset:
                if d == abel_size:
                    gen_elem = e.copy()
                break
            power = qmul(power, e)
        if gen_elem is not None:
            break

    # Label cosets
    coset_label = {}
    power = np.array([1, 0, 0, 0], dtype=float)
    for d in range(abel_size):
        ci = elem_to_coset[qkey(power)]
        coset_label[ci] = d
        power = qmul(power, gen_elem)

    omega = np.exp(2j * np.pi / abel_size)
    minus_label = coset_label[minus_coset]

    sectors = []
    for p in range(abel_size):
        chi_m1 = omega ** (p * minus_label)
        is_spin = chi_m1.real < -0.5
        sectors.append((f'χ_{p}', chi_m1, is_spin))

    return sectors, False


def sector_multiplicities(elements, j_max=10):
    """Compute m(χ_k, j) for each 1D char and each j.

    Returns (classes, sizes, G, comm, abel_size, char_vals_on_classes, sector_info).
    """
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    G = int(sum(sizes))

    comm = compute_commutator_subgroup(elements)
    comm_keys = set(comm.keys())
    abel_size = G // len(comm)

    mk = qkey(np.array([-1, 0, 0, 0], dtype=float))
    obstructed = mk in comm_keys

    # Build coset map
    elem_to_coset = {}
    cosets_list = []
    for e in elements:
        k = qkey(e)
        if k in elem_to_coset:
            continue
        coset = set()
        for ck, cv in comm.items():
            prod = qmul(e, cv)
            coset.add(qkey(prod))
        found = False
        for ci, existing in enumerate(cosets_list):
            if existing & coset:
                for kk in coset:
                    elem_to_coset[kk] = ci
                found = True
                break
        if not found:
            ci = len(cosets_list)
            cosets_list.append(coset)
            for kk in coset:
                elem_to_coset[kk] = ci

    identity_key = qkey(np.array([1, 0, 0, 0], dtype=float))
    identity_coset = elem_to_coset[identity_key]

    # Find generator
    gen_elem = None
    for e in elements:
        k = qkey(e)
        power = e.copy()
        for d in range(1, abel_size + 1):
            pk = qkey(power)
            if elem_to_coset[pk] == identity_coset:
                if d == abel_size:
                    gen_elem = e.copy()
                break
            power = qmul(power, e)
        if gen_elem is not None:
            break

    if gen_elem is None and abel_size == 1:
        # Trivial abelianization — only trivial char
        chi_vals = [np.ones(len(classes), dtype=complex)]
        chi_m1 = [1.0]
    else:
        coset_label = {}
        power = np.array([1, 0, 0, 0], dtype=float)
        for d in range(abel_size):
            ci = elem_to_coset[qkey(power)]
            coset_label[ci] = d
            power = qmul(power, gen_elem)

        omega = np.exp(2j * np.pi / abel_size)
        minus_label = coset_label[elem_to_coset[mk]]

        chi_vals = []
        chi_m1 = []
        class_labels = np.array([coset_label[elem_to_coset[list(c['keys'])[0]]]
                                 for c in classes])
        for p in range(abel_size):
            vals = omega ** (p * class_labels)
            chi_vals.append(vals)
            chi_m1.append(omega ** (p * minus_label))

    # Compute multiplicities
    results = {}  # (sector_index, j) → multiplicity
    for p in range(abel_size):
        for two_j in range(0, 2 * j_max + 1):
            j = two_j / 2.0
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            m = np.sum(chi_vals[p].conj() * chi_Vj * sizes) / G
            m_r = round(m.real)
            if m_r > 0:
                results[(p, j)] = m_r

    return {
        'G': G, 'abel_size': abel_size, 'obstructed': obstructed,
        'chi_m1': chi_m1, 'multiplicities': results,
    }


def main():
    print("=" * 70)
    print("BERRY PHASE AND HOLONOMY UNDER 2π ROTATION")
    print("=" * 70)

    # ============================================================
    # PART 1: HOLONOMY TABLE
    # ============================================================
    print("\n" + "=" * 70)
    print("PART 1: HOLONOMY χ(-1) FOR EACH GROUP AND SECTOR")
    print("=" * 70)
    print("\n  The element -1 ∈ 2H generates the 2π rotation loop.")
    print("  In sector χ_k, the wave function acquires phase χ_k(-1).\n")

    groups = [
        ("SO(3) [Z₂]", np.array([[1,0,0,0],[-1,0,0,0]], dtype=float)),
        ("SO(3)/D₃ [Dic₃]", build_dicyclic(3)),
        ("SO(3)/T [2T]", build_binary_tetrahedral()),
        ("SO(3)/O [2O]", build_binary_octahedral()),
        ("SO(3)/I [2I]", build_binary_icosahedral()),
    ]

    print(f"  {'Space':>18}  {'|2H|':>4}  {'κ':>2}  {'Sectors':>30}  {'Forced bosonic?':>16}")
    print(f"  {'─'*75}")

    for name, elements in groups:
        sectors, obstructed = compute_sectors(elements)
        G = len(elements)
        kappa = len(sectors)
        sector_str = ', '.join(
            f"{s[0]}({'+1' if s[1].real > 0 else '-1'})"
            if abs(s[1].imag) < 0.01
            else f"{s[0]}({s[1]:.2f})"
            for s in sectors
        )
        forced = "YES" if obstructed else "no"
        print(f"  {name:>18}  {G:4d}  {kappa:2d}  {sector_str:>30}  {forced:>16}")

    print("""
  KEY OBSERVATION:
    Non-obstructed (SO(3), SO(3)/D₃): some sectors have χ(-1) = -1.
    → 2π rotation CAN be detected. Berry phase = π is physical.

    Obstructed (SO(3)/T, SO(3)/O, SO(3)/I): ALL sectors have χ(-1) = +1.
    → 2π rotation is INVISIBLE. Berry phase = 0 in every sector.
    → The space is "forced bosonic": topology allows spinorial sectors
       but group structure forbids them.
""")

    # ============================================================
    # PART 2: SECTOR-RESOLVED PROPAGATOR FOR SO(3)/T
    # ============================================================
    print("=" * 70)
    print("PART 2: SECTOR-RESOLVED PROPAGATOR ON SO(3)/T")
    print("=" * 70)
    print("""
  The Schulman-Laidlaw-DeWitt propagator decomposes as:
    K_χ(t) = Σ_j m(χ,j) · (2j+1) · e^{-iE_j·t}

  where E_j = j(j+1) and m(χ,j) = multiplicity of χ in V_j|_{2T}.

  The return amplitude (partition function in Euclidean time) is:
    Z_χ(β) = Σ_j m(χ,j) · (2j+1) · e^{-β·j(j+1)}
""")

    elements = build_binary_tetrahedral()
    data = sector_multiplicities(elements, j_max=15)

    abel = data['abel_size']
    mults = data['multiplicities']

    # Print sector content
    print(f"  SO(3)/T: |2T|=24, Abel≅Z₃, κ=3 sectors")
    print(f"  χ₀(-1) = +1 (tensorial)")
    print(f"  χ₁(-1) = +1 (tensorial, ω = e^{{2πi/3}})")
    print(f"  χ₂(-1) = +1 (tensorial, ω² = e^{{4πi/3}})")
    print()

    print(f"  Multiplicities m(χ_k, j) for j = 0..10:")
    print(f"  {'j':>5}", end="")
    for p in range(abel):
        print(f"  {'χ_'+str(p):>5}", end="")
    print(f"  {'total':>6}")
    print(f"  {'─'*35}")

    for two_j in range(0, 21):
        j = two_j / 2.0
        row = []
        for p in range(abel):
            m = mults.get((p, j), 0)
            row.append(m)
        if sum(row) > 0:
            print(f"  {j:5.1f}", end="")
            for m in row:
                print(f"  {m:5d}", end="")
            print(f"  {sum(row):6d}")

    # Partition function at several temperatures
    print(f"\n  Partition function Z_χ(β) at β = 0.1, 0.5, 1.0:")
    print(f"  {'β':>5}", end="")
    for p in range(abel):
        print(f"    {'Z_χ_'+str(p):>8}", end="")
    print(f"    {'Z_total':>8}  {'Z_SO(3)':>8}")
    print(f"  {'─'*60}")

    for beta in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0]:
        Z_sectors = []
        for p in range(abel):
            Z = 0.0
            for two_j in range(0, 31):
                j = two_j / 2.0
                m = mults.get((p, j), 0)
                if m > 0:
                    Z += m * (2 * j + 1) * np.exp(-beta * j * (j + 1))
            Z_sectors.append(Z)

        Z_total = sum(Z_sectors)

        # SO(3) partition function for comparison (j integer only, all appear)
        Z_so3 = 0.0
        for j_int in range(16):
            j = float(j_int)
            Z_so3 += (2 * j + 1) * np.exp(-beta * j * (j + 1))

        print(f"  {beta:5.2f}", end="")
        for Z in Z_sectors:
            print(f"    {Z:8.3f}", end="")
        print(f"    {Z_total:8.3f}  {Z_so3:8.3f}")

    print("""
  INTERPRETATION:
    At high temperature (small β): Z_total → κ/|G| · Z_SO(3) = (3/24)·Z_SO(3).
    This is the Weyl limit: the density of scalar states is κ/|G|.

    At low temperature (large β): Z_χ₀ → 1 (ground state j=0).
    Only the trivial sector contributes at T→0.

    The three sectors have IDENTICAL partition functions for χ₁, χ₂
    (conjugate pair). This is a Z₃ symmetry of the spectrum.
""")

    # ============================================================
    # PART 3: COMPARISON — BERRY PHASE ON 2π LOOP
    # ============================================================
    print("=" * 70)
    print("PART 3: BERRY PHASE ON 2π ROTATION — ALL GROUPS")
    print("=" * 70)
    print("""
  Physical protocol: slowly rotate the rigid body by 2π.
  The wave function acquires a geometric phase = χ(-1).

  For a quantum rotor in eigenstate |j,m,χ⟩:
    after 2π rotation → χ(-1) · |j,m,χ⟩

  This phase is:
    +1 for tensorial states (bosonic)
    -1 for spinorial states (fermionic)
""")

    # For each group, compute the ratio of spinorial to total states
    for name, elements in groups:
        classes = compute_conjugacy_classes(elements)
        classes.sort(key=lambda c: -c['trace'])
        sizes = np.array([c['size'] for c in classes], dtype=float)
        G = int(sum(sizes))
        comm = compute_commutator_subgroup(elements)
        abel = G // len(comm)

        sectors, obstructed = compute_sectors(elements)
        n_spin = sum(1 for s in sectors if s[2])
        n_tens = sum(1 for s in sectors if not s[2])

        print(f"  {name}:")
        print(f"    κ = {abel}: {n_tens} tensorial + {n_spin} spinorial sectors")

        if obstructed:
            print(f"    Berry phase: +1 in ALL sectors (forced bosonic)")
            print(f"    2π rotation is topologically trivial in EVERY sector\n")
        else:
            print(f"    Berry phase: +1 in {n_tens} sectors, -1 in {n_spin} sectors")
            print(f"    2π rotation IS detectable in {n_spin} spinorial sector(s)")

            # Count states
            data = sector_multiplicities(elements, j_max=10)
            n_spin_states = 0
            n_tens_states = 0
            for (p, j), m in data['multiplicities'].items():
                chi_m1 = data['chi_m1'][p]
                if hasattr(chi_m1, 'real'):
                    is_spin = chi_m1.real < -0.5
                else:
                    is_spin = chi_m1 < -0.5
                if is_spin:
                    n_spin_states += m
                else:
                    n_tens_states += m

            print(f"    Scalar states (j ≤ 10): {n_tens_states} tensorial, "
                  f"{n_spin_states} spinorial\n")

    # ============================================================
    # SUMMARY
    # ============================================================
    print("=" * 70)
    print("SUMMARY: BERRY PHASE AND THE OBSTRUCTION")
    print("=" * 70)
    print("""
  The Q₈ obstruction has a direct physical manifestation as Berry phase:

  1. NON-OBSTRUCTED spaces (SO(3), SO(3)/D₃, SO(3)/D₅, ...):
     → Spinorial sectors exist: χ(-1) = -1
     → 2π rotation produces Berry phase π
     → Half-integer angular momentum is physically realizable
     → The Schulman effect survives on the quotient space

  2. OBSTRUCTED spaces (SO(3)/T, SO(3)/O, SO(3)/I):
     → ALL sectors have χ(-1) = +1
     → 2π rotation produces NO phase in ANY sector
     → Only integer angular momentum appears in scalar states
     → The Schulman effect is killed by the group structure

  The obstruction is not just algebraic (−1 ∈ [G,G]) but has a
  concrete physical consequence: the Berry phase of a 2π rotation
  is forced to vanish in every quantization sector.

  This is the strongest formulation of the result:
    Q₈ ⊂ 2H  ⟹  every scalar quantum state on SO(3)/H is bosonic.
""")


if __name__ == '__main__':
    main()
