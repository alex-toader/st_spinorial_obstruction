#!/usr/bin/env python3
"""
Sector-resolved spectrum of the rigid rotor on SO(3)/D₄.

2D₄ = Dic₄ (order 16) is OBSTRUCTED: Q₈ ⊂ 2D₄.
Compare with 2D₃ (non-obstructed) to show obstruction toggle.

Dic₄ character analysis:
  Relations: a⁸=1, x²=a⁴=-1, xax⁻¹=a⁻¹
  → χ(a) ∈ {±1}, χ(x) ∈ {±1}  (4 characters, Abel ≅ Z₂×Z₂)
  → χ(-1) = χ(a⁴) = χ(a)⁴ = +1 for ALL characters
  → ALL sectors tensorial (integer j only)

D₃ → D₄ transition: spinorial sectors disappear.
This is the "dihedral parity toggle" — obstruction alternates with N.

RAW OUTPUT:
===========
D₃ (N=3, odd): NON-obstructed, 2T+2S, 40 states, spinorial at j=1.5+
D₄ (N=4, even): OBSTRUCTED, 4T+0S, 32 states, ALL integer j
D₅ (N=5, odd): NON-obstructed, 2T+2S, 24 states, spinorial at j=2.5+

Toggle: N odd ↔ spinorial exists, N even ↔ forbidden
j_min confirms: D₃→3/2, D₅→5/2 = |[G,G]|/2 ✓
D₄ has Abel≅Z₂×Z₂ (not Z₄), all 4 chars tensorial

Key: D₃→D₄ transition kills ALL half-integer scalar states.
     D₄→D₅ transition restores them.
"""
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))

import numpy as np
from src.quaternion import qkey, qmul
from src.group import compute_conjugacy_classes, compute_commutator_subgroup


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


def build_dicyclic_characters(N, elements, classes):
    """Build all 1D characters of Dic_N analytically.

    For Dic_N with relations a^{2N}=1, x²=a^N=-1, xax⁻¹=a⁻¹:
      - χ(a) must be real (since χ(a)=χ(a⁻¹)=χ(a)*) → χ(a) ∈ {±1}
      - x²=a^N → χ(x)²=χ(a)^N

    N odd:  χ(a)=+1 → χ(x)²=1; χ(a)=-1 → χ(x)²=-1 → χ(x)=±i
            4 characters, 2 tensorial + 2 spinorial
    N even: χ(a)=+1 → χ(x)²=1; χ(a)=-1 → χ(x)²=1
            4 characters, ALL tensorial (χ(-1)=χ(a)^N=(±1)^N=+1)
    """
    # Build element quaternions for a^k and x·a^k
    a_powers = {}
    for k in range(2 * N):
        angle = k * np.pi / N
        a_powers[k] = np.array([np.cos(angle), 0, 0, np.sin(angle)])

    x_elem = np.array([0, 0, 1, 0], dtype=float)
    xa_powers = {}
    for k in range(2 * N):
        xa_powers[k] = qmul(x_elem, a_powers[k])

    # Character definitions
    if N % 2 == 1:  # N odd
        char_defs = [
            ('χ₀', 1.0, 1.0+0j),      # a→+1, x→+1
            ('χ₁', 1.0, -1.0+0j),     # a→+1, x→-1
            ('χ₂', -1.0, 1j),         # a→-1, x→+i
            ('χ₃', -1.0, -1j),        # a→-1, x→-i
        ]
    else:  # N even
        char_defs = [
            ('χ₀', 1.0, 1.0+0j),      # a→+1, x→+1
            ('χ₁', 1.0, -1.0+0j),     # a→+1, x→-1
            ('χ₂', -1.0, 1.0+0j),     # a→-1, x→+1
            ('χ₃', -1.0, -1.0+0j),    # a→-1, x→-1
        ]

    chars = {}
    for cname, chi_a, chi_x in char_defs:
        elem_vals = {}
        for k in range(2 * N):
            elem_vals[qkey(a_powers[k])] = chi_a ** k
        for k in range(2 * N):
            elem_vals[qkey(xa_powers[k])] = chi_x * (chi_a ** k)

        chi_on_classes = []
        for c in classes:
            rep_key = list(c['keys'])[0]
            chi_on_classes.append(elem_vals[rep_key])
        chars[cname] = np.array(chi_on_classes)

    return chars


def analyze_group(N):
    """Full spectral analysis for Dic_N."""
    elements = build_dicyclic(N)
    G_order = len(elements)

    # Verify closure
    key_set = {qkey(e) for e in elements}
    for i in range(G_order):
        for j_idx in range(G_order):
            assert qkey(qmul(elements[i], elements[j_idx])) in key_set

    # Commutator
    comm = compute_commutator_subgroup(elements)
    mk = qkey(np.array([-1, 0, 0, 0], dtype=float))
    obstructed = mk in set(comm.keys())

    # Classes
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    G = int(sum(sizes))

    label = f"D_{N}"
    print(f"\n{'═'*60}")
    print(f"  2{label} = Dic_{N}, |G| = {G}")
    print(f"{'═'*60}")
    print(f"  |[G,G]| = {len(comm)}, |Abel| = {G // len(comm)}")
    print(f"  -1 ∈ [G,G]: {obstructed} → {'OBSTRUCTED' if obstructed else 'non-obstructed'}")

    # Characters
    chars = build_dicyclic_characters(N, elements, classes)

    # Check orthogonality
    print(f"\n  Orthogonality check:")
    names = list(chars.keys())
    all_ok = True
    for i, n1 in enumerate(names):
        for j_idx, n2 in enumerate(names):
            if j_idx >= i:
                inner = np.sum(chars[n1].conj() * chars[n2] * sizes) / G
                expected = 1.0 if i == j_idx else 0.0
                ok = abs(inner - expected) < 1e-10
                if not ok:
                    all_ok = False
                    print(f"    ⟨{n1},{n2}⟩ = {inner.real:.4f} ✗")
    if all_ok:
        print(f"    All ⟨χᵢ,χⱼ⟩ = δᵢⱼ ✓")

    # Types
    print(f"\n  Character types:")
    for name, chi in chars.items():
        for i, c in enumerate(classes):
            if abs(c['trace'] + 2) < 1e-10:
                val = chi[i]
                typ = "tensorial" if val.real > 0.5 else "SPINORIAL"
                print(f"    {name}: χ(-1) = {val.real:+.0f} → {typ}")
                break

    # Multiplicity table
    print(f"\n  {'j':>5}  {'dim':>4}", end="")
    for name in names:
        print(f"  {name:>6}", end="")
    print(f"  {'E':>9}")
    print(f"  {'-'*55}")

    for two_j in range(0, 15):
        j = two_j / 2.0
        chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
        line = f"  {j:5.1f}  {int(2*j+1):4d}"
        for name in names:
            m = round((np.sum(chars[name].conj() * chi_Vj * sizes) / G).real)
            line += f"  {m:6d}"
        line += f"  {j*(j+1):9.2f}"
        print(line)

    # Sector summary
    print(f"\n  Sector summary:")
    total_all = 0
    for name, chi in chars.items():
        for i, c in enumerate(classes):
            if abs(c['trace'] + 2) < 1e-10:
                is_spin = chi[i].real < -0.5
                break

        j_range = ([d / 2.0 for d in range(1, 15, 2)] if is_spin
                   else [float(d) for d in range(8)])

        states = []
        for j in j_range:
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            m = round((np.sum(chi.conj() * chi_Vj * sizes) / G).real)
            if m > 0:
                states.append((j, m))

        total = sum(m for _, m in states)
        total_all += total
        typ = 'S' if is_spin else 'T'
        if states:
            j_list = ", ".join(f"{j:.1f}" for j, _ in states[:6])
            print(f"    {name}({typ}): j = {j_list}  [{total} states]")
        else:
            print(f"    {name}({typ}): no states up to j=7")

    print(f"    Total: {total_all} scalar states")
    return total_all


def main():
    print("=" * 70)
    print("DIHEDRAL PARITY TOGGLE: D₃ vs D₄ SPECTRAL COMPARISON")
    print("=" * 70)

    t3 = analyze_group(3)
    t4 = analyze_group(4)

    # Also do D₅ for the pattern
    t5 = analyze_group(5)

    print("\n" + "=" * 70)
    print("PARITY TOGGLE SUMMARY")
    print("=" * 70)
    print(f"""
Dihedral index N and scalar spectrum:

  N=3 (odd):  NON-obstructed, 2T+2S sectors, {t3} scalar states
  N=4 (even): OBSTRUCTED,     4T+0S sectors, {t4} scalar states
  N=5 (odd):  NON-obstructed, 2T+2S sectors, {t5} scalar states

Pattern: N odd ↔ spinorial scalars exist
         N even ↔ all sectors tensorial

Physical meaning:
  D₃ → D₄: adding one rotation axis KILLS spinorial sectors
  D₄ → D₅: adding one more RESTORES them
  The obstruction toggles with the parity of the dihedral index.

This is the spectral manifestation of the Q₈ criterion:
  Q₈ ⊂ Dic_N ⟺ N even ⟺ no spinorial scalars
""")


if __name__ == '__main__':
    main()
