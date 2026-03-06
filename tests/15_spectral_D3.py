#!/usr/bin/env python3
"""
Sector-resolved spectrum of the rigid rotor on SO(3)/D₃.

D₃ is NOT obstructed (Q₈ ⊄ 2D₃), so scalar spinorial sectors exist.
Compare with SO(3)/O where spinorial sectors require vector bundles.

2D₃ = Dic₃ = binary dihedral group of order 12.
Abel(2D₃) ≅ Z₄ (|[G,G]| = 3), so 4 scalar sectors.

Character table of Z₄ quotient:
  Generator relations: a⁶=1, x²=a³=-1, xax⁻¹=a⁻¹
  Since xax⁻¹=a⁻¹ → χ(a)=χ(a)* → χ(a) ∈ {+1,-1}
  Since x²=a³ → χ(x)²=χ(a)³

  χ₀: a→+1, x→+1   χ₀(-1)=+1 (tensorial)
  χ₁: a→+1, x→-1   χ₁(-1)=+1 (tensorial)
  χ₂: a→-1, x→+i   χ₂(-1)=-1 (SPINORIAL)
  χ₃: a→-1, x→-i   χ₃(-1)=-1 (SPINORIAL)

Key result: 2 tensorial (integer j) + 2 spinorial (half-integer j) sectors.
Spinorial states appear as ORDINARY SCALARS — no vector bundle needed.

RAW OUTPUT:
===========
======================================================================
SECTOR-RESOLVED SPECTRUM: RIGID ROTOR ON SO(3)/D₃
======================================================================

Closure verified ✓
-1 ∈ [G,G]: False → NOT obstructed
|[G,G]| = 3, |Abel(2D₃)| = 4

2D₃: 6 classes, |G| = 12
  C0: size=1, trace=2.0000
  C1: size=2, trace=1.0000
  C2: size=3, trace=0.0000
  C3: size=3, trace=0.0000
  C4: size=2, trace=-1.0000
  C5: size=1, trace=-2.0000

Character orthogonality:
  ⟨χ₀, χ₀⟩ = 1.0000 ✓
  ⟨χ₁, χ₁⟩ = 1.0000 ✓
  ⟨χ₂, χ₂⟩ = 1.0000 ✓
  ⟨χ₃, χ₃⟩ = 1.0000 ✓

Character values at -1:
  χ₀: χ(-1) = +1+0j → tensorial
  χ₁: χ(-1) = +1+0j → tensorial
  χ₂: χ(-1) = -1+0j → SPINORIAL
  χ₃: χ(-1) = -1+0j → SPINORIAL

======================================================================
FULL MULTIPLICITY TABLE
======================================================================

    j   dim   type      χ₀      χ₁      χ₂      χ₃          E
-------------------------------------------------------------
   0.0     1    (T)       1       0       0       0       0.00
   0.5     2    (S)       0       0       0       0       0.75
   1.0     3    (T)       0       1       0       0       2.00
   1.5     4    (S)       0       0       1       1       3.75
   2.0     5    (T)       1       0       0       0       6.00
   2.5     6    (S)       0       0       1       1       8.75
   3.0     7    (T)       1       2       0       0      12.00
   3.5     8    (S)       0       0       1       1      15.75
   4.0     9    (T)       2       1       0       0      20.00
   4.5    10    (S)       0       0       2       2      24.75
   5.0    11    (T)       1       2       0       0      30.00
   5.5    12    (S)       0       0       2       2      35.75
   6.0    13    (T)       3       2       0       0      42.00
   6.5    14    (S)       0       0       2       2      48.75

======================================================================
PHYSICAL SPECTRUM PER SECTOR
======================================================================

Sector χ₀ (tensorial):
  → only integer j contribute
      j   E=j(j+1)   mult
  -------------------------
    0.0       0.00      1
    2.0       6.00      1
    3.0      12.00      1
    4.0      20.00      2
    5.0      30.00      1
    6.0      42.00      3
  Total states (up to j=6.5): 9

Sector χ₁ (tensorial):
  → only integer j contribute
      j   E=j(j+1)   mult
  -------------------------
    1.0       2.00      1
    3.0      12.00      2
    4.0      20.00      1
    5.0      30.00      2
    6.0      42.00      2
  Total states (up to j=6.5): 8

Sector χ₂ (SPINORIAL):
  → only half-integer j contribute
      j   E=j(j+1)   mult
  -------------------------
    1.5       3.75      1
    2.5       8.75      1
    3.5      15.75      1
    4.5      24.75      2
    5.5      35.75      2
    6.5      48.75      2
  Total states (up to j=6.5): 9

Sector χ₃ (SPINORIAL):
  → only half-integer j contribute
      j   E=j(j+1)   mult
  -------------------------
    1.5       3.75      1
    2.5       8.75      1
    3.5      15.75      1
    4.5      24.75      2
    5.5      35.75      2
    6.5      48.75      2
  Total states (up to j=6.5): 9

======================================================================
COMPARISON: SO(3)/D₃ vs SO(3)/O
======================================================================

SO(3)/O (obstructed, Q₈ ⊂ 2O):
  Abel(2O) ≅ Z₂ → 2 scalar sectors, BOTH tensorial
  All scalar states: integer j only
  Spinorial states FORBIDDEN as scalars → require rank-2 vector bundle H₁
  Scalar spectrum: A₀ at j=0,4,6; A₁ at j=3,6 (extremely sparse)

SO(3)/D₃ (non-obstructed, Q₈ ⊄ 2D₃):
  Abel(2D₃) ≅ Z₄ → 4 scalar sectors: 2 tensorial + 2 spinorial
  Spinorial scalar states at half-integer j!
  No vector bundle needed — spinorial wavefunctions exist as ordinary scalars

SPECTRAL MANIFESTATION OF Q₈ OBSTRUCTION:
  The obstruction doesn't just forbid spinorial sectors —
  it forces ALL scalar sectors to be tensorial,
  decimating the spectrum to integer j only.
"""
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))

import numpy as np
from src.group import compute_conjugacy_classes, compute_commutator_subgroup
from src.quaternion import qkey, qmul


def build_binary_dihedral_3():
    """Construct 2D₃ = Dic₃ (order 12) as unit quaternions.

    Generated by a = e^{iπ/3} and x = j, with relations:
      a⁶ = 1, x² = a³ = -1, xax⁻¹ = a⁻¹

    Elements: {aᵏ, x·aᵏ} for k = 0,...,5
    """
    elements = []
    # a^k = (cos(kπ/3), 0, 0, sin(kπ/3))
    for k in range(6):
        angle = k * np.pi / 3
        q = np.array([np.cos(angle), 0, 0, np.sin(angle)])
        elements.append(q)
    # x·a^k where x = j = (0,0,1,0)
    x = np.array([0, 0, 1, 0], dtype=float)
    for k in range(6):
        angle = k * np.pi / 3
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


def assign_character_values(elements, classes):
    """Build the 4 characters of Z₄ = G/[G,G] for 2D₃.

    χ₀: a→+1, x→+1   (trivial)
    χ₁: a→+1, x→-1
    χ₂: a→-1, x→+i
    χ₃: a→-1, x→-i

    Returns dict mapping character name to array of values on classes.
    """
    # Identify which element is which: a^k vs x·a^k
    # a^k have form (cos, 0, 0, sin) — component 1,2 are 0 or just comp 2
    # x·a^k have nonzero components 1 or 2

    # Classify each element: "rotation" (a^k) or "reflection" (x·a^k)
    # a^k: q = (cos(kπ/3), 0, 0, sin(kπ/3)) → q[1]=q[2]=0
    # x·a^k: j·(cos,0,0,sin) = (0, sin, cos, 0) rotated → q[0]=q[3]=0... let me check
    # Actually qmul(j, (w,x,y,z)) = (0,0,1,0)·(w,x,y,z) = (-y, z, w, -x)
    # So x·a^k = (-sin(kπ/3), 0, cos(kπ/3), 0) when a^k = (cos, 0, 0, sin)
    # Wait, j = (0,0,1,0), a^k = (cos(kπ/3), 0, 0, sin(kπ/3))
    # qmul((0,0,1,0), (c,0,0,s)) = (0·c - 0·0 - 1·0 - 0·s,
    #                                 0·0 + 0·c + 1·s - 0·0,  ... )
    # Let me just compute it numerically

    # For character assignment, we need: for each class, determine its coset in G/[G,G]
    # [G,G] = {1, a², a⁴}
    # Cosets of [G,G]:
    #   C0 = {1, a², a⁴}          → a^0 mod [G,G] → χ(a)^0 · χ(x)^0
    #   C1 = {a, a³, a⁵}          → a^1 mod [G,G] → χ(a)^1
    #   C2 = {x, xa², xa⁴}        → x mod [G,G]   → χ(x)^1
    #   C3 = {xa, xa³, xa⁵}       → xa mod [G,G]  → χ(a)·χ(x)

    # We need to identify which conjugacy class corresponds to which coset.
    # First, classify elements

    a_powers = {}  # k -> quaternion for a^k
    for k in range(6):
        angle = k * np.pi / 3
        a_powers[k] = np.array([np.cos(angle), 0, 0, np.sin(angle)])

    x_elem = np.array([0, 0, 1, 0], dtype=float)
    xa_powers = {}  # k -> quaternion for x·a^k
    for k in range(6):
        xa_powers[k] = qmul(x_elem, a_powers[k])

    # Character values:
    # χ(a^k) = χ(a)^k
    # χ(x·a^k) = χ(x)·χ(a)^k

    chars = {}
    char_defs = {
        'χ₀': (1.0, 1.0+0j),      # a→+1, x→+1
        'χ₁': (1.0, -1.0+0j),     # a→+1, x→-1
        'χ₂': (-1.0, 1j),         # a→-1, x→+i
        'χ₃': (-1.0, -1j),        # a→-1, x→-i
    }

    for cname, (chi_a, chi_x) in char_defs.items():
        # Build value for each GROUP ELEMENT, then extract per class
        elem_vals = {}  # qkey -> chi value
        for k in range(6):
            key = qkey(a_powers[k])
            elem_vals[key] = chi_a ** k
        for k in range(6):
            key = qkey(xa_powers[k])
            elem_vals[key] = chi_x * (chi_a ** k)

        # Now assign to conjugacy classes
        chi_on_classes = []
        for c in classes:
            rep_key = list(c['keys'])[0]
            chi_on_classes.append(elem_vals[rep_key])
        chars[cname] = np.array(chi_on_classes)

    return chars


def main():
    print("=" * 70)
    print("SECTOR-RESOLVED SPECTRUM: RIGID ROTOR ON SO(3)/D₃")
    print("=" * 70)

    # Build 2D₃ and verify
    elements = build_binary_dihedral_3()
    assert len(elements) == 12, f"Expected 12 elements, got {len(elements)}"

    # Verify closure
    key_set = {qkey(e) for e in elements}
    for i in range(len(elements)):
        for j_idx in range(len(elements)):
            prod = qmul(elements[i], elements[j_idx])
            assert qkey(prod) in key_set, f"Closure failed: {i}×{j_idx}"
    print("\nClosure verified ✓")

    # Compute -1 ∈ [G,G]?
    minus_one = np.array([-1, 0, 0, 0], dtype=float)
    mk = qkey(minus_one)
    comm = compute_commutator_subgroup(elements)
    has_minus_one = mk in comm
    print(f"-1 ∈ [G,G]: {has_minus_one} → {'obstructed' if has_minus_one else 'NOT obstructed'}")
    print(f"|[G,G]| = {len(comm)}, |Abel(2D₃)| = {len(elements) // len(comm)}")

    # Conjugacy classes
    classes = compute_conjugacy_classes(elements)
    classes.sort(key=lambda c: -c['trace'])
    sizes = [c['size'] for c in classes]
    G = sum(sizes)

    print(f"\n2D₃: {len(classes)} classes, |G| = {G}")
    for i, c in enumerate(classes):
        print(f"  C{i}: size={c['size']}, trace={c['trace']:.4f}")

    # Build 4 characters of Z₄ quotient
    chars = assign_character_values(elements, classes)
    sizes_arr = np.array(sizes, dtype=float)

    # Verify orthogonality
    print("\nCharacter orthogonality:")
    names = list(chars.keys())
    for i, n1 in enumerate(names):
        for j_idx, n2 in enumerate(names):
            if j_idx >= i:
                inner = np.sum(chars[n1].conj() * chars[n2] * sizes_arr) / G
                if abs(inner) > 1e-10 or i == j_idx:
                    expected = 1.0 if i == j_idx else 0.0
                    status = "✓" if abs(inner - expected) < 1e-10 else "✗"
                    print(f"  ⟨{n1}, {n2}⟩ = {inner.real:.4f} {status}")

    # Classify tensorial vs spinorial
    print("\nCharacter values at -1:")
    for name, chi in chars.items():
        for i, c in enumerate(classes):
            if abs(c['trace'] + 2) < 1e-10:
                val = chi[i]
                typ = "tensorial" if val.real > 0 else "SPINORIAL"
                print(f"  {name}: χ(-1) = {val:+.0f} → {typ}")
                break

    # Full multiplicity table
    print("\n" + "=" * 70)
    print("FULL MULTIPLICITY TABLE")
    print("=" * 70)

    header = f"{'j':>5}  {'dim':>4}  {'type':>5}"
    for name in names:
        header += f"  {name:>6}"
    header += f"  {'E':>9}"
    print(f"\n{header}")
    print("-" * len(header))

    for two_j in range(0, 14):
        j = two_j / 2.0
        chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
        typ = "(T)" if two_j % 2 == 0 else "(S)"
        E = j * (j + 1)

        line = f"  {j:4.1f}  {int(2*j+1):4d}  {typ:>5}"
        for name in names:
            m = np.sum(chars[name].conj() * chi_Vj * sizes_arr) / G
            line += f"  {round(m.real):6d}"
        line += f"  {E:9.2f}"
        print(line)

    # Physical spectrum per sector
    print("\n" + "=" * 70)
    print("PHYSICAL SPECTRUM PER SECTOR")
    print("=" * 70)

    for name, chi in chars.items():
        # Determine if tensorial or spinorial
        for i, c in enumerate(classes):
            if abs(c['trace'] + 2) < 1e-10:
                is_spinorial = chi[i].real < 0
                break

        j_type = "half-integer" if is_spinorial else "integer"
        label = "SPINORIAL" if is_spinorial else "tensorial"
        print(f"\nSector {name} ({label}):")
        print(f"  → only {j_type} j contribute")
        print(f"  {'j':>5}  {'E=j(j+1)':>9}  {'mult':>5}")
        print(f"  {'-'*25}")
        total = 0
        if is_spinorial:
            j_range = [k / 2.0 for k in range(1, 14, 2)]
        else:
            j_range = [float(k) for k in range(7)]
        for j in j_range:
            chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
            m = np.sum(chi.conj() * chi_Vj * sizes_arr) / G
            m_r = round(m.real)
            if m_r > 0:
                E = j * (j + 1)
                print(f"  {j:5.1f}  {E:9.2f}  {m_r:5d}")
                total += m_r
        print(f"  Total states (up to j=6.5): {total}")

    # Comparison summary
    print("\n" + "=" * 70)
    print("COMPARISON: SO(3)/D₃ vs SO(3)/O")
    print("=" * 70)

    print("""
SO(3)/O (obstructed, Q₈ ⊂ 2O):
  Abel(2O) ≅ Z₂ → 2 scalar sectors, BOTH tensorial
  All scalar states: integer j only
  Spinorial states FORBIDDEN as scalars → require rank-2 vector bundle H₁
  Scalar spectrum: A₀ at j=0,4,6; A₁ at j=3,6 (extremely sparse)

SO(3)/D₃ (non-obstructed, Q₈ ⊄ 2D₃):
  Abel(2D₃) ≅ Z₄ → 4 scalar sectors: 2 tensorial + 2 spinorial
  Spinorial scalar states at half-integer j!
  No vector bundle needed — spinorial wavefunctions exist as ordinary scalars

SPECTRAL MANIFESTATION OF Q₈ OBSTRUCTION:
  The obstruction doesn't just forbid spinorial sectors —
  it forces ALL scalar sectors to be tensorial,
  decimating the spectrum to integer j only.""")


if __name__ == '__main__':
    main()
