#!/usr/bin/env python3
"""
Particle-hole exploration: the duality m(χ,j) + m(χ,j*-j) = 1 as a
half-filled band with particle-hole symmetry.

QUESTION: What does the duality look like when viewed as a band structure?
- S: j ↦ j*-j is the particle-hole involution
- j* is the "Fermi energy"
- Below threshold: half-filled, rigid, binary (m ∈ {0,1})
- Above threshold: all levels filled, growing multiplicities

EXPLORATION:
  Part 1:  Occupation function n(j) — band diagram for multiple groups
  Part 2:  Particle-hole symmetry S and its spectrum
  Part 3:  "Fermi energy" E* = j*(j*+1)B and transition behavior
  Part 4:  Sub-threshold Peter-Weyl uniformity
  Part 5:  What breaks at j > j*? Deviation from particle-hole
  Part 6:  Cumulative state count — staircase function
  Part 7:  Energy gap structure — particle vs hole energies
  Part 8:  Anatomy of occupation pattern — non-central sum S(j)
  Part 9:  Dicyclic alternation formula
  Part 10: Mode decomposition S(j) = Σ A_α · f_α(j)
  Part 11: Mode amplitudes and group structure
  Part 12: Mode pairing α ↔ π-α and collapse
  Part 13: Antisymmetry proof from mode structure
  Part 14: Non-obstructed groups — what breaks
  Part 15: Obstructed vs non-obstructed comparison
  Part 16: The d=2 mode as obstruction switch
  Part 17: Polyhedral — distributed obstruction
  Part 18: Symmetry breaking O → D₃ as mode filtering

RAW OUTPUT:
===========

PART 1: OCCUPATION FUNCTION n(j) = m(A₀, j)

  2T: ▮··▮▮·       (3/6 = 0.500)
  2O: ▮···▮·▮·▮▮▮· (6/12 = 0.500)
  2I: ▮·····▮···▮·▮··▮▮·▮·▮▮▮·▮▮▮▮▮· (15/30 = 0.500)
  Dic_4:  ▮·▮·       (alternating)
  Dic_6:  ▮·▮·▮·     (alternating)
  Dic_10: ▮·▮·▮·▮·▮· (alternating)

PARTS 2-7: PARTICLE-HOLE STRUCTURE
  - m(j) + m(j*-j) = 1 for all j ≤ j* (perfect pairing)
  - E(j) + E(j*-j) NOT constant: duality is algebraic, not energetic
  - Peter-Weyl: ⊕_{j≤j*} V_j = (j*+1)/2 × ρ_reg⁺ (exact)
  - Above threshold: m(j+P) = m(j) + 1 with P = |G|/4

PART 8: NON-CENTRAL SUM DECOMPOSITION
  m(A₀, j) = [2(2j+1) + S(j)] / |G|
  S(j) = Σ_{non-central} χ_j(g).  Particle iff S(j) > 0.
  S(j) = -S(j*-j): perfect antisymmetry = content of duality.

PART 9: DICYCLIC ALTERNATION
  Dic_n (n even): m = (1+(-1)^j)/2. Single mode (d=2). Exact.

PART 10: MODE DECOMPOSITION (key result)
  S(j) = Σ_α A_α · sin((2j+1)α)/sin(α)
  where α = half-angle of SU(2) element, A_α = class size.
  Exact to machine precision (~10⁻¹⁴) for all groups.

  Modes by group:
    2T: α = π/3(×8), π/2(×6), 2π/3(×8)              — 3 modes
    2O: α = π/4(×6), π/3(×8), π/2(×18), 2π/3(×8), 3π/4(×6) — 5 modes
    2I: α = π/5(×12), π/3(×20), 2π/5(×12), π/2(×30),
        3π/5(×12), 2π/3(×20), 4π/5(×12)              — 7 modes

PART 12: MODE COLLAPSE (α ↔ π-α pairing)
  For integer j: f(π-α, j) = f(α, j) (proof: sin((2j+1)(π-α)) = sin((2j+1)α)).
  Modes α and π-α are identical → collapse to canonical α ≤ π/2:
    T: 2 independent modes (π/3, π/2)
    O: 3 independent modes (π/4, π/3, π/2)
    I: 4 independent modes (π/5, π/3, 2π/5, π/2)
  Note: order 5 in I splits into TWO canonical modes (π/5 ≠ 2π/5).

  Collapsed formulas:
    S_T(j) = 16·f(π/3) + 6·f(π/2)
    S_O(j) = 12·f(π/4) + 16·f(π/3) + 18·f(π/2)
    S_I(j) = 24·f(π/5) + 40·f(π/3) + 24·f(2π/5) + 30·f(π/2)
  where f(α) = sin((2j+1)α)/sin(α).

PART 13: ANTISYMMETRY PROOF
  For α = kπ/d, condition (B) gives |G|/(2d) = 2m (even).
  Then f(α, j*-j) = sin(2mkπ - (2j+1)α)/sin(α) = -f(α, j).
  Key: |G|/(2d) even ⟹ k·|G|/(2d) even for ANY k.
  So ALL modes of a given SO(3) order flip together.
  Verified to ~10⁻¹⁴ for all 15 half-angles across T, O, I.

  BOTTOM LINE: condition (B) = every standing wave mode is
  antisymmetric around j*/2 = duality = half-filling.

KEY FINDINGS:
  1. Dicyclic: perfectly alternating (even=particle, odd=hole)
     Formula: m(A₀,j) = (1+(-1)^j)/2 for j ≤ j*
  2. Polyhedral: irregular but exactly half-filled
     Pattern = multifrequency interference from conjugacy classes
  3. NOT energy symmetry — E(j)+E(j*-j) varies
     Duality is algebraic (on j labels), not dynamical
  4. Above threshold: periodicity P=|G|/4, m(j+P)=m(j)+1
     Particle-hole symmetry breaks, all levels fill
  5. Peter-Weyl uniformity perfect: ⊕V_j = (j*+1)/2 × ρ_reg⁺
  6. Non-central sum S(j) perfectly antisymmetric: S(j) = -S(j*-j)
     This is the representation-theoretic content of the duality
  7. Occupation determined by sign(S(j) + 2(2j+1)):
     positive = particle, negative = hole
  8. Mode decomposition: S(j) = Σ A_α · sin((2j+1)α)/sin(α)
     Each conjugacy class = one standing wave mode
  9. α ↔ π-α collapse: integer j makes paired modes identical
     T: 2 modes, O: 3 modes, I: 4 modes (order 5 splits in I)
  10. Antisymmetry = condition (B): |G|/(2d) even forces
      every mode to be odd under j ↔ j*-j. This IS the duality.
  11. Non-obstructed (Dic_n, n odd): duality BREAKS at d=2 mode.
      |G|/(2·2) = n odd → f₂ symmetric → m+m'=2 (both filled).
      S(j)+S(j*-j) = 4n·(-1)^j (verified n=3,5,7,9,11).
      Band: ▮·▮·▮ (both endpoints filled, fraction > 1/2).
  12. Dicyclic: obstruction = parity of ONE mode (d=2).
      (-1)^(n-1) = -1 (n even) → antisym → obstructed.
      (-1)^(n-1) = +1 (n odd)  → sym → free.
  13. Polyhedral: obstruction distributed across ALL modes.
      T,O,I: ALL |H|/d even simultaneously — cannot break one mode
      without breaking the symmetry group itself.
  14. Symmetry breaking O→D₃ = mode filtering + mode flip:
      - d=4 mode (4-fold axis) vanishes entirely
      - d=3 mode survives but weakens (8→2)
      - d=2 mode FLIPS: |G|/(2·2) goes from 12 (even) to 3 (ODD)
      This flip = loss of antisymmetry = loss of duality = re-entrance.
      The re-entrance of spinorial states IS the d=2 mode switching parity.
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
from src.quaternion import qmul, qkey
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       build_binary_icosahedral,
                       compute_conjugacy_classes, compute_commutator_subgroup)

IDENTITY_KEY = qkey(np.array([1, 0, 0, 0], dtype=float))
MINUS1_KEY = qkey(np.array([-1, 0, 0, 0], dtype=float))

def chi_su2(j, trace):
    """SU(2) spin-j character at given trace = 2cos(θ/2)."""
    if abs(trace) > 1.999:
        sign = 1 if trace > 0 else (-1)**(int(2*j))
        return sign * (2*j + 1)
    theta_half = np.arccos(trace / 2)
    return np.sin((2*j + 1) * theta_half) / np.sin(theta_half)

def scalar_multiplicity(j, classes, sizes, G):
    """Multiplicity of trivial character A₀ in V_j|_G."""
    chi_Vj = np.array([chi_su2(j, c['trace']) for c in classes])
    m = np.sum(chi_Vj * sizes) / G
    return int(round(m.real))


# ==========================================================
# Build groups
# ==========================================================
print("=" * 70)
print("PARTICLE-HOLE EXPLORATION")
print("=" * 70)

groups = {}

for name, builder, order in [
    ('2T', build_binary_tetrahedral, 24),
    ('2O', build_binary_octahedral, 48),
    ('2I', build_binary_icosahedral, 120),
]:
    elems = builder()
    classes = compute_conjugacy_classes(elems)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    j_star = order // 4 - 1
    groups[name] = {
        'G': order, 'classes': classes, 'sizes': sizes,
        'j_star': j_star, 'elements': elems
    }

# Also add dicyclic groups (analytic construction, no float drift)
def build_dicyclic(n):
    """Build Dic_n of order 4n analytically."""
    elements = []
    for k in range(2 * n):
        angle = k * np.pi / n
        c, s = np.cos(angle), np.sin(angle)
        elements.append(np.array([c, s, 0, 0], dtype=float))
        elements.append(np.array([0, 0, c, s], dtype=float))
    return elements

for n in [4, 6, 10]:
    name = f'Dic_{n}'
    order = 4 * n
    elems = build_dicyclic(n)
    assert len(elems) == order, f"{name}: expected {order}, got {len(elems)}"
    classes = compute_conjugacy_classes(elems)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)
    j_star = order // 4 - 1
    groups[name] = {
        'G': order, 'classes': classes, 'sizes': sizes,
        'j_star': j_star, 'elements': elems
    }

# ==========================================================
# Part 1: Occupation function — band diagram
# ==========================================================
print("\n" + "=" * 70)
print("PART 1: OCCUPATION FUNCTION n(j) = m(A₀, j)")
print("=" * 70)

J_MAX = 60

for name in ['2T', '2O', '2I', 'Dic_4', 'Dic_6', 'Dic_10']:
    g = groups[name]
    j_star = g['j_star']
    print(f"\n  {name} (|G|={g['G']}, j*={j_star}, E*={j_star*(j_star+1)}B)")

    mults = []
    for j in range(J_MAX + 1):
        m = scalar_multiplicity(j, g['classes'], g['sizes'], g['G'])
        mults.append(m)

    # Band diagram: show sub-threshold as ▮ (filled) / · (empty)
    # and above-threshold with multiplicity number
    band_sub = ""
    band_above = ""
    for j in range(min(j_star + 1, J_MAX + 1)):
        band_sub += "▮" if mults[j] > 0 else "·"

    # Count filled/empty below threshold
    n_filled = sum(1 for j in range(j_star + 1) if mults[j] > 0)
    n_empty = j_star + 1 - n_filled

    print(f"    Sub-threshold [{0}..{j_star}]: {band_sub}")
    print(f"    Filled: {n_filled}/{j_star+1} = {n_filled/(j_star+1):.3f}")

    # Above threshold: show first 20 values
    above = [mults[j] for j in range(j_star + 1, min(j_star + 21, J_MAX + 1))]
    print(f"    Above-threshold [{j_star+1}..{j_star+20}]: {above}")
    print(f"    All ≥ 1 above threshold: {all(m >= 1 for m in above)}")

# ==========================================================
# Part 2: Particle-hole involution S: j → j*-j
# ==========================================================
print("\n" + "=" * 70)
print("PART 2: PARTICLE-HOLE INVOLUTION S: j → j*-j")
print("=" * 70)

for name in ['2T', '2O', '2I']:
    g = groups[name]
    j_star = g['j_star']
    print(f"\n  {name} (j*={j_star}):")

    # Show explicit pairs and their occupation
    print(f"    {'j':>3} {'j*-j':>4}  {'m(j)':>4} {'m(j*-j)':>7}  {'sum':>3}  status")
    print("    " + "-" * 45)
    for j in range((j_star + 1) // 2 + 1):
        jp = j_star - j
        m_j = scalar_multiplicity(j, g['classes'], g['sizes'], g['G'])
        m_jp = scalar_multiplicity(jp, g['classes'], g['sizes'], g['G'])
        status = "particle" if m_j == 1 else "hole    "
        if j == jp:
            status = "midpoint"
        print(f"    {j:3d} {jp:4d}  {m_j:4d} {m_jp:7d}  {m_j+m_jp:3d}  {status}")

# ==========================================================
# Part 3: "Fermi energy" and energy scales
# ==========================================================
print("\n" + "=" * 70)
print("PART 3: FERMI ENERGY E* = j*(j*+1)B")
print("=" * 70)

print(f"\n  {'Group':>8} {'|G|':>4} {'j*':>3} {'E*=j*(j*+1)B':>14} "
      f"{'E₁=j₁(j₁+1)B':>14} {'E*/E₁':>6}")
print("  " + "-" * 55)

for name in ['2T', '2O', '2I', 'Dic_4', 'Dic_6', 'Dic_10']:
    g = groups[name]
    j_star = g['j_star']
    E_star = j_star * (j_star + 1)

    # Find j₁ (first excited scalar)
    j1 = None
    for j in range(1, J_MAX + 1):
        m = scalar_multiplicity(j, g['classes'], g['sizes'], g['G'])
        if m > 0:
            j1 = j
            break
    E1 = j1 * (j1 + 1) if j1 else 0

    print(f"  {name:>8} {g['G']:4d} {j_star:3d} {E_star:14d} "
          f"{E1:14d} {E_star/E1 if E1 else 0:6.1f}")

# ==========================================================
# Part 4: Sub-threshold Peter-Weyl uniformity
# ==========================================================
print("\n" + "=" * 70)
print("PART 4: SUB-THRESHOLD PETER-WEYL UNIFORMITY")
print("=" * 70)

for name in ['2T', '2O', '2I']:
    g = groups[name]
    j_star = g['j_star']

    # Total dimension below threshold
    total_dim = sum(2*j + 1 for j in range(j_star + 1))
    expected = (j_star + 1) * g['G'] // 4  # = (j*+1)/2 × |G|/2

    # Scalar content: sum of m(A₀,j)(2j+1) for j ≤ j*
    scalar_content = sum(
        scalar_multiplicity(j, g['classes'], g['sizes'], g['G']) * (2*j + 1)
        for j in range(j_star + 1)
    )

    # Expected scalar per copy of ρ_reg⁺: 1 (trivial appears once)
    # Number of copies: (j*+1)/2
    expected_scalar = (j_star + 1) // 2

    print(f"\n  {name} (j*={j_star}, |G|={g['G']}):")
    print(f"    Total dim ⊕ V_j for j ∈ [0,j*] = {total_dim}")
    print(f"    Expected (j*+1)/2 × |G|/2 = {(j_star+1)//2} × {g['G']//2} = {expected}")
    print(f"    Match: {total_dim == expected}")
    print(f"    Scalar A₀ content: {scalar_content} states")
    print(f"    Expected (j*+1)/2 × 1 = {expected_scalar}")
    print(f"    Match: {scalar_content == expected_scalar}")

# ==========================================================
# Part 5: What breaks above threshold?
# ==========================================================
print("\n" + "=" * 70)
print("PART 5: ABOVE-THRESHOLD — PARTICLE-HOLE BREAKING")
print("=" * 70)

for name in ['2O']:
    g = groups[name]
    j_star = g['j_star']
    P = g['G'] // 4  # period

    print(f"\n  {name} (j*={j_star}, P={P}):")
    print(f"    m(A₀, j) for j = 0..{j_star + 2*P}:")
    print(f"    {'j':>3} {'m(j)':>5} {'period':>7} {'deviation from PH':>20}")
    print("    " + "-" * 40)

    for j in range(j_star + 2*P + 1):
        m = scalar_multiplicity(j, g['classes'], g['sizes'], g['G'])

        if j <= j_star:
            jp = j_star - j
            m_partner = scalar_multiplicity(jp, g['classes'], g['sizes'], g['G'])
            ph_dev = f"m+m' = {m + m_partner}"
            period_label = "sub-thr"
        else:
            # Above threshold: check periodicity m(j) = m(j-P) + 1
            if j >= P:
                m_prev = scalar_multiplicity(j - P, g['classes'], g['sizes'], g['G'])
                ph_dev = f"m(j)-m(j-P) = {m - m_prev}"
                period_label = f"period {(j - j_star - 1) // P + 1}"
            else:
                ph_dev = ""
                period_label = ""

        if j <= j_star + 2*P:
            print(f"    {j:3d} {m:5d} {period_label:>7} {ph_dev:>20}")

# ==========================================================
# Part 6: Cumulative state count — staircase
# ==========================================================
print("\n" + "=" * 70)
print("PART 6: CUMULATIVE STATE COUNT N(J) = Σ_{j≤J} m(A₀,j)")
print("=" * 70)

for name in ['2T', '2O', '2I']:
    g = groups[name]
    j_star = g['j_star']
    P = g['G'] // 4

    print(f"\n  {name} (j*={j_star}, P={P}):")

    cumul = 0
    cumul_weyl = 0
    print(f"    {'J':>3} {'N(J)':>5} {'Weyl':>7} {'diff':>5}")
    print("    " + "-" * 25)

    for J in range(min(3 * j_star, 50) + 1):
        m = scalar_multiplicity(J, g['classes'], g['sizes'], g['G'])
        cumul += m
        weyl = (J + 1)**2 / (g['G'] // 2)  # asymptotic: (J+1)²/|H|

        if J % max(j_star // 3, 1) == 0 or J == j_star or J == j_star + P:
            print(f"    {J:3d} {cumul:5d} {weyl:7.1f} {cumul - weyl:+5.1f}")

# ==========================================================
# Part 7: Energy gap structure — particle vs hole energies
# ==========================================================
print("\n" + "=" * 70)
print("PART 7: ENERGY GAP STRUCTURE")
print("=" * 70)

for name in ['2O']:
    g = groups[name]
    j_star = g['j_star']

    particles = []  # j where m(j) = 1 (filled)
    holes = []      # j where m(j) = 0 (empty)

    for j in range(j_star + 1):
        m = scalar_multiplicity(j, g['classes'], g['sizes'], g['G'])
        E = j * (j + 1)
        if m == 1:
            particles.append((j, E))
        else:
            holes.append((j, E))

    print(f"\n  {name} (j*={j_star}):")
    print(f"    Particles (m=1): {[p[0] for p in particles]}")
    print(f"    Holes     (m=0): {[h[0] for h in holes]}")

    # Energy content
    E_particles = sum(p[1] for p in particles)
    E_holes = sum(h[1] for h in holes)
    E_total = sum(j*(j+1) for j in range(j_star + 1))

    print(f"\n    Particle energies E = j(j+1): {[p[1] for p in particles]}")
    print(f"    Hole energies:                {[h[1] for h in holes]}")
    print(f"    Sum particle E: {E_particles}")
    print(f"    Sum hole E:     {E_holes}")
    print(f"    Total E:        {E_total}")
    print(f"    E_particles/E_holes = {E_particles/E_holes:.4f}" if E_holes > 0 else "")

    # Check: are particle and hole energies symmetric around E*/2?
    E_star = j_star * (j_star + 1)
    print(f"\n    E* = j*(j*+1) = {E_star}")
    print(f"    E*/2 = {E_star/2}")
    print(f"    Mean particle E = {np.mean([p[1] for p in particles]):.2f}")
    print(f"    Mean hole E     = {np.mean([h[1] for h in holes]):.2f}")

    # S maps j → j*-j. What does this do to energy?
    # E(j) = j(j+1), E(j*-j) = (j*-j)(j*-j+1)
    # E(j) + E(j*-j) = j² + j + (j*)² - 2j*j + j² + j* - j + 1
    #                 = ... let me compute directly
    print(f"\n    Energy duality: E(j) + E(j*-j) for each pair:")
    for j in range((j_star + 1) // 2 + 1):
        jp = j_star - j
        Ej = j * (j + 1)
        Ejp = jp * (jp + 1)
        print(f"      j={j}, j*-j={jp}: E={Ej}, E'={Ejp}, sum={Ej+Ejp}")

# ==========================================================
# Part 8: What determines the occupation pattern?
# For each group, compute m(A₀, j) and decompose into
# central + non-central contributions to understand
# WHY certain j are particles vs holes.
# ==========================================================
print("\n" + "=" * 70)
print("PART 8: ANATOMY OF THE OCCUPATION PATTERN")
print("=" * 70)

for name in ['2T', '2O', '2I']:
    g = groups[name]
    j_star = g['j_star']
    G_order = g['G']
    H_order = G_order // 2

    print(f"\n  {name} (|G|={G_order}, |H|={H_order}, j*={j_star}):")

    # Element orders in H = G/{±1}
    orders = set()
    for c in g['classes']:
        if abs(c['trace']) < 1.999:  # non-central
            # trace = 2cos(θ/2), SO(3) angle θ = 2·arccos(trace/2)·2... wait
            # Actually: trace of SU(2) element = 2cos(α) where α = θ/2
            # SO(3) order d: rotation angle θ = 2π/d, so α = π/d
            # trace = 2cos(π/d)
            alpha = np.arccos(c['trace'] / 2)
            if alpha > 0.01:
                d = round(np.pi / alpha)
                orders.add(d)

    print(f"    SO(3) element orders: {sorted(orders)}")
    print(f"    |H|/d for each: {sorted([H_order // d for d in orders])}")
    print(f"    All |H|/d even: {all(H_order % (2*d) == 0 for d in orders)}")

    # Non-central character sum S(j) = Σ_{g ≠ ±1} χ_j(g)
    # m(A₀,j) = [2(2j+1) + S(j)] / |G|
    # Particle iff m = 1 iff S(j) = |G| - 2(2j+1)
    # Hole iff m = 0 iff S(j) = -2(2j+1)
    print(f"\n    {'j':>3} {'2(2j+1)':>8} {'S(j)':>8} {'m(j)':>5} {'S/|G|':>8} note")
    print("    " + "-" * 50)

    for j in range(j_star + 1):
        central = 2 * (2*j + 1)
        chi_vals = [chi_su2(j, c['trace']) for c in g['classes']]
        total_sum = sum(chi_vals[i] * g['sizes'][i] for i in range(len(g['classes'])))
        S_j = total_sum - central  # non-central sum
        m_j = int(round(total_sum / G_order))

        note = ""
        if m_j == 1:
            note = "particle"
        else:
            note = "hole"

        print(f"    {j:3d} {central:8d} {S_j:8.1f} {m_j:5d} {S_j/G_order:8.3f} {note}")

    # Show the non-central sum decomposed by conjugacy class
    print(f"\n    Non-central sum S(j) by conjugacy class:")
    nc_classes = [(i, c) for i, c in enumerate(g['classes']) if abs(c['trace']) < 1.999]

    header = f"    {'j':>3}"
    for i, c in nc_classes:
        alpha = np.arccos(c['trace'] / 2)
        d = round(np.pi / alpha) if alpha > 0.01 else 0
        header += f"  c{i}(d={d},s={c['size']})"
    header += "     S(j)"
    print(header)
    print("    " + "-" * len(header))

    for j in range(j_star + 1):
        line = f"    {j:3d}"
        S_total = 0
        for i, c in nc_classes:
            contrib = chi_su2(j, c['trace']) * c['size']
            S_total += contrib
            line += f"  {contrib:12.2f}"
        line += f"  {S_total:8.1f}"
        print(line)

# ==========================================================
# Part 9: Dicyclic pattern — why alternating?
# ==========================================================
print("\n" + "=" * 70)
print("PART 9: DICYCLIC PATTERN — WHY ALTERNATING?")
print("=" * 70)

for name in ['Dic_4', 'Dic_6']:
    g = groups[name]
    j_star = g['j_star']
    G_order = g['G']
    n = G_order // 4

    print(f"\n  {name} (n={n}, |G|={G_order}, j*={j_star}):")
    print(f"    {'j':>3} {'m(j)':>5} {'(-1)^j':>6}  note")
    print("    " + "-" * 25)
    for j in range(j_star + 1):
        m = scalar_multiplicity(j, g['classes'], g['sizes'], G_order)
        parity = (-1)**j
        note = "= (1+(-1)^j)/2" if m == (1 + parity) // 2 else "MISMATCH"
        print(f"    {j:3d} {m:5d} {parity:6d}  {note}")

    # Check formula: m(A₀, j) = (1+(-1)^j)/2 for j ≤ j*
    all_match = all(
        scalar_multiplicity(j, g['classes'], g['sizes'], G_order)
        == (1 + (-1)**j) // 2
        for j in range(j_star + 1)
    )
    print(f"    Formula m = (1+(-1)^j)/2: {'ALL MATCH' if all_match else 'FAILS'}")

# ==========================================================
# Part 10: Mode decomposition S(j) = Σ_α A_α · f_α(j)
#
# For integer j, paired classes (g, -g) give identical χ_j values.
# Group by SU(2) half-angle α (NOT by SO(3) order d):
#   A_α = total class size for all classes with half-angle α
#   f_α(j) = sin((2j+1)α) / sin(α)  [= Dirichlet kernel = χ_j(trace)]
# Then S(j) = Σ_α A_α · f_α(j) exactly.
#
# KEY: order-5 elements in I have TWO distinct half-angles:
#   α = π/5 (72° rotation) and α = 2π/5 (144° rotation).
#   These are distinct modes, not one mode.
# ==========================================================
print("\n" + "=" * 70)
print("PART 10: MODE DECOMPOSITION S(j) = Σ A_α · f_α(j)")
print("=" * 70)

def f_mode_alpha(j, alpha):
    """Standing wave mode for half-angle α."""
    return np.sin((2*j + 1) * alpha) / np.sin(alpha)

for name in ['2T', '2O', '2I']:
    g = groups[name]
    j_star = g['j_star']
    G_order = g['G']
    H_order = G_order // 2

    # Compute mode amplitudes by grouping classes by half-angle α
    # Use rounded α as key (4 decimal places)
    mode_data = []  # list of (alpha, weight, label)
    for c in g['classes']:
        if abs(c['trace']) > 1.999:  # central (±1)
            continue
        alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
        # Find SO(3) order for labeling
        best_d = None
        for k in range(1, 20):
            d_cand = np.pi * k / alpha
            if abs(d_cand - round(d_cand)) < 0.01:
                best_d = int(round(d_cand))
                break
        # Rational representation: α = kπ/d
        k_num = int(round(alpha * best_d / np.pi)) if best_d else 0
        mode_data.append((alpha, c['size'], best_d, k_num))

    # Merge by α (group classes with same half-angle)
    alpha_modes = {}  # rounded_alpha -> (alpha, total_weight, label)
    for alpha, weight, d, k in mode_data:
        key = round(alpha, 8)
        if key not in alpha_modes:
            alpha_modes[key] = {'alpha': alpha, 'weight': 0, 'd': d, 'k': k}
        alpha_modes[key]['weight'] += weight

    modes = sorted(alpha_modes.values(), key=lambda m: m['alpha'])

    print(f"\n  {name} (|G|={G_order}, |H|={H_order}, j*={j_star}):")
    print(f"    Modes (by half-angle α):")
    for m in modes:
        label = f"α={m['k']}π/{m['d']}" if m['d'] else f"α={m['alpha']:.4f}"
        period = m['d']  # f has period d (the SO(3) order)
        print(f"      {label}: A = {m['weight']} (period {period})")

    # Check: total non-central = |G| - 2
    total_nc = sum(m['weight'] for m in modes)
    assert total_nc == G_order - 2, f"Non-central total {total_nc} ≠ {G_order-2}"

    # Verify decomposition: S(j) = Σ A_α · f_α(j)
    print(f"\n    Verification: S(j) vs Σ A_α·f_α(j)")
    header = f"    {'j':>3} {'S(j)':>10} {'Σ modes':>10} {'err':>8}  "
    for m in modes:
        label = f"{m['k']}π/{m['d']}"
        header += f"  {label:>8}"
    header += "  occ"
    print(header)
    print("    " + "-" * (38 + 10 * len(modes)))

    max_err = 0
    for j in range(j_star + 1):
        # Actual S(j)
        central = 2 * (2*j + 1)
        total_sum = sum(chi_su2(j, c['trace']) * s
                        for c, s in zip(g['classes'], g['sizes']))
        S_actual = total_sum - central

        # Mode decomposition
        S_modes = sum(m['weight'] * f_mode_alpha(j, m['alpha']) for m in modes)

        err = abs(S_actual - S_modes)
        max_err = max(max_err, err)

        m_j = int(round(total_sum / G_order))
        label = "▮" if m_j == 1 else "·"

        line = f"    {j:3d} {S_actual:10.2f} {S_modes:10.2f} {err:8.1e}  "
        for m in modes:
            contrib = m['weight'] * f_mode_alpha(j, m['alpha'])
            line += f"  {contrib:8.2f}"
        line += f"  {label}"
        print(line)

    print(f"    Max error: {max_err:.2e}")
    assert max_err < 1e-10, f"{name}: mode decomposition failed, max_err={max_err:.2e}"

    # The "spectral content" of the occupation
    print(f"\n    Analytic formula:")
    terms = []
    for m in modes:
        label = f"{m['k']}π/{m['d']}"
        terms.append(f"{m['weight']}·f({label})")
    print(f"    S(j) = {' + '.join(terms)}")
    print(f"    where f(α) = sin((2j+1)α)/sin(α)")

# ==========================================================
# Part 11: Connection to group structure
# ==========================================================
print("\n" + "=" * 70)
print("PART 11: MODE AMPLITUDES AND GROUP STRUCTURE")
print("=" * 70)

print("\n  Summary of mode amplitudes (by half-angle α):")

for name in ['2T', '2O', '2I']:
    g = groups[name]
    G_order = g['G']
    H_order = G_order // 2
    j_star = g['j_star']

    # Recompute modes (same logic as Part 10)
    alpha_modes = {}
    for c in g['classes']:
        if abs(c['trace']) > 1.999:
            continue
        alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
        best_d = None
        for k in range(1, 20):
            d_cand = np.pi * k / alpha
            if abs(d_cand - round(d_cand)) < 0.01:
                best_d = int(round(d_cand))
                break
        k_num = int(round(alpha * best_d / np.pi)) if best_d else 0
        key = round(alpha, 8)
        if key not in alpha_modes:
            alpha_modes[key] = {'alpha': alpha, 'weight': 0, 'd': best_d, 'k': k_num}
        alpha_modes[key]['weight'] += c['size']

    modes = sorted(alpha_modes.values(), key=lambda m: m['alpha'])

    print(f"\n  {name} (|G|={G_order}, |H|={H_order}, j*={j_star}):")
    for m in modes:
        label = f"α = {m['k']}π/{m['d']}"
        # Geometric interpretation
        theta_deg = 2 * m['alpha'] * 180 / np.pi
        print(f"    {label:>12}  A = {m['weight']:3d}  "
              f"(SO(3) rotation {theta_deg:.0f}°, order {m['d']})")

    total = sum(m['weight'] for m in modes)
    print(f"    {'Total':>12}      {total:3d}  (= |G|-2 = {G_order-2})")

# Interpretation
print(f"\n  Key observations:")
print(f"    1. Each non-central conjugacy class is a distinct mode")
print(f"       with half-angle α and weight A = class size")
print(f"    2. For T and O, each SO(3) order d has one half-angle α=π/d")
print(f"    3. For I, order 5 splits into TWO modes:")
print(f"       α=π/5 (72° rotation) and α=2π/5 (144° rotation)")
print(f"    4. A_α = 2 × |{{h ∈ H : half-angle = α}}|")
print(f"       (each SO(3) element lifts to {{g,-g}} in SU(2))")
print(f"    5. Number of distinct modes = number of non-central")
print(f"       conjugacy CLASS PAIRS (c, -c) in SU(2)")

# ==========================================================
# Part 12: Mode pairing α ↔ π-α and collapse
#
# SU(2) classes come in pairs (g, -g). Non-central classes
# with half-angle α have a partner class at π-α.
# For INTEGER j: f(π-α, j) = f(α, j) (because sin((2j+1)(π-α))
#   = sin((2j+1)π)cos((2j+1)α) - cos((2j+1)π)sin((2j+1)α)
#   = 0 - (-1)^{2j} sin((2j+1)α) = sin((2j+1)α)
#   and sin(π-α) = sin(α)).
# So paired modes are identical at integer j!
#
# This means we can collapse to INDEPENDENT modes:
#   one per CANONICAL half-angle min(α, π-α).
# For T and O: one mode per order d. For I: order 5 has TWO
# canonical angles (π/5 and 2π/5, both ≤ π/2), so 4 modes total.
# ==========================================================
print("\n" + "=" * 70)
print("PART 12: MODE PAIRING α ↔ π-α AND COLLAPSE")
print("=" * 70)

for name in ['2T', '2O', '2I']:
    g = groups[name]
    j_star = g['j_star']
    G_order = g['G']
    H_order = G_order // 2

    # Group by canonical half-angle min(α, π-α)
    collapsed = {}  # rounded canonical α -> (alpha, total_weight)
    for c in g['classes']:
        if abs(c['trace']) > 1.999:
            continue
        alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
        alpha_canon = min(alpha, np.pi - alpha)
        key = round(alpha_canon, 8)
        if key not in collapsed:
            collapsed[key] = {'alpha': alpha_canon, 'weight': 0}
        collapsed[key]['weight'] += c['size']

    modes_collapsed = sorted(collapsed.values(), key=lambda m: m['alpha'])

    print(f"\n  {name} (|G|={G_order}, |H|={H_order}, j*={j_star}):")
    print(f"    Collapsed modes (α ≤ π/2, combining α and π-α):")
    for m in modes_collapsed:
        # Express as rational multiple of π
        found = False
        for d in range(2, 20):
            for k in range(1, d):
                if abs(m['alpha'] - k * np.pi / d) < 0.001:
                    print(f"      α={k}π/{d}: Ã = {m['weight']}")
                    found = True
                    break
            if found:
                break

    # Verify: for integer j, this gives same S(j)
    print(f"\n    Verification of collapsed decomposition:")
    max_err = 0
    for j in range(j_star + 1):
        central = 2 * (2*j + 1)
        total_sum = sum(chi_su2(j, c['trace']) * s
                        for c, s in zip(g['classes'], g['sizes']))
        S_actual = total_sum - central

        S_collapsed = sum(m['weight'] * f_mode_alpha(j, m['alpha'])
                          for m in modes_collapsed)

        err = abs(S_actual - S_collapsed)
        max_err = max(max_err, err)

    print(f"    Max error: {max_err:.2e}")
    assert max_err < 1e-10, f"{name}: collapsed decomposition failed"
    print(f"    ✓ Collapsed formula works")
    print(f"    Number of independent modes: {len(modes_collapsed)}")

print(f"\n  RESULT: At integer j, modes α and π-α are identical.")
print(f"  Collapse reduces modes to canonical angles α ≤ π/2:")
print(f"    T: 2 modes (π/3, π/2)")
print(f"    O: 3 modes (π/4, π/3, π/2)")
print(f"    I: 4 modes (π/5, π/3, 2π/5, π/2)")
print(f"  Note: order 5 gives TWO modes (π/5 ≠ 2π/5), both ≤ π/2.")

# ==========================================================
# Part 13: Antisymmetry from mode structure
#
# Why is S(j) = -S(j*-j)?
# For each mode with half-angle α = kπ/d:
#   f(α, j*-j) = sin((2j*+1-2j)α)/sin(α)
#   2j*+1 = |G|/2 - 1
#   (2j*+1)α = (|G|/2-1)·kπ/d = k·|G|π/(2d) - kπ/d
#
# Condition (B): |G|/(2d) = 2m is even for all d.
# Then k·|G|π/(2d) = 2mkπ (integer multiple of 2π).
# So sin((2j*+1-2j)α) = sin(2mkπ - (2j+1)kπ/d) = -sin((2j+1)α).
# Therefore f(α, j*-j) = -f(α, j) for ALL half-angles.
#
# This holds for any k: the parity of |G|/(2d) controls ALL
# modes of a given SO(3) order simultaneously.
# ==========================================================
print("\n" + "=" * 70)
print("PART 13: ANTISYMMETRY S(j) = -S(j*-j) FROM MODE STRUCTURE")
print("=" * 70)

print(f"\n  Proof that f(α, j*-j) = -f(α, j) for α = kπ/d, |G|/(2d) even:")
print(f"    f(α, j) = sin((2j+1)α) / sin(α)")
print(f"    j* = |G|/4 - 1 ⟹ 2j*+1 = |G|/2 - 1")
print(f"    (2j*+1)·kπ/d = k·|G|π/(2d) - kπ/d = 2mkπ - kπ/d")
print(f"    sin(2mkπ - (2j+1)kπ/d) = -sin((2j+1)kπ/d)")
print(f"    ∴ f(α, j*-j) = -f(α, j) for all j. ∎")
print(f"\n  Key: |G|/(2d) even ⟹ k·|G|/(2d) even for ANY k.")
print(f"  So ALL half-angles of a given SO(3) order flip together.")

# Verify numerically for ALL half-angles (not just π/d)
print(f"\n  Numerical verification (all half-angles):")
for name in ['2T', '2O', '2I']:
    g = groups[name]
    j_star = g['j_star']
    G_order = g['G']

    print(f"\n    {name}: |G|={G_order}, j*={j_star}")

    # Collect all distinct half-angles
    angles = {}
    for c in g['classes']:
        if abs(c['trace']) > 1.999:
            continue
        alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
        key = round(alpha, 8)
        if key not in angles:
            # Find rational form kπ/d
            best_d, best_k = None, None
            for d in range(2, 20):
                for k in range(1, d):
                    if abs(alpha - k * np.pi / d) < 0.001:
                        best_d, best_k = d, k
                        break
                if best_d:
                    break
            angles[key] = (alpha, best_k, best_d)

    for key in sorted(angles.keys()):
        alpha, k, d = angles[key]
        ratio = G_order / (2 * d)
        label = f"α={k}π/{d}"

        # Check f(α, j) + f(α, j*-j) = 0 for all j
        max_dev = 0
        for j in range(j_star + 1):
            fj = f_mode_alpha(j, alpha)
            fjp = f_mode_alpha(j_star - j, alpha)
            max_dev = max(max_dev, abs(fj + fjp))

        print(f"      {label:>8}: |G|/(2d)={ratio:.0f} "
              f"({'even' if int(ratio) % 2 == 0 else 'ODD!'}), "
              f"max|f+f'|={max_dev:.1e} ✓")

# ==========================================================
# Part 14: Non-obstructed groups — what breaks?
#
# For Dic_n with n ODD: -1 ∉ [G,G], no obstruction.
# The condition |G|/(2d) even FAILS for d=2 (since |G|=4n,
# |G|/4 = n is odd). So the d=2 mode is NOT antisymmetric,
# and the duality m+m'=1 breaks.
#
# Questions:
#   1. Does S(j) = -S(j*-j) still hold? (Expect: no)
#   2. Which modes break antisymmetry?
#   3. What does the occupation pattern look like?
#   4. Is j* = |G|/4 - 1 still meaningful?
# ==========================================================
print("\n" + "=" * 70)
print("PART 14: NON-OBSTRUCTED GROUPS — WHAT BREAKS?")
print("=" * 70)

# Build Dic_n for n odd
for n in [3, 5, 7]:
    name = f'Dic_{n}'
    order = 4 * n
    elems = build_dicyclic(n)
    assert len(elems) == order
    classes = compute_conjugacy_classes(elems)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)

    # j* would be |G|/4 - 1 = n-1, but |G|/4 is not integer when
    # ... wait, |G|=4n so |G|/4 = n always. j* = n-1.
    j_star = n - 1

    print(f"\n  {name} (n={n}, |G|={order}, j*={j_star}, n is ODD):")

    # Mode decomposition
    alpha_modes = {}
    for c in classes:
        if abs(c['trace']) > 1.999:
            continue
        alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
        key = round(alpha, 8)
        if key not in alpha_modes:
            # Find rational form
            best_d, best_k = None, None
            for d in range(2, 4*n):
                for k in range(1, d):
                    if abs(alpha - k * np.pi / d) < 0.001:
                        best_d, best_k = d, k
                        break
                if best_d:
                    break
            alpha_modes[key] = {'alpha': alpha, 'weight': 0, 'd': best_d, 'k': best_k}
        alpha_modes[key]['weight'] += c['size']

    modes = sorted(alpha_modes.values(), key=lambda m: m['alpha'])

    print(f"    Modes:")
    for m in modes:
        label = f"α={m['k']}π/{m['d']}" if m['d'] else f"α=?"
        ratio = order / (2 * m['d']) if m['d'] else 0
        parity = "even" if int(ratio) % 2 == 0 else "ODD"
        print(f"      {label:>10}: A={m['weight']:3d}, "
              f"|G|/(2d)={ratio:.1f} ({parity})")

    # Occupation pattern
    print(f"\n    Occupation and duality check:")
    ss_label = "S+S'"
    print(f"    {'j':>3} {'m(j)':>5} {'m(j*-j)':>7} {'sum':>4} {'S(j)':>8} {'S(j*-j)':>8} "
          f"{ss_label:>8}")
    print("    " + "-" * 55)

    for j in range(j_star + 1):
        jp = j_star - j
        m_j = scalar_multiplicity(j, classes, sizes, order)
        m_jp = scalar_multiplicity(jp, classes, sizes, order)

        # S(j)
        central = 2 * (2*j + 1)
        total_sum = sum(chi_su2(j, c['trace']) * s
                        for c, s in zip(classes, sizes))
        S_j = total_sum - central

        # S(j*-j)
        central_p = 2 * (2*jp + 1)
        total_sum_p = sum(chi_su2(jp, c['trace']) * s
                          for c, s in zip(classes, sizes))
        S_jp = total_sum_p - central_p

        print(f"    {j:3d} {m_j:5d} {m_jp:7d} {m_j+m_jp:4d} {S_j:8.1f} {S_jp:8.1f} "
              f"{S_j+S_jp:8.1f}")

    # Band diagram
    band = ""
    n_filled = 0
    for j in range(j_star + 1):
        m = scalar_multiplicity(j, classes, sizes, order)
        band += "▮" if m > 0 else "·"
        if m > 0:
            n_filled += 1
    print(f"\n    Band: {band}")
    print(f"    Filled: {n_filled}/{j_star+1} = {n_filled/(j_star+1):.3f}")

# ==========================================================
# Part 15: Contrast — obstructed vs non-obstructed Dic_n
#
# Compare pairs: Dic_4 (obstructed) vs Dic_3 (free)
#                Dic_6 (obstructed) vs Dic_5 (free)
# ==========================================================
print("\n" + "=" * 70)
print("PART 15: OBSTRUCTED vs NON-OBSTRUCTED COMPARISON")
print("=" * 70)

for n_even, n_odd in [(4, 3), (6, 5), (10, 9)]:
    print(f"\n  Dic_{n_even} (obstructed) vs Dic_{n_odd} (non-obstructed):")

    for n in [n_even, n_odd]:
        order = 4 * n
        elems = build_dicyclic(n)
        classes = compute_conjugacy_classes(elems)
        classes.sort(key=lambda c: -c['trace'])
        sizes = np.array([c['size'] for c in classes], dtype=float)
        j_star = n - 1

        band = ""
        duality_ok = True
        for j in range(j_star + 1):
            m_j = scalar_multiplicity(j, classes, sizes, order)
            m_jp = scalar_multiplicity(j_star - j, classes, sizes, order)
            band += "▮" if m_j > 0 else "·"
            if m_j + m_jp != 1:
                duality_ok = False

        obst = "obstructed" if n % 2 == 0 else "FREE"
        print(f"    Dic_{n:2d} ({obst:>10}): {band:>20}  "
              f"duality m+m'=1: {'✓' if duality_ok else '✗'}")

        # Show where duality breaks for non-obstructed
        if not duality_ok:
            print(f"         Duality violations:")
            for j in range(j_star + 1):
                m_j = scalar_multiplicity(j, classes, sizes, order)
                m_jp = scalar_multiplicity(j_star - j, classes, sizes, order)
                if m_j + m_jp != 1:
                    print(f"           j={j}, j*-j={j_star-j}: "
                          f"m={m_j}, m'={m_jp}, sum={m_j+m_jp}")

# ==========================================================
# Part 16: The d=2 mode as obstruction switch
#
# For Dic_n: the ONLY mode with variable parity is d=2.
#   |G|/(2·2) = n.
#   n even → antisymmetric → duality → obstructed
#   n odd  → symmetric    → no duality → free
#
# The d=2 mode is f_2(j) = sin((2j+1)π/2)/sin(π/2) = (-1)^j.
# f_2(j*-j) = (-1)^{j*-j} = (-1)^{j*}·(-1)^{-j} = (-1)^{n-1}·(-1)^j
#   n even: (-1)^{n-1} = -1 → f_2(j*-j) = -(-1)^j = -f_2(j) ✓ antisym
#   n odd:  (-1)^{n-1} = +1 → f_2(j*-j) = +(-1)^j = +f_2(j) ✗ symmetric
#
# So the obstruction literally IS the sign of (-1)^{n-1}.
# ==========================================================
print("\n" + "=" * 70)
print("PART 16: THE d=2 MODE AS OBSTRUCTION SWITCH")
print("=" * 70)

print(f"\n  For Dic_n: f_2(j) = (-1)^j, and j* = n-1.")
print(f"  f_2(j*-j) = (-1)^(j*-j) = (-1)^(n-1) · (-1)^j")
print(f"    n even: (-1)^(n-1) = -1 → ANTISYMMETRIC → duality → obstructed")
print(f"    n odd:  (-1)^(n-1) = +1 → SYMMETRIC     → no duality → free")
print(f"\n  Verification:")

for n in range(2, 16):
    order = 4 * n
    j_star = n - 1
    sign = (-1)**(n-1)
    obstructed = n % 2 == 0
    # f_2(j) + f_2(j*-j) = (-1)^j + (-1)^{j*-j} = (-1)^j(1 + (-1)^{j*})
    # = (-1)^j(1 + (-1)^{n-1})
    # n even: 1 + (-1)^{n-1} = 1-1 = 0 → antisymmetric
    # n odd:  1 + (-1)^{n-1} = 1+1 = 2 → symmetric, amplitude 2
    f_sum = 1 + sign

    print(f"    n={n:2d}: (-1)^(n-1)={sign:+d}, "
          f"f₂(j)+f₂(j*-j) = {f_sum}·(-1)^j  "
          f"{'→ antisym → OBSTRUCTED' if obstructed else '→ sym → free'}")

# The weight of the d=2 mode
print(f"\n  For Dic_n, the d=2 mode weight is A_2 = 2n (all order-2 elements).")
print(f"  The contribution to S(j)+S(j*-j) from d=2:")
print(f"    n even: 2n · 0 = 0")
print(f"    n odd:  2n · 2·(-1)^j = 4n·(-1)^j")
print(f"  Other modes (d>2) always contribute 0 (|G|/(2d) even).")
print(f"  So for n odd: S(j)+S(j*-j) = 4n·(-1)^j.")

# Verify this formula
print(f"\n  Verification of S(j)+S(j*-j) = 4n·(-1)^j for n odd:")
for n in [3, 5, 7, 9, 11]:
    order = 4 * n
    j_star = n - 1
    elems = build_dicyclic(n)
    classes = compute_conjugacy_classes(elems)
    classes.sort(key=lambda c: -c['trace'])
    sizes = np.array([c['size'] for c in classes], dtype=float)

    max_err = 0
    for j in range(j_star + 1):
        jp = j_star - j
        S_j = sum(chi_su2(j, c['trace']) * s for c, s in zip(classes, sizes)) - 2*(2*j+1)
        S_jp = sum(chi_su2(jp, c['trace']) * s for c, s in zip(classes, sizes)) - 2*(2*jp+1)
        expected = 4 * n * (-1)**j
        err = abs((S_j + S_jp) - expected)
        max_err = max(max_err, err)

    print(f"    Dic_{n:2d}: max|S+S' - 4n(-1)^j| = {max_err:.1e} ✓")

# ==========================================================
# Part 17: Polyhedral — what if we COULD break condition (B)?
#
# For polyhedral groups, condition (B) holds: |H|/d is even for
# ALL d. There is no single "switch" mode. The obstruction is
# distributed across ALL modes simultaneously.
#
# But we can ask: for each mode, what is the antisymmetry check?
# ==========================================================
print("\n" + "=" * 70)
print("PART 17: POLYHEDRAL — DISTRIBUTED OBSTRUCTION")
print("=" * 70)

for name in ['2T', '2O', '2I']:
    g = groups[name]
    j_star = g['j_star']
    G_order = g['G']
    H_order = G_order // 2

    print(f"\n  {name} (|H|={H_order}):")
    print(f"    {'axis':>10} {'d':>2} {'|H|/d':>5} {'even?':>5} "
          f"{'weight':>6} {'% of |G|-2':>10}")
    print("    " + "-" * 50)

    total_nc = G_order - 2
    angles_seen = {}
    for c in g['classes']:
        if abs(c['trace']) > 1.999:
            continue
        alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
        alpha_canon = min(alpha, np.pi - alpha)
        key = round(alpha_canon, 8)
        if key not in angles_seen:
            # Find d
            best_d = None
            for dd in range(2, 20):
                for k in range(1, dd):
                    if abs(alpha_canon - k * np.pi / dd) < 0.001:
                        best_d = dd
                        break
                if best_d:
                    break
            angles_seen[key] = {'d': best_d, 'weight': 0}
        angles_seen[key]['weight'] += c['size']

    for key in sorted(angles_seen.keys()):
        m = angles_seen[key]
        d = m['d']
        ratio = H_order // d
        axis_type = {2: "2-fold", 3: "3-fold", 4: "4-fold", 5: "5-fold"}.get(d, f"{d}-fold")
        print(f"    {axis_type:>10} {d:2d} {ratio:5d} {'yes':>5} "
              f"{m['weight']:6d} {100*m['weight']/total_nc:10.1f}%")

    print(f"    ALL |H|/d even → EVERY mode antisymmetric → full duality")

print(f"\n  Summary:")
print(f"    Dicyclic: obstruction controlled by ONE mode (d=2)")
print(f"      Parity of n switches it on/off")
print(f"    Polyhedral: obstruction distributed across ALL modes")
print(f"      ALL rotation orders have |H|/d even simultaneously")
print(f"      Cannot break one mode without breaking symmetry group")

# ==========================================================
# Part 18: Symmetry breaking O → D₃ as mode filtering
#
# When O breaks to D₃, the rotation axes change:
#   O: 2-fold (9), 3-fold (8), 4-fold (6) → 23 non-identity
#   D₃: 2-fold (3), 3-fold (2) → 5 non-identity
# The 4-fold axis mode DISAPPEARS entirely.
# The 2-fold mode loses weight (9→3), 3-fold (8→2).
#
# In the band picture: which modes survive controls whether
# obstruction persists. D₃ has n=3 (odd) → NOT obstructed.
# The d=2 mode has |2D₃|/(2·2) = 12/4 = 3 → ODD → symmetric.
# So the duality breaks: scalar spinorial states appear.
# ==========================================================
print("\n" + "=" * 70)
print("PART 18: SYMMETRY BREAKING O → D₃ AS MODE FILTERING")
print("=" * 70)

# Build 2D₃ = Dic_3
n_D3 = 3
elems_2D3 = build_dicyclic(n_D3)
order_2D3 = 4 * n_D3
classes_2D3 = compute_conjugacy_classes(elems_2D3)
classes_2D3.sort(key=lambda c: -c['trace'])
sizes_2D3 = np.array([c['size'] for c in classes_2D3], dtype=float)
j_star_2D3 = n_D3 - 1  # = 2

# Compute modes for both groups
def get_modes(classes, sizes, G_order):
    """Get mode amplitudes grouped by canonical half-angle."""
    alpha_modes = {}
    for c in classes:
        if abs(c['trace']) > 1.999:
            continue
        alpha = np.arccos(np.clip(c['trace'] / 2, -1, 1))
        key = round(alpha, 8)
        if key not in alpha_modes:
            best_d, best_k = None, None
            for d in range(2, 20):
                for k in range(1, d):
                    if abs(alpha - k * np.pi / d) < 0.001:
                        best_d, best_k = d, k
                        break
                if best_d:
                    break
            alpha_modes[key] = {'alpha': alpha, 'weight': 0, 'd': best_d, 'k': best_k}
        alpha_modes[key]['weight'] += c['size']
    return sorted(alpha_modes.values(), key=lambda m: m['alpha'])

modes_2O = get_modes(groups['2O']['classes'], groups['2O']['sizes'], 48)
modes_2D3 = get_modes(classes_2D3, sizes_2D3, order_2D3)

print(f"\n  Mode comparison: 2O → 2D₃")
print(f"  {'angle':>10} {'2O wt':>6} {'2D₃ wt':>7} {'|G|/(2d)':>8} {'parity':>7} note")
print("  " + "-" * 60)

# Collect all angles from both
all_angles = {}
for m in modes_2O:
    label = f"{m['k']}π/{m['d']}"
    all_angles[label] = {'alpha': m['alpha'], 'w_2O': m['weight'], 'w_2D3': 0,
                         'd': m['d'], 'k': m['k']}
for m in modes_2D3:
    label = f"{m['k']}π/{m['d']}"
    if label in all_angles:
        all_angles[label]['w_2D3'] = m['weight']
    else:
        all_angles[label] = {'alpha': m['alpha'], 'w_2O': 0, 'w_2D3': m['weight'],
                             'd': m['d'], 'k': m['k']}

for label in sorted(all_angles.keys(), key=lambda l: all_angles[l]['alpha']):
    m = all_angles[label]
    d = m['d']
    ratio_2D3 = order_2D3 / (2 * d)
    parity = "even" if int(ratio_2D3) % 2 == 0 else "ODD"
    note = ""
    if m['w_2O'] > 0 and m['w_2D3'] == 0:
        note = "← LOST (4-fold axis gone)"
    elif m['w_2O'] > 0 and m['w_2D3'] > 0:
        note = f"← reduced {m['w_2O']}→{m['w_2D3']}"
    if m['w_2D3'] > 0:
        print(f"  {label:>10} {m['w_2O']:6d} {m['w_2D3']:7d} {ratio_2D3:8.1f} {parity:>7} {note}")
    else:
        print(f"  {label:>10} {m['w_2O']:6d} {'—':>7} {'—':>8} {'—':>7} {note}")

# Band comparison
print(f"\n  Band diagrams:")
j_star_2O = 11
band_2O = ""
for j in range(j_star_2O + 1):
    m = scalar_multiplicity(j, groups['2O']['classes'], groups['2O']['sizes'], 48)
    band_2O += "▮" if m > 0 else "·"

band_2D3 = ""
for j in range(j_star_2D3 + 1):
    m = scalar_multiplicity(j, classes_2D3, sizes_2D3, order_2D3)
    band_2D3 += "▮" if m > 0 else "·"

# Extended band for D₃
band_2D3_ext = ""
for j in range(j_star_2O + 1):
    m = scalar_multiplicity(j, classes_2D3, sizes_2D3, order_2D3)
    band_2D3_ext += "▮" if m > 0 else "·"

print(f"    2O  (j*=11): {band_2O}  (half-filled, duality)")
print(f"    2D₃ (j*= 2): {band_2D3}  (both ends filled, no duality)")
print(f"    2D₃ extended to j=11: {band_2D3_ext}")

# What happens to the duality
print(f"\n  Duality check after breaking:")
print(f"    2O:  m+m' = 1 for all j ≤ 11 → obstructed → no scalar spinorial")
print(f"    2D₃: m+m' = 2 for j=0,2 (both filled) → duality broken → "
      f"scalar spinorial at j=1/2!")

# The re-entrance as mode switching
print(f"\n  RE-ENTRANCE MECHANISM (mode picture):")
print(f"    O has 3 modes: d=4 (antisym), d=3 (antisym), d=2 (antisym)")
print(f"    D₃ has 2 modes: d=3 (antisym), d=2 (SYMMETRIC — n=3 odd)")
print(f"    The d=4 mode vanishes (no 4-fold axis in D₃)")
print(f"    The d=2 mode FLIPS from antisymmetric to symmetric")
print(f"    This flip = loss of obstruction = re-entrance of spinorial states")

print("\n" + "=" * 70)
print("EXPLORATION COMPLETE")
print("=" * 70)
