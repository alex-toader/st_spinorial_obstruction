#!/usr/bin/env python3
"""
Branching rules: restriction of 2O irreps to 2D₄ ⊂ 2O.

2D₄ embedded in 2O as stabilizer of [100] (tetragonal axis):
  a = (1+k)/√2  (order 8, 90° about z-axis in SO(3))
  x = i          (order 4, 180° about x-axis)
Relations: a⁸ = 1, x² = a⁴ = -1, xax⁻¹ = a⁻¹.

CONTRAST with trigonal (2O → 2D₃):
  Trigonal: L(S) → σ₁(S) ⊕ σ₃(S) ⊕ φ₁(S) — two 1D spinorial sectors restored
  Tetragonal: ALL spinorial 2O irreps → rank-2 irreps of 2D₄ only

Mechanism: Q₈ ⊂ 2D₄, so 2D₄ has NO 1D spinorial characters.
By central parity, spinorial 2O irreps (ρ(-1)=-I) can only restrict to
spinorial 2D₄ irreps, which are all 2-dimensional.

RAW OUTPUT:
===========

2D₄ ⊂ 2O verified: |2D₄| = 16, all elements in 2O
Relations verified: a⁴ = x² = -1, a⁸ = 1, xax⁻¹ = a⁻¹
Q₈ ⊂ 2D₄: True

2D₄: 7 conjugacy classes, 7 irreps (4×dim1 + 3×dim2)
  All 4 one-dimensional chars: tensorial (χ(-1) = +1)
  ψ₁, ψ₃: spinorial (dim 2)
  ψ₂: tensorial (dim 2)

BRANCHING RULES: 2O → 2D₄ (cubic → tetragonal)
  A₀(T) dim 1:  χ₀(T)
  A₁(T) dim 1:  χ₂(T)
   E(T) dim 2:  χ₀(T) ⊕ χ₂(T)
  T₁(T) dim 3:  χ₃(T) ⊕ ψ₂(T)
  T₂(T) dim 3:  χ₁(T) ⊕ ψ₂(T)
  H₁(S) dim 2:  ψ₁(S)
  H₂(S) dim 2:  ψ₃(S)
   L(S) dim 4:  ψ₁(S) ⊕ ψ₃(S)

CONTRAST:
  Trigonal L → σ₁(S) ⊕ σ₃(S) ⊕ ψ₁(S) → 2 NEW scalar spinorial sectors
  Tetragonal L → ψ₁(S) ⊕ ψ₃(S) → rank-2 ONLY → NO new scalar sectors

5 assertions passed.
"""
import sys
sys.path.insert(0, '/Users/alextoader/Sites/st_spinorial_obstruction')

import numpy as np
from src.quaternion import qmul, qconj, qkey
from src.group import build_binary_octahedral

SQ2 = np.sqrt(2)

# Note: qconj computes quaternion conjugate (negate vector part).
# For unit quaternions |q|=1, the conjugate equals the inverse: q⁻¹ = q*.
# All quaternions in this script are unit quaternions.


# =============================================================
# 1. Build 2D₄ inside 2O
# =============================================================

a = np.array([1/SQ2, 0, 0, 1/SQ2])    # (1+k)/√2, order 8
x = np.array([0, 1, 0, 0])              # i, order 4

def generate_group(generators, max_order=200):
    elements = [np.array([1., 0., 0., 0.])]
    seen = {qkey(elements[0])}
    queue = list(elements)
    idx = 0
    while idx < len(queue):
        g = queue[idx]; idx += 1
        for gen in generators:
            for h in [qmul(g, gen), qmul(gen, g)]:
                k = qkey(h)
                if k not in seen:
                    seen.add(k)
                    queue.append(h)
                    elements.append(h)
        if len(elements) > max_order:
            raise RuntimeError(f"Group too large: {len(elements)}")
    return elements

elements_2D4 = generate_group([a, x, qconj(a)])
assert len(elements_2D4) == 16, f"|2D₄| = {len(elements_2D4)}, expected 16"

# Verify embedding in 2O
elements_2O = build_binary_octahedral()
keys_2O = {qkey(g) for g in elements_2O}
for g in elements_2D4:
    assert qkey(g) in keys_2O, f"{g} not in 2O"
print(f"2D₄ ⊂ 2O verified: |2D₄| = 16, all elements in 2O")

# Verify relations
a4 = qmul(qmul(a, a), qmul(a, a))
assert np.allclose(a4, [-1, 0, 0, 0]), f"a⁴ ≠ -1, got {a4}"
a8 = qmul(a4, a4)
assert np.allclose(a8, [1, 0, 0, 0]), f"a⁸ ≠ 1"
assert np.allclose(qmul(x, x), [-1, 0, 0, 0]), "x² ≠ -1"
xax_inv = qmul(qmul(x, a), qconj(x))  # x has x⁻¹ = -x for order 4
# Actually x⁻¹ = conj(x)/|x|² = -x since |x|=1 and x is pure imaginary
assert np.allclose(xax_inv, qconj(a)), f"xax⁻¹ ≠ a⁻¹"
print("Relations verified: a⁴ = x² = -1, a⁸ = 1, xax⁻¹ = a⁻¹")

# Verify Q₈ ⊂ 2D₄
q8_elements = [
    np.array([1, 0, 0, 0]),
    np.array([-1, 0, 0, 0]),
    np.array([0, 1, 0, 0]),   # i = x
    np.array([0, -1, 0, 0]),
    np.array([0, 0, 1, 0]),   # j
    np.array([0, 0, -1, 0]),
    np.array([0, 0, 0, 1]),   # k = a²
    np.array([0, 0, 0, -1]),
]
keys_2D4 = {qkey(g) for g in elements_2D4}
q8_in = all(qkey(q) in keys_2D4 for q in q8_elements)
print(f"Q₈ ⊂ 2D₄: {q8_in}")
assert q8_in, "Q₈ not contained in 2D₄!"

# =============================================================
# 2. Conjugacy classes of 2D₄
# =============================================================

def conjugacy_classes(elements):
    remaining = list(range(len(elements)))
    classes = []
    while remaining:
        rep_idx = remaining[0]
        rep = elements[rep_idx]
        cls, new_remaining = [], []
        for idx in remaining:
            g = elements[idx]
            is_conj = any(np.allclose(qmul(qmul(h, g), qconj(h)), rep, atol=1e-10)
                         for h in elements)
            (cls if is_conj else new_remaining).append(idx)
        classes.append(cls)
        remaining = new_remaining
    return classes

classes_2D4 = conjugacy_classes(elements_2D4)
assert len(classes_2D4) == 7, f"Expected 7 classes, got {len(classes_2D4)}"
class_sizes = [len(c) for c in classes_2D4]

print(f"\n2D₄: 7 conjugacy classes, sizes = {class_sizes}")
for i, cls in enumerate(classes_2D4):
    rep = elements_2D4[cls[0]]
    print(f"  C{i}: size {len(cls)}, rep = {np.round(rep, 4)}, tr = {2*rep[0]:.4f}")

# =============================================================
# 3. Character table of 2D₄ (Dic₄)
# =============================================================
#
# Dic₄: order 16, abelianization Z₂ × Z₂
# Relations: a⁸=1, x²=a⁴=-1, xax⁻¹=a⁻¹
#
# 4 one-dimensional characters (N=4 even):
#   χ(a) ∈ {+1,-1}, χ(x) ∈ {+1,-1}
#   χ(-1) = χ(a⁴) = χ(a)⁴ = 1 for ALL → ALL tensorial
#
# 3 two-dimensional irreps:
#   ψ_m(a) = diag(ω^m, ω^{-m}) where ω = e^{2πi/8}
#   ψ_m(x) = [[0,1],[1,0]] (or variant)
#   m = 1, 2, 3
#   ψ_m(-1) = ψ_m(a⁴) = diag(ω^{4m}, ω^{-4m}) = diag((-1)^m, (-1)^m)
#   m=1: spinorial (-I), m=2: tensorial (+I), m=3: spinorial (-I)

# Sort classes by trace (descending) for standard ordering
class_order = sorted(range(7), key=lambda i: -2*elements_2D4[classes_2D4[i][0]][0])
classes_2D4 = [classes_2D4[i] for i in class_order]
class_sizes = [len(c) for c in classes_2D4]

print(f"\n2D₄ classes (sorted by trace):")
class_traces = []
for i, cls in enumerate(classes_2D4):
    rep = elements_2D4[cls[0]]
    tr = 2 * rep[0]
    class_traces.append(tr)
    print(f"  C{i}: size {len(cls)}, tr = {tr:.4f}")

# Build character values on classes
# We need to identify which element represents each class
# The classes of Dic₄ are:
#   {1}: tr=2, size 1
#   {a, a⁷}: tr=√2, size 2
#   {a², a⁶}: tr=0, size 2  (these are ±k)
#   {a³, a⁵}: tr=-√2, size 2
#   {a⁴=-1}: tr=-2, size 1
#   {x, a²x, a⁴x, a⁶x}: tr=0, size 4  (type 1)
#   {ax, a³x, a⁵x, a⁷x}: tr=0, size 4  (type 2)

# Identify x-type classes (tr=0, size 4) vs a²-type (tr=0, size 2)
x_key = qkey(x)
a2 = qmul(a, a)  # = k
a2_key = qkey(a2)

# Map our sorted classes to standard Dic₄ classes
print("\nIdentifying classes:")
std_map = {}
for i, cls in enumerate(classes_2D4):
    size = len(cls)
    tr = class_traces[i]
    keys_in_cls = {qkey(elements_2D4[idx]) for idx in cls}

    if size == 1 and abs(tr - 2) < 0.01:
        std_map[i] = 'e'
        print(f"  C{i} → {{1}}")
    elif size == 2 and abs(tr - SQ2) < 0.01:
        std_map[i] = 'a'
        print(f"  C{i} → {{a, a⁷}}")
    elif size == 2 and abs(tr) < 0.01:
        std_map[i] = 'a2'
        print(f"  C{i} → {{a², a⁶}} = {{k, -k}}")
    elif size == 2 and abs(tr + SQ2) < 0.01:
        std_map[i] = 'a3'
        print(f"  C{i} → {{a³, a⁵}}")
    elif size == 1 and abs(tr + 2) < 0.01:
        std_map[i] = 'minus1'
        print(f"  C{i} → {{-1}}")
    elif size == 4 and x_key in keys_in_cls:
        std_map[i] = 'x'
        print(f"  C{i} → {{x, a²x, a⁴x, a⁶x}}")
    elif size == 4:
        std_map[i] = 'ax'
        print(f"  C{i} → {{ax, a³x, a⁵x, a⁷x}}")
    else:
        raise ValueError(f"Unidentified class C{i}: size={size}, tr={tr}")

# Character table values in standard order: e, a, a², a³, -1, x, ax
omega = np.exp(2j * np.pi / 8)

char_table_std = {
    # 1D characters: (N=4 even, all tensorial)
    'χ₀': {'e': 1, 'a': 1, 'a2': 1, 'a3': 1, 'minus1': 1, 'x': 1, 'ax': 1},
    'χ₁': {'e': 1, 'a': 1, 'a2': 1, 'a3': 1, 'minus1': 1, 'x': -1, 'ax': -1},
    'χ₂': {'e': 1, 'a': -1, 'a2': 1, 'a3': -1, 'minus1': 1, 'x': 1, 'ax': -1},
    'χ₃': {'e': 1, 'a': -1, 'a2': 1, 'a3': -1, 'minus1': 1, 'x': -1, 'ax': 1},
    # 2D characters: tr(ψ_m) on each class
    # ψ₁: a → ω¹+ω⁻¹=√2, a² → ω²+ω⁻²=0, a³ → ω³+ω⁻³=-√2, -1 → ω⁴+ω⁻⁴=-2
    'ψ₁': {'e': 2, 'a': SQ2, 'a2': 0, 'a3': -SQ2, 'minus1': -2, 'x': 0, 'ax': 0},
    # ψ₂: a → ω²+ω⁻²=0, a² → ω⁴+ω⁻⁴=-2, a³ → ω⁶+ω⁻⁶=0, -1 → ω⁸+ω⁻⁸=2
    'ψ₂': {'e': 2, 'a': 0, 'a2': -2, 'a3': 0, 'minus1': 2, 'x': 0, 'ax': 0},
    # ψ₃: a → ω³+ω⁻³=-√2, a² → ω⁶+ω⁻⁶=0, a³ → ω⁹+ω⁻⁹=√2, -1 → ω¹²+ω⁻¹²=-2
    'ψ₃': {'e': 2, 'a': -SQ2, 'a2': 0, 'a3': SQ2, 'minus1': -2, 'x': 0, 'ax': 0},
}

# Convert to array indexed by our class ordering
CHAR_2D4 = {}
for name, vals in char_table_std.items():
    CHAR_2D4[name] = np.array([vals[std_map[i]] for i in range(7)], dtype=complex)

# Print character table
print("\n2D₄ character table:")
header = "".join(f"  C{i}({class_sizes[i]})" for i in range(7))
print(f"  {'':>8}{header}")
for name, chi in CHAR_2D4.items():
    vals = "".join(f"  {v.real:6.2f}" for v in chi)
    # Spinorial status
    minus1_idx = [i for i in range(7) if std_map[i] == 'minus1'][0]
    v_m1 = chi[minus1_idx].real
    if name.startswith('χ'):
        dim = 1
        spin = 'S' if v_m1 < -0.5 else 'T'
    else:
        dim = 2
        spin = 'S' if v_m1 < -0.5 else 'T'
    print(f"  {name:>8}{vals}  dim={dim} {spin}")

# Verify orthogonality
print("\nOrthogonality check:")
names = list(CHAR_2D4.keys())
sizes_arr = np.array(class_sizes, dtype=float)
all_ok = True
for i, n1 in enumerate(names):
    for j, n2 in enumerate(names):
        if j >= i:
            ip = np.sum(np.conj(CHAR_2D4[n1]) * CHAR_2D4[n2] * sizes_arr) / 16
            expected = 1.0 if i == j else 0.0
            if abs(ip - expected) > 0.01:
                print(f"  FAIL: ⟨{n1},{n2}⟩ = {ip}")
                all_ok = False
if all_ok:
    print("  All ⟨χᵢ,χⱼ⟩ = δᵢⱼ ✓")

# Verify sum of dim² = |G|
dim_sq = sum(int(CHAR_2D4[n][0].real)**2 for n in names)
assert dim_sq == 16, f"Σ dim² = {dim_sq}, expected 16"
print(f"  Σ dim² = {dim_sq} = |2D₄| ✓")

# Spinorial status summary
print("\nSpinorial status:")
minus1_idx = [i for i in range(7) if std_map[i] == 'minus1'][0]
for name, chi in CHAR_2D4.items():
    v = chi[minus1_idx].real
    dim = int(chi[0].real)
    if v < -0.5:
        print(f"  {name}(-1) = {v:+.0f}  SPINORIAL (dim {dim})")
    else:
        print(f"  {name}(-1) = {v:+.0f}  tensorial (dim {dim})")

n_1d_spin = sum(1 for n in ['χ₀','χ₁','χ₂','χ₃']
                if CHAR_2D4[n][minus1_idx].real < -0.5)
print(f"\n  1D spinorial characters: {n_1d_spin}")
assert n_1d_spin == 0, "Expected NO 1D spinorial characters for even dihedral!"
print(f"  → NO scalar spinorial sectors available (Q₈ ⊂ 2D₄)")

# =============================================================
# 4. 2O character table
# =============================================================

classes_2O = conjugacy_classes(elements_2O)
class_info_2O = []
for cls in classes_2O:
    rep = elements_2O[cls[0]]
    class_info_2O.append({
        'size': len(cls), 'trace': round(2 * rep[0], 6),
        'keys': {qkey(elements_2O[i]) for i in cls}
    })

CHAR_2O = {
    'A₀': {(1,2.0):1,(12,0.0):1,(8,-1.0):1,(6,0.0):1,(1,-2.0):1,(6,-SQ2):1,(8,1.0):1,(6,SQ2):1},
    'A₁': {(1,2.0):1,(12,0.0):-1,(8,-1.0):1,(6,0.0):1,(1,-2.0):1,(6,-SQ2):-1,(8,1.0):1,(6,SQ2):-1},
    'E':  {(1,2.0):2,(12,0.0):0,(8,-1.0):-1,(6,0.0):2,(1,-2.0):2,(6,-SQ2):0,(8,1.0):-1,(6,SQ2):0},
    'T₁': {(1,2.0):3,(12,0.0):1,(8,-1.0):0,(6,0.0):-1,(1,-2.0):3,(6,-SQ2):-1,(8,1.0):0,(6,SQ2):-1},
    'T₂': {(1,2.0):3,(12,0.0):-1,(8,-1.0):0,(6,0.0):-1,(1,-2.0):3,(6,-SQ2):1,(8,1.0):0,(6,SQ2):1},
    'H₁': {(1,2.0):2,(12,0.0):0,(8,-1.0):-1,(6,0.0):0,(1,-2.0):-2,(6,-SQ2):-SQ2,(8,1.0):1,(6,SQ2):SQ2},
    'H₂': {(1,2.0):2,(12,0.0):0,(8,-1.0):-1,(6,0.0):0,(1,-2.0):-2,(6,-SQ2):SQ2,(8,1.0):1,(6,SQ2):-SQ2},
    'L':  {(1,2.0):4,(12,0.0):0,(8,-1.0):1,(6,0.0):0,(1,-2.0):-4,(6,-SQ2):0,(8,1.0):-1,(6,SQ2):0},
}

# Verify (size, trace) keys are unique across 2O classes
st_keys = [(ci['size'], round(ci['trace'], 4)) for ci in class_info_2O]
assert len(set(st_keys)) == len(class_info_2O), \
    f"Non-unique (size,trace) keys in 2O: {st_keys}"

def char_2O_val(irrep, q):
    k = qkey(q)
    for ci in class_info_2O:
        if k in ci['keys']:
            key = (ci['size'], round(ci['trace'], 4))
            for (s, t), v in CHAR_2O[irrep].items():
                if s == key[0] and abs(t - key[1]) < 0.01:
                    return v
    raise ValueError(f"Element {q} not found in 2O classes")

# =============================================================
# 5. Branching rules: 2O → 2D₄ (tetragonal)
# =============================================================

print("\n" + "=" * 70)
print("BRANCHING RULES: 2O → 2D₄ (cubic → tetragonal)")
print("=" * 70)
print()

irrep_names_2O = ['A₀', 'A₁', 'E', 'T₁', 'T₂', 'H₁', 'H₂', 'L']
irrep_names_2D4 = list(CHAR_2D4.keys())

labels_2D4 = {
    'χ₀': 'χ₀(T)', 'χ₁': 'χ₁(T)', 'χ₂': 'χ₂(T)', 'χ₃': 'χ₃(T)',
    'ψ₁': 'ψ₁(S)', 'ψ₂': 'ψ₂(T)', 'ψ₃': 'ψ₃(S)'
}

for name_2O in irrep_names_2O:
    # Restricted character on 2D₄ classes
    chi_res = np.array([char_2O_val(name_2O, elements_2D4[classes_2D4[i][0]])
                        for i in range(7)], dtype=complex)
    dim = int(chi_res[0].real)
    rho_m1_idx = [i for i in range(7) if std_map[i] == 'minus1'][0]
    rho_m1 = chi_res[rho_m1_idx].real
    spin = 'S' if rho_m1 < -0.5 else 'T'

    # Decompose
    decomp = []
    dim_check = 0
    for irr in irrep_names_2D4:
        ip = np.sum(np.conj(CHAR_2D4[irr]) * chi_res * sizes_arr) / 16
        m = int(round(ip.real))
        assert abs(ip.imag) < 0.1 and abs(ip.real - m) < 0.1, \
            f"{name_2O}, {irr}: multiplicity not integer: {ip}"
        if m > 0:
            d = int(CHAR_2D4[irr][0].real)
            dim_check += m * d
            label = labels_2D4[irr]
            decomp.append(f"{'%d' % m if m > 1 else ''}{label}")

    assert dim_check == dim, f"{name_2O}: dim mismatch {dim_check} vs {dim}"
    decomp_str = ' ⊕ '.join(decomp)
    print(f"  {name_2O:>3}({spin}) dim {dim}:  {decomp_str}")

# =============================================================
# 6. Contrast: trigonal vs tetragonal
# =============================================================

print()
print("=" * 70)
print("CONTRAST: TRIGONAL vs TETRAGONAL")
print("=" * 70)

print("""
Spinorial 2O irreps and their branching:

                    Trigonal (2O → 2D₃)        Tetragonal (2O → 2D₄)
  ─────────────────────────────────────────────────────────────────────""")

# Trigonal branching results from script 3 (3_branching_2O_to_2D3.py).
# Independently verified: 2D₃ = Dic₃ with generators a=(1+i+j+k)/2, b=(i-j)/√2.
# L(S) → σ₁(S) ⊕ σ₃(S) ⊕ ψ₁(S) is the key: two 1D spinorial sectors emerge.
trigonal = {
    'H₁': 'ψ₁(S)',
    'H₂': 'ψ₁(S)',
    'L':  'σ₁(S) ⊕ σ₃(S) ⊕ ψ₁(S)',
}

for name_2O in ['H₁', 'H₂', 'L']:
    # Tetragonal decomposition
    chi_res = np.array([char_2O_val(name_2O, elements_2D4[classes_2D4[i][0]])
                        for i in range(7)], dtype=complex)
    decomp = []
    for irr in irrep_names_2D4:
        ip = np.sum(np.conj(CHAR_2D4[irr]) * chi_res * sizes_arr) / 16
        m = int(round(ip.real))
        if m > 0:
            decomp.append(f"{'%d' % m if m > 1 else ''}{labels_2D4[irr]}")
    tet_str = ' ⊕ '.join(decomp)

    print(f"  {name_2O}(S) dim {int(chi_res[0].real)}:  {trigonal[name_2O]:30s}  {tet_str}")

print("""
Key result:
  Trigonal:   L → σ₁(S) ⊕ σ₃(S) ⊕ ψ₁(S)  → 2 NEW scalar spinorial sectors
  Tetragonal: ALL spinorial irreps → rank-2 irreps ONLY → NO new scalar sectors

  The direction of symmetry breaking determines re-entrance.
  Trigonal (O→D₃): Q₈ ⊄ 2D₃ → obstruction lifted
  Tetragonal (O→D₄): Q₈ ⊂ 2D₄ → obstruction persists
""")

# =============================================================
# 7. Assertions
# =============================================================

print("=" * 70)
print("VERIFICATION ASSERTIONS")
print("=" * 70)

# Check: no 1D spinorial in 2D₄
for name in ['χ₀', 'χ₁', 'χ₂', 'χ₃']:
    v = CHAR_2D4[name][minus1_idx].real
    assert v > 0.5, f"{name}(-1) = {v}, expected +1"
print("  [1] All 1D chars of 2D₄ are tensorial ✓")

# Check: all spinorial 2O irreps → only spinorial 2D₄ irreps
for name_2O in ['H₁', 'H₂', 'L']:
    chi_res = np.array([char_2O_val(name_2O, elements_2D4[classes_2D4[i][0]])
                        for i in range(7)], dtype=complex)
    for irr in ['χ₀', 'χ₁', 'χ₂', 'χ₃']:  # 1D chars
        ip = np.sum(np.conj(CHAR_2D4[irr]) * chi_res * sizes_arr) / 16
        m = int(round(ip.real))
        assert m == 0, f"{name_2O}|_2D₄ contains 1D char {irr} with mult {m}!"
print("  [2] No spinorial 2O irrep branches to any 1D char of 2D₄ ✓")

# Check: tensorial 2O irreps → only tensorial 2D₄ irreps
for name_2O in ['A₀', 'A₁', 'E', 'T₁', 'T₂']:
    chi_res = np.array([char_2O_val(name_2O, elements_2D4[classes_2D4[i][0]])
                        for i in range(7)], dtype=complex)
    for irr in ['ψ₁', 'ψ₃']:  # spinorial 2D irreps
        ip = np.sum(np.conj(CHAR_2D4[irr]) * chi_res * sizes_arr) / 16
        m = int(round(ip.real))
        assert m == 0, f"{name_2O}|_2D₄ contains spinorial {irr} with mult {m}!"
print("  [3] No tensorial 2O irrep branches to spinorial 2D₄ irrep ✓")

# Dimension checks
print(f"  [4] |2D₄| = 16, Σ dim² = {dim_sq} ✓")

# Q₈ containment
print(f"  [5] Q₈ ⊂ 2D₄: True ✓")

print(f"\nAll assertions passed.")


if __name__ == '__main__':
    pass  # All computation runs at module level (research script)
