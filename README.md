# st_spinorial_obstruction

Spectral reflection symmetry and periodic multiplicity structure for quantum rotors on SO(3)/H.

## Problem

A quantum rigid rotor with orientational symmetry H ⊂ SO(3) lives on the configuration space M = SO(3)/H. Because M is multiply connected (π₁ = 2H ⊂ SU(2)), quantization splits into superselection sectors labelled by one-dimensional characters χ: 2H → U(1). Sectors with χ(−1) = +1 are tensorial, those with χ(−1) = −1 are spinorial. For polyhedral groups (T, O, I), −1 lies in the commutator subgroup of 2H, forcing all sectors tensorial — no spinorial sector exists.

## Main result

Triple equivalence (A) ⟺ (B) ⟺ (C) for all finite H ⊂ SO(3) with |H| even:

- **(A)** Commutator obstruction: −1 ∈ [2H, 2H]
- **(B)** Divisibility condition: |H|/ord(h) is even for every h ∈ H \ {e}
- **(C)** Spectral reflection: m(χ, j) + m(χ, j\*−j) = 1 for every scalar character χ and every j ≤ j\* = |H|/2 − 1 of matching parity

Under (B), the multiplicity decomposes as m(χ, j) = ⌊j/P⌋ + ε(χ, j) where P = |H|/2, ε is periodic with period P, binary ({0,1}), and reflecting (ε + ε' = 1). Binary occupation, regular-representation tiling, staircase periodicity, and threshold sharpness follow as corollaries.

Under symmetry breaking H → K with Q₈ ⊄ 2K, spinorial sectors reappear; for odd dihedral groups, the reflection fails globally but persists in spinorial sectors (sector-selective duality).

## Results

Classification of sub-threshold spectral structure for all finite rotation groups:

| Group H | \|H\| | j\* | (B) holds | Reflection | Spinorial sectors |
|---------|-------|-----|-----------|------------|-------------------|
| D_n (n even) | 2n | n−1 | yes | all χ | blocked |
| T | 12 | 5 | yes | all χ | blocked |
| O | 24 | 11 | yes | all χ | blocked |
| I | 60 | 29 | yes | all χ | blocked |
| D_n (n odd) | 2n | n−1 | no | spinorial only | accessible (j = 1/2) |
| C_n | n | — | |H| odd | — | accessible (j = 1/2) |

Verified computationally on all binary polyhedral groups (2T, 2O, 2I), all dicyclic groups Dic_n (n = 2..101), and all dihedral groups D_n (n = 3..15).

## Requirements

- Python 3.9+
- NumPy, SciPy
- pytest

## Running tests

```bash
# All tests (70 files, 1740 tests)
.venv/bin/python -m pytest tests/ -v

# Single file
.venv/bin/python -m pytest tests/51_table_verification.py -v

# Quick smoke test (core obstruction + table verification)
.venv/bin/python -m pytest tests/1_test_spinor_obstruction.py tests/51_table_verification.py -v
```

## Test suite → Paper claims

Each test file maps to specific claims in the paper (`paper/main_v10.tex`). The test suite contains 70 files with 1740 pytest tests and ~3200 assertions.

### Core theorem verification

| Test file | Paper ref | What it verifies |
|-----------|-----------|------------------|
| `51_table_verification.py` (197 tests) | Tables 2, 5, 6 | Every m(χ,j) entry in paper tables; dual column cross-check |
| `52_converse_formula.py` (112 tests) | Thm 3.3 | Converse: S_odd from group elements; (B) direct verification |
| `53_above_threshold.py` (266 tests) | Cor 3.7, 3.8 | m(j\*+1) = 2; periodicity m(j+P) = m(j)+1; extended duality |
| `35_stress_test_duality.py` (833 tests) | Thm 3.3, Cor 3.4 | Reflection identity on Dic_n, n = 2..101; all scalar characters |
| `34_analytic_equivalence.py` | Thm 3.3 | (B)⟺(C) classification-free via sum-to-product identity |

### Topological obstruction (A) ⟺ (B)

| Test file | Paper ref | What it verifies |
|-----------|-----------|------------------|
| `58_vector3_sylow.py` (47 tests) | Prop 4.1 | Sylow chain: (B) ⟹ Syl₂ non-cyclic ⟹ V₄ ⊂ H ⟹ Q₈ ⊂ 2H ⟹ (A) |
| `59_gen_quaternion_no_so3.py` (77 tests) | Prop 4.1 | Q_{2^n} ⊄ SO(3); V₄ ⊂ H ⟺ (B) on 16 groups; Burnside orbit-counting |
| `50_structural_chain.py` (12 tests) | Prop 4.1 | V₄ ⟺ Q₈ ⟺ single commutator ⟺ s > 0; ADE-free |

### Symmetry breaking and selective duality

| Test file | Paper ref | What it verifies |
|-----------|-----------|------------------|
| `54_reentrance_D_odd.py` (81 tests) | §4.4 | Spinorial duality holds, tensorial fails for D₃..D₁₁ |
| `66_selective_duality_general.py` (39 tests) | §4.4 | Duality ⟺ immunity to (B)-violating classes |
| `30_symmetry_breaking_phase_diagram.py` | §4.2, Table 6 | O → H' phase diagram: 4/7 channels restore spinorial |

### Supporting verifications

| Test file | Paper ref | What it verifies |
|-----------|-----------|------------------|
| `56_bridge_B2_B3.py` (94 tests) | §3.2, §3.3 | exp(H) = \|H\|/2 (integrality argument); d=2 uniqueness |
| `49_2Dn_analytic.py` (12 tests) | §3.1, §3.3 | Analytic proof for 2D_n; mode antisymmetry |
| `65_molien_gorenstein.py` (14 tests) | §5 | Molien series; Gorenstein palindrome = A₀ reflection |
| `37_frobenius_reformulation.py` (23 tests) | Cor 3.6 | V_j ⊕ V_{j\*−j} ≅ ρ_reg; Peter-Weyl uniformity |

### Negative results and scope

| Test file | What it verifies |
|-----------|------------------|
| `55_su3_counterexample.py` (318 tests) | SU(3)/Δ(27): no complementary pairing (rank-1 specific) |
| `73_rigidity_and_spin4.py` (20 tests) | Spin(4) = SU(2)×SU(2): reflection fails; phenomenon SU(2)-specific |
| `68_spectral_zeta_heat.py` (26 tests) | Heat kernel: no Seeley-DeWitt coefficient distinguishes (B) |
| `72_twisted_hilbert.py` (30 tests) | Twisted Hilbert functions: Gorenstein universal, not (B)-specific |

## Structure

```
src/
├── __init__.py
├── quaternion.py    # Hamilton product, conjugate, inverse, SU(2)/SO(3) conversions
└── group.py         # Build 2T (24), 2O (48), 2I (120); closure, conjugacy classes

tests/
├── 1_test_spinor_obstruction.py    # Core obstruction proof (44 tests)
├── 2_test_general_obstruction.py   # All binary polyhedral + dicyclic (97 tests)
├── ...                             # 68 more test files (see tests/tests_map.md)
├── 74_blue_phase.py                # Blue phase application (22 tests)
└── tests_map.md                    # Complete test inventory with paper mappings

paper/
└── main_v10.tex     # LaTeX manuscript (J. Phys. A target)
```

## Key source files

`src/quaternion.py` implements quaternion arithmetic for unit quaternions in SU(2): Hamilton product (`qmul`), inverse (`qinv`), and the maps to SU(2) matrices (`quat_to_su2`) and SO(3) rotation matrices (`quat_to_so3`).

`src/group.py` constructs the binary polyhedral groups 2T (24 elements), 2O (48 elements), and 2I (120 elements) as unit quaternion arrays, and provides group algorithms: closure verification, commutator subgroup computation, and conjugacy class identification.

## License

MIT
