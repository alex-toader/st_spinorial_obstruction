# Test Map — st_spinorial_obstruction

**Date:** Feb 2026

---

## Overview

55 files in `tests/`. Pytest suites (1404 test functions) plus verification/demonstration scripts with `assert` or print-based checks. All run as `/usr/bin/python3 -m pytest tests/XX_name.py -v`.

**Topic:** Spinorial obstruction on SO(3)/H quotient spaces. Central theorem: -1 ∈ [G,G] ⟺ Q₈ ⊂ G ⟺ no scalar (rank-1) spinorial sector exists.

---

## Source Code

| File | Description |
|------|-------------|
| `src/quaternion.py` | Hamilton product, conjugate, inverse, SU(2)/SO(3) conversion |
| `src/group.py` | Build 2T (24), 2O (48), 2I (120); closure, commutator subgroup, conjugacy classes |

---

## Test Files

### A. Core Obstruction (files 1-2) — 141 pytest tests

| # | File | Tests | Asserts | What it does | Key result |
|---|------|-------|---------|--------------|------------|
| 1 | 1_test_spinor_obstruction.py | 44 | 42 | 2O group structure, conjugacy classes, commutator subgroup, obstruction proof | -1 ∈ [2O,2O]=2T; abelianization Z₂; κ=2 sectors |
| 2 | 2_test_general_obstruction.py | 97 | 129 | All binary polyhedral (2T,2O,2I) + dicyclic + cyclic; comparative analysis | Monotonic: κ(2T)=3, κ(2O)=2, κ(2I)=1; Q₈ ⊂ [G,G] for all obstructed |

### B. Branching Rules (files 3-5) — scripts with asserts

| # | File | Asserts | What it does | Key result |
|---|------|---------|--------------|------------|
| 3 | 3_branching_2O_to_2D3.py | 8 | 2O→2D₃ restriction (cubic→trigonal) | Only L (rank-4 spinorial) produces 1D spinorial; H₁,H₂ stay rank-2 |
| 4 | 4_branching_2T_to_Z6.py | 9 | 2T→Z₆ restriction (tetrahedral→trigonal) | All 3 spinorial irreps produce 1D sectors; abelian Z₆ allows full restoration |
| 5 | 5_branching_2I_to_2D5.py | 12 | 2I→2D₅ restriction (icosahedral→pentagonal) | Only D (rank-6, maximal) produces 1D spinorial; pattern: highest-rank only |

### C. Structural Proofs (files 6-8) — scripts, no asserts

| # | File | What it does | Key result |
|---|------|--------------|------------|
| 6 | 6_restoration_threshold.py | Threshold: which spinorial irreps give 1D sectors upon restriction to Dic_N | dim ≥ N+1 required; b-type vanishing lemma + weight decomposition |
| 7 | 7_unification_Q8_branching.py | Unified criterion: Q₈ ⊄ Dic_N ⟺ N odd ⟺ restoration possible | Parity unifies obstruction and branching threshold |
| 8 | 8_Q8_universal_criterion.py | Complete proof: -1 ∈ [G,G] ⟺ Q₈ ⊂ G for all finite SU(2) subgroups | Uses ADE classification; Q₈ is minimal obstructed group |

### D. Spectral & Topology (files 9-11) — scripts

| # | File | What it does | Key result |
|---|------|--------------|------------|
| 9 | 9_spectral_and_topology.py | V_j decomposition tables, H¹ counting, McKay irreducibility pattern | 2T: 3 consecutive irred, 2O: 4, 2I: 6; dihedral parity toggle |
| 10 | 10_classification_free_proof.py | Classification-free proof: [a,b]=-1 ⟹ Q₈ ⊂ G via anticommutation | Pure imaginary pair generates Q₈; verified for all obstructed groups |
| 11 | 11_mckay_branch_point.py | McKay graph theorem: consecutive irred count = dist(extending node, branch)+1 | 2T(Ê₆): 3, 2O(Ê₇): 4, 2I(Ê₈): 6; SU(2) recursion walks McKay graph |

### E. General Threshold & Scope (files 12-13) — scripts

| # | File | What it does | Key result |
|---|------|--------------|------------|
| 12 | 12_general_threshold.py | Verifies dim_min = \|[G,G]\|+1 across cyclic, dicyclic, obstructed | All predictions match; obstructed groups: threshold = ∞ |
| 13 | 13_anomaly_and_scope.py | Z₂ center screening interpretation; D₄ counterexample for SO(3) groups | Q₈ criterion specific to SU(2) subgroups, not a 't Hooft anomaly |

### F. Spectral Sectors — Group-by-Group (files 14-20) — scripts

| # | File | Group | What it does | Key result |
|---|------|-------|--------------|------------|
| 14 | 14_spectral_sectors.py | 2O | Sector-resolved rigid rotor spectrum on SO(3)/O | 2 tensorial sectors (A₀,A₁); spinorial needs rank-2 bundle H₁ |
| 15 | 15_spectral_D3.py | 2D₃ | Sector-resolved spectrum on SO(3)/D₃ (non-obstructed) | 4 scalar sectors (2 tensorial + 2 spinorial); j=1/2 accessible |
| 16 | 16_jmin_spinorial.py | Cyclic+Dic | Verifies j_min^spinorial = \|[G,G]\|/2 for 12 groups | Formula confirmed: cyclic j_min=1/2, Dic_N j_min=N/2 |
| 17 | 17_spectral_2T.py | 2T | Sector-resolved spectrum on SO(3)/T (obstructed, κ=3) | All 3 sectors tensorial; spinorial blocked |
| 18 | 18_spectral_summary.py | All | Comparison table SO(3), SO(3)/T, SO(3)/O, SO(3)/D₃; Weyl check | Fraction → κ/\|G\|; obstruction decimates spectrum |
| 19 | 19_spectral_D4.py | D₃,D₄,D₅ | Dihedral parity toggle: N odd (free) vs N even (obstructed) | Parity of N controls obstruction |
| 20 | 20_spectral_2I.py | 2I | Sector-resolved spectrum on SO(3)/I (perfect group, κ=1) | Most sparse: first excited at j=6 (gap=42B) |

### G. Asymptotics & Analytic (files 21-22) — scripts

| # | File | What it does | Key result |
|---|------|--------------|------------|
| 21 | 21_weyl_convergence.py | Weyl law: N_scalar/N_total → κ/\|G\| as J→∞ | 0–2.5% error at J=50; range 1.0 (SO(3)) to 0.008 (SO(3)/I) |
| 22 | 22_onset_analytic.py | Analytic proof j_min = N/2 for Dic_N (N odd) | Confirmed N=1,3,5,7,9,11,13; from b-type vanishing |

### H. Universal Verification & Physics (files 23-27) — scripts

| # | File | Asserts | What it does | Key result |
|---|------|---------|--------------|------------|
| 23 | 23_universal_obstruction.py | 0 | -1 ∈ [G,G] test on 33 finite SU(2) subgroups | 33/33 pass: theorem universal |
| 24 | 24_berry_holonomy.py | 0 | Berry phase χ(-1) on 5 quotient spaces | Obstructed → forced bosonic; spinorial π only in non-obstructed |
| 25 | 25_O_h_extension.py | 0 | O_h vs O: does spatial inversion change obstruction? | O(3)/O_h ≅ SO(3)/O; inversion irrelevant; obstruction unchanged |
| 26 | 26_analytic_proof_and_H2.py | 0 | Three-line proof: -1 ∈ [G,G] ⟹ no spinorial 1D char; H² classification | χ([a,b])=1 for U(1); spinorial only in rank≥2 bundles |
| 27 | 27_thermodynamic_observables.py | 0 | Partition function, spectral gaps, specific heat for rigid rotors | Gap enhancement 56× (SO(3)/I vs SO(3)); rotational freezing from obstruction |

### I. Paper Verification & Challenges (files 28-31) — scripts with asserts

| # | File | Asserts | What it does | Key result |
|---|------|---------|--------------|------------|
| 28 | 28_challenge_verification.py | 29 | 116 physics assertions (SF₆, Peter-Weyl, Schottky, etc.) | All 6 challenges pass; SLD selection rule weaker than nuclear spin |
| 29 | 29_branching_2O_to_2D4.py | 16 | 2O→2D₄ restriction (cubic→tetragonal vs trigonal) | Tetragonal: obstruction persists (Q₈⊂2D₄); trigonal: re-entrance |
| 30 | 30_symmetry_breaking_phase_diagram.py | 16 | Complete O→H' phase diagram (7 subgroups) | Re-entrance ⟺ Q₈ ⊄ 2H'; 4/7 channels restore spinorial |
| 31 | 31_v3_claims_verification.py | 13 | Paper v3 claim verification (H₁ bundle, Weyl bound, odd dihedrals) | All 20 assertions pass; m(H₁,j) ~ (2j+1)/12 asymptotic |

### J. Post-Paper Results (files 32-33) — scripts with asserts

| # | File | Asserts | What it does | Key result |
|---|------|---------|--------------|------------|
| 32 | 32_weyl_remainder_bound.py | 9 | N_sc(J) = (2κ/\|G\|)(J+1)² + R(J) with explicit C(G) bound | 0 ≤ R(J) ≤ C(G); Dirichlet kernel identity; tightness 80.5% |
| 33 | 33_multiplicity_duality.py | 13 | Spectral duality: m(χ,j) + m(χ,j*-j) = 1 below threshold j*=\|G\|/4-1 | Duality ⟺ obstruction ⟺ divisibility; all ADE verified; zero fraction = 1/2 |
| 34 | 34_analytic_equivalence.py | 11 | Analytic proof: (B)⟺(C) classification-free via sum-to-product; (A)⟺(B) on ADE | χ_j+χ_{j*-j} = 2sin((j*+1)α)cos((2j-j*)α)/sin(α); vanishes iff \|H\|/d even |
| 35 | 35_stress_test_duality.py | 833 | Stress test duality on Dic_n for n=2..101 | A₀ duality, binary multiplicities, all scalar chars, j* last zero, triple equivalence |
| 36 | 36_reentrance_thermodynamics.py | 22 | Sector-resolved spectral comparison O→D₃ with physical units | Sector correspondence A₀→χ₀, A₁→χ₂; densification 3.3×/6×; Schottky peak; convergence |
| 37 | 37_frobenius_reformulation.py | 23 | V_j⊕V_{j*-j}\|_G ≅ ρ_reg⁺ ⟺ obstruction; full irrep duality | Peter-Weyl uniformity; m(ρ,j)+m(ρ,j*-j)=dim(ρ); threshold V_{j*}=ρ_reg⁺−A₀ |
| 38 | 38_last_zero_proof.py | 43 | Proof: m(χ,j) ≥ 1 above last zero, for all scalar χ | S(j) periodic P=\|G\|/4; m(j+P)=m(j)+1; Dic_n fully analytic (all 4 chars); polyhedral 48-value; m(χ,0)=0→m(χ,j*)=1; cosets from group multiplication |
| 39 | 39_H1_bundle_symmetry_breaking.py | 296 | H₁ bundle spectrum under O→D₃; branching consistency; vector bundle duality | Gap closes 15B→3B; scalar spinorial from L only; density 6×; m(ρ,j)+m(ρ,j*-j)=dim(ρ) for all irreps; ρ_reg⁻ verified |
| 40 | 40_particle_hole_exploration.py | 0 | Mode decomposition of non-central sum; particle-hole structure; re-entrance as mode flip | S(j) = Σ A_α · f_α(j) exact to 10⁻¹⁴; α↔π-α collapse (T:2, O:3, I:4 modes); cond(B)⟹ mode antisymmetry; O→D₃ = d=2 mode parity flip |

### K. Open Problem Investigations (files 41-43) — pytest suites

| # | File | Tests | Asserts | What it does | Key result |
|---|------|-------|---------|--------------|------------|
| 41 | 41_commutator_width.py | 8 | 18 | Commutator width of -1: single [a,b] or product needed? | Width=1 always; witness [a^{n/2},b]=-1 for 2D_n; all obstructed groups |
| 42 | 42_frobenius_commutator.py | 7 | 10 | Frobenius formula N(-1)=\|G\|(n₊-n₋) verified computationally | n₊≠n₋ ⟺ obstructed for all ADE; s=3 (Q₈,2D_n even), s=2 (2O), s=1 (2T,2I) |
| 43 | 43_class_counting.py | 4 | 3 | Class counting under G→H=G/{±1}: s = #{self-conjugate classes} | s = 2·#classes(H)-#classes(G); self-conj ⟺ half-angle π/2 (order 2 in SO(3)) |

### L. Paper Finalization (files 48-54) — pytest suites

| # | File | Tests | Asserts | What it does | Key result |
|---|------|-------|---------|--------------|------------|
| 48 | 48_stress_tests.py | 9 | 47 | Comprehensive stress tests: Q₈ minimal, brute force all chars, cl(-1), defect, vector bundle, Hasse | All scalar chars duality ✓; δ(j)=(-1)^j; binary m∈{0,1} rank-1 specific; Q₈ subset ⟺ commutator |
| 49 | 49_2Dn_analytic.py | 12 | 19 | Analytic proof: (B) for 2D_n ⟺ n even; K=2·gcd(k,n) formula; mode antisymmetry | Reflection K=n is sole obstruction; polyhedral K all even; k'K sign correct |
| 50 | 50_structural_chain.py | 12 | 25 | Structural chain V₄⟺Q₈⟺single commutator⟺s>0; ADE-free analysis | Trace argument; self-conjugacy=commutator; V₄⟹perp axes⟹Q₈; regular rep; (A)⟹(B) structural |
| 51 | 51_table_verification.py | 197 | 197 | Entry-level verification of paper tables: 2O, 2T, 2I + obstruction consistency | All m(χ,j) pinned; dual column cross-check; 2O corrected table; 2T κ=3 A₁=A₂ conjugacy; 2I 30-row; -1∈[G,G] check |
| 52 | 52_converse_formula.py | 112 | 112 | Converse formula m(A₀,j*)=2\|S_odd\|/\|H\|; S_odd from group elements; condition (B) direct | C_n n=2..20, D_n n=3..15 odd: \|S_odd\|=\|H\|/2, m(j*)=1; D_n even + T,O,I: S_odd=∅, cond(B) verified |
| 53 | 53_above_threshold.py | 266 | 266 | Above-threshold breakdown: m(j*+1)=2 (j*+1=P), periodicity, extended duality sum=2, pinned values | 2T,2O,2I: binary breaks at j*+1; duality sum shifts 1→2; self-consistency cross-check; 2O A₁ |
| 54 | 54_reentrance_D_odd.py | 81 | 81 | Re-entrance duality in D_n (n odd): spinorial duality holds, tensorial fails | D₃..D₁₁: spinorial m(j)+m(j*-j)=1; tensorial sums ∈ {0,2}; d=2 mode vanishing; periodicity m(j+P)=m(j)+1 |
| 55 | 55_su3_counterexample.py | 318 | 318 | SU(3) negative result: Δ(27) has no complementary pairing | Δ(27)⊂SU(3): 9 chars, N-ality selection (p-q≡0 mod 3); no constant m(χ,(p,q))+m(χ,(P-p,Q-q)) at any (P,Q); no affine involution on support (deduped pairs); diagonal slice; dim degree 2; all 8 non-trivial chars identical; class function property; (p,q)↔(q,p) symmetry; imaginary-part guard on multiplicities |
| 56 | 56_bridge_B2_B3.py | 94 | 94 | Bridge tests: d=2 uniqueness (B2) + periodicity chain (B3) + integrality argument | B2: f(π/2,j)=0 identic at half-int j, f(π/d,j)≠0 for d≥3; analytic reason (2j+1 even); **arithmetic uniqueness proof** (d\|2k ∀k ⇒ d\|2 ⇒ d=2, no sampling). B3: exp(H)=\|H\|/2 computed from group elements (non-tautological); SO(3) orders cross-checked; mode period=d exact; m(A₀,j+P)=m(A₀,j)+1; threshold sharpness m(A₀,j*+1)=2. **Integrality argument** (classification-free proof exp(H)=\|H\|/2): step 1 P\|\|H\|/2 from (B); step 2 noncentral periodic; step 3 central shift=2P/\|H\|; step 4 shift∈Z; step 5 0<2P/\|H\|≤1 ∧ ∈Z ⇒ =1 ⇒ P=\|H\|/2 |
| 57 | 56_perturbative_stability.py | 19 | 19 | Perturbative stability H1+H2: duality under SO(3)-invariant perturbation + lifting criterion | H1: j good quantum number ⟹ m(χ,j) invariant; H2: (B) fails for K ⟺ Q₈⊄2K, verified all O-subchannels |
| 58 | 58_vector3_sylow.py | 47 | 47 | Sylow-based classification-free proof: (B) ⟹ Syl₂ non-cyclic ⟹ V₄ ⊂ H ⟹ Q₈ ⊂ 2H ⟹ (A) | Full chain verified on all ADE groups; v₂(ord) < v₂(\|H\|) arithmetic; involution count ≥ 3 for (B)-groups |
| 59 | 59_gen_quaternion_no_so3.py | 77 | 77 | Q_{2^n} ⊄ SO(3) + (A)⟹(B) gap closure via Burnside orbit-counting | Axis forcing (1440-dir scan); commutation contradiction; V₄ orthogonality (2θ formula); end-to-end on D₂,D₄,D₆,T,O,I; **V₄⊂H ⟺ (B) on 16 groups**; **Burnside eq. for 2-groups: only cyclic/dihedral (a=1..20)**; dihedral exponent bound (a=2..7); **Sylow transfer: v₂(ord h)<v₂(\|H\|) on full H** (7 groups); **contrapositive: C₄,D₃,D₅,D₇ have witness violating (B₃)** |

### M. Investigation Phase (files 60-63) — pytest suites

| # | File | Tests | What it does | Key result |
|---|------|-------|--------------|------------|
| 60 | 60_nuclear_spin_isomers.py | 16 | I6: Nuclear spin weights on SO(3)/O and O→D₃ | g_spinorial=0 always (ρ_spin(-1)=+Id); statistical ratios 5:1→2:1 at O→D₃ |
| 61 | 61_thermodynamic_reflection.py | 5 | I7: Partition function reflection identity Z_χ+Z_refl=Z_trunc | Identity verified to 10⁻¹⁰; NOT a T-duality; restatement of m+m'=1 |
| 62 | 62_gap_universality.py | 42 | I3: Gap universality for non-obstructed groups | Fails for D_n odd: j_min=n/2, gap=n(n+2)B/4; cyclic universal at 1/2 |
| 63 | 63_berry_phase.py | 8 | I17: Berry phase on SO(3)/H configuration space | Berry phase = χ(-1) = sector label; reformulation, no new content |
| 64 | 64_palindromic_raman.py | 17 | I22: Raman complementarity from reflection identity | Paired Raman lines (j, j*-j-2) cannot both be allowed; verified 2T, 2O, 2I; SF₆ A₀ has 3 lines |
| 65 | 65_molien_gorenstein.py | 14 | I26: Molien series and Gorenstein property | m(A₀,j)=[t^{2j}]M(t); reflection = Gorenstein palindrome; rational forms verified 2T,2O,2I |
| 66 | 66_selective_duality_general.py | 39 | I30: Generalized sector-selective duality | Duality ⟺ immunity to (B)-violating classes; D_n odd spinorial only; polyhedral all chars |
| 67 | 67_sf6_matrix_intensities.py | 21 | I29: SF₆ nuclear spin weights and stick spectra | g(A₀)=10, g(A₁)=2 for O; g(χ₀)=16, g(χ₂)=8 for D₃; reflection fails under D₃; line density 1/2→2/3 |
| 68 | 68_spectral_zeta_heat.py | 26 | D2: Heat kernel / spectral zeta investigation | a₋₁=1/\|H\| (standard); no Seeley-DeWitt coeff distinguishes (B); Weyl ratio monotone; NEGATIVE |
| 69 | 69_s3_quadratic_rigidity.py | 16 | D4: S³/Γ quadratic rigidity | Δ²=P ⟺ (B); m_S³(nP−1)=Pn²/2 exact; non-(B) oscillates {P−1,P+1} |
| 70 | 70_chebyshev_mckay_spectral.py | 22 | D5: Chebyshev–McKay spectral interpretation | χ_j+χ_{j*-j}=2·T·U_{P-1}; U_{P-1}=0 at McKay eigenvalues ⟺ (B); μ palindrome universal |
| 71 | 71_descent_and_mckay.py | 38 | C1/C2/C4: Descent theorem, McKay eigenvalues, Coxeter check | C1 POSITIVE: reflection = universal palindrome + staircase, (B) = descent condition; C2 NEGATIVE: palindrome = Lagrange, not Poincaré; C4 PARTIAL NEG: affine Dynkin eigs = class traces but ≠ Coxeter eigs |
| 72 | 72_twisted_hilbert.py | 30 | RP3: Twisted Hilbert functions for semi-invariant rings | NEGATIVE as discriminator: binary palindromic numerator universal (Gorenstein), independent of (B). All scalar chars of 2T, 2O, 2I + D_n have rational forms with binary palindromic numerators and shared hsop |
| 73 | 73_rigidity_and_spin4.py | 20 | §4/§5/§6: Rank-1 rigidity, McKay involution, Spin(4) test | §0 NEG: McKay bipartiteness ≠ (B) (D̂_n trees → trivially bipartite). §1 POS: rank-1 rigidity (staircase slope=1 on real 2T/2O/2I data + SU(3) dim non-constant). §2 NEG: Spin(4) product 50% fail (2T,2O,2I,mixed). §3 NEG: diagonal values>1 (2T,2O,2I) |
| 74 | 74_blue_phase.py | 22 | §3: Blue phase application | No new math. BP I/II both point group O. Defect stabilizers C₄,C₃,C₂ all cyclic → spinorial. j=2 nematic: forbidden in O, allowed at defects (C₂ has 3 modes). All matches paper Tables 2,5,6 |

---

## Summary Statistics

| Category | Files | pytest tests | assert count |
|----------|-------|-------------|--------------|
| Core obstruction (A) | 2 | 141 | 171 |
| Branching rules (B) | 3 | 0 | 29 |
| Structural proofs (C) | 3 | 0 | 0 |
| Spectral & topology (D) | 3 | 0 | 1 |
| Threshold & scope (E) | 2 | 0 | 0 |
| Spectral per-group (F) | 7 | 0 | 3 |
| Asymptotics (G) | 2 | 0 | 0 |
| Universal & physics (H) | 5 | 0 | 0 |
| Paper verification (I) | 4 | 0 | 74 |
| Post-paper (J) | 9 | 0 | 1244 |
| Open problem (K) | 3 | 19 | 31 |
| Paper finalization (L) | 12 | 1244 | 1302 |
| Investigation (M) | 13 | 294 | ~294 |
| **TOTAL** | **70** | **1740** | **~3236** |

Note: files 51-53 added Feb 2026 (tracker_v4 items 0a, 0b, 0e). File 54 added Feb 2026 (item 0h). File 55 added Feb 2026 (tracker_v5 item 0f). File 56 added Feb 2026 (tracker_v5 items B2, B3). Files 57-59 added Feb 2026 (A1 classification-free proof). Files 60-63 added Feb 2026 (investigation phase: I6, I7, I3, I17). Files 64-66 added Feb 2026 (I22, I26, I30). Files 67-70 added Feb 2026 (I29 SF₆, D2 heat/zeta, D4 S³ rigidity, D5 Chebyshev-McKay). File 71 added Feb 2026 (C1 descent, C2 Poincaré neg, C4 Coxeter neg). File 72 added Feb 2026 (RP3 twisted Hilbert neg).

---

## Open Directions (from tracker_v2.md)

- ~~Classification-free proof of Q₈ ⊂ G ⟺ -1 ∈ [G,G] without ADE~~ **DONE** (files 58-59: Sylow + Q_{2^n} ⊄ SO(3))
- Spectral gap lower bound λ₁ ≥ C(G) (medium)
- Extension to SU(3) (hard)
- BaTiO₃ re-entrance temperature (low-medium)
- Forbidden Raman/IR lines for C₆₀ (medium)
- Multiplicity duality: Chebyshev reformulation (7b), particle-hole operator (7c), McKay connection (7d)
- ~~Classification-free (A)⟺(B)~~ **DONE** (= item above)
- Mode decomposition for paper: S(j)=Σ A_α·f_α(j) gives constructive proof of duality (each mode antisymmetric under cond(B)); collapse to T:2, O:3, I:4 independent modes; re-entrance = d=2 mode parity flip (file 40)
- Invariant ring structure (speculative, wip/w_1_mass_spectrum/): Molien series (1+t⁹)/((1-t⁴)(1-t⁶)) for O on S²; generators K₄,K₆,K₉; ratio 3/2 from Coxeter exponents B₃; connection to constituent quark mass ratio m_s/m_ud = 3/2 (NOT for paper)
