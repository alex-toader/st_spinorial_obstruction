#!/usr/bin/env python3
"""
Frobenius reformulation of spectral duality.

MAIN RESULT:
  For obstructed G ⊂ SU(2) with 4 | |G|, and j* = |G|/4 - 1:

    V_j ⊕ V_{j*-j}|_G  ≅  ρ_reg⁺     for all j ∈ [0, j*]

  where ρ_reg⁺ = bosonic half of the regular representation:
    χ_{ρ_reg⁺}(g) = |G|/2  for g = ±1  (central)
    χ_{ρ_reg⁺}(g) = 0      for g ≠ ±1  (non-central)

CONSEQUENCES:
  1. m(ρ,j) + m(ρ,j*-j) = dim(ρ) for every bosonic irrep ρ of G.
     For 1D characters: m(χ,j) + m(χ,j*-j) = 1 (= spectral duality).
  2. The sub-threshold Peter-Weyl decomposition is perfectly uniform:
       ⊕_{j=0}^{j*} V_j|_G  =  (j*+1)/2  copies of ρ_reg⁺
     Each bosonic irrep ρ appears with total multiplicity (j*+1)/2 · dim(ρ).
  3. V_{j*}|_G = ρ_reg⁺ - A₀ (threshold representation).

EQUIVALENCE:
  V_j ⊕ V_{j*-j}|_G = ρ_reg⁺  ⟺  obstruction (A): -1 ∈ [G,G]
  Verified on all obstructed groups (2T, 2O, 2I, Dic_{n even})
  and shown to FAIL on all non-obstructed groups (Dic_{n odd}).

PROOF STRUCTURE:
  The character identity χ_j(g) + χ_{j*-j}(g) = 0 on non-central g
  follows from the sum-to-product identity:
    χ_j + χ_{j*-j} = 2 sin((j*+1)α) cos((2j-j*)α) / sin(α)
  Under divisibility (B), (j*+1)α = integer·π, so sin = 0.
  The trigonometry is minimal (sin(nπ) = 0), but unavoidable:
  the vanishing is arithmetic (divisibility of |G|/4 by the element order),
  not purely representation-theoretic.

  The conceptual gain is in the FORMULATION, not the proof technique:
  duality is the statement that paired SU(2) representations restrict
  to the regular representation. This connects to Frobenius reciprocity
  and makes the duality a structural fact about G ⊂ SU(2).

Asserts: 23
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from src.quaternion import qkey, qmul
from src.group import (build_binary_tetrahedral, build_binary_octahedral,
                       build_binary_icosahedral,
                       compute_conjugacy_classes, compute_commutator_subgroup)


def chi_su2(j, trace):
    """SU(2) spin-j character at group element with given trace."""
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
    """Build Dic_N = 2D_N of order 4N."""
    elements = []
    for k in range(2 * N):
        angle = k * np.pi / N
        elements.append(np.array([np.cos(angle), 0, 0, np.sin(angle)]))
    b = np.array([0, 0, 1, 0], dtype=float)
    for k in range(2 * N):
        angle = k * np.pi / N
        ak = np.array([np.cos(angle), 0, 0, np.sin(angle)])
        elements.append(qmul(b, ak))
    return np.array(elements)


def rho_reg_plus(trace, G):
    """Character of bosonic regular representation."""
    if abs(abs(trace) - 2) < 0.01:
        return G / 2
    return 0.0


def main():
    n_asserts = 0

    # ===========================================================
    # PART 1: CHARACTER IDENTITY ON POLYHEDRAL GROUPS
    # ===========================================================
    print("=" * 70)
    print("PART 1: V_j + V_{j*-j} = rho_reg+ (polyhedral groups)")
    print("=" * 70)

    groups = [
        ('2T', build_binary_tetrahedral()),
        ('2O', np.array(build_binary_octahedral())),
        ('2I', build_binary_icosahedral()),
    ]

    for name, elements in groups:
        cls = compute_conjugacy_classes(elements)
        cls.sort(key=lambda c: -c['trace'])
        G = len(elements)
        j_star = G // 4 - 1

        # Check character identity on every class, every j
        max_err = 0.0
        for j in range(j_star + 1):
            jp = j_star - j
            for c in cls:
                val = chi_su2(j, c['trace']) + chi_su2(jp, c['trace'])
                expected = rho_reg_plus(c['trace'], G)
                err = abs(val - expected)
                max_err = max(max_err, err)

        n_j = j_star + 1
        n_cls = len(cls)
        print(f"\n  {name}: |G|={G}, j*={j_star}, "
              f"checked {n_j} × {n_cls} = {n_j*n_cls} values, "
              f"max error = {max_err:.2e}")
        assert max_err < 1e-10, f"{name}: max_err = {max_err}"; n_asserts += 1

    # ===========================================================
    # PART 2: DICYCLIC FAMILY — EQUIVALENCE WITH OBSTRUCTION
    # ===========================================================
    print("\n" + "=" * 70)
    print("PART 2: DICYCLIC FAMILY (Dic_n, n=2..30)")
    print("=" * 70)

    obstructed_pass = 0
    obstructed_total = 0
    non_obstructed_fail = 0
    non_obstructed_total = 0

    for n in range(2, 31):
        elements = build_dicyclic(n)
        G = len(elements)  # 4n
        j_star = G // 4 - 1  # n-1
        cls = compute_conjugacy_classes(elements)
        cls.sort(key=lambda c: -c['trace'])

        ok = True
        for j in range(j_star + 1):
            jp = j_star - j
            for c in cls:
                val = chi_su2(j, c['trace']) + chi_su2(jp, c['trace'])
                expected = rho_reg_plus(c['trace'], G)
                if abs(val - expected) > 0.01:
                    ok = False
                    break
            if not ok:
                break

        obstructed = (n % 2 == 0)
        if obstructed:
            obstructed_total += 1
            if ok:
                obstructed_pass += 1
        else:
            non_obstructed_total += 1
            if not ok:
                non_obstructed_fail += 1

    print(f"\n  Obstructed (n even): {obstructed_pass}/{obstructed_total} PASS")
    print(f"  Non-obstructed (n odd): {non_obstructed_fail}/{non_obstructed_total} FAIL")
    assert obstructed_pass == obstructed_total; n_asserts += 1
    assert non_obstructed_fail == non_obstructed_total; n_asserts += 1
    print(f"  Equivalence: identity holds ⟺ obstruction (verified n=2..30)")

    # ===========================================================
    # PART 3: PETER-WEYL UNIFORMITY
    # ===========================================================
    print("\n" + "=" * 70)
    print("PART 3: SUB-THRESHOLD PETER-WEYL UNIFORMITY")
    print("=" * 70)

    print(f"\n  Claim: ⊕_{{j=0}}^{{j*}} V_j|_G = (j*+1)/2 copies of ρ_reg⁺")

    for name, elements in groups:
        G = len(elements)
        j_star = G // 4 - 1
        n_pairs = (j_star + 1) // 2
        pw_dim = sum(2 * j + 1 for j in range(j_star + 1))
        reg_dim = n_pairs * G // 2
        print(f"\n  {name}: j*={j_star}, (j*+1)/2={n_pairs}")
        print(f"    Σ(2j+1) = {pw_dim},  (j*+1)/2 × |G|/2 = {n_pairs}×{G//2} = {reg_dim}")
        assert pw_dim == reg_dim; n_asserts += 1

    # Also verify for dicyclic
    for n in [2, 4, 6, 10, 20]:
        G = 4 * n
        j_star = n - 1
        n_pairs = (j_star + 1) // 2
        pw_dim = sum(2 * j + 1 for j in range(j_star + 1))
        reg_dim = n_pairs * G // 2
        assert pw_dim == reg_dim; n_asserts += 1

    print(f"\n  Dimension check passed for 2T, 2O, 2I + Dic_{{2,4,6,10,20}}")

    # ===========================================================
    # PART 4: THRESHOLD REPRESENTATION
    # ===========================================================
    print("\n" + "=" * 70)
    print("PART 4: THRESHOLD REPRESENTATION V_{j*}|_G = ρ_reg⁺ - A₀")
    print("=" * 70)

    for name, elements in groups:
        cls = compute_conjugacy_classes(elements)
        cls.sort(key=lambda c: -c['trace'])
        G = len(elements)
        j_star = G // 4 - 1

        max_err = 0.0
        for c in cls:
            chi_jstar = chi_su2(j_star, c['trace'])
            # rho_reg+ - A0: G/2 - 1 on central, -1 on non-central
            if abs(abs(c['trace']) - 2) < 0.01:
                expected = G / 2 - 1
            else:
                expected = -1.0
            max_err = max(max_err, abs(chi_jstar - expected))

        print(f"  {name}: max error = {max_err:.2e}")
        assert max_err < 1e-10; n_asserts += 1

    # ===========================================================
    # PART 5: FULL IRREP DUALITY ON 2O
    # ===========================================================
    print("\n" + "=" * 70)
    print("PART 5: FULL IRREP DUALITY m(ρ,j)+m(ρ,j*-j) = dim(ρ) [2O]")
    print("=" * 70)

    elements = np.array(build_binary_octahedral())
    cls = compute_conjugacy_classes(elements)
    cls.sort(key=lambda c: -c['trace'])
    G = 48
    j_star = 11
    sizes = np.array([c['size'] for c in cls])
    traces = [c['trace'] for c in cls]

    # Extract bosonic irreps from V_0, V_1, V_2, V_3, V_4
    chi_on_cls = {}
    for j in range(j_star + 1):
        chi_on_cls[j] = np.array([chi_su2(j, t) for t in traces])

    # Known bosonic irreps of 2O: A0(1), A1(1), E(2), T1(3), T2(3)
    # V_0 = A0, V_1 = T1
    A0 = chi_on_cls[0]
    T1 = chi_on_cls[1]
    # V_2 (dim 5): m(A0,2) and m(T1,2) = 0, so V_2 must be E+T2
    # V_3 (dim 7): contains A1+T1+T2 (1+3+3=7)
    # Extract A1 from V_3:
    m_A0_3 = round(np.sum(A0.conj() * chi_on_cls[3] * sizes).real / G)
    m_T1_3 = round(np.sum(T1.conj() * chi_on_cls[3] * sizes).real / G)
    # V_3 - m_A0_3*A0 - m_T1_3*T1 should give A1+T2 or similar
    # Actually, for the duality check we don't need to extract irreps.
    # We have T1 from V_1. Check m(T1,j)+m(T1,j*-j) = 3 for all j.

    print(f"\n  T1 (dim 3): m(T1,j) + m(T1,j*-j) for j=0..{j_star}")
    all_ok = True
    for j in range(j_star + 1):
        jp = j_star - j
        m_j = round(np.sum(T1.conj() * chi_on_cls[j] * sizes).real / G)
        m_jp = round(np.sum(T1.conj() * chi_on_cls[jp] * sizes).real / G)
        s = m_j + m_jp
        if s != 3:
            all_ok = False
        if j <= 3 or j >= j_star - 1:
            print(f"    j={j:2d}: m(T1,{j})={m_j}, m(T1,{jp})={m_jp}, sum={s}")
    print(f"  Result: {'PASS' if all_ok else 'FAIL'}")
    assert all_ok; n_asserts += 1

    # Extract E (dim 2) from V_2: V_2 = E + T2, so E+T2 = chi_on_cls[2]
    # T2 = V_2 - E. But we don't have E separately.
    # Use: E = V_2 - T2, T2 = V_3 - A1 - T1. This needs A1.
    # Easier: A1 = sign character, +1 on [2O,2O]=2T, -1 outside.
    comm = compute_commutator_subgroup(elements)
    comm_keys = set(comm.keys())
    A1 = np.array([1.0 if list(c['keys'])[0] in comm_keys else -1.0
                    for c in cls])
    # Verify A1 is correct: dim 1, A1(±1) = 1
    assert abs(A1[0] - 1) < 0.01 and abs(A1[-1] - 1) < 0.01; n_asserts += 1

    print(f"\n  A1 (dim 1): m(A1,j) + m(A1,j*-j) for j=0..{j_star}")
    all_ok = True
    for j in range(j_star + 1):
        jp = j_star - j
        m_j = round(np.sum(A1.conj() * chi_on_cls[j] * sizes).real / G)
        m_jp = round(np.sum(A1.conj() * chi_on_cls[jp] * sizes).real / G)
        s = m_j + m_jp
        if s != 1:
            all_ok = False
        if j <= 2 or j >= j_star - 1:
            print(f"    j={j:2d}: m(A1,{j})={m_j}, m(A1,{jp})={m_jp}, sum={s}")
    print(f"  Result: {'PASS' if all_ok else 'FAIL'}")
    assert all_ok; n_asserts += 1

    # T2 from V_3: V_3(7) = A1(1) + T1(3) + T2(3)
    # But m(A1,3) might not be 1. Let me check via the identity instead.
    # T2 = V_3 - m(A1,3)*A1 - m(T1,3)*T1 - m(A0,3)*A0
    m_A0_3 = round(np.sum(A0.conj() * chi_on_cls[3] * sizes).real / G)
    m_A1_3 = round(np.sum(A1.conj() * chi_on_cls[3] * sizes).real / G)
    m_T1_3 = round(np.sum(T1.conj() * chi_on_cls[3] * sizes).real / G)
    remaining_dim = 7 - m_A0_3 - m_A1_3 - 3 * m_T1_3
    # V_3 should have T2: dim 3
    T2 = chi_on_cls[3] - m_A0_3 * A0 - m_A1_3 * A1 - m_T1_3 * T1
    # Check dim
    assert abs(T2[0] - 3) < 0.01; n_asserts += 1  # T2 has dim 3

    print(f"\n  T2 (dim 3): m(T2,j) + m(T2,j*-j) for j=0..{j_star}")
    all_ok = True
    for j in range(j_star + 1):
        jp = j_star - j
        m_j = round(np.sum(T2.conj() * chi_on_cls[j] * sizes).real / G)
        m_jp = round(np.sum(T2.conj() * chi_on_cls[jp] * sizes).real / G)
        s = m_j + m_jp
        if s != 3:
            all_ok = False
        if j <= 2 or j >= j_star - 1:
            print(f"    j={j:2d}: m(T2,{j})={m_j}, m(T2,{jp})={m_jp}, sum={s}")
    print(f"  Result: {'PASS' if all_ok else 'FAIL'}")
    assert all_ok; n_asserts += 1

    # E from V_2: V_2(5) = E(2) + T2(3)
    m_T2_2 = round(np.sum(T2.conj() * chi_on_cls[2] * sizes).real / G)
    E = chi_on_cls[2] - m_T2_2 * T2 - round(np.sum(A0.conj() * chi_on_cls[2] * sizes).real / G) * A0
    # Subtract any A1, T1 if present
    m_A1_2 = round(np.sum(A1.conj() * chi_on_cls[2] * sizes).real / G)
    m_T1_2 = round(np.sum(T1.conj() * chi_on_cls[2] * sizes).real / G)
    E = E - m_A1_2 * A1 - m_T1_2 * T1
    assert abs(E[0] - 2) < 0.01; n_asserts += 1  # E has dim 2

    print(f"\n  E (dim 2): m(E,j) + m(E,j*-j) for j=0..{j_star}")
    all_ok = True
    for j in range(j_star + 1):
        jp = j_star - j
        m_j = round(np.sum(E.conj() * chi_on_cls[j] * sizes).real / G)
        m_jp = round(np.sum(E.conj() * chi_on_cls[jp] * sizes).real / G)
        s = m_j + m_jp
        if s != 2:
            all_ok = False
        if j <= 2 or j >= j_star - 1:
            print(f"    j={j:2d}: m(E,{j})={m_j}, m(E,{jp})={m_jp}, sum={s}")
    print(f"  Result: {'PASS' if all_ok else 'FAIL'}")
    assert all_ok; n_asserts += 1

    # ===========================================================
    # PART 6: PHYSICAL SUMMARY
    # ===========================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"""
FROBENIUS REFORMULATION OF SPECTRAL DUALITY

1. MAIN IDENTITY (verified on 2T, 2O, 2I + Dic_{{n=2..30}}):
   V_j ⊕ V_{{j*-j}}|_G  ≅  ρ_reg⁺    for all j ∈ [0, j*]
   This holds ⟺ G is obstructed (-1 ∈ [G,G]).

2. CONSEQUENCES:
   - 1D characters: m(χ,j) + m(χ,j*-j) = 1  [= spectral duality]
   - All bosonic irreps: m(ρ,j) + m(ρ,j*-j) = dim(ρ)
   - Peter-Weyl uniformity: ⊕_{{j=0}}^{{j*}} V_j|_G = (j*+1)/2 × ρ_reg⁺

3. THRESHOLD REPRESENTATION:
   V_{{j*}}|_G = ρ_reg⁺ - A₀
   The threshold level contains every bosonic irrep with multiplicity = dim,
   except the trivial representation (which is absent).

4. CONCEPTUAL MEANING:
   Duality is not a trigonometric accident. It is the statement that
   paired SU(2) representations restrict to the regular representation
   of G. The obstruction forces this uniformity below the threshold.

5. WHAT REMAINS TRIGONOMETRIC:
   The vanishing χ_j(g) + χ_{{j*-j}}(g) = 0 on non-central g
   reduces to sin((|G|/4)α) = 0, which holds because divisibility (B)
   forces (|G|/4)α to be an integer multiple of π. This is arithmetic
   (sin(nπ) = 0), not deep representation theory.
""")

    print(f"All checks passed ({n_asserts} asserts).")


if __name__ == '__main__':
    main()
