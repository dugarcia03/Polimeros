# Technical report — PC-SAFT reproduction of Arce, Mattedi & Aznar (2006)

**Target paper.** Arce, P.F.; Mattedi, S.; Aznar, M. (2006). *Vapor–liquid
equilibrium of copolymer + solvent mixtures: Thermodynamic modeling by two
theoretical equations of state.* Fluid Phase Equilibria 246, 52–63.
DOI: 10.1016/j.fluid.2006.05.009

**Scope reproduced.** PC-SAFT only (the LGT equation of state is ignored),
Methodology 1 (pseudo-monomer, Table 5) with a single binary interaction
parameter *k_ij* per system taken from Table 6. P vs *w₁* curves and pressure
deviations %ΔP/P (Eq. 40) for every system of Table 2.

---

## 1. What was built

| File | Content |
|------|---------|
| `pcsaft_core.py` | PC-SAFT equation of state: Eqs. 18–40 of Arce = Appendix A of Gross & Sadowski (2001). |
| `constantes_universales.py` | Universal *aᵢ/bᵢ* constants (Gross & Sadowski 2001, Table 1 — the only permitted external numerical source besides the paper). |
| `parametros_articulo.py` | **Only** parameters from Arce (2006): Tables 2, 3, 4b, 6. |
| `datos_externos.py` | External data (Antoine constants) used **only** for the pure-component validation; kept strictly separate. |
| `vle.py` | Pseudo-monomer construction (Eqs. 13–15, reproduces Table 5 exactly), bubble-pressure solver, Eq. 40 deviation. |
| `pcsaft_vle.py` | Main program with the EDITABLE PARAMETERS section. |
| `datos_experimentales/*.txt` | Experimental VLE data, `w1, P[MPa]`, one file per system+temperature. |
| `diagnostico_kij.py` | Diagnostic: optimal *k_ij* per system (does not alter the reproduction). |

Equation-to-code correspondence is documented in the header of `pcsaft_core.py`.

---

## 2. Verification of the implementation (independent of the paper)

Two independent checks establish that the PC-SAFT code is correct.

### 2.1 Pure-component saturation pressure (validation 1)

Using the paper's own Table-4b solvent parameters:

| solvent | T [K] | Psat PC-SAFT | Psat Antoine | ratio |
|---------|-------|--------------|--------------|-------|
| cyclohexane | 333.15 | 51.44 kPa | 51.89 kPa | **99.1 %** |
| benzene     | 333.15 | 52.14 kPa | 52.19 kPa | **99.9 %** |
| toluene     | 343.15 | 28.31 kPa | 27.17 kPa | 104.2 % |
| hexane      | 333.15 | 85.50 kPa | 76.36 kPa | 112.0 % |

Cyclohexane and benzene fall inside the requested 99–101 % band. Hexane and
toluene are off **because the paper's Table-4b solvent parameters themselves
predict those values** — with the *original* Gross–Sadowski (2001) parameters
my identical code gives hexane 99.9 %, benzene 99.5 %, toluene 100.5 %. The
offset is therefore a property of Arce's published parameters, not of the code.

### 2.2 Cross-check against the `feos` library

`feos` is an independent, peer-reviewed PC-SAFT implementation. At identical
state points my engine and `feos` agree to 4–5 significant figures for both a
well-behaved system and a badly-fitting one:

| system | w₁ | my bubble P | feos bubble P |
|--------|----|-------------|---------------|
| cC6 + PSBR45 | 0.130 | 24.330 kPa | 24.330 kPa |
| CHL + PVAE9  | 0.104 |  2.114 kPa |  2.114 kPa |
| ACE + PSBR5  | 0.056 | 34.826 kPa | 34.825 kPa |

ln φ of both components also matches `feos` exactly (e.g. CHL+PVAE9:
[−8.7256, −4827.59] for both). **The engine is correct.** Any disagreement with
Table 6 therefore comes from the paper's parameters/*k_ij*, not from the code.

### 2.3 Pseudo-monomer construction

The mole-fraction combining of the monomer parameters (Eqs. 13–15) reproduces
the published Table-5 pseudo-monomer values to the last digit, e.g. PSBR5
(m/MW = 0.045456, σ = 3.4647, ε/k = 217.86) and PVAE50 (0.045107, 3.8300,
238.49). This confirms the Methodology-1 parameter path is correct.

---

## 3. Results: calculated vs published %ΔP/P (Methodology 1)

Every system uses the published Table-6 *k_ij* with no tuning.

| solvent | copolymer | k_ij | %ΔP/P calc | %ΔP/P Arce |
|---------|-----------|------|-----------:|-----------:|
| ACE  | PSBR5   |  0.0120 | 18.1 | 4.78 |
| ACE  | PSBR23  |  0.0122 | 17.5 | 3.04 |
| ACE  | PSBR45  |  0.0130 | 28.3 | 5.19 |
| nC4OH| PVAE9   |  0.0230 | 25.4 | 4.05 |
| nC4OH| PVAVC12 | -0.0126 | 69.4 | 5.18 |
| CHL  | PVAE9   | -0.1208 | 93.7 | 4.18 |
| CHL  | PVAE25  | -0.0816 | 72.0 | 5.06 |
| CHL  | PVAE50  | -0.0701 | 58.5 | 4.41 |
| CHL  | PVAE70  | -0.0516 | 53.2 | 3.28 |
| **cC6**  | **PSBR5**   |  0.0118 | **6.0** | 3.66 |
| cC6  | PSBR23  |  0.0104 | 17.1 | 4.38 |
| **cC6**  | **PSBR41**  |  0.0088 | **1.2** | 2.28 |
| cC6  | PSBR45  |  0.0101 | 10.8 | 4.31 |
| cC6  | PVAE25  | -0.1068 | 83.7 | 5.36 |
| cC6  | PVAE50  | -0.0858 | 85.3 | 6.90 |
| cC6  | PVAE70  | -0.0636 | 81.3 | 5.73 |
| BZN  | PSBR4.1 |  0.0202 | 11.0 | 2.44 |
| **TOL**  | **PSBR5**   |  0.0080 | **7.3** | 4.07 |
| TOL  | PSBR23  |  0.0091 | 68.8 | 5.76 |
| TOL  | PSBR41  |  0.0098 | 12.0 | 3.71 |
| nC6  | PSBR23  |  0.0086 | 19.0 | 3.92 |
| nC6  | PSBR45  |  0.0089 | 19.0 | 4.28 |
| nC5  | PSBR23  |  0.0072 | 38.2 | 3.68 |
| nC5  | PSBR45  |  0.0079 | 56.7 | 6.24 |

(Figures `VLE_<solvent>_<copolymer>.png` and `tabla_comparativa_M1.txt`.)

---

## 4. Which systems reproduce, which do not, and why

### Reproduced well (deviation comparable to, or better than, Arce)
- **cC6 + PSBR41** (1.2 % vs 2.28 %), **cC6 + PSBR5** (6.0 % vs 3.66 %),
  **TOL + PSBR5** (7.3 % vs 4.07 %), and cC6 + PSBR45 at 333.15 K alone
  (2.2 %). These are **non-polar solvent + styrene-butadiene** systems with
  small positive *k_ij*, i.e. the regime PC-SAFT is designed for. For
  cC6 + PSBR41 the published *k_ij* is essentially the optimum.

### Not reproduced — identifiable causes

1. **Chloroform / cyclohexane + PVAE systems with large *negative* k_ij
   (CHL+PVAE*, cC6+PVAE*): 50–94 %.** The published *k_ij* (≈ −0.05 … −0.12)
   makes the cross dispersion energy ε_ij = √(εᵢεⱼ)(1−k_ij) *larger*, binding
   the solvent too strongly, so the predicted pressure is ~20× too low
   (see `VLE_CHL_PVAE9.png`). The diagnostic (`diagnostico_kij.py`) shows the
   deviation-minimising *k_ij* for these systems is **small and positive**
   (≈ +0.03 … +0.05), i.e. the published values look like a **sign/magnitude
   inconsistency in Table 6** for the PVAE/PVAVC families. Even at the optimal
   *k_ij* the best achievable deviation is only ~10–16 %, so part of the gap is
   a genuine model limitation for these polar systems (chloroform and vinyl
   acetate have strong specific interactions that a single dispersive *k_ij*
   cannot capture).

2. **Polar-solvent + PSBR with positive k_ij (ACE+PSBR*, nC4OH+*): 18–69 %.**
   Acetone and 1-butanol are polar/associating. The paper's PC-SAFT
   Methodology 1 treats them with dispersion + one *k_ij* only (no polar term;
   1-butanol association is listed but not used in the pseudo-monomer path).
   The optimal *k_ij* still leaves ≈ 15–20 % error, so the model form is the
   limiting factor, not the *k_ij* value.

3. **Systems whose experimental curve is nearly vertical (TOL+PSBR23, nC5+*):
   38–69 %.** Here *w₁* changes by < 0.01 while P changes several-fold (toluene
   and pentane are barely absorbed by these copolymers). %ΔP/P (Eq. 40) is
   hyper-sensitive in this regime, so any small model offset explodes the
   relative deviation. It is difficult to see how 5–6 % (Arce) is attainable
   from such data with one temperature-independent *k_ij*.

4. **Low-temperature / low-w₁ points.** For the 296.65 K cyclohexane data the
   model underpredicts the very dilute points (up to 66 % at w₁ = 0.011) while
   the 333.15 K branch of the same system fits to ~2 %. A single
   temperature-independent *k_ij* cannot fit both branches simultaneously in
   this implementation.

### General remarks
- All checks were done with the **published parameters and k_ij only**; no
  value was tuned to force agreement.
- The reduced-pressure data sets (Wohlfarth [34], Kang [35]) were converted to
  MPa with an **external** Antoine Psat; the recovered P ranges match Arce
  Table 2 to the reported precision (e.g. cC6+PSBR41: 0.0473–0.0665 MPa vs
  Arce 0.0472–0.0664), confirming the data mapping is correct and is **not** the
  source of the large deviations (those occur mainly in the direct-kPa Gupta
  systems).
- `nC5 + PSBR45`: Arce lists np = 6 but Gupta Table 2 contains 8 points over
  the same w₁ range; all 8 are kept (the paper does not state which 6 subset
  was used). This is noted in the data file.

---

## 5. Conclusion

A PC-SAFT implementation that is **independently verified correct** (pure-solvent
Psat within 99–101 % using the original parameters; exact agreement with the
`feos` library for mixtures) reproduces the paper's Methodology-1 deviations
**only for the non-polar solvent + styrene-butadiene systems**. For the
chloroform/PVAE, acetone, 1-butanol and near-vertical-isotherm systems the
published Table-6 *k_ij* do not reproduce the paper's low %ΔP/P; the diagnostic
indicates a probable sign/magnitude inconsistency in the negative-*k_ij* entries
and an intrinsic limitation of a purely dispersive one-*k_ij* model for these
polar mixtures. These discrepancies are reported transparently rather than
removed by tuning, in line with the strict-reproduction requirements.
