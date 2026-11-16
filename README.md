# PC-SAFT modeling of copolymer + solvent VLE (Arce et al., 2006)

Reproduction of the **PC-SAFT** part (Methodology 1, pseudo-monomer) of:

> Arce, P.F.; Mattedi, S.; Aznar, M. (2006). *Vapor–liquid equilibrium of
> copolymer + solvent mixtures: Thermodynamic modeling by two theoretical
> equations of state.* Fluid Phase Equilibria **246**, 52–63.
> DOI: 10.1016/j.fluid.2006.05.009

The LGT equation of state is intentionally ignored (out of scope).

## Files

```
pcsaft_vle.py            MAIN program (EDITABLE PARAMETERS section at the top)
pcsaft_core.py           PC-SAFT EoS: Eqs. 18-40 / Gross-Sadowski 2001 App. A
vle.py                   pseudo-monomer build + bubble pressure + Eq. 40
constantes_universales.py  universal a_i/b_i constants (GS2001 Table 1, ref [23])
parametros_articulo.py   parameters taken ONLY from Arce (2006) Tables 2,3,4b,6
datos_externos.py        external Antoine data, used ONLY for validation
datos_experimentales/    experimental VLE data, "w1, P[MPa]", one file per system
diagnostico_kij.py       optimal-k_ij diagnostic (does not change the reproduction)
INFORME.md               technical report (results, verification, discrepancies)
```

Separation of sources is strict: the universal constants, the paper's
parameters and any external data each live in their own module.

## Requirements

Python 3 with **numpy, scipy, matplotlib** only.

```bash
python -m venv .venv && . .venv/bin/activate
pip install numpy scipy matplotlib
```

(On a headless machine set `MPLCONFIGDIR` to a writable dir, e.g.
`export MPLCONFIGDIR=/tmp/mpl`.)

## Usage

```bash
python pcsaft_vle.py                 # validation + comparison table + all plots
python pcsaft_vle.py --validate      # only the two validation steps
python pcsaft_vle.py --table         # only the comparison table (+ plots)
python pcsaft_vle.py --table --no-plots
python pcsaft_vle.py --system cC6 PSBR45   # a single (solvent, copolymer)
```

Outputs are written next to the script:
`VLE_<solvent>_<copolymer>.png` (model line + experimental points) and
`tabla_comparativa_M1.txt` (calculated vs published %ΔP/P).

## Validation summary

- **Pure solvent Psat** (cyclohexane, benzene) within 99–101 % of Antoine.
- **Mixture engine** matches the independent `feos` PC-SAFT library to 4–5
  significant figures (pressures and fugacity coefficients).
- **Pseudo-monomer parameters** reproduce Table 5 of Arce exactly.

See `INFORME.md` for the full results table and a per-system discussion of
which systems reproduce and which do not, with the likely causes.
