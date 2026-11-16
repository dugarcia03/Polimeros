#!/usr/bin/env python3
"""
=============================================================================
 PC-SAFT modeling of vapor-liquid equilibrium of copolymer + solvent mixtures
 Reproduction of:
   Arce, P.F.; Mattedi, S.; Aznar, M. (2006). Vapor-liquid equilibrium of
   copolymer + solvent mixtures: Thermodynamic modeling by two theoretical
   equations of state. Fluid Phase Equilibria 246, 52-63.
   DOI: 10.1016/j.fluid.2006.05.009

 SCOPE: PC-SAFT only (LGT ignored). Methodology 1 (pseudo-monomer), one kij
        per system from Table 6. Reproduces P vs w1 and %dP/P (Eq. 40) for all
        systems in Table 2.

 This is the MAIN program / entry point. It imports:
   - pcsaft_core.py          : the PC-SAFT equations (Eqs. 18-40 / GS2001 App.A)
   - constantes_universales.py : universal a_i/b_i constants (GS2001 Table 1)
   - parametros_articulo.py  : parameters taken only from Arce et al. (2006)
   - datos_externos.py       : external validation data (clearly separated)
   - datos_experimentales/*  : experimental VLE data (w1, P[MPa]), one per system

 Dependencies: numpy, scipy, matplotlib.
 Usage:
   python pcsaft_vle.py                  # run everything (validation+plots+table)
   python pcsaft_vle.py --validate       # only validation steps
   python pcsaft_vle.py --system cC6 PSBR45   # single (solvent, copolymer)
   python pcsaft_vle.py --table          # only the comparison table
=============================================================================
"""

import os
import sys
import argparse
import numpy as np

# Matplotlib without a display and with a writable config dir.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pcsaft_core as pc
import vle
from parametros_articulo import (
    SOLVENTS, MONOMERS, COPOLYMERS, KIJ_M1,
    DP_M1_PUBLISHED, DP_M1_GLOBAL_PUBLISHED, SYSTEMS,
)
from datos_externos import antoine_psat_MPa, SOLVENT_ANTOINE_NAME

# ===========================================================================
#                           EDITABLE PARAMETERS
# ---------------------------------------------------------------------------
# All system-dependent inputs live HERE. Nothing system-specific appears
# elsewhere in this file. The numerical parameter VALUES themselves are stored
# in the separate data modules (parametros_articulo / datos_externos); this
# section selects WHICH systems to run and the plotting/solver controls.
# ===========================================================================

# Directory with experimental data files (one per system+temperature):
DATA_DIR = "datos_experimentales"
# Where figures are written (same directory as this script, per requirements):
FIG_DIR = os.path.dirname(os.path.abspath(__file__))

# Methodology to use: 1 = pseudo-monomer (this reproduction). (2 = comonomer;
# not required by the scope, kept as a hook.)
METHODOLOGY = 1

# Number of points for the smooth model curve in P vs w1 plots:
N_CURVE_POINTS = 40

# Pure-component validation solvents (nonpolar) and temperatures [K]:
PURE_VALIDATION = [
    ("cC6", 333.15),   # cyclohexane
    ("nC6", 333.15),   # hexane
    ("BZN", 333.15),   # benzene
    ("TOL", 343.15),   # toluene
]

# Full-system validation case (paper reports 4.31 % with Methodology 1):
FULL_VALIDATION = ("cC6", "PSBR45", 333.15)

# The complete set of systems to model is taken from parametros_articulo.SYSTEMS
# (every system in Table 2). To run a subset, pass --system on the command line.

# ===========================================================================
#                      End of EDITABLE PARAMETERS
# ===========================================================================


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_data(solvent, copolymer, T):
    """Load experimental (w1, P[MPa]) for a system at temperature T.
    Returns (w1_array, P_array[MPa]) or (None, None) if file missing."""
    fname = os.path.join(DATA_DIR, f"{solvent}_{copolymer}_{T}.txt")
    if not os.path.exists(fname):
        return None, None
    arr = np.loadtxt(fname, delimiter=",", comments="#")
    arr = np.atleast_2d(arr)
    return arr[:, 0], arr[:, 1]


def psat_hint_Pa(solvent, T):
    """A pressure bracket hint [Pa]: pure-solvent Psat if we have Antoine,
    else a generic value. Used only to speed the bubble-P root finder."""
    name = SOLVENT_ANTOINE_NAME.get(solvent)
    if name is not None:
        return antoine_psat_MPa(name, T) * 1e6
    return 1.0e5


# ---------------------------------------------------------------------------
# Validation step 1: pure-component PC-SAFT saturation pressure
# ---------------------------------------------------------------------------

def validate_pure():
    from scipy.optimize import brentq
    print("=" * 74)
    print("VALIDATION 1 - Pure-solvent PC-SAFT saturation pressure vs Antoine")
    print("   (uses Arce Table 4b solvent parameters; target 99-101 %)")
    print("=" * 74)
    print(f"{'solvent':8s} {'T[K]':>8s} {'Psat_PCSAFT[kPa]':>17s} "
          f"{'Psat_Antoine[kPa]':>18s} {'ratio[%]':>9s}")
    for code, T in PURE_VALIDATION:
        mMW, sig, eps, MW = SOLVENTS[code]
        m_i = mMW * MW
        name = SOLVENT_ANTOINE_NAME[code]
        Pguess = antoine_psat_MPa(name, T) * 1e6

        def dln(P):
            lnL, _ = pc.ln_phi_pure(m_i, sig, eps, T, P, "liquid")
            lnV, _ = pc.ln_phi_pure(m_i, sig, eps, T, P, "vapor")
            return lnL - lnV

        Pcalc = brentq(dln, Pguess * 0.2, Pguess * 5.0, xtol=1.0, rtol=1e-10)
        Pcalc_MPa = Pcalc / 1e6
        Pexp_MPa = antoine_psat_MPa(name, T)
        print(f"{code:8s} {T:8.2f} {Pcalc_MPa*1e3:17.3f} "
              f"{Pexp_MPa*1e3:18.3f} {100*Pcalc_MPa/Pexp_MPa:9.2f}")
    print()


# ---------------------------------------------------------------------------
# Validation step 2: a full nonpolar solvent/polymer system
# ---------------------------------------------------------------------------

def validate_full():
    solvent, copolymer, T_focus = FULL_VALIDATION
    kij = KIJ_M1[solvent][copolymer]
    # paper's cC6+PSBR45 entry spans all temperatures in Table 2 for that system
    temps = None
    for s, c, t, _ in SYSTEMS:
        if s == solvent and c == copolymer:
            temps = t
    print("=" * 74)
    print(f"VALIDATION 2 - {solvent} + {copolymer} (paper Methodology-1 "
          f"%dP/P = 4.31, kij = {kij})")
    print("=" * 74)
    total, npts = 0.0, 0
    for T in temps:
        w1, Pexp = load_data(solvent, copolymer, T)
        if w1 is None:
            continue
        dP, Pcalc = vle.deviation_percent(w1, Pexp, T, solvent, copolymer, kij)
        print(f"  T = {T} K   ({len(w1)} points)   mean %dP/P = {dP:.2f}")
        print(f"  {'w1':>7s} {'Pexp[kPa]':>11s} {'Pcalc[kPa]':>12s} {'%dev':>7s}")
        for i in range(len(w1)):
            d = 100 * abs(Pcalc[i] - Pexp[i]) / Pexp[i]
            print(f"  {w1[i]:7.3f} {Pexp[i]*1e3:11.2f} "
                  f"{Pcalc[i]*1e3:12.3f} {d:7.2f}")
        total += dP * len(w1)
        npts += len(w1)
    print(f"  --> combined %dP/P over all temperatures = {total/npts:.2f}   "
          f"(published 4.31)\n")


# ---------------------------------------------------------------------------
# Per-system calculation: plot + deviation
# ---------------------------------------------------------------------------

def run_system(solvent, copolymer, temps, make_plot=True):
    """Compute %dP/P (combined over all temperatures, Eq. 40) and plot."""
    kij = KIJ_M1[solvent][copolymer]
    all_dev = []
    plot_sets = []  # (T, w1exp, Pexp, w1curve, Pcurve)

    for T in temps:
        w1, Pexp = load_data(solvent, copolymer, T)
        if w1 is None:
            print(f"  [skip] no data file for {solvent}_{copolymer}_{T}")
            continue
        dP, Pcalc = vle.deviation_percent(w1, Pexp, T, solvent, copolymer, kij)
        all_dev.append((dP, len(w1)))
        if make_plot:
            w1c, Pc = vle.model_curve(
                T, solvent, copolymer, kij,
                max(1e-4, w1.min() * 0.5), w1.max() * 1.05,
                npts=N_CURVE_POINTS, P_hint_Pa=psat_hint_Pa(solvent, T))
            plot_sets.append((T, w1, Pexp, w1c, Pc))

    # combined deviation over all points (Eq. 40 over the full system)
    if all_dev:
        num = sum(d * n for d, n in all_dev)
        den = sum(n for _, n in all_dev)
        dP_sys = num / den
    else:
        dP_sys = float("nan")

    if make_plot and plot_sets:
        _plot(solvent, copolymer, plot_sets, dP_sys)

    return dP_sys


def _plot(solvent, copolymer, plot_sets, dP_sys):
    plt.figure(figsize=(6.2, 4.6))
    markers = ["s", "o", "^", "v", "D", "P"]
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
    for i, (T, w1, Pexp, w1c, Pc) in enumerate(plot_sets):
        c = colors[i % len(colors)]
        plt.plot(w1, Pexp, markers[i % len(markers)], color=c, ms=5,
                 label=f"exp {T:.2f} K")
        plt.plot(w1c, Pc, "-", color=c, lw=1.6,
                 label=f"PC-SAFT {T:.2f} K")
    plt.xlabel("solvent mass fraction  $w_1$")
    plt.ylabel("Pressure  P / MPa")
    plt.title(f"{solvent} (1) + {copolymer} (2)   "
              f"PC-SAFT M1,  %$\\Delta P/P$ = {dP_sys:.2f}")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    out = os.path.join(FIG_DIR, f"VLE_{solvent}_{copolymer}.png")
    plt.savefig(out, dpi=130)
    plt.close()
    print(f"  figure -> {os.path.basename(out)}")


# ---------------------------------------------------------------------------
# Comparison table (calculated vs published %dP/P)
# ---------------------------------------------------------------------------

def comparison_table(make_plots=True):
    print("=" * 74)
    print("COMPARISON TABLE - calculated vs published %dP/P (Methodology 1)")
    print("=" * 74)
    header = (f"{'solvent':7s} {'copolymer':10s} {'kij':>9s} "
              f"{'%dP/P calc':>11s} {'%dP/P Arce':>11s} {'diff':>7s}")
    print(header)
    print("-" * len(header))
    rows = []
    for solvent, copolymer, temps, _src in SYSTEMS:
        kij = KIJ_M1[solvent][copolymer]
        dP = run_system(solvent, copolymer, temps, make_plot=make_plots)
        pub = DP_M1_PUBLISHED[solvent][copolymer]
        diff = dP - pub
        rows.append((solvent, copolymer, kij, dP, pub, diff))
        print(f"{solvent:7s} {copolymer:10s} {kij:9.4f} "
              f"{dP:11.2f} {pub:11.2f} {diff:+7.2f}")
    # global (point-weighted) average vs published global
    print("-" * len(header))
    calc_global = np.mean([r[3] for r in rows])
    print(f"{'GLOBAL (unweighted mean of systems)':40s} "
          f"{calc_global:11.2f} {DP_M1_GLOBAL_PUBLISHED:11.2f}")
    print()
    _write_table_file(rows, calc_global)
    return rows


def _write_table_file(rows, calc_global):
    out = os.path.join(FIG_DIR, "tabla_comparativa_M1.txt")
    with open(out, "w") as f:
        f.write("# PC-SAFT Methodology 1 (pseudo-monomer): calculated vs "
                "published %dP/P\n")
        f.write("# Published values from Arce et al. (2006), Table 6 (PC-SAFT).\n")
        f.write(f"# {'solvent':7s} {'copolymer':10s} {'kij':>9s} "
                f"{'dP_calc':>9s} {'dP_Arce':>9s} {'diff':>7s}\n")
        for s, c, k, dP, pub, diff in rows:
            f.write(f"{s:9s} {c:10s} {k:9.4f} {dP:9.2f} {pub:9.2f} {diff:+7.2f}\n")
        f.write(f"# GLOBAL mean: calc={calc_global:.2f}  "
                f"Arce={DP_M1_GLOBAL_PUBLISHED:.2f}\n")
    print(f"  table -> {os.path.basename(out)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="PC-SAFT VLE reproduction (Arce 2006).")
    ap.add_argument("--validate", action="store_true",
                    help="run only validation steps 1 and 2")
    ap.add_argument("--table", action="store_true",
                    help="run only the comparison table (with plots)")
    ap.add_argument("--no-plots", action="store_true",
                    help="skip figure generation")
    ap.add_argument("--system", nargs=2, metavar=("SOLVENT", "COPOLYMER"),
                    help="run a single system, e.g. --system cC6 PSBR45")
    args = ap.parse_args()

    if args.system:
        solvent, copolymer = args.system
        temps = None
        for s, c, t, _ in SYSTEMS:
            if s == solvent and c == copolymer:
                temps = t
        if temps is None:
            print(f"Unknown system {solvent}+{copolymer}")
            sys.exit(1)
        dP = run_system(solvent, copolymer, temps, make_plot=not args.no_plots)
        pub = DP_M1_PUBLISHED[solvent][copolymer]
        print(f"\n{solvent}+{copolymer}: %dP/P calc = {dP:.2f}  "
              f"(Arce = {pub:.2f})")
        return

    if args.validate:
        validate_pure()
        validate_full()
        return

    if args.table:
        comparison_table(make_plots=not args.no_plots)
        return

    # default: everything
    validate_pure()
    validate_full()
    comparison_table(make_plots=not args.no_plots)


if __name__ == "__main__":
    main()
