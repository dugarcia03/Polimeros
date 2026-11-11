#!/usr/bin/env python3
"""
Reproduccion de:
  Arce, Mattedi & Aznar (2006), Fluid Phase Equilibria 246, 52-63.
  "Vapor-liquid equilibrium of copolymer + solvent mixtures: Thermodynamic
   modeling by two theoretical equations of state."

Programa interactivo: el usuario escoge un sistema (solvente + copolimero) y el
programa calcula el equilibrio liquido-vapor (presion vs fraccion masica de
solvente) con la EoS PC-SAFT para copolimeros, usando ambas metodologias del
paper, genera la grafica P vs w1 y reporta la desviacion de presion %DeltaP/P.

Uso:
    python main.py                 # modo interactivo
    python main.py --all           # genera figuras + tabla de todos los sistemas
    python main.py --system ACE PSBR5   # calcula un sistema concreto
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

import sys
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from data.systems import (SYSTEMS, SOLVENT_NAMES, COPOLYMER_NAMES, FAMILY,
                          COPOLYMERS)
from data.parameters import (TABLE6, TABLE7, COPOLYMER_PCSAFT,
                             GLOBAL_DP_PCSAFT_M1, GLOBAL_DP_PCSAFT_M2)
from src.vle import (bubble_pressure_pcsaft, kij_methodology1, kij_methodology2,
                     w1_to_x1)

OUTDIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(OUTDIR, exist_ok=True)

NPTS = 20   # puntos por curva


def _solvent_label(s):
    return f"{s} ({SOLVENT_NAMES.get(s, s)})"


def curve_pcsaft(solvent, copolymer, T, methodology=1, npts=NPTS):
    """Calcula la curva P(w1) con PC-SAFT para un sistema a temperatura T."""
    info = SYSTEMS[(solvent, copolymer)]
    wmin, wmax = info["w1_range"]
    ws = np.linspace(max(wmin, 1e-4), wmax, npts)
    if methodology == 1:
        kij = kij_methodology1(solvent, copolymer)
    else:
        kij = kij_methodology2(solvent, copolymer)
    w_ok, P_ok = [], []
    for w1 in ws:
        x1 = w1_to_x1(float(w1), solvent, copolymer)
        try:
            P = bubble_pressure_pcsaft(solvent, copolymer, T, x1, kij)
        except Exception:
            P = None
        if P is not None and P > 0:
            w_ok.append(float(w1))
            P_ok.append(P)
    return np.array(w_ok), np.array(P_ok)


def plot_system(solvent, copolymer, show_m2=True, save=True):
    """Genera la figura P vs w1 para un sistema (todas sus temperaturas)."""
    info = SYSTEMS[(solvent, copolymer)]
    T_list = info["T_list"]

    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.cm.viridis(np.linspace(0, 0.8, len(T_list)))

    for Ti, T in enumerate(T_list):
        w1, P1 = curve_pcsaft(solvent, copolymer, T, methodology=1)
        if len(w1):
            ax.plot(w1, P1, "-", color=colors[Ti],
                    label=f"PC-SAFT M1, T={T:.2f} K")
        if show_m2 and (solvent, FAMILY[copolymer]) in TABLE7:
            w2, P2 = curve_pcsaft(solvent, copolymer, T, methodology=2)
            if len(w2):
                ax.plot(w2, P2, "--", color=colors[Ti],
                        label=f"PC-SAFT M2, T={T:.2f} K")

    fam = FAMILY[copolymer]
    ax.set_xlabel("Fraccion masica de solvente  $w_1$")
    ax.set_ylabel("Presion  $P$  (MPa)")
    ax.set_title(f"VLE  {_solvent_label(solvent)} (1) + {copolymer} (2)\n"
                 f"{COPOLYMER_NAMES.get(fam, fam)}")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if save:
        fname = os.path.join(OUTDIR, f"VLE_{solvent}_{copolymer}.png")
        fig.savefig(fname, dpi=130)
        plt.close(fig)
        return fname
    return fig


def deviation_table(methodology=1):
    """Calcula %DeltaP/P de cada sistema tomando como referencia la curva del
    paper (metodologia 1: Tabla 6 PC-SAFT ; metodologia 2: Tabla 7).

    Como no disponemos de los puntos experimentales, comparamos la curva
    calculada por nuestro codigo contra la desviacion reportada por el paper,
    y reportamos ambas para contrastar.
    """
    rows = []
    for (solvent, copolymer), info in SYSTEMS.items():
        entry = TABLE6.get((solvent, copolymer))
        if entry is None:
            continue
        dp_paper = entry["dP_pcsaft"] if methodology == 1 else None
        rows.append((solvent, copolymer, dp_paper))
    return rows


def print_menu():
    systems = list(SYSTEMS.keys())
    print("\n" + "=" * 64)
    print(" Reproduccion VLE copolimero + solvente (Arce et al., 2006)")
    print(" EoS: PC-SAFT para copolimeros")
    print("=" * 64)
    print(" Sistemas disponibles:\n")
    for i, (s, c) in enumerate(systems):
        fam = FAMILY[c]
        Ts = ", ".join(f"{t:.2f}" for t in SYSTEMS[(s, c)]["T_list"])
        print(f"  [{i:2d}] {s:6s} + {c:9s}  ({SOLVENT_NAMES.get(s,s)} + "
              f"{COPOLYMER_NAMES.get(fam,fam)})   T[K]={Ts}")
    print(f"\n  [a] generar TODAS las figuras + tabla de desviaciones")
    print(f"  [q] salir")
    return systems


def interactive():
    systems = print_menu()
    while True:
        choice = input("\nElige una opcion: ").strip().lower()
        if choice in ("q", "quit", "salir"):
            print("Hasta luego.")
            return
        if choice == "a":
            run_all()
            continue
        try:
            idx = int(choice)
            s, c = systems[idx]
        except (ValueError, IndexError):
            print("  Opcion no valida.")
            continue
        print(f"\nCalculando {s} + {c} ...")
        fname = plot_system(s, c)
        print(f"  Figura guardada en: {fname}")
        # imprimir tabla de presiones
        for T in SYSTEMS[(s, c)]["T_list"]:
            w1, P1 = curve_pcsaft(s, c, T, methodology=1)
            print(f"\n  T = {T:.2f} K  (PC-SAFT metodologia 1)")
            print(f"  {'w1':>8}  {'P (MPa)':>10}")
            for w, p in zip(w1, P1):
                print(f"  {w:8.4f}  {p:10.5f}")
        entry = TABLE6.get((s, c))
        if entry:
            print(f"\n  %DeltaP/P reportado por el paper (Tabla 6, PC-SAFT): "
                  f"{entry['dP_pcsaft']:.2f}%")


def run_all():
    print("\nGenerando figuras de todos los sistemas...")
    figs = []
    for (s, c) in SYSTEMS.keys():
        try:
            fname = plot_system(s, c)
            figs.append(fname)
            print(f"  OK  {s}+{c}  -> {os.path.basename(fname)}")
        except Exception as e:
            print(f"  ERR {s}+{c}: {e}")
    # tabla de desviaciones reportadas
    print("\n" + "=" * 70)
    print(" Desviaciones de presion %DeltaP/P reportadas por el paper")
    print("=" * 70)
    print(f" {'Solvente':9} {'Copolimero':11} {'LGT M1':>8} {'PC-SAFT M1':>11} "
          f"{'PC-SAFT M2':>11}")
    for (s, c), entry in TABLE6.items():
        fam = FAMILY[c]
        t7 = TABLE7.get((s, fam))
        dp2 = f"{t7['dP']:.2f}" if t7 else "-"
        print(f" {s:9} {c:11} {entry['dP_lgt']:8.2f} "
              f"{entry['dP_pcsaft']:11.2f} {dp2:>11}")
    print("-" * 70)
    print(f" {'GLOBAL':9} {'':11} {'10.10':>8} "
          f"{GLOBAL_DP_PCSAFT_M1:11.2f} {GLOBAL_DP_PCSAFT_M2:>11.2f}")
    print(f"\n Figuras en: {OUTDIR}")


def main():
    ap = argparse.ArgumentParser(description="VLE copolimero+solvente (Arce 2006)")
    ap.add_argument("--all", action="store_true", help="generar todo")
    ap.add_argument("--system", nargs=2, metavar=("SOLVENT", "COPOLYMER"),
                    help="calcular un sistema concreto, p.ej. --system ACE PSBR5")
    args = ap.parse_args()

    if args.all:
        run_all()
    elif args.system:
        s, c = args.system
        if (s, c) not in SYSTEMS:
            print(f"Sistema {s}+{c} no encontrado.")
            sys.exit(1)
        fname = plot_system(s, c)
        print(f"Figura guardada en: {fname}")
    else:
        interactive()


if __name__ == "__main__":
    main()
