#!/usr/bin/env python3
"""
Diagnostic: for each system, find the kij that MINIMISES %dP/P and compare it to
the published Table-6 kij. This does NOT change the reproduction (which always
uses the published kij); it only characterises WHY the published kij do not
reproduce the paper's deviations.

Run:  python diagnostico_kij.py            (all systems)
      python diagnostico_kij.py cC6 PSBR45 (single system)
"""
import sys
import numpy as np
import vle
import pcsaft_vle as M
from parametros_articulo import SYSTEMS, KIJ_M1, DP_M1_PUBLISHED


def dev_system(kij, solvent, copolymer, temps):
    tot, n = 0.0, 0
    for T in temps:
        w1, Pexp = M.load_data(solvent, copolymer, T)
        if w1 is None:
            continue
        d, _ = vle.deviation_percent(w1, Pexp, T, solvent, copolymer, kij)
        tot += d * len(w1)
        n += len(w1)
    return tot / n


def best_kij(solvent, copolymer, temps):
    """Single coarse grid search (diagnostic resolution ~0.025 in kij)."""
    grid = np.linspace(-0.15, 0.15, 13)
    devs = [dev_system(k, solvent, copolymer, temps) for k in grid]
    kb = grid[int(np.argmin(devs))]
    return kb, min(devs)


def run(systems):
    print(f"{'solvent':7s} {'copolymer':10s} {'kij_pub':>8s} {'dP_pub':>7s} "
          f"{'kij_opt':>8s} {'dP_opt':>7s} {'dP_Arce':>8s}")
    print("-" * 60)
    for solvent, copolymer, temps, _ in systems:
        kpub = KIJ_M1[solvent][copolymer]
        dpub = dev_system(kpub, solvent, copolymer, temps)
        kopt, dopt = best_kij(solvent, copolymer, temps)
        print(f"{solvent:7s} {copolymer:10s} {kpub:8.4f} {dpub:7.2f} "
              f"{kopt:8.4f} {dopt:7.2f} {DP_M1_PUBLISHED[solvent][copolymer]:8.2f}")


if __name__ == "__main__":
    if len(sys.argv) == 3:
        s, c = sys.argv[1], sys.argv[2]
        systems = [x for x in SYSTEMS if x[0] == s and x[1] == c]
    else:
        systems = SYSTEMS
    run(systems)
