#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grafica UNA figura por SOLVENTE (como las Figs. 2-10 del paper): en cada figura
aparecen todas las composiciones de copolimero de ese solvente, cada una con su
curva PC-SAFT (linea) y sus puntos experimentales (simbolos).

Cuando un solvente tiene dos familias de copolimero a T distinta (p.ej. cC6 con
PSBR a 333.15 K y PVAE a 353.15 K), se generan figuras separadas por familia,
igual que hace el paper.

Reutiliza parametros/datos de correr_todos.py y el motor de pcsaft_vle.py.
Figuras en ./figuras_por_solvente/
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import importlib.util
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# cargar parametros/datos
src = open("correr_todos.py").read().split("if __name__")[0]
ns = {}
exec(compile(src, "correr_todos.py", "exec"), ns)
SOLV, COP, DATA, T6 = ns["SOLV"], ns["COP"], ns["DATA"], ns["KIJ_DP"]

# motor PC-SAFT
spec = importlib.util.spec_from_file_location("pv", "pcsaft_vle.py")
pv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras_por_solvente")
os.makedirs(OUT, exist_ok=True)

NOMBRES = {"ACE": "acetona", "nC4OH": "1-butanol", "CHL": "cloroformo",
           "cC6": "ciclohexano", "nC6": "n-hexano", "nC5": "n-pentano",
           "TOL": "tolueno", "BZN": "benceno"}
# familia del copolimero (para separar figuras cuando un solvente tiene 2)
def familia(cop):
    if cop.startswith("PSBR"): return "PSBR"
    if cop.startswith("PVAE"): return "PVAE"
    if cop.startswith("PVAVC"): return "PVAVC"
    return cop

MARKERS = ["o", "s", "^", "v", "D", "P", "*", "X"]
COLORS = plt.cm.tab10.colors


def curva(solv_k, cop_k, T, kij, npts=35):
    ms = SOLV[solv_k]; ps = COP[cop_k]
    pv.SOLVENTE = {"nombre": solv_k, "MW": ms[3], "m_MW": ms[0], "sigma": ms[1], "eps_k": ms[2]}
    pv.COPOLIMERO = {"nombre": cop_k, "MW": float(ps[3]), "m_MW": ps[0], "sigma": ps[1], "eps_k": ps[2]}
    # nC4OH es auto-asociante (Tabla 4b): activa el termino de asociacion (Ecs.36-38)
    if solv_k == "nC4OH":
        solv = pv.Componente(solv_k, ms[0]*ms[3], ms[1], ms[2], 0.008752, 2489.23)
    else:
        solv = pv.construir_solvente()
    cop = pv.construir_copolimero()
    ml = pv.PCSAFT([solv, cop], [[0, kij], [kij, 0]]); mv = pv.PCSAFT([solv], [[0.0]])
    pts = DATA[(solv_k, cop_k, T)]
    wmin = min(p[0] for p in pts); wmax = max(p[0] for p in pts)
    ws = np.linspace(max(wmin * 0.5, 1e-4), wmax * 1.05, npts)
    wok, Pok = [], []
    for w1 in ws:
        x1 = pv.w1_a_x1(float(w1))
        try:
            P = pv.presion_de_burbuja(ml, mv, T, x1)
        except Exception:
            P = None
        if P and P > 0:
            wok.append(float(w1)); Pok.append(P)
    return wok, Pok, pts


def main():
    # agrupar sistemas por (solvente, familia)
    grupos = defaultdict(list)
    for (s, c, T) in DATA:
        grupos[(s, familia(c))].append((c, T))

    print("Generando figuras en:", OUT)
    for (solv_k, fam), sistemas in grupos.items():
        sistemas.sort()
        fig, ax = plt.subplots(figsize=(7, 5))
        for idx, (cop_k, T) in enumerate(sistemas):
            kij, dp_pap = T6[(solv_k, cop_k)]
            col = COLORS[idx % len(COLORS)]
            mk = MARKERS[idx % len(MARKERS)]
            try:
                wok, Pok, pts = curva(solv_k, cop_k, T, kij)
            except Exception as e:
                print("  ERR %s+%s: %s" % (solv_k, cop_k, e)); continue
            # curva calculada
            if wok:
                ax.plot(wok, Pok, "-", color=col, lw=1.8)
            # puntos experimentales
            ax.plot([p[0] for p in pts], [p[1] for p in pts], mk, color=col,
                    ms=6, markerfacecolor="white", markeredgewidth=1.3,
                    label="%s (T=%.2f K)" % (cop_k, T))
        ax.set_xlabel("Fraccion masica de solvente  $w_1$")
        ax.set_ylabel("Presion  $P$  (MPa)")
        nom = NOMBRES.get(solv_k, solv_k)
        ax.set_title("VLE PC-SAFT (metod. 1)\n%s (%s) + %s\n(linea = PC-SAFT, simbolos = experimental)"
                     % (solv_k, nom, fam))
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, title="Copolimero")
        fig.tight_layout()
        fname = os.path.join(OUT, "VLE_%s_%s.png" % (solv_k, fam))
        fig.savefig(fname, dpi=130); plt.close(fig)
        print("  OK  %s + %s  (%d copolimeros)" % (solv_k, fam, len(sistemas)))

    n = len(os.listdir(OUT))
    print("\nListo: %d figuras (una por solvente/familia) en %s" % (n, OUT))


if __name__ == "__main__":
    main()
