#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grafica TODOS los sistemas (metodologia 1) en figuras temporales.
Reutiliza parametros y datos de correr_todos.py y el motor de pcsaft_vle.py.
Cada figura: curva PC-SAFT (linea) + puntos experimentales.
Las figuras se guardan en la carpeta ./figuras_tmp/
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import importlib.util
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# cargar parametros/datos de correr_todos.py (sin ejecutar su main)
src = open("correr_todos.py").read().split("if __name__")[0]
ns = {}
exec(compile(src, "correr_todos.py", "exec"), ns)
SOLV, COP, DATA, T6 = ns["SOLV"], ns["COP"], ns["DATA"], ns["KIJ_DP"]

# motor PC-SAFT
spec = importlib.util.spec_from_file_location("pv", "pcsaft_vle.py")
pv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figuras_tmp")
os.makedirs(OUT, exist_ok=True)

NOMBRES = {"ACE": "acetona", "nC4OH": "1-butanol", "CHL": "cloroformo",
           "cC6": "ciclohexano", "nC6": "n-hexano", "nC5": "n-pentano",
           "TOL": "tolueno", "BZN": "benceno"}


def curva(solv_k, cop_k, T, kij, npts=40):
    ms = SOLV[solv_k]; ps = COP[cop_k]
    pv.SOLVENTE = {"nombre": solv_k, "MW": ms[3], "m_MW": ms[0], "sigma": ms[1], "eps_k": ms[2]}
    pv.COPOLIMERO = {"nombre": cop_k, "MW": float(ps[3]), "m_MW": ps[0], "sigma": ps[1], "eps_k": ps[2]}
    solv = pv.construir_solvente(); cop = pv.construir_copolimero()
    ml = pv.PCSAFT([solv, cop], [[0, kij], [kij, 0]]); mv = pv.PCSAFT([solv], [[0.0]])
    # rango de w1 a partir de los datos experimentales del sistema
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
    print("Generando figuras en:", OUT)
    n = 0
    for (solv_k, cop_k, T), pts in DATA.items():
        kij, dp_pap = T6[(solv_k, cop_k)]
        try:
            wok, Pok, pts = curva(solv_k, cop_k, T, kij)
        except Exception as e:
            print("  ERR %s+%s: %s" % (solv_k, cop_k, e)); continue
        fig, ax = plt.subplots(figsize=(6.5, 4.6))
        if wok:
            ax.plot(wok, Pok, "-", color="#3b0a70", lw=2, label="PC-SAFT (metod. 1)")
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", color="#d62728",
                ms=6, markerfacecolor="white", markeredgewidth=1.3, label="Experimental")
        ax.set_xlabel("Fraccion masica de solvente  $w_1$")
        ax.set_ylabel("Presion  $P$  (MPa)")
        ax.set_title("%s (%s) + %s   T=%.2f K\nkij=%.4f  |  paper %%dP=%.2f%%"
                     % (solv_k, NOMBRES.get(solv_k, solv_k), cop_k, T, kij, dp_pap))
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
        fig.tight_layout()
        fname = os.path.join(OUT, "VLE_%s_%s_%.2f.png" % (solv_k, cop_k, T))
        fig.savefig(fname, dpi=120); plt.close(fig)
        n += 1
        print("  OK  %s+%s (%.2f K)" % (solv_k, cop_k, T))
    print("\nListo: %d figuras en %s" % (n, OUT))


if __name__ == "__main__":
    main()
