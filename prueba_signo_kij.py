#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRUEBA DEFINITIVA: re-ejecutar los 8 sistemas con kij NEGATIVO (Tabla 6),
comparando el %DeltaP con:
  (a) el kij tal como esta en la Tabla 6 (negativo)
  (b) el kij con el signo INVERTIDO (hipotesis de errata de imprenta)
frente al %DeltaP reportado por el paper.

Si con el signo invertido todos caen a ~4-5%, confirma la errata de signo.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import importlib.util

# cargar parametros/datos de correr_todos.py
src = open("correr_todos.py").read().split("if __name__")[0]
ns = {}; exec(compile(src, "correr_todos.py", "exec"), ns)
SOLV, COP, DATA, T6 = ns["SOLV"], ns["COP"], ns["DATA"], ns["KIJ_DP"]

spec = importlib.util.spec_from_file_location("pv", "pcsaft_vle.py")
pv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)


def dp(solv_k, cop_k, T, pts, kij):
    ms = SOLV[solv_k]; ps = COP[cop_k]
    pv.SOLVENTE = {"nombre": solv_k, "MW": ms[3], "m_MW": ms[0], "sigma": ms[1], "eps_k": ms[2]}
    pv.COPOLIMERO = {"nombre": cop_k, "MW": float(ps[3]), "m_MW": ps[0], "sigma": ps[1], "eps_k": ps[2]}
    solv = pv.construir_solvente(); cop = pv.construir_copolimero()
    ml = pv.PCSAFT([solv, cop], [[0, kij], [kij, 0]]); mv = pv.PCSAFT([solv], [[0.0]])
    s = 0.0; n = 0
    for w1, Pe in pts:
        x1 = pv.w1_a_x1(w1)
        try:
            Pc = pv.presion_de_burbuja(ml, mv, T, x1)
        except Exception:
            Pc = None
        if Pc and Pc > 0:
            s += abs(Pc - Pe) / Pe; n += 1
    return (100 * s / n) if n else None


if __name__ == "__main__":
    # sistemas con kij negativo en Tabla 6
    neg = [(s, c, T) for (s, c, T) in DATA if T6[(s, c)][0] < 0]
    print("PRUEBA DE ERRATA DE SIGNO EN kij (sistemas con kij negativo en Tabla 6)")
    print("=" * 74)
    print("%-6s %-9s %9s | %10s %10s | %9s" % (
        "Solv", "Copol", "kij_T6", "dP(kij_T6)", "dP(-kij)", "dP_paper"))
    print("-" * 74)
    mejoras = 0
    for (s, c, T) in neg:
        kij_t6, dp_pap = T6[(s, c)]
        pts = DATA[(s, c, T)]
        d_neg = dp(s, c, T, pts, kij_t6)      # tal cual (negativo)
        d_pos = dp(s, c, T, pts, -kij_t6)     # signo invertido
        mejor = "<-- +signo mejor" if (d_pos is not None and d_neg is not None and d_pos < d_neg) else ""
        if d_pos is not None and d_pos < d_neg:
            mejoras += 1
        print("%-6s %-9s %+9.4f | %10.1f %10.1f | %9.2f  %s" % (
            s, c, kij_t6, d_neg, d_pos, dp_pap, mejor))
    print("-" * 74)
    print("Sistemas donde invertir el signo MEJORA: %d de %d" % (mejoras, len(neg)))
