# -*- coding: utf-8 -*-
"""
Tabla comparativa de los TRES metodos sobre los mismos sistemas:
  M1         : pseudo-monomero (1 segmento)           -- del paper
  M2         : copolimero de 2 segmentos (alpha/beta) -- del paper
  M1 + polar : M1 mas termino dipolar Gross-Vrabec    -- EXTENSION (no del paper)

Genera stdout + imagen PNG (tabla_tres_metodos.png).
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import metodologia_polar as M


def calcular():
    filas = []
    s1 = n1 = s2 = n2 = sp = npp = 0.0
    for (solv, cop, T), pts in M.DATA.items():
        if (solv, cop) not in M.KIJ_M1:
            continue
        # M1 (sin polar)
        M.USAR_POLAR = False
        d1, _ = M.dp_system(M.modelos_m1(solv, cop), solv, T, pts)
        # M2
        mm2 = M.modelos_m2(solv, cop)
        d2, _ = M.dp_system(mm2, solv, T, pts) if mm2 else (None, 0)
        if cop not in M.COP_M2:
            d2 = None
        # M1 + polar
        M.USAR_POLAR = True
        dp, _ = M.dp_system(M.modelos_m1(solv, cop), solv, T, pts)
        filas.append((solv, cop, T, d1, d2, dp))
        if d1 is not None: s1 += d1; n1 += 1
        if d2 is not None: s2 += d2; n2 += 1
        if dp is not None: sp += dp; npp += 1
    g = (s1 / n1, s2 / n2, sp / npp, int(n1), int(n2), int(npp))
    return filas, g


def imprimir(filas, g):
    print("%-6s %-9s %7s | %8s | %8s | %10s" %
          ("Solv", "Copol", "T[K]", "M1", "M2", "M1+polar"))
    print("-" * 60)
    for (solv, cop, T, d1, d2, dp) in filas:
        f = lambda d: ("%8.2f" % d) if d is not None else "    None"
        print("%-6s %-9s %7.2f | %s | %s | %10s" %
              (solv, cop, T, f(d1), f(d2), ("%.2f" % dp) if dp is not None else "None"))
    print("-" * 60)
    print("GLOBAL (prom. por sistema):")
    print("  M1         : %.2f%% (%d sist)  [paper 4.41%%]" % (g[0], g[3]))
    print("  M2         : %.2f%% (%d sist)  [paper 2.04%%]" % (g[1], g[4]))
    print("  M1 + polar : %.2f%% (%d sist)  [extension]" % (g[2], g[5]))


def imagen(filas, g, out="tabla_tres_metodos.png"):
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), out)
    enc = ["Solvente", "Copolimero", "T [K]", "M1  %", "M2  %", "M1+polar  %"]
    celdas = []
    for (solv, cop, T, d1, d2, dp) in filas:
        celdas.append([solv, cop, "%.2f" % T,
                       "%.2f" % d1 if d1 is not None else "None",
                       "%.2f" % d2 if d2 is not None else "None",
                       "%.2f" % dp if dp is not None else "None"])
    celdas.append(["GLOBAL", "(prom.)", "",
                   "%.2f" % g[0], "%.2f" % g[1], "%.2f" % g[2]])
    celdas.append(["paper", "", "", "4.41", "2.04", "(no aplica)"])

    fig_h = 0.42 * (len(celdas) + 2)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    ax.axis("off")
    ax.set_title("Desviacion media en presion %dP (Ec. 40)\n"
                 "Comparacion de tres metodos - Arce et al. (2006) + extension polar",
                 fontsize=13, fontweight="bold", pad=14)
    tabla = ax.table(cellText=celdas, colLabels=enc, cellLoc="center", loc="center")
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(9)
    tabla.scale(1.0, 1.35)
    for j in range(len(enc)):
        c = tabla[0, j]
        c.set_facecolor("#2f5d62"); c.set_text_props(color="white", fontweight="bold")
    # resaltar las dos filas resumen
    for rr in (len(celdas) - 1, len(celdas)):
        for j in range(len(enc)):
            c = tabla[rr, j]
            c.set_facecolor("#dbe9e9"); c.set_text_props(fontweight="bold")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Imagen guardada en:", outpath)


if __name__ == "__main__":
    filas, g = calcular()
    imprimir(filas, g)
    imagen(filas, g)
