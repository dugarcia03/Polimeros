# Comparativa Metodología 1 — mi cálculo vs paper

Reproducción PC-SAFT (metodología 1, pseudo-monómero) del artículo:

> Arce, Mattedi & Aznar (2006), *Fluid Phase Equilibria* 246, 52–63.
> DOI: 10.1016/j.fluid.2006.05.009

Desviación de presión %ΔP/P calculada según la Ec. 40 del paper, comparada con
los valores reportados en la Tabla 6 (columna PC-SAFT).

| Solvente | Copolímero | T [K] | kij | mi %ΔP | paper %ΔP | diferencia |
|----------|------------|-------|---------|--------|-----------|------------|
| ACE      | PSBR5      | 323.15 | 0.0120  | 22.32  | 4.78  | ↑ |
| ACE      | PSBR23     | 333.15 | 0.0122  | 25.01  | 3.04  | ↑ |
| ACE      | PSBR45     | 333.15 | 0.0130  | 33.40  | 5.19  | ↑ |
| nC4OH    | PVAE9      | 353.15 | 0.0230  | 25.35  | 4.05  | ↑ |
| nC4OH    | PVAVC12    | 353.15 | −0.0126 | 69.40  | 5.18  | ↑↑ |
| CHL      | PVAE9      | 333.15 | −0.1208 | 93.74  | 4.18  | ↑↑ |
| CHL      | PVAE25     | 333.15 | −0.0816 | 71.95  | 5.06  | ↑↑ |
| CHL      | PVAE50     | 333.15 | −0.0701 | 58.48  | 4.41  | ↑↑ |
| CHL      | PVAE70     | 333.15 | −0.0516 | 53.21  | 3.28  | ↑↑ |
| cC6      | PSBR45     | 333.15 | 0.0101  | 2.17   | 4.31  | ✓ mejor |
| cC6      | PSBR5      | 333.15 | 0.0118  | 6.02   | 3.66  | ✓ cerca |
| cC6      | PSBR23     | 333.15 | 0.0104  | 11.99  | 4.38  | ~ |
| cC6      | PVAE25     | 353.15 | −0.1068 | 83.73  | 5.36  | ↑↑ |
| cC6      | PVAE50     | 353.15 | −0.0858 | 85.26  | 6.90  | ↑↑ |
| cC6      | PVAE70     | 353.15 | −0.0636 | 81.28  | 5.73  | ↑↑ |
| TOL      | PSBR5      | 308.15 | 0.0080  | 7.28   | 4.07  | ✓ cerca |
| TOL      | PSBR23     | 308.15 | 0.0091  | 68.84  | 5.76  | ↑↑ |
| nC6      | PSBR23     | 343.15 | 0.0086  | 19.04  | 3.92  | ↑ |
| nC6      | PSBR45     | 343.15 | 0.0089  | 19.01  | 4.28  | ↑ |
| nC5      | PSBR23     | 333.15 | 0.0072  | 38.23  | 3.68  | ↑ |
| nC5      | PSBR45     | 333.15 | 0.0079  | 56.69  | 6.24  | ↑↑ |
| **GLOBAL** |          |        |         | —      | **4.41** | |

**Leyenda:** ✓ = reproduce bien (cerca o mejor que el paper) · ~ = aceptable ·
↑ = mi error mayor · ↑↑ = mi error mucho mayor.

## Notas

- **Reproducen bien** (solventes no polares + PSBR): cC6+PSBR45 (2.17 %, mejor que
  el paper), cC6+PSBR5, TOL+PSBR5.
- **Reproducen mal**: solventes polares/asociantes (ACE, nC4OH) y todos los
  copolímeros PVAE.
- Causas identificadas en el diagnóstico:
  1. Acetona y 1-butanol son polares/asociantes → PC-SAFT no polar reproduce mal
     su presión de vapor.
  2. Los kij negativos grandes de la Tabla 6 (PVAE/PVAVC con CHL y cC6) parecen
     tener posible errata de signo: con el kij de signo corregido, varios de esos
     sistemas bajan a ~4–5 %.
- **GLOBAL 4.41 %** es la desviación global PC-SAFT (metodología 1) reportada por
  el paper en la Tabla 6.
