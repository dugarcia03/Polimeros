# Tabla de errores %ΔP: Metodología 1 y 2 (calculado) vs Paper

Comparación sistema por sistema del error relativo de presión (%ΔP, Ec. 40)
entre nuestros cálculos y los valores reportados por Arce et al. (2006).

- **M1** = metodología 1 (pseudo-monómero, 1 segmento). Paper = Tabla 6.
- **M2** = metodología 2 (copolímero como 2 segmentos α/β). Paper = Tabla 7 (por familia).
- M2 solo aplica a los 6 copolímeros con fracciones de enlace en Tabla 3
  (PSBR23/45, PVAE9/25/50/70); por eso PSBR5, PVAVC12 y familias sin Tabla 3
  aparecen con "–" en M2.

| Solvente | Copolímero | M1 calc | M1 paper | M2 calc | M2 paper | Mejor M2 vs M1 |
|---|---|---:|---:|---:|---:|:---:|
| ACE | PSBR5 | 22.32 | 4.78 | – | 2.83 | – |
| ACE | PSBR23 | 25.01 | 3.04 | 23.15 | 2.83 | ✓ |
| ACE | PSBR45 | 33.40 | 5.19 | 30.43 | 2.83 | ✓ |
| nC4OH | PVAE9 | 10.81 | 4.05 | 17.30 | 1.88 | ✗ |
| nC4OH | PVAVC12 | 82.10 | 5.18 | – | – | – |
| CHL | PVAE9 | 93.74 | 4.18 | 58.62 | 2.59 | ✓ |
| CHL | PVAE25 | 71.95 | 5.06 | 19.45 | 2.59 | ✓ |
| CHL | PVAE50 | 58.48 | 4.41 | 7.37 | 2.59 | ✓ |
| CHL | PVAE70 | 53.21 | 3.28 | 24.35 | 2.59 | ✓ |
| cC6 | PSBR5 | 6.02 | 3.66 | – | 2.11 | – |
| cC6 | PSBR23 | 11.99 | 4.38 | 20.75 | 2.11 | ✗ |
| cC6 | PSBR45 | 2.17 | 4.31 | 5.46 | 2.11 | ✗ |
| cC6 | PVAE25 | 83.73 | 5.36 | 29.04 | 2.09 | ✓ |
| cC6 | PVAE50 | 85.26 | 6.90 | 55.22 | 2.09 | ✓ |
| cC6 | PVAE70 | 81.28 | 5.73 | 70.13 | 2.09 | ✓ |
| TOL | PSBR5 | 7.28 | 4.07 | – | 2.13 | – |
| TOL | PSBR23 | 68.84 | 5.76 | 64.38 | 2.13 | ✓ |
| nC6 | PSBR23 | 19.04 | 3.92 | 18.92 | 1.84 | ✓ |
| nC6 | PSBR45 | 19.01 | 4.28 | 24.40 | 1.84 | ✗ |
| nC5 | PSBR23 | 38.23 | 3.68 | 30.16 | 2.07 | ✓ |
| nC5 | PSBR45 | 56.69 | 6.24 | 49.77 | 2.07 | ✓ |

## Promedios

| Metodología | Calc (nuestro) | Paper | n sistemas |
|---|---:|---:|---:|
| M1 | 44.31 % | 4.64 % | 21 |
| M2 | 32.29 % | 2.26 % | 17 |

Global reportado por el paper (sobre sus 24 sistemas): **M1 = 4.41 %, M2 = 2.04 %**.

(Los promedios "paper" de la tabla difieren un poco del global del paper porque se
calculan solo sobre los sistemas que tenemos cargados, y el paper promedia sobre 24
sistemas, incluidos los de Wohlfarth/Kang que no están cargados y las segundas
isotermas de los sistemas de doble temperatura.)

## Observaciones

- **M2 mejora respecto a M1 en 11 de 15 sistemas** donde aplica, especialmente en
  los PVAE: CHL+PVAE50 (58→7 %), CHL+PVAE25 (72→19 %), cC6+PVAE25 (84→29 %).
- M2 empeora en 4 sistemas (nC4OH+PVAE9, cC6+PSBR23, cC6+PSBR45, nC6+PSBR45),
  todos de bajo error en M1.
- Ambas metodologías quedan por encima del error del paper. El motor y las
  ecuaciones están verificados (ver `hallazgos_*.md`); la brecha se atribuye a una
  decisión metodológica no documentada del paper sobre la ponderación de segmentos.
