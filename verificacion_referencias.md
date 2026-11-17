# Verificación de las referencias de los datos experimentales

Se verificó que cada sistema use los datos de la referencia correcta, cotejando
la Tabla 2 del paper de Arce (que indica la fuente de cada sistema) contra los
datos cargados y contra los PDFs de las fuentes primarias disponibles en la rama.

## Referencias del paper y PDFs disponibles

| Ref | Fuente | PDF en repo | Sistemas que la usan |
|---|---|---|---|
| [1] | Gupta & Prausnitz, J. Chem. Eng. Data 40 (1995) 784 | `gupta1995.pdf` ✓ | casi todos |
| [33] | Wohlfarth, Int. Electron. J. Phys. Chem. Data 2 (1996) | — (no está) | cC6+PSBR41, TOL+PSBR41 |
| [34] | Wohlfarth, Macromol. Chem. Phys. 198 (1997) 2689 | `wohlfarth1997.pdf` ✓ | cC6+PSBR41, TOL+PSBR41 |
| [35] | Kang et al., J. Chem. Eng. Data 47 (2002) 788 | `kang2002.pdf` ✓ | BZN+PSBR4.1 |

PDFs confirmados por su título/autores:
- `gupta1995.pdf` = Gupta & Prausnitz 1995 (ref [1]). ✓
- `wohlfarth1997.pdf` = "VLE of concentrated poly(styrene-co-butadiene)...", Macromol.
  Chem. Phys. 198 (1997) = ref [34]. ✓
- `kang2002.pdf` = "VLE of Several Copolymer + Solvent Systems", Kang/Huang/Fu/Liu/Hu
  (benceno y tolueno; incluye poly(styrene-co-butadiene)) = ref [35]. ✓

## Resultado del cotejo sistema por sistema

### Sistemas cargados con referencia CORRECTA ([1] Gupta), verificados

Todos los sistemas de Gupta cargados coinciden con la Tabla 2 en **referencia,
rango de w1, rango de P y temperatura**. Verificados además celda por celda contra
el PDF (CHL+PVAE y TOL, en revisiones anteriores): datos fieles a la fuente.

ACE+PSBR5/23/45, nC4OH+PVAE9, nC4OH+PVAVC12, CHL+PVAE9/25/50/70, cC6+PSBR5,
cC6+PVAE25/50/70, TOL+PSBR5/23, nC6+PSBR23/45, nC5+PSBR23/45 → todos [1] Gupta. ✓

### Matiz 1: sistemas con DOBLE temperatura (datos parciales)

| Sistema | Tabla 2 | np Tabla 2 | np cargado | Nota |
|---|---|---|---|---|
| cC6+PSBR23 | 296.65–333.15 K | 15 | 8 | solo isoterma 333.15 K |
| cC6+PSBR45 | 296.65–333.15 K | 16 | 9 | solo isoterma 333.15 K |
| ACE+PSBR23 | 323.15–333.15 K | 22 | 9 | solo isoterma 333.15 K |
| ACE+PSBR45 | 323.15–333.15 K | 24 | 9 | solo isoterma 333.15 K |

Verificado en el PDF de Gupta: cC6+PSBR23 y PSBR45 TIENEN datos a 23.5 °C (296.65 K)
Y a 60 °C (333.15 K). La referencia es correcta (Gupta), pero solo cargamos la
isoterma de 333.15 K (coherente con la T del modelo). El paper promedia ambas
isotermas, por lo que su np es mayor y su %ΔP NO es directamente comparable con el
nuestro en estos 4 sistemas.

### Matiz 2: sistemas de Wohlfarth/Kang NO cargados

| Sistema | Ref | Estado |
|---|---|---|
| cC6+PSBR41 | [33,34] Wohlfarth | NO cargado |
| TOL+PSBR41 | [33,34] Wohlfarth | NO cargado |
| BZN+PSBR4.1 | [35] Kang | NO cargado |

Estos 3 sistemas (de las referencias [33,34] y [35], no de Gupta) no están en
nuestro `DATA`, así que no entran en los cálculos actuales. El paper sí los incluye
en su promedio global. Sus PDFs (`wohlfarth1997.pdf`, `kang2002.pdf`) están en la
rama y podrían cargarse para completar la reproducción.

## Conclusión

- **Las referencias asignadas son correctas.** Cada sistema cargado usa la fuente
  que indica la Tabla 2, confirmado por rangos de w1/P/T y por cotejo directo con
  los PDFs primarios.
- No se detectó ningún sistema usando datos de una fuente equivocada.
- Hallazgos menores: (1) los 4 sistemas de doble temperatura solo tienen cargada
  una isoterma (np parcial, %ΔP no comparable 1:1 con el paper); (2) faltan por
  cargar los 3 sistemas de Wohlfarth/Kang.
