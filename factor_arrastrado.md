# El error se arrastra desde M1: es el término de energía, no la metodología

Hipótesis del usuario (confirmada): si M2 no corrige el error de M1 de forma
consistente, el problema NO está en cómo se combinan los segmentos (M1 vs M2),
sino en un **factor común a ambas**. Este documento lo identifica.

## 1. La diferencia M1−M2 NO coincide entre cálculo y paper

| | D_calc = M1c−M2c | D_paper = M1p−M2p |
|---|---|---|
| Promedio | +15.5 | +2.4 |
| Rango | −8.8 a +54.7 | +0.2 a +4.8 |

En el paper, pasar de M1 a M2 mejora ~2 puntos (modesto, porque M1 ya daba ~4%).
En nuestro cálculo, M2 cambia el resultado de forma enorme y errática. Si
capturáramos la física, D_calc ≈ D_paper. No lo hace → **el error viene de antes**,
de algo que M1 y M2 comparten.

## 2. El factor común: el ε (energía) del copolímero en la dispersión

Experimento (M1, cC6+PVAE25, error base 83.7%):

| Cambio | %ΔP |
|---|---|
| ε_copol = 202.57 (Tabla 5) | 83.7 |
| ε_copol = 140 | 67.6 |
| ε_copol = 100 | **7.1** |
| σ_copol variado 3.0–3.6 | ~84 (sin efecto) |

Control cC6+PSBR45 (error base 2.2%): subir ε de 246 a 300 lo lleva a 76%.

→ El %ΔP es **extremadamente sensible al ε del copolímero** (σ casi no importa).
  Los PVAE necesitan un ε mucho menor del que da la Tabla 5.

## 3. El ε de Tabla 5 está bien calculado, pero es "demasiado alto" para el motor

El ε del pseudo-monómero de Tabla 5 = media aritmética por fracción molar de los
monómeros (verificado exacto). El del PVAE sale alto porque el **monómero VA tiene
ε=421 K** (anómalamente alto; el resto de monómeros/solventes están en 160–320 K).

## 4. El kij negativo AGRAVA el problema (no lo causa)

El ε cruzado solvente-copolímero es ε_cruz = √(ε_s·ε_c)·(1−kij). Con kij<0,
(1−kij)>1 → ε_cruz aún mayor → copolímero más afín → presión más baja.

Prueba (M1), sistemas de kij negativo (PVAE) con kij del paper vs kij=0:

| Sistema | kij paper | %ΔP (kij paper) | %ΔP (kij=0) |
|---|---|---|---|
| CHL+PVAE50 | −0.0701 | 58.5 | **6.0** |
| CHL+PVAE25 | −0.0816 | 72.0 | **18.1** |
| CHL+PVAE70 | −0.0516 | 53.2 | **10.5** |
| cC6+PVAE25 | −0.1068 | 83.7 | **29.0** |
| CHL+PVAE9 | −0.1208 | 93.7 | **58.7** |

Los sistemas de kij positivo (PSBR, buenos) empeoran con kij=0, como debe ser.

- Invertir la convención (1+kij): empeora todo (promedio 77%). Descartado: la
  convención estándar (1−kij) es la correcta.
- Con kij=0 los PVAE ya mejoran mucho, pero no llegan a ~5%: el ε base del
  copolímero (202 K) sigue siendo demasiado alto incluso sin kij.

## 5. Conclusión

El usuario tiene razón: **el error se arrastra desde M1 y no es la forma de calcular
la combinación de segmentos (M1↔M2).** El factor común y dominante es el **término
de energía de dispersión del copolímero**: el ε del pseudo-monómero/segmentos (que
hereda el ε anómalo del VA = 421 K) produce, en nuestro motor, un ε cruzado
solvente-copolímero demasiado alto, y los kij negativos del paper lo agravan.

El motor y la regla de combinación están verificados; el ε de Tabla 5 está bien
calculado según la regla del paper. Por tanto, la pieza que falla es **cómo el ε
del copolímero (y su ε cruzado con el solvente) se usa en la dispersión** — una
diferencia respecto al paper que no se recupera de las tablas/ecuaciones publicadas,
y que afecta por igual a M1 y M2.
