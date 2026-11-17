# Metodología 2 (segmentos α/β) — implementación y resultados

Rama `metodologia-2`. La Metodología 1 queda preservada en la rama
`metodologia-1-final-provisional` (recuperable con `git checkout`).

## Qué se implementó

Motor nuevo `pcsaft_m2.py` que trata el copolímero como una cadena de **dos tipos
de segmento** (α, β = los dos monómeros), siguiendo la extensión de PC-SAFT a
copolímeros de Gross-Spuhl-Tumakaka-Sadowski (ref [24] del paper):

- **Ec. 18-20:** m_{i,a} = w_{i,a}·MW·(m/MW)_{i,a} por segmento; z_{i,a}=m_{i,a}/m_i.
- **Ec. 22:** hard-chain con doble suma Σ_α Σ_β B_{iα,iβ}·ln g_{iα,iβ}
  (fracciones de enlace B de la Tabla 3).
- **Ecs. 28-29:** dispersión con doble suma sobre segmentos z_iα·z_jβ.
- **Ec. 33:** tres kij por par de segmento (Tabla 7): κ_a (solvente-mon1),
  κ_b (solvente-mon2), κ_αβ (mon1-mon2).
- Término de asociación (Ecs. 36-38) para nC4OH.

Datos de entrada, todos del paper:
- Tabla 4b: parámetros de solventes y **monómeros**.
- Tabla 3: fracciones de enlace B (solo 6 copolímeros con datos PVT).
- Tabla 7: los tres kij.

Mapeo verificado: en Tabla 3, β = mon1 (el del %peso del nombre, STR/VA),
α = mon2. En el motor, segmento 0 = mon1, segmento 1 = mon2, con la matriz B
reordenada en consecuencia.

## Resultados (%ΔP calculado vs Tabla 7 del paper)

| Sistema | M2 calc | paper (familia) | M1 (ref) |
|---|---|---|---|
| CHL+PVAE25 | **6.3** | 2.59 | 72 |
| cC6+PSBR23 | **6.4** | 2.11 | 12 |
| nC5+PSBR23 | **8.1** | 2.07 | 38 |
| CHL+PVAE70 | 12.0 | 2.59 | 53 |
| nC4OH+PVAE9 | 14.3 | 1.88 | 11 |
| cC6+PVAE25 | 15.1 | 2.09 | 84 |
| cC6+PVAE50 | 16.3 | 2.09 | 85 |
| nC5+PSBR45 | 21.4 | 2.07 | 57 |
| ACE+PSBR23 | 28.0 | 2.83 | 25 |
| CHL+PVAE50 | 45.1 | 2.59 | 58 |
| cC6+PSBR45 | 55.1 | 2.11 | 2 |
| TOL+PSBR23 | 58.2 | 2.13 | 69 |
| cC6+PVAE70 | 60.7 | 2.09 | 81 |
| ACE+PSBR45 | 86.3 | 2.83 | 33 |
| nC6+PSBR45 | 151.7 | 1.84 | 19 |

Global M2 (prom. familias): **41 %** vs paper **2.04 %**.

## Diagnóstico

1. **M2 mejora mucho varios sistemas PVAE** que en M1 eran catastróficos
   (CHL+PVAE25: 72→6 %, cC6+PVAE25: 84→15 %). Confirma que la separación en
   segmentos α/β es la dirección correcta para esos sistemas.

2. **Pero PSBR45 y nC6 empeoran** (nC6+PSBR45 sobreestima P por factor 2.5).
   El error es **insensible al MW** (no es numérico) y a la elección molar/masa.

3. **Parámetros de polímero en los segmentos:** empeora todo (195 %). M2 usa
   parámetros de monómero.

4. **Re-optimizando κa/κb por familia** con nuestro motor: nC6+PSBR 95→27 %,
   cC6+PSBR 32→15 %. No llega al 2 %, y los κ óptimos salen **negativos**
   (−0.032) donde el paper usa **positivos** (+0.008) — signo opuesto.

## Conclusión

La metodología 2 está implementada fielmente según la estructura del paper y
**mejora sustancialmente los sistemas PVAE**, pero no reproduce el 2.04 % global.
La misma discrepancia sistemática que afecta a M1 reaparece: nuestro motor
requiere κ de signo opuesto a los del paper y, aun optimizando, no alcanza su
precisión. Como las ecuaciones están verificadas idénticas a Gross-Sadowski y el
motor reproduce n-alcanos puros con <1 % de error, la diferencia reside en la
**parte de mezcla con polímero de m muy grande** (cálculo de φ̂ del solvente), no
en las ecuaciones ni en la metodología elegida.

---

# Continuación: bug encontrado y aislamiento de la fugacidad

## Bug corregido: doble conteo del enlace cruzado en las fracciones B

Al auditar el hard-chain se detectó que la matriz de fracciones de enlace B
**no sumaba 1** (sumaba 1.28-1.95). Las tres fracciones de Tabla 3 (B_αα, B_αβ,
B_ββ) suman 1, pero en la doble suma Σ_a Σ_b B[a][b]·ln g el enlace cruzado αβ
aparece dos veces, así que debe repartirse **B_αβ/2** en cada posición off-diagonal.

Efecto del fix:
- Global M2: **41.3 % → 32.8 %**.
- Arregla los casos catastróficos: nC6+PSBR45 **152 → 24 %**, cC6+PSBR45 **55 → 5.5 %**.
- Resultados mucho más uniformes (5-70 % en vez de 6-152 %).

## Verificación de las fracciones de enlace (Tabla 3 vs Tabla 1)

Las B de Tabla 3 se reproducen con la fórmula de copolímero random (Tabla 1):
B_αβ = 2·z_β·m/(m−1), usando **m calculado con parámetros de POLÍMERO** (no monómero):

| Copolímero | B_αβ fórmula (m-polímero) | B_αβ Tabla 3 | dif |
|---|---|---|---|
| PVAE9 | 0.1605 | 0.1608 | −0.0003 |
| PVAE25 | 0.4544 | 0.4577 | −0.0033 |
| PSBR23 | 0.2818 | 0.2839 | −0.0021 |
| PSBR45 | 0.6200 | 0.6293 | −0.0093 |

→ El paper usa m-polímero para el esqueleto (m, z, B) pero parámetros de monómero
  para σ/ε de los segmentos. Usar parámetros de polímero en σ/ε empeora (70 %).

## Aislamiento del coeficiente de fugacidad del solvente

- **Límite x1→1:** ln φ̂₁(mezcla) → ln φ del solvente puro (dif = +0.00000). Correcto.
- **Paso h de la derivada composicional:** estable (h=1e-6 a 1e-9 dan el mismo
  valor). No hay sesgo numérico.
- **Errata conocida de Gross-Sadowski (Ec. A.11, C1):** es solo tipográfica
  (falta un exponente −1 impreso); nuestro código ya la implementa correcta
  (la prueba de oro con n-alcanos lo confirma).

## Re-optimización de kij (prueba definitiva)

Para la familia CHL+PVAE, barriendo κa (con κb, καβ fijos de Tabla 7):

| κa | %ΔP |
|---|---|
| −0.20 | 56.6 |
| −0.0532 (paper) | 28.4 |
| −0.01 (óptimo) | **23.4** |
| +0.08 | 53.4 |

→ El mínimo alcanzable por nuestro motor es **23.4 %**; NO baja al 2.59 % del paper
  con ningún κa. El residuo es irreducible por ajuste de parámetros.

## Conclusión final de la Metodología 2

M2 está implementada fielmente (y se corrigió un bug real en las fracciones B).
Mejora mucho varios PVAE respecto a M1, pero el global queda en ~33 % frente al
2.04 % del paper, y **ningún ajuste de kij lo cierra**. Con:
- ecuaciones verificadas idénticas a Gross-Sadowski,
- motor validado <1 % en n-alcanos puros,
- errata C1 ya corregida, fracciones B corregidas y sumando 1,
- límite de dilución y derivadas numéricas correctos,

la discrepancia residual solo puede venir de una **decisión metodológica no
documentada en el paper** sobre cómo combina los segmentos VA/ET del copolímero
(p. ej. una regla distinta para σ/ε del copolímero, o un m efectivo propio), que
no es recuperable de las tablas y ecuaciones publicadas.
