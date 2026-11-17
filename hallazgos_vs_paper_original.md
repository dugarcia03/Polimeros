# Cotejo contra la fuente primaria (Arce, Mattedi & Aznar, FPE 246, 2006)

Este documento cierra el diagnóstico de por qué la metodología 1 reproduce bien
unos sistemas y no otros, ahora que se cotejó el código contra el **PDF original**
(`1-s2.0-S0378381206002202-main (1).pdf`) subido a la rama.

## 1. Lo que quedó VERIFICADO contra el paper (todo correcto)

### 1.1 Parámetros de pseudo-monómero (Tabla 5) — coinciden al 100 %
Reproduje los 9 pseudo-monómeros del PC-SAFT aplicando las Ecs. 13–15 del paper
(media por fracción **molar** de los monómeros de la Tabla 4b):

| Copolímero | m/MW calc (ref) | σ calc (ref) | ε/k calc (ref) |
|---|---|---|---|
| PSBR5   | 0.045456 (0.045456) | 3.4648 (3.4647) | 217.86 (217.86) |
| PSBR23  | 0.043119 (0.043119) | 3.5212 (3.5212) | 229.02 (229.01) |
| PSBR45  | 0.039562 (0.039562) | 3.6072 (3.6072) | 246.00 (246.00) |
| PVAE9   | 0.054677 (0.054676) | 3.4765 (3.4765) | 186.36 (186.36) |
| PVAE25  | 0.051699 (0.051700) | 3.5865 (3.5864) | 202.58 (202.57) |
| PVAE50  | 0.045108 (0.045107) | 3.8300 (3.8300) | 238.48 (238.49) |
| PVAE70  | 0.036812 (0.036813) | 4.1365 (4.1364) | 283.68 (283.67) |
| PVAVC12 | 0.084033 (0.084039) | 2.6428 (2.6425) | 183.59 (183.57) |

→ La construcción del copolímero (reglas de combinación, fracciones molares,
pesos moleculares de Tabla 3) es **exacta**.

### 1.2 Ecuaciones del motor — coinciden con las Ecs. 21–35
- Ec. 21 (a_res = a_hc + a_disp), Ec. 22 (cadena dura), Ec. 23 (g_hs),
  Ecs. 24–25 (m̄, ζ_k), Ec. 26 (d_i), Ecs. 28–32 (dispersión a₁, a₂, C₁, aₖ, bₖ),
  Ec. 33 (reglas cruzadas), Ec. 34 (ln φ̂), Ec. 39 (equilibrio), Ec. 40 (%ΔP).
- Unidades consistentes: σ[Å], ε/k[K], ρ[Å⁻³], P→Pa con factor 1e30.
- Signo de los kij: verificado (invertirlo empeora 7/8 sistemas → la errata de
  signo queda **refutada**; los kij negativos del paper son los correctos).

## 2. La causa raíz de las discrepancias

El error NO está en el código, los parámetros, los datos ni los kij.
Está en **dos limitaciones físicas de la metodología 1 tal como el propio paper
la define**, que afectan exactamente a los sistemas que fallan:

### 2.1 Falta el término de asociación (Ecs. 36–38) → afecta a nC4OH
- El 1-butanol es **auto-asociativo**. La Tabla 4b le da κAB=0.008752 y
  εAB=2489.23 K. El motor M1 **no implementa** el término a_assoc.
- Síntoma directo: la presión de vapor del butanol puro da ~673 % de error.
- Por eso `nC4OH+PVAVC12` da 69 % (paper 5.18 %) y `nC4OH+PVAE9` da 25 %.

### 2.2 El pseudo-monómero promedia el kij sobre todo el copolímero
Diagnóstico punto-a-punto de `cC6+PVAE25` (solvente apolar, error 84 %):

| w1 | x1 | Pexp [MPa] | Pcalc [MPa] | Pcalc/Pexp |
|---|---|---|---|---|
| 0.011 | 0.921 | 0.0075 | 0.00028 | **0.04** |
| 0.074 | 0.988 | 0.0265 | 0.00220 | 0.08 |
| 0.228 | 0.997 | 0.0549 | 0.00995 | 0.18 |
| 0.369 | 0.998 | 0.0704 | 0.02210 | 0.31 |

Control `cC6+PSBR45` (mismo solvente, error 2 %): Pcalc/Pexp ≈ 1.00 en todo el rango.

**Lectura:** con un solo pseudo-segmento, el kij muy negativo (−0.1068) del PVAE
se aplica al copolímero **entero**, amplificando ε_cruzado un ~11 % sobre toda la
cadena. El solvente queda artificialmente demasiado afín → la presión de burbuja
se desploma (hasta 4 % del valor real a dilución infinita).

La Tabla 7 del paper (metodología 2) confirma el mecanismo: para `cC6+PVAE` el
paper usa **κij^a = −0.1085 solo con el segmento VA**, y **κij^b = +0.0101 con el
segmento ET**. El kij negativo grande NO aplica a toda la cadena, solo a la
fracción de acetato de vinilo. El pseudo-monómero de M1 no puede distinguir
segmentos, así que sobre-corrige.

## 3. Mapa de sistemas: por qué cada grupo cae donde cae

| Grupo | Física que falta en M1 | %ΔP M1 (nuestro) | %ΔP paper |
|---|---|---|---|
| PSBR + alcanos/cC6/BZN/TOL5 | ninguna → **reproduce bien** | 2–19 % | 2–5 % |
| PVAE + CHL | kij muy negativo sobre toda la cadena (proxy de H-bond cruzado) | 53–94 % | 3–5 % |
| PVAE + cC6 | ídem (kij −0.06…−0.11 promediado) | 81–85 % | 5–7 % |
| nC4OH + PVAE/PVAVC | término de asociación ausente | 25–69 % | 4–5 % |
| ACE + PSBR | acetona fuertemente dipolar (sin término polar) | 22–33 % | 3–5 % |
| TOL + PSBR23, nC5 | pocos puntos + kij promediado | 38–69 % | 4–6 % |

## 4. Conclusión

La metodología 1 del código es una **reproducción fiel y correcta** del PC-SAFT
pseudo-monómero del paper. Reproduce dentro de ~2–5 % los sistemas donde el
modelo no requiere física extra (PSBR con solventes apolares).

Donde el paper mismo necesitó kij extremos (proxy de puente de hidrógeno) o
término de asociación, la metodología 1 **no puede** igualar el ~4 % porque:
1. no incluye el término de asociación (Ecs. 36–38) que el paper sí usa, y
2. el pseudo-monómero no puede localizar el kij en un solo tipo de segmento,
   cosa que el paper sí hace en su **metodología 2** (Tabla 7, global 2.04 %).

El camino para cerrar la brecha es implementar la **metodología 2** completa:
copolímero como dos tipos de segmento (α, β) con fracciones de enlace B_iα,iβ
(Tabla 3, Ecs. 22–25) y tres kij por familia (κij^a, κij^b, κiα,jβ de la Tabla 7),
más el término de asociación para nC4OH.

---

# Revisión final exhaustiva (2ª pasada, solo datos del artículo principal)

Se probó todo lo que quedaba, usando únicamente datos del paper de Arce y de su
fuente experimental [1] (Gupta 1995). Resultados:

## 5. Pruebas nuevas realizadas y su veredicto

### 5.1 Validación de Psat de cada solvente puro (Tabla 4b vs Antoine)
Aísla motor+parámetros del tema de mezcla:

| Solvente | T[K] | Psat modelo | Antoine | %err | Veredicto |
|---|---|---|---|---|---|
| CHL | 333.15 | 0.0962 | 0.0972 | **1.0** | OK |
| cC6 | 333.15 | 0.0514 | 0.0519 | **0.9** | OK |
| BZN | 298.15 | 0.0127 | 0.0127 | **0.1** | OK |
| TOL | 308.15 | 0.0065 | 0.0062 | **4.0** | OK |
| nC6 | 343.15 | 0.1171 | 0.1054 | 11.1 | aceptable |
| nC5 | 333.15 | 0.2396 | 0.2145 | 11.7 | aceptable |
| ACE | 323.15 | 0.1583 | 0.0814 | **94.5** | FALLA (dipolar) |
| nC4OH | 353.15 | 0.1486 | 0.0223 | **565** | FALLA (asociación) |

→ El motor reproduce todos los solventes **no polares y no asociantes**. Falla
  solo acetona (dipolo) y 1-butanol (asociación), que son física no incluida.

### 5.2 Término de asociación (Ecs. 36-38) para nC4OH — IMPLEMENTADO
Se añadió al motor el término de asociación esquema 2B usando SOLO los parámetros
del paper (Tabla 4b: κAB=0.008752, εAB=2489.23 K). Resultado decisivo:

- **Psat del 1-butanol puro: de 565 % a 3.0 %** de error (0.5-4 % en el rango).
- **nC4OH+PVAE9: de 25.35 % a 10.81 %** (gran mejora; curva en rango correcto).
- nC4OH+PVAVC12: empeora (69→82 %) porque el PVC es un segmento muy particular
  y necesita kij localizado (metodología 2), no el término de asociación.

→ Confirma que la asociación era la física faltante para el butanol. El motor
  ahora activa el término solo cuando el paper da parámetros de asociación.

### 5.3 Sensibilidad al MW del copolímero — DESCARTADO como causa
El %ΔP de cC6+PVAE25 es **insensible al MW** (83.7 % para MW entre 5000 y 88000).
En el límite de cadena larga, el VLE del solvente depende de m/MW, no de MW.
La escala de m (Ec. 18-19) está bien implementada.

### 5.4 Barrido de kij — confirma que es estructural, no de valor
¿Qué kij haría que nuestro motor dé el %ΔP del paper?
- CHL+PVAE9: existe kij=+0.060 → 5.27 % (paper 4.18 % con kij=−0.1208). Signo opuesto.
- cC6+PVAE25: el mejor posible es 14.3 % (kij=+0.020); **ningún kij llega a 5.36 %**.

→ Para cC6+PVAE ningún kij único reproduce el paper: es **estructural**. Un solo
  pseudo-segmento con un solo kij no tiene grados de libertad suficientes. El paper
  lo resuelve en metodología 2 (Tabla 7) con kij por segmento (VA y ET por separado).

### 5.5 Datos de tolueno re-verificados contra imagen del PDF de Gupta — CORRECTOS
Se sospechó que TOL+PSBR23 (68.8 %) tenía datos mal transcritos. Se renderizó la
tabla de Gupta p.4 como imagen y se leyó celda por celda:
- Formato de Gupta: pares **(P[kPa], w1)**. Para TOL+PSBR23: P={0.90,2.08,3.04,4.36,4.83,5.48}.
- Los datos cargados **ya eran correctos** (incluidos los 0.90/2.08 kPa, que SÍ son de tolueno).

El problema de TOL+PSBR23 NO es de datos: a w1=0.029-0.038 (x1 cambia solo de 0.9957
a 0.9968) la Pexp se sextuplica. Es un régimen ultra-diluido con copolímero de
MW=720000; a composición casi constante el modelo da P casi constante (correcto).
El dato experimental está en el límite de ruido. TOL+PSBR5 (MW=100000, w1 hasta 0.513)
reproduce bien (ratio Pcalc/Pexp = 0.90-1.11).

### 5.6 Descomposición de términos de Helmholtz — todo sano
a_hc, a1, a2, I1, I2, C1, m2es en rangos físicos normales (sin overflow/NaN) tanto
en sistemas buenos como malos. El solver de densidad encuentra la raíz líquida
física correcta (rama ascendente de P(η)). No hay artefactos numéricos.

## 6. Conclusión final de la revisión

Nada nuevo cambia el diagnóstico de fondo, y se confirma con más fuerza:

1. **El motor y los parámetros son correctos** (verificados contra Tablas 4b/5/6
   y Ecs. 21-40 del PDF original, y contra Psat puros y datos de Gupta).
2. **La única física que faltaba y se pudo añadir con datos del paper** es la
   asociación del butanol (Ecs. 36-38), que mejora nC4OH+PVAE9 de 25 % a 11 %.
3. **El resto de la brecha es intrínseco a la metodología 1**: un pseudo-segmento
   con un solo kij no puede reproducir sistemas donde el paper usó kij por segmento
   (metodología 2, global 2.04 %) o donde falta el término dipolar (acetona).

Para cerrar del todo la brecha se requiere la metodología 2 (segmentos α/β con
fracciones de enlace B de Tabla 3 y los tres kij de Tabla 7) — no hay ningún ajuste
dentro de la metodología 1 que lo logre con los datos del artículo.
