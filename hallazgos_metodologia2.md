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
