# Experimento: añadir el término dipolar a PC-SAFT

Rama: `experimento-termino-polar`. Archivo: `metodologia_polar.py`.

## Qué se probó

El paper de Arce et al. (2006) **no usa término polar** en PC-SAFT: trata acetona,
cloroformo y 1-butanol como no polares (al butanol solo le pone asociación). Este
experimento añade el **término dipolar de Gross & Vrabec (2006)** —la formulación
PCP-SAFT estándar— usando los parámetros m, σ, ε/k del paper **tal cual, sin
reajustar**, y el momento dipolar μ como dato físico tabulado (no ajustado):

| Solvente | μ (Debye) | fuente |
|----------|-----------|--------|
| acetona (ACE) | 2.88 | NIST/CRC |
| cloroformo (CHL) | 1.04 | NIST/CRC |
| 1-butanol (nC4OH) | 1.66 | NIST/CRC |

Implementación: término a_polar = a₂/(1 − a₃/a₂) con las integrales J₂ (coef. a_n,
b_n) y J₃ (coef. c_n) y el momento dipolar reducido μ̃² = μ²/(m ε k_B σ³ · 4πε₀).
Verificado: μ̃²(acetona) = 2.35 (valor físico correcto, orden 1-3).

## Resultados

### Presión de saturación de solventes puros (MPa)

| Solvente | sin polar | con polar | referencia | 
|----------|-----------|-----------|------------|
| acetona (323 K)   | 0.1583 (94%) | **0.0569 (30%)** | 0.0814 |
| cloroformo (333 K)| 0.0962       | 0.0942           | 0.0828 |
| 1-butanol (353 K) | 0.0217       | 0.0203           | — |

El término polar **mejora** la acetona (de 94% a 30% de error), aunque la
sobrecorrige (se pasa por debajo). En cloroformo y butanol el efecto es pequeño.

### Desviación %ΔP en los sistemas VLE (metodología 1)

| Sistema | sin polar | con polar | paper |
|---------|-----------|-----------|-------|
| ACE + PSBR5   | 22.33 | 36.32 | 4.78 |
| ACE + PSBR23  | 25.01 | **20.98** | 3.04 |
| ACE + PSBR45  | 33.40 | **9.57**  | 5.19 |
| CHL + PVAE9   | 93.74 | 93.76 | 4.18 |
| CHL + PVAE25  | 71.95 | 72.14 | 5.06 |
| nC4OH + PVAE9 | 10.81 | 10.99 | 4.05 |
| **GLOBAL (26 sist.)** | **37.62%** | **37.13%** | 4.41% |

Los 20 sistemas con solvente **no polar** (cC6, TOL, nC5, nC6, BZN) dan valores
idénticos con y sin el término (control: el término solo actúa donde μ > 0).

## Interpretación

1. **El término está bien implementado**: tiene efecto físico correcto (mejora la
   Psat de la acetona, μ̃² razonable, no afecta los no polares).

2. **Con los parámetros "tal cual" el resultado es mixto, no una mejora clara**:
   - Acetona: mejora unos sistemas (PSBR45: 33→10%) y empeora otros (PSBR5: 22→36%).
     La razón es el **doble conteo parcial**: los m, σ, ε/k de Arce se ajustaron sin
     término polar, así que ya absorbían algo de la polaridad dentro de la
     dispersión; al sumar el término dipolar explícito se sobrecorrige, y además los
     kij (ajustados sin polaridad) quedan descuadrados.
   - Cloroformo: μ=1.04 D es pequeño → efecto casi nulo (93.74→93.76%). Confirma que
     el gran error de los PVAE **no es la polaridad del solvente** sino el
     pseudo-monómero de baja densidad (diagnóstico de la auditoría).
   - Butanol: la asociación ya captura su no-idealidad → el término dipolar aporta poco.

3. **Global prácticamente igual** (37.62% → 37.13%): solo 4 de 26 sistemas tienen
   solvente polar, y el efecto neto es casi nulo.

## Conclusión

Añadir el término polar **con los parámetros que ya hay** es técnicamente viable
(solo requiere μ, que es dato tabulado) pero **no mejora los resultados de forma
consistente**, por el doble conteo con parámetros ajustados sin polaridad. Para que
el término polar ayude de verdad habría que **re-regresar m, σ, ε/k** de los
solventes polares con μ fijo — eso ya genera parámetros propios, distintos del paper.

Este experimento confirma el diagnóstico de la auditoría: el error de los sistemas
PVAE viene del pseudo-monómero, no de la falta de término polar.
