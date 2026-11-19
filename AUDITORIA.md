# Auditoría de `metodologia_1_y_2.py` (PC-SAFT, Arce et al. 2006)

Revisión a fondo de ambas metodologías buscando errores de unidades, escritura de
ecuaciones, valores de constantes y cualquier otra fuente de desviación de la precisión.

## Resumen ejecutivo

**El motor termodinámico PC-SAFT está correctamente implementado.** Las constantes,
las ecuaciones (hard-sphere, hard-chain, dispersión, C1, presión, fugacidad) y las
unidades se verificaron una por una contra el paper y contra la formulación estándar
de Gross & Sadowski (2001). La validación empírica más potente —reproducir la presión
de saturación de solventes puros— confirma el motor:

| Solvente | Psat calc (MPa) | Psat Antoine/NIST | error |
|----------|-----------------|-------------------|-------|
| ciclohexano (333 K) | 0.05144 | 0.05189 | **0.9 %** |
| benceno (308 K)     | 0.01975 | 0.01976 | **0.1 %** |
| tolueno (343 K)     | 0.02831 | 0.02717 | **4.2 %** |
| n-pentano (333 K)   | 0.23959 | 0.21450 | 11.7 % |
| n-hexano (343 K)    | 0.11712 | 0.10539 | 11.1 % |
| acetona (323 K)     | 0.15829 | 0.08137 | 94.5 % |

La densidad líquida también sale bien (acetona: 0.7594 g/cm³ calc vs 0.765 real).

**Conclusión:** las desviaciones grandes NO vienen de errores de código. Vienen de:
1. Los parámetros de solvente del paper (acetona, n-alcanos) no reproducen bien su
   propia Psat porque se ajustaron sin término polar y, para los n-alcanos, no son
   los parámetros canónicos de Gross-Sadowski (con los canónicos, el motor da n-pentano
   a 0.2 % de error — prueba de que el motor es correcto).
2. La física del pseudo-monómero en los sistemas polares (PVAE/PVAVC), ya documentada.

## Verificación ítem por ítem (motor)

| Componente | Estado | Nota |
|------------|--------|------|
| Constantes A_COEF, B_COEF | ✅ correcto | fila por fila vs Gross-Sadowski 2001 |
| `a_hs` (Boublík-Mansoori) | ✅ correcto | 3 términos y prefactor 1/ζ0 |
| `g_seg` (Ec. 23) | ✅ correcto | `d_i d_j/(d_i+d_j)`, exponentes 1 y 2 |
| `a_hc` (Ec. 22) | ✅ correcto | se reduce bien a 1 segmento; doble suma B para 2 |
| `_C1` (Ec. 30) | ✅ correcto | signos y exponentes exactos |
| `a_disp` (Ec. 28-29) | ✅ correcto | factores 2π y π; I1 usa a_k, I2 usa b_k |
| `_m2es` (Ec. 33 cruzados) | ✅ correcto | σ=½(σi+σj), ε=√(εiεj)(1-kij) |
| `pressure` (unidades) | ✅ correcto | ρ·1e30 → 1/m³, KB en J/K → P en Pa |
| `ln_phi` (Ec. 34-35) | ✅ correcto | se reduce a a_res+(Z-1)-lnZ para puro |
| `a_assoc` (Ec. 36-38) | ✅ fiel al paper | usa σ³ (Ec. 38 del paper usa σ³, no d³); ρ_i·Δ adimensional OK |
| `Z` (derivada numérica) | ✅ aceptable | paso dρ=1e-6ρ; sin impacto en precisión |
| `resolver_densidad` | ✅ correcto en la práctica | ver abajo |

### Sobre el solver de densidad
Se verificó en todos los sistemas reales que las raíces encontradas son las físicas:
η_liq ≈ 0.42-0.43, η_vap ≈ 0.0001-0.0005, y dP/dρ > 0 (ramas estables) en todos los
puntos. No toma raíces metaestables. El techo η=0.5 es suficiente porque los líquidos
de estos sistemas (P < 0.07 MPa) están en η ≈ 0.42.

### Sobre la asociación (σ³ vs d³)
El paper de Arce (Ec. 38) escribe la fuerza de asociación con **σ³**. La convención
estándar de Gross & Sadowski usa **d(T)³** (diferencia ~4 % en Δ a 353 K). El código
usa σ³ para ser **fiel al paper**. Solo afecta a los sistemas con 1-butanol (nC4OH).
Si se quisiera la convención Gross-Sadowski, cambiar `sa.sigma**3` por
`self.d_seg(i,a,T)**3` en `a_assoc`.

## Hallazgo en la metodología 2 (inconsistencia del paper)

Al verificar la fracción de segmento ziβ contra la Tabla 3 del paper:

| Copolímero | ziβ con monómeros | ziβ con polímeros | Tabla 3 |
|------------|-------------------|-------------------|---------|
| PSBR23 | 0.1364 | **0.1409** | 0.1409 |
| PSBR45 | 0.3020 | **0.3100** | 0.3100 |
| PVAE9  | 0.0199 | 0.0949 | 0.0802 |
| PVAE50 | 0.1705 | 0.5147 | 0.4686 |

El paper (sec. 3) dice que M2 usa los parámetros de **poly(α) y poly(β)** (polímeros
PS, PB, PVA, PE), y con ellos ziβ reproduce EXACTAMENTE la Tabla 3 para PSBR.

**Pero** al usar los parámetros de polímero con los kij publicados (Tabla 7), el error
global de M2 sube a ~70 %, mientras que con los parámetros de monómero el error es
~33 %. Esto significa que los kij de la Tabla 7 fueron ajustados por los autores con
los parámetros de monómero (o con una construcción no documentada). No es posible
reproducir el 2.04 % del paper sin reajustar los kij.

**Decisión:** se mantienen los parámetros de MONÓMERO porque dan el mejor acuerdo con
los kij tal como están publicados. El código documenta esta inconsistencia del paper.

## Veredicto global

- **Unidades:** correctas (verificado dimensionalmente y con Psat/densidad).
- **Ecuaciones:** fieles al paper y a Gross-Sadowski.
- **Constantes:** exactas.
- **No hay bugs que desvíen la precisión.** Las desviaciones vienen de los parámetros
  y la física del paper, no de la implementación.
