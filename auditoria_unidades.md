# Auditoría dimensional completa de los motores PC-SAFT (M1 y M2)

Se verifica que cada ecuación sea dimensionalmente congruente. Sistema de unidades
de Gross-Sadowski: longitudes en **Å**, densidad como **densidad numérica** ρ en
**Å⁻³** (moléculas por Å³), energías reducidas (ε/k) en **K**, Helmholtz reducida
ã = A/(NkT) **adimensional**.

## 0. Unidades de las variables base

| Variable | Código | Unidad | Nota |
|---|---|---|---|
| σ (sigma) | `sigma` | Å | diámetro de segmento |
| ε/k (eps_k) | `eps_k` | K | energía / Boltzmann |
| m | `m` | – (adimensional) | nº de segmentos por molécula |
| ρ (rho) | `rho` | Å⁻³ | densidad numérica de moléculas |
| T | `T` | K | temperatura |
| k_B (KB) | `KB` | J/K = 1.380649e-23 | Boltzmann |
| d_i(T) | `d_seg` | Å | diámetro dependiente de T |
| x_i | `x` | – | fracción molar |
| z_{i,α} | `z` | – | fracción de segmento (Σ=1) |
| ζ_k | `zeta[k]` | Å^(k−3) | ver abajo |
| η = ζ₃ | – | – (adimensional) | fracción de empaquetamiento |

## 1. Diámetro d_i(T) — Ec. 26 / A.9

```
d = sigma * (1 - 0.12*exp(-3*eps_k/T))
```
- `eps_k/T` = K/K = **adimensional** ✓ (por eso ε/k se guarda en K, no en J)
- `exp(...)` adimensional; `0.12` adimensional
- `d = sigma·(adim)` → **[Å]** ✓

## 2. ζ_k (zeta) — Ec. 25 / A.8

```
zeta[k] = (pi/6) * rho * sum_i x_i m_i sum_a z_ia * d_ia^k
```
- `rho` [Å⁻³], `x_i` [–], `m_i` [–], `z_ia` [–], `d^k` [Å^k]
- `zeta[k]` = Å⁻³ · Å^k = **[Å^(k−3)]**
  - ζ₀ [Å⁻³], ζ₁ [Å⁻²], ζ₂ [Å⁻¹], ζ₃ [Å⁰ = adimensional] ✓
- ζ₃ = η debe ser adimensional y <1. Verificado numéricamente: η=0.365 ✓

## 3. a_hs (ã^hs) — Ec. A.6 (Boublík-Mansoori)

```
a_hs = (1/z0)*( 3*z1*z2/(1-z3) + z2^3/(z3*(1-z3)^2) + (z2^3/z3^2 - z0)*ln(1-z3) )
```
Chequeo de cada término (debe dar adimensional, pues ã = A/NkT):
- `3*z1*z2/(1-z3)` = [Å⁻²·Å⁻¹]/[–] = Å⁻³ ; dividido por z0 [Å⁻³] → **adim** ✓
- `z2^3/(z3*(1-z3)^2)` = [Å⁻¹]³/[–] = Å⁻³ ; /z0 → **adim** ✓
- `z2^3/z3^2` = Å⁻³/[–] = Å⁻³ ; menos `z0` [Å⁻³] → Å⁻³ (homogéneo ✓); ×ln(1−z3)[–]; /z0 → **adim** ✓
- `1/z0` = [Å³]
→ a_hs **adimensional** ✓. (La resta `z2^3/z3^2 − z0` es dimensionalmente homogénea: ambos Å⁻³.)

## 4. g_seg (g^hs) — Ec. 23 / A.7

```
dd = di*dj/(di+dj)                       # [Å·Å/Å] = [Å]
g  = 1/(1-z3) + dd*3*z2/(1-z3)^2 + dd^2*2*z2^2/(1-z3)^3
```
- `1/(1-z3)` [–]
- `dd*z2` = Å·Å⁻¹ = **[–]** ✓
- `dd^2*z2^2` = Å²·Å⁻² = **[–]** ✓
→ g^hs **adimensional** ✓ (como debe ser, es una función de distribución radial)

## 5. a_hc (ã^hc) — Ec. 22

```
a_hc = m_bar * a_hs - sum_i x_i (m_i-1) * sum_a sum_b B_ab * ln g_ab
```
- `m_bar*a_hs` = [–]·[–] = **adim** ✓
- `x_i (m_i-1)` [–]; `B_ab` [–]; `ln g` [–] → **adim** ✓
→ a_hc **adimensional** ✓

## 6. Dispersión: a_k, b_k, I1, I2, C1 — Ecs. 28-32 / A.11,A.16-19

```
a_k, b_k = polinomios en (m_bar-1)/m_bar     # [–], constantes aik/bik adimensionales
I1 = sum a_k * eta^k ; I2 = sum b_k * eta^k  # eta[–] -> I1,I2 [–]
C1 = 1/(1 + m_bar*f(eta) + (1-m_bar)*g(eta)) # todo [–] -> C1 [–]
```
→ I1, I2, C1 **adimensionales** ✓

## 7. m2es (m²εσ³) — Ec. 28/29, A.12/A.13

```
m2es = sum_i sum_j x_i x_j m_i m_j sum_a sum_b z_ia z_jb (eps_ab/T)^p * sig_ab^3
```
- `eps_ab/T` = K/K = **adim** ✓ (clave: ε/k [K] dividido por T [K]); elevado a p [–]
- `sig_ab^3` = **[Å³]**
- `x·x·m·m·z·z` todos [–]
→ m2es = **[Å³]** (para p=1 y p=2 igual, porque (ε/T)^p es adimensional) ✓

## 8. a_disp (ã^disp) — Ecs. 28-29

```
a1 = -2*pi * rho * I1 * m2es(p=1)
a2 = -pi * rho * m_bar * C1 * I2 * m2es(p=2)
```
- `rho` [Å⁻³] · `m2es` [Å³] = **[–]** ✓ ; I1,I2,C1,m_bar,2π [–]
→ a1, a2, a_disp **adimensionales** ✓

## 9. a_assoc (ã^assoc) — Ecs. 36-38

```
delta = sigma^3 * g * kAB * (exp(epsAB/T) - 1)
rho_i = rho * x_i * z_ia
rd    = rho_i * delta
XA    = (-1+sqrt(1+4*rd))/(2*rd)
a_assoc = sum x_i z_ia * 2*(ln XA - XA/2 + 1/2)
```
- `epsAB/T` = K/K = **adim** ✓ ; `exp(..)-1` [–]
- `kAB` [–] (volumen de asociación reducido, adimensional); `sigma^3` [Å³]
- `delta` = **[Å³]**
- `rho_i` [Å⁻³] · `delta` [Å³] = **rd [–]** ✓ (argumento de sqrt debe ser adim)
- `XA` [–]; `a_assoc` [–] ✓
→ a_assoc **adimensional** ✓

## 10. Z (factor de compresibilidad)

```
Z = 1 + rho * d(a_res)/d(rho)
```
- `a_res` [–], derivada respecto a `rho` [Å⁻³] → d(a_res)/drho [Å³]
- `rho * [Å³]` = Å⁻³·Å³ = **[–]** → Z **adimensional** ✓

## 11. Presión — Ec. A.23

```
pressure = Z * (rho * 1e30) * KB * T
```
- `rho*1e30`: convierte Å⁻³ → m⁻³ (1 Å⁻³ = 1e30 m⁻³) → [m⁻³]
- `Z` [–] · [m⁻³] · `KB` [J/K] · `T` [K] = J/m³ = N·m/m³ = N/m² = **[Pa]** ✓
→ P en Pascales. Factor 1e30 verificado: (10¹⁰ Å/m)³ = 1e30. Coincide con Ec. A.23 de GS.

## 12. ln_phi (coef. de fugacidad) — Ecs. 34-35

```
ln_phi_i = a_res + (Z-1) + da/dx_i - sum_j x_j da/dx_j - ln Z
```
- Todos los términos [–] (a_res, Z, derivadas ∂ã/∂x_i con x adim, ln Z) → **adim** ✓
→ ln φ adimensional, como debe ser ✓

## 13. Conversión w1 → x1

```
n_s = w1 / MW_s ; n_c = (1-w1)/MW_c ; x1 = n_s/(n_s+n_c)
```
- `w1/MW` = (g/g)/(g/mol) = mol/g ·(base 1 g) → proporcional a moles
- x1 = cociente de moles → **adim** ✓. Las unidades de MW (g/mol) se cancelan.

## 14. Equilibrio — Ec. 39

```
log(x1) + ln_phi_L1 - ln_phi_V1 = 0
```
- `x1` [–], `ln_phi` [–] → argumento del log adimensional ✓

## 15. %ΔP — Ec. 40

```
100/np * sum |P_calc - P_exp| / P_exp
```
- `|Pc-Pe|/Pe` = Pa/Pa = **adim** ; ×100 → porcentaje ✓

## Conclusión

**Todas las ecuaciones son dimensionalmente congruentes.** Puntos clave que hacen
que funcione:
1. ε/k se guarda en **K**, de modo que ε/(kT) aparece como `eps_k/T` [K/K] sin
   necesidad de multiplicar por k_B. Correcto y consistente en d_i, m2es y a_assoc.
2. ρ es densidad **numérica** en Å⁻³; con σ en Å, todos los productos ρ·σ³, ρ·d^k·...
   salen adimensionales donde deben.
3. El factor **1e30** en la presión convierte Å⁻³ → m⁻³ para dar Pa con k_B en J/K.
4. ã = A/NkT es adimensional en todos sus términos (hs, hc, disp, assoc).

Verificación numérica (n-hexano, 330 K): η=0.365 (<1, adim), densidad másica
0.62 g/cm³ (física), P=0.070 MPa (objetivo). Dimensionalmente sano.

---

# Verificación: los datos de las tablas se pasan en las unidades correctas

Revisión de TODOS los archivos con datos (no solo los motores), comprobando que
los valores de cada tabla entren a las ecuaciones en la unidad que estas esperan.

## Tabla 4b (solventes, monómeros, polímeros)

| Parámetro | Rótulo paper | Unidad real de uso | Transformación en código | ¿OK? |
|---|---|---|---|---|
| m/MW | "(mol)⁻¹" | **[mol/g]** | `m = (m/MW)·MW` → m **adim** | ✓ |
| σ | Å | Å | directo en d_i, σ_ij | ✓ |
| ε/k | K | K | `eps_k/T` [K/K] adim | ✓ |
| kAB | – | adim | directo | ✓ |
| εAB | K | K | `epsAB/T` [K/K] adim | ✓ |

**Punto clave (m/MW):** el rótulo "(mol)⁻¹" de la Tabla 4b es una simplificación;
dimensionalmente m/MW es **[mol/g]** (segmentos·mol/g). Al multiplicar por MW [g/mol]
las unidades se cancelan y m sale adimensional en rango físico:
- solventes: m ≈ 2.2–2.9 (verificado: nC6 m=2.94, cercano al GS 3.06).
- copolímeros: m = (m/MW)·MW_copol ≈ 4000–30000 segmentos.

**Consistencia MW_copol:** el mismo MW_copol [g/mol] se usa en m (Ec.18) y en la
conversión w1→x1. Coherente.

## Tabla 5 (pseudo-monómeros): m/MW [mol/g], σ [Å], ε/k [K] — mismas unidades. ✓

## Tabla 3 (fracciones): z, B adimensionales; Σz=1 y ΣB=1 verificados. ✓

## Tabla 6 / 7 (kij): adimensionales, entran en (1−kij). ✓

## Tabla 2 (datos experimentales)

| Dato | Unidad | Nota | ¿OK? |
|---|---|---|---|
| w1 | – (frac. másica) | directo en w1→x1 | ✓ |
| T | K | directo | ✓ |
| P | MPa (paper) | datos de Gupta en **kPa** convertidos a **MPa** (÷1000) | ✓ |

**Conversión de presión verificada:** los rangos de P cargados (en MPa) coinciden
exactamente con la Tabla 2 en los 7 sistemas chequeados (ACE+PSBR5 0.0137–0.0613,
CHL+PVAE9 0.0140–0.0871, etc.). El motor devuelve P en Pa y divide por 1e6 → MPa;
se compara MPa contra MPa. Congruente.

## MW de copolímeros (código vs Tabla 3)

PSBR5=100000, PSBR23=720000, PSBR41=100000, PSBR45=600000, PVAE9=100000,
PVAE25=88000, PVAE50=250000, PVAE70=280000, PVAVC12=30000. **Todos coinciden
con la columna "MW copolymer" de la Tabla 3.** ✓

## Constantes universales (Gross-Sadowski)

A_COEF, B_COEF: 7×3 adimensionales (polinomios en η y m̄, ambos adim); a00=0.9105631445. ✓

## Conclusión

**Todos los datos de las tablas se pasan a las ecuaciones en las unidades correctas,
sin conversiones incorrectas ni implícitas erróneas.** Las dos conversiones no
triviales están bien hechas:
1. m/MW [mol/g] × MW [g/mol] → m adimensional.
2. Presiones de Gupta kPa → MPa (÷1000), consistente con Tabla 2 y con la salida
   del motor (Pa → MPa, ÷1e6).
No se detectó ningún dato en unidad incorrecta en ningún archivo.
