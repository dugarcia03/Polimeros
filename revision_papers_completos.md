# Revisión con los papers completos (Arce 2006 + Wohlfarth 1997)

El usuario adjuntó el PDF completo de Arce (pp. 52-63) y el de Wohlfarth 1997
(ref [34], pp. 2689-2702). Se revisaron ambos a fondo. Hallazgos:

## 1. Ecuaciones de Arce confirmadas (texto legible completo)

- **Ec. 22 (a_hc):** ã^hc = m̄·ã^hs − Σ x_i(m_i−1)·Σ_α Σ_β B_{iα,iβ}·ln g. ✓ (código OK)
- **Ec. 28 (a1):** −2πρ(Σ aₖηᵏ)·ΣΣ x x m m ΣΣ z z (ε/kT)σ³. ✓
- **Ec. 29 (a2):** −πρm̄·(1+Z^hc+ρ∂Z^hc/∂ρ)⁻¹·(Σ bₖηᵏ)·ΣΣ... (ε/kT)²σ³. ✓
  El término (1+Z^hc+...)⁻¹ es exactamente nuestro C1. ✓
- **Ec. 30, 31, 32, 33:** todas coinciden con el código.

## 2. Fracciones de enlace B — relación correcta confirmada (Tabla 1)

La Tabla 1 de Arce: para random con z_iβ<z_iα:
  B_{iα,iα} = 1 − B_{iα,iβ} − B_{iβ,iα}  ;  B_{iβ,iβ} = 0

En la Tabla 3, para PVAE9: B_aa=0.8392, B_ab=0.1608, B_bb=0 → **suman 1**.
Esto significa que el 0.1608 de la Tabla 3 es el **cruzado total** (B_iα,iβ + B_iβ,iα).
En la doble suma ΣΣ se reparte 0.0804 en [α][β] y 0.0804 en [β][α].

**Verificado que el código M2 lo hace correcto:** matriz B = [[0,0.0804],[0.0804,0.8392]],
coef cruzado total = 0.1608, diagonal = 0.8392, ΣΣB = 1.0000. El fix /2 era correcto.

## 3. Wohlfarth (ref [34]): número de segmentos con otra definición

Wohlfarth define r_B = M_B/M_A (MW polímero / MW solvente), dando r_B ≈ 1000-1700.
Arce usa m = (m/MW)·MW_copol, dando m ≈ 23000-31000 para los mismos SBR. Son
definiciones de "segmento" distintas (cada EoS tiene la suya); no comparables.
El observable común es P vs w1(másico).

**Dato nuevo:** vía r_B, Wohlfarth implica MW_copol ~86000-100000 para SBR23/41,
mientras Arce usa MW=720000 para PSBR23. Posible diferencia de MW en Arce, pero
el %ΔP es insensible al MW (ya verificado), así que no afecta.

## 4. Wohlfarth Tabla 1a: datos de validación independientes

La Tabla 1a de Wohlfarth da G^E y γ₁/γ₀ para SBR+solvente (acetona, cloroformo,
cC6, hexano, tolueno...). Confirma el signo de las no-idealidades:
- SBR+cloroformo: γ₁/γ₀ **negativo** (−1.07 a −1.09) → desviación negativa fuerte
  (el cloroformo es muy afín al copolímero). Consistente con los kij negativos que
  Arce usa para CHL+PVAE.
- SBR+cyclohexane: γ cercano a 0 (sistema casi atérmico).
- SBR+acetona: γ positivo moderado.

Esto respalda físicamente los signos de kij del paper de Arce.

## 5. Conclusión de la revisión

Tener ambos papers completos **confirma** todo lo verificado antes:
- Las ecuaciones del código son transcripción fiel de Arce (y de Gross-Sadowski).
- Las fracciones de enlace B están bien interpretadas (reparto /2 correcto).
- Los signos de las no-idealidades (y por tanto de los kij) son físicamente
  correctos según los datos independientes de Wohlfarth.

No aparece ningún error nuevo. La discrepancia de magnitud en los PVAE sigue
atribuida a la sensibilidad de nuestro motor al pseudo-monómero de baja densidad,
un detalle de implementación de Arce no documentado, no a un error de replicación.
