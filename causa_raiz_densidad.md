# Causa raíz: los parámetros de monómero dan densidad de copolímero irreal

Siguiendo la pista del ε del VA, se revisó la **coherencia de los parámetros**
calculando la densidad líquida que predice cada "componente puro". Esto reveló la
causa física raíz del error que se arrastra desde M1 (y afecta también a M2).

## Densidad líquida predicha (g/cm³, 1 bar)

### Polímeros puros (Tabla 4b) — CORRECTOS
| Componente | ρ modelo | ρ real | ¿OK? |
|---|---|---|---|
| PVA | 1.144 | ~1.19 | ✓ |
| PE | 0.818 | ~0.92 | ✓ |
| PS | 1.059 | ~1.05 | ✓ |
| PB | 0.873 | ~0.90 | ✓ |

Los polímeros se calibraron con datos PVT (ref [44]) → densidades correctas.

### Pseudo-monómeros (Tabla 5, derivados de MONÓMEROS) — IRREALES
| Pseudo-monómero | ρ modelo | ρ real esperada | error |
|---|---|---|---|
| PVAE25 | **0.552** | ~1.0–1.1 | ~50 % bajo |
| PVAE50 | **0.545** | ~1.0–1.1 | ~50 % bajo |
| PVAE70 | **0.555** | ~1.0–1.1 | ~50 % bajo |
| PSBR5 | 0.723 | ~0.92 | ~22 % bajo |
| PSBR23 | 0.737 | ~0.92 | ~20 % bajo |
| PSBR45 | 0.761 | ~0.92 | ~17 % bajo |

### Copolímero con parámetros de POLÍMERO (fracción másica) — RAZONABLES
| Copolímero | ρ modelo |
|---|---|
| PVAE25 | 0.885 |
| PVAE50 | 0.961 |
| PSBR45 | 0.906 |

## Diagnóstico

1. **Los parámetros de monómero (Tabla 4b) producen densidades de copolímero
   demasiado bajas.** El defecto es grave en PVAE (0.55 vs ~1.0, la mitad) y
   moderado en PSBR (0.74 vs 0.92).

2. **Causa:** el monómero VA tiene σ=5.07 Å (anómalamente grande; su polímero PVA
   tiene σ=3.51) y m/MW=0.0115 (muy bajo). Segmentos grandes y escasos → mucho
   volumen libre → densidad irreal → termodinámica de mezcla distorsionada.

3. **Esto explica todo lo observado y es común a M1 y M2:**
   - PSBR funciona razonablemente (densidad algo baja, pero cercana).
   - PVAE falla catastróficamente (densidad = mitad de lo físico).
   - Bajar el ε del copolímero "arregla" numéricamente el VLE porque compensa el
     exceso de volumen libre.
   - El error se arrastra desde M1 porque M1 y M2 usan los mismos parámetros de
     monómero para el esqueleto del copolímero.

## Por qué el paper sí funciona (hipótesis consistente)

El paper calibró **sus kij con los parámetros de monómero tal cual** (densidad
"mala" incluida). Su conjunto {parámetros de monómero + kij} es internamente
consistente: los kij negativos de los PVAE compensan, dentro de SU implementación,
la distorsión del volumen libre. En nuestro motor (verificado contra Gross-Sadowski
y n-alcanos), esa misma combinación no se reproduce porque el acoplamiento exacto
entre el término de volumen libre del pseudo-monómero y el kij depende de detalles
de cómo el paper evalúa la dispersión del copolímero, que no están publicados.

Usar parámetros de polímero da densidades correctas pero rompe la consistencia con
los kij del paper (que fueron ajustados con parámetros de monómero), por eso empeora
el VLE (70 %).

## Conclusión

La causa raíz del error NO es la metodología (M1 vs M2) ni las ecuaciones ni los
datos, sino la **falta de autoconsistencia de los parámetros de monómero de la
Tabla 4b**: no reproducen la densidad del polímero (crítico en VA). El paper absorbe
esa inconsistencia en sus kij ajustados; nosotros no podemos replicar ese
acoplamiento sin la información de cómo el paper calibró los kij contra ese esqueleto
de baja densidad.

---

# Corrección y matiz importante (revisión posterior)

Al verificar si los parámetros de monómero están "mal", se encontró un matiz que
corrige la interpretación anterior:

## Los parámetros de monómero SÍ reproducen el monómero puro

Densidad del MONÓMERO puro líquido (Tabla 4b) vs real:

| Monómero | ρ modelo | ρ real | ¿OK? |
|---|---|---|---|
| VA (298 K) | 0.906 | ~0.93 (acetato de vinilo) | ✓ |
| STR (298 K) | 0.894 | ~0.906 (estireno) | ✓ |
| ET (150 K) | 0.594 | ~0.57 | ✓ |
| VC (260 K) | 0.968 | ~0.91 | ✓ |
| BTD (268 K) | 0.650 | ~0.62 | ✓ |

→ Los parámetros de monómero (incluido VA con σ=5.07) **NO son una errata**: están
  bien calibrados al monómero puro. El σ grande del VA es correcto para la molécula
  libre de acetato de vinilo.

## El problema es el pseudo-monómero, no los parámetros base

El pseudo-monómero (combinación aritmética Ec.13-15 de monómeros) da densidad de
copolímero irreal (~0.56 vs ~0.92-0.97 del EVA fundido real), porque hereda el σ
grande del monómero VA libre (5.07), mientras el polímero PVA real es compacto
(σ=3.51). La combinación de monómeros "esponjosos" produce un pseudo-segmento con
demasiado volumen por masa.

## Lo que esto significa para la replicación

- Esto es **inherente al método del paper** (pseudo-monómero de monómeros), no un
  error nuestro ni de los datos.
- El paper usó exactamente este pseudo-monómero de baja densidad y obtuvo ~4 %.
- Nuestro motor (verificado contra Gross-Sadowski y n-alcanos) es muy sensible a
  esa baja densidad en el régimen de los PVAE; el del paper, aparentemente, no.
- **No corresponde reajustar kij** (sería fabricar datos nuevos, no replicar).

## Estado del diagnóstico (honesto)

Verificado y correcto, usando SOLO datos del artículo:
- Ecuaciones = Gross-Sadowski (idénticas); motor validado <1 % en n-alcanos.
- Parámetros de monómero, solvente y polímero = Tabla 4b (reproducen puros).
- Pseudo-monómeros = Tabla 5 (media aritmética molar, exacta).
- Fracciones B = Tabla 3 (fórmula Tabla 1, exacta); datos = Gupta (verificados).

La discrepancia persiste y su origen está localizado en la **sensibilidad del VLE
al pseudo-monómero de baja densidad de los PVAE**, un acoplamiento que nuestro motor
PC-SAFT estándar maneja distinto a como lo hizo el paper. No es recuperable de la
información publicada sin conocer el detalle exacto de la implementación del paper,
y no se puede cerrar reajustando parámetros sin dejar de replicar el artículo.
