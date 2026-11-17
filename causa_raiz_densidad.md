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
