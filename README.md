# Reproducción VLE de mezclas copolímero + solvente

Implementación en Python que reproduce los resultados del artículo:

> **P.F. Arce, S. Mattedi, M. Aznar** (2006).
> *Vapor–liquid equilibrium of copolymer + solvent mixtures: Thermodynamic
> modeling by two theoretical equations of state.*
> **Fluid Phase Equilibria 246, 52–63.** doi:10.1016/j.fluid.2006.05.009

El programa modela el **equilibrio líquido–vapor (VLE)** de 24 mezclas
copolímero + solvente mediante dos ecuaciones de estado teóricas:

- **PC-SAFT** (Perturbed-Chain SAFT) extendida a copolímeros — *modelo principal,
  totalmente funcional*.
- **LGT** (Lattice Gas Theory, contribución de grupos de Mattedi-Tavares-Castier)
  — *parcialmente implementada* (ver sección «Estado de la LGT»).

---

## 1. Qué hace

El programa es **interactivo**: eliges un sistema (solvente + copolímero) y
calcula la presión de equilibrio frente a la fracción másica de solvente,
`P` vs `w1`, generando:

- la **gráfica** `P` (MPa) vs `w1` para cada temperatura del sistema,
  con las dos metodologías de PC-SAFT del paper;
- la **tabla de presiones** calculadas punto a punto;
- la **desviación de presión** `%ΔP/P` reportada por el paper (Ec. 40) para
  contrastar.

Siguiendo lo acordado, **sólo se grafican las curvas calculadas por los
modelos**; los puntos experimentales (que provienen de las referencias
[1,33,34,35] del paper y no se incluyen en él) se omiten.

## 2. Instalación

```bash
cd Polimeros
python -m venv .venv && source .venv/bin/activate   # opcional
pip install -r requirements.txt
```

## 3. Uso

```bash
# modo interactivo (menú con los 24 sistemas)
python main.py

# calcular un sistema concreto
python main.py --system ACE PSBR5

# generar TODAS las figuras + tabla de desviaciones
python main.py --all
```

Las figuras se guardan en `results/`.

> Nota: si tu entorno tiene el `$HOME` de sólo lectura, exporta
> `MPLCONFIGDIR=/tmp/mpl` antes de ejecutar (matplotlib lo necesita para su
> caché).

## 4. Metodologías del paper

El artículo usa dos estrategias de parametrización (ambas implementadas en
PC-SAFT):

- **Metodología 1 (pseudo-monómero, Tabla 6).** El copolímero se trata con un
  único parámetro de interacción binaria `kij` solvente–copolímero.
  Desviación global reportada: **PC-SAFT 4.41 %**, LGT 10.10 %.
- **Metodología 2 (co-monómeros separados, Tabla 7).** Se ajustan
  `kij_a` (solvente–monómero 1), `kij_b` (solvente–monómero 2) y
  `kiαjβ` (monómero 1–monómero 2). Desviación global reportada:
  **PC-SAFT 2.04 %** (la mejor del paper).

En las gráficas, las curvas **sólidas** son la metodología 1 y las
**discontinuas** la metodología 2.

## 5. Estructura del código

```
Polimeros/
├── main.py                 # programa interactivo + generación batch
├── requirements.txt
├── data/
│   ├── parameters.py       # Tablas 3, 4a, 4b, 5, 6, 7 del paper (transcritas)
│   └── systems.py          # Tabla 2 (24 sistemas) + composición de copolímeros
├── src/
│   ├── pcsaft.py           # EoS PC-SAFT para copolímeros (Ecs. 18–38)
│   ├── lgt.py              # EoS LGT contribución de grupos (Ecs. 1–12)
│   ├── vle.py              # equilibrio PC-SAFT (Ec. 39) y %ΔP (Ec. 40)
│   └── vle_lgt.py          # equilibrio LGT (pseudo-monómero)
└── results/                # figuras generadas
```

## 6. Base teórica (resumen)

Como el copolímero no es volátil, la condición de equilibrio (Ec. 39) es:

```
x1 · φ̂1^L(T, P, x) = φ1^V(T, P, y1 = 1)
```

donde «1» es el solvente y el vapor es solvente puro. A `T` y composición de
líquido dadas, se busca la presión `P` que satisface esa igualdad.

**PC-SAFT para copolímeros** (Gross et al., 2003; Arce & Aznar): cada
copolímero contiene dos tipos de segmento α y β (= los dos monómeros). El
número de segmentos de cada tipo se obtiene de los parámetros del homopolímero
correspondiente (Tabla 4b) vía `m_{i,α} = w_{i,α}·MW_copol·(m/MW)_α` (Ec. 18).
Las fracciones de enlace `B` (Tabla 1 / Tabla 3) ponderan el término de cadena
dura (Ec. 22). La contribución residual = cadena dura + dispersión (Ecs. 21–32).

## 7. Validación

- El motor PC-SAFT produce isotermas `P(ρ)` físicamente correctas (gas ideal al
  diluir, lazo de van der Waals a densidades líquidas).
- Las fracciones de segmento calculadas coinciden con la Tabla 3 del paper
  (p.ej. PSBR23: `z_β(PS) = 0.1409`, exacto).
- Las curvas `P` vs `w1` reproducen la forma y la magnitud de las Figs. 2–10
  (rangos de presión coincidentes con la Tabla 2).

## 8. Estado de la LGT

La EoS LGT de contribución de grupos (Mattedi-Tavares-Castier) está
implementada en `src/lgt.py` y `src/vle_lgt.py`, pero **no converge de forma
fiable**. El motivo es que las ecuaciones (1) y (10) del PDF del paper están
corrompidas por el OCR (términos partidos, el término combinatorio `l` y el
término energético aparecen incompletos), de modo que el factor de
compresibilidad `Z` no tiende correctamente a 1 en el límite de gas ideal.

Reproducir la LGT con fidelidad requiere las ecuaciones originales de las
referencias [19] (tesis de S.M. Silva) o [20] (Mattedi et al., *Fluid Phase
Equilib.* 142, 1998). Por eso el entregable principal es **PC-SAFT**, que es
además el modelo que el propio paper recomienda por su mejor desempeño.

Las desviaciones **LGT** que el paper reporta se incluyen en la tabla de
`main.py --all` a modo de referencia bibliográfica.
