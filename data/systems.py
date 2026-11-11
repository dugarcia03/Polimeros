"""
Definicion de las mezclas copolimero + solvente estudiadas en el paper
(Tabla 2) y la composicion de cada copolimero.

Nomenclatura de copolimeros (seccion 5 del paper):
  "las letras mayusculas representan el nombre del copolimero (compuesto por
   dos unidades-monomero) y el numero es el porcentaje en peso de la primera
   unidad-monomero".

  PSBR<n>  : poly(styrene-co-butadiene)        -> n = %peso de estireno (STR)
                                                   resto = butadieno (BTD)
  PVAE<n>  : poly(vinyl acetate-co-ethylene)   -> n = %peso de acetato de vinilo (VA)
                                                   resto = etileno (ET)
  PVAVC<n> : poly(vinyl acetate-co-vinyl chlor)-> n = %peso de acetato de vinilo (VA)
                                                   resto = cloruro de vinilo (VC)

Para cada copolimero se indica:
  mon1, mon2 : claves de los dos monomeros (mon1 = el del % del nombre)
  w1_mon     : fraccion MASICA del monomero 1 en el copolimero
  homopoly1/2: homopolimero PC-SAFT correspondiente a cada monomero (Tabla 4b)

La fraccion molar de cada monomero (para Eqs. 13-17) se calcula a partir de
las fracciones masicas y los pesos moleculares de los monomeros.
"""

from data.parameters import LGT_PURE


# --- definicion de copolimeros ------------------------------------------------
# w1_mon = porcentaje en peso del primer monomero / 100
COPOLYMERS = {
    "PSBR4.1": {"mon1": "STR", "mon2": "BTD", "w1_mon": 0.041, "homo1": "PS", "homo2": "PB"},
    "PSBR5":   {"mon1": "STR", "mon2": "BTD", "w1_mon": 0.050, "homo1": "PS", "homo2": "PB"},
    "PSBR23":  {"mon1": "STR", "mon2": "BTD", "w1_mon": 0.230, "homo1": "PS", "homo2": "PB"},
    "PSBR41":  {"mon1": "STR", "mon2": "BTD", "w1_mon": 0.410, "homo1": "PS", "homo2": "PB"},
    "PSBR45":  {"mon1": "STR", "mon2": "BTD", "w1_mon": 0.450, "homo1": "PS", "homo2": "PB"},
    "PVAE9":   {"mon1": "VA",  "mon2": "ET",  "w1_mon": 0.090, "homo1": "PVA", "homo2": "PE"},
    "PVAE25":  {"mon1": "VA",  "mon2": "ET",  "w1_mon": 0.250, "homo1": "PVA", "homo2": "PE"},
    "PVAE50":  {"mon1": "VA",  "mon2": "ET",  "w1_mon": 0.500, "homo1": "PVA", "homo2": "PE"},
    "PVAE70":  {"mon1": "VA",  "mon2": "ET",  "w1_mon": 0.700, "homo1": "PVA", "homo2": "PE"},
    "PVAVC12": {"mon1": "VA",  "mon2": "VC",  "w1_mon": 0.120, "homo1": "PVA", "homo2": "PVC"},
}

# --- familias (para metodologia 2, Tabla 7) ----------------------------------
FAMILY = {
    "PSBR4.1": "PSBR", "PSBR5": "PSBR", "PSBR23": "PSBR", "PSBR41": "PSBR", "PSBR45": "PSBR",
    "PVAE9": "PVAE", "PVAE25": "PVAE", "PVAE50": "PVAE", "PVAE70": "PVAE",
    "PVAVC12": "PVAVC",
}

# --- Tabla 2: sistemas solvente+copolimero, rango de w1 (solvente) y T -------
# w1_range = (min, max) fraccion masica de solvente en fase copolimero
# T_list   = temperaturas [K] a las que se estudio el sistema
SYSTEMS = {
    ("ACE",   "PSBR5"):   {"np": 15, "w1_range": (0.0003, 0.1200), "T_list": [323.15]},
    ("ACE",   "PSBR23"):  {"np": 22, "w1_range": (0.0090, 0.2310), "T_list": [323.15, 333.15]},
    ("ACE",   "PSBR45"):  {"np": 24, "w1_range": (0.0040, 0.2580), "T_list": [323.15, 333.15]},

    ("nC4OH", "PVAE9"):   {"np": 7,  "w1_range": (0.0060, 0.0360), "T_list": [353.15]},
    ("nC4OH", "PVAVC12"): {"np": 7,  "w1_range": (0.0070, 0.0460), "T_list": [353.15]},

    ("CHL",   "PVAE9"):   {"np": 9,  "w1_range": (0.0310, 0.2530), "T_list": [333.15]},
    ("CHL",   "PVAE25"):  {"np": 8,  "w1_range": (0.0620, 0.4810), "T_list": [333.15]},
    ("CHL",   "PVAE50"):  {"np": 8,  "w1_range": (0.0980, 0.5520), "T_list": [333.15]},
    ("CHL",   "PVAE70"):  {"np": 8,  "w1_range": (0.0590, 0.5730), "T_list": [333.15]},

    ("cC6",   "PSBR5"):   {"np": 9,  "w1_range": (0.0270, 0.3230), "T_list": [333.15]},
    ("cC6",   "PSBR23"):  {"np": 15, "w1_range": (0.0150, 0.3280), "T_list": [296.65, 333.15]},
    ("cC6",   "PSBR41"):  {"np": 8,  "w1_range": (0.2169, 0.4766), "T_list": [343.15]},
    ("cC6",   "PSBR45"):  {"np": 16, "w1_range": (0.0110, 0.3010), "T_list": [296.65, 333.15]},

    ("cC6",   "PVAE25"):  {"np": 9,  "w1_range": (0.0110, 0.3690), "T_list": [353.15]},
    ("cC6",   "PVAE50"):  {"np": 9,  "w1_range": (0.0050, 0.2870), "T_list": [353.15]},
    ("cC6",   "PVAE70"):  {"np": 9,  "w1_range": (0.0040, 0.2420), "T_list": [353.15]},

    ("BZN",   "PSBR4.1"): {"np": 20, "w1_range": (0.0459, 0.4297), "T_list": [298.15, 308.15]},

    ("TOL",   "PSBR5"):   {"np": 6,  "w1_range": (0.0430, 0.5130), "T_list": [308.15]},
    ("TOL",   "PSBR23"):  {"np": 6,  "w1_range": (0.0290, 0.0380), "T_list": [308.15]},
    ("TOL",   "PSBR41"):  {"np": 16, "w1_range": (0.1832, 0.5692), "T_list": [343.15, 373.15]},

    ("nC6",   "PSBR23"):  {"np": 6,  "w1_range": (0.0140, 0.2620), "T_list": [343.15]},
    ("nC6",   "PSBR45"):  {"np": 6,  "w1_range": (0.0160, 0.2490), "T_list": [343.15]},

    ("nC5",   "PSBR23"):  {"np": 9,  "w1_range": (0.0080, 0.0890), "T_list": [333.15]},
    ("nC5",   "PSBR45"):  {"np": 6,  "w1_range": (0.0100, 0.0340), "T_list": [333.15]},
}

# nombres legibles de los solventes
SOLVENT_NAMES = {
    "ACE": "acetona", "nC4OH": "1-butanol", "CHL": "cloroformo",
    "cC6": "ciclohexano", "BZN": "benceno", "TOL": "tolueno",
    "nC6": "n-hexano", "nC5": "n-pentano",
}

COPOLYMER_NAMES = {
    "PSBR": "poli(estireno-co-butadieno)",
    "PVAE": "poli(acetato de vinilo-co-etileno)",
    "PVAVC": "poli(acetato de vinilo-co-cloruro de vinilo)",
}


def monomer_mole_fractions(copolymer):
    """Fraccion MOLAR de cada monomero en el copolimero, a partir de la
    fraccion masica (nombre) y los pesos moleculares de monomero (Tabla 4a).

    Devuelve (x_mon1, x_mon2).
    """
    cp = COPOLYMERS[copolymer]
    w1 = cp["w1_mon"]
    w2 = 1.0 - w1
    mw1 = LGT_PURE[cp["mon1"]]["MW"]
    mw2 = LGT_PURE[cp["mon2"]]["MW"]
    n1 = w1 / mw1
    n2 = w2 / mw2
    x1 = n1 / (n1 + n2)
    x2 = n2 / (n1 + n2)
    return x1, x2


def monomer_weight_fractions(copolymer):
    """Fraccion MASICA de cada monomero (w_mon1, w_mon2)."""
    cp = COPOLYMERS[copolymer]
    return cp["w1_mon"], 1.0 - cp["w1_mon"]
