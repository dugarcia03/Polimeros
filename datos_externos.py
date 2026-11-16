"""
EXTERNAL data and assumptions NOT contained in Arce et al. (2006).

Everything here is kept strictly separate from parametros_articulo.py so that
the reproduction of the paper never relies on tuned or borrowed numbers.
These values are used ONLY for independent validation (e.g. checking that the
PC-SAFT pure-solvent saturation pressure is physically correct), never to fit
or correct any parameter of the paper.

Antoine constants (NIST/standard form):  log10(P[bar]) = A - B/(T[K] + C)
Source: NIST WebBook / standard Antoine compilations.  These are used only in
the pure-component PC-SAFT validation (validation step 1).
"""

# Antoine (P in bar, T in K):  log10(P) = A - B/(T + C)
ANTOINE = {
    # name    A         B          C           valid T range (K)   source note
    "cyclohexane": (3.96988, 1203.526, -50.287, (293, 354)),  # NIST (Willingham)
    "hexane":      (4.00266, 1171.530, -48.784, (286, 343)),  # NIST (Carruth/Kobayashi)
    "benzene":     (4.01814, 1203.835, -53.226, (287, 354)),  # NIST (Willingham)
    "toluene":     (4.07827, 1343.943, -53.773, (308, 385)),  # NIST (Willingham)
}

# Map the paper's solvent codes to the Antoine names above (nonpolar solvents).
SOLVENT_ANTOINE_NAME = {
    "cC6": "cyclohexane",
    "nC6": "hexane",
    "BZN": "benzene",
    "TOL": "toluene",
}


def antoine_psat_MPa(name, T):
    """Saturation pressure [MPa] from Antoine constants (P in bar internally)."""
    A, B, C, _ = ANTOINE[name]
    P_bar = 10.0 ** (A - B / (T + C))
    return P_bar * 0.1  # bar -> MPa
