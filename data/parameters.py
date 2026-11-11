"""
Parametros del articulo:

    Arce, P.F.; Mattedi, S.; Aznar, M.
    "Vapor-liquid equilibrium of copolymer + solvent mixtures:
     Thermodynamic modeling by two theoretical equations of state"
    Fluid Phase Equilibria 246 (2006) 52-63.

Todas las tablas de parametros estan transcritas aqui tal como aparecen en el
paper. Las referencias a Tablas/Ecuaciones (Eq.) son las del articulo.

Unidades:
  LGT:
    MW      : g/mol
    v_inf   : 1e-6 m3/mol   (va_inf, Tabla 4a)
    Q       : adimensional  (parametro de area)
    Baa     : K             (Eq. 12)
    u0_R    : K             (uaa0/R, Eq. 11/12)  -> siempre negativo
  PC-SAFT:
    m_MW    : 1/(mol) = (m/MW)   (segmentos por masa)
    sigma   : Angstrom
    eps_k   : K  (epsilon/k)
    kAB     : volumen de asociacion (adim)
    epsAB_k : K  (energia de asociacion)
"""

# ---------------------------------------------------------------------------
# Tabla 4a: parametros de componente puro para la LGT EoS
#   (solventes y monomeros)
#   columnas: MW, v_inf [1e-6 m3/mol], Q, Baa [K], u0_R [K]
# ---------------------------------------------------------------------------
LGT_PURE = {
    # --- solventes ---
    "nC5":    {"MW": 72.15,  "v_inf": 25.613, "Q": 3.0708,  "Baa": 96.331,  "u0_R": -604.796},
    "nC6":    {"MW": 86.18,  "v_inf": 27.284, "Q": 3.4492,  "Baa": 80.985,  "u0_R": -679.584},
    "cC6":    {"MW": 84.16,  "v_inf": 22.923, "Q": 2.5248,  "Baa": 168.80,  "u0_R": -705.286},
    "BZN":    {"MW": 78.12,  "v_inf": 37.975, "Q": 4.5067,  "Baa": 106.84,  "u0_R": -507.720},
    "TOL":    {"MW": 92.14,  "v_inf": 21.601, "Q": 3.6697,  "Baa": 108.28,  "u0_R": -742.905},
    "CHL":    {"MW": 119.38, "v_inf": 22.048, "Q": 4.1747,  "Baa": 50.722,  "u0_R": -642.003},
    "ACE":    {"MW": 58.08,  "v_inf": 27.542, "Q": 3.5557,  "Baa": 111.73,  "u0_R": -584.812},
    "nC4OH":  {"MW": 74.12,  "v_inf": 61.697, "Q": 10.192,  "Baa": 108.14,  "u0_R": -385.077},
    # --- monomeros ---
    "STR":    {"MW": 104.14, "v_inf": 82.074, "Q": 8.5932,  "Baa": 0.7607,  "u0_R": -323.224},
    "VA":     {"MW": 86.09,  "v_inf": 57.695, "Q": 6.1134,  "Baa": 5.3522,  "u0_R": -249.269},
    "BTD":    {"MW": 54.10,  "v_inf": 48.334, "Q": 4.9983,  "Baa": 6.4662,  "u0_R": -231.991},
    "ET":     {"MW": 28.10,  "v_inf": 30.520, "Q": 3.1701,  "Baa": 6.7838,  "u0_R": -243.983},
    "VC":     {"MW": 62.50,  "v_inf": 34.485, "Q": 5.6125,  "Baa": -1.0051, "u0_R": -761.811},
}

# ---------------------------------------------------------------------------
# Tabla 4b: parametros de componente puro para la PC-SAFT EoS
#   columnas: m/MW [1/mol], sigma [A], eps/k [K], kAB, epsAB [K]
# ---------------------------------------------------------------------------
PCSAFT_PURE = {
    # --- solventes ---
    "nC5":    {"m_MW": 0.034090, "sigma": 3.8729, "eps_k": 240.11, "kAB": None,     "epsAB_k": None},
    "nC6":    {"m_MW": 0.034088, "sigma": 3.8383, "eps_k": 239.57, "kAB": None,     "epsAB_k": None},
    "cC6":    {"m_MW": 0.029946, "sigma": 3.8399, "eps_k": 279.21, "kAB": None,     "epsAB_k": None},
    "BZN":    {"m_MW": 0.032712, "sigma": 3.5878, "eps_k": 281.55, "kAB": None,     "epsAB_k": None},
    "TOL":    {"m_MW": 0.030984, "sigma": 3.6869, "eps_k": 282.69, "kAB": None,     "epsAB_k": None},
    "CHL":    {"m_MW": 0.021826, "sigma": 3.4085, "eps_k": 265.94, "kAB": None,     "epsAB_k": None},
    "ACE":    {"m_MW": 0.037879, "sigma": 3.4956, "eps_k": 272.26, "kAB": None,     "epsAB_k": None},
    "nC4OH":  {"m_MW": 0.035126, "sigma": 3.7563, "eps_k": 263.18, "kAB": 0.008752, "epsAB_k": 2489.23},
    # --- monomeros ---
    "STR":    {"m_MW": 0.024339, "sigma": 3.9750, "eps_k": 318.70, "kAB": None, "epsAB_k": None},
    "VA":     {"m_MW": 0.011525, "sigma": 5.0706, "eps_k": 421.42, "kAB": None, "epsAB_k": None},
    "BTD":    {"m_MW": 0.046033, "sigma": 3.4508, "eps_k": 215.10, "kAB": None, "epsAB_k": None},
    "ET":     {"m_MW": 0.056070, "sigma": 3.4250, "eps_k": 178.77, "kAB": None, "epsAB_k": None},
    "VC":     {"m_MW": 0.091211, "sigma": 2.4024, "eps_k": 160.05, "kAB": None, "epsAB_k": None},
    # --- polimeros (homopolimeros) ---
    "PS":     {"m_MW": 0.033238, "sigma": 3.5022, "eps_k": 320.14, "kAB": None, "epsAB_k": None},
    "PVA":    {"m_MW": 0.029905, "sigma": 3.5086, "eps_k": 310.14, "kAB": None, "epsAB_k": None},
    "PB":     {"m_MW": 0.060552, "sigma": 2.9845, "eps_k": 240.21, "kAB": None, "epsAB_k": None},
    "PE":     {"m_MW": 0.028193, "sigma": 4.0125, "eps_k": 320.24, "kAB": None, "epsAB_k": None},
    "PVC":    {"m_MW": 0.021369, "sigma": 3.7486, "eps_k": 370.42, "kAB": None, "epsAB_k": None},
}

# ---------------------------------------------------------------------------
# Tabla 5: parametros de pseudo-monomero (metodologia 1) ya calculados
#   por el paper via Eqs. (13)-(16).  Se incluyen para validacion cruzada
#   con el calculo propio a partir de las reglas de combinacion.
#   LGT: v_inf, Q, Baa, u0_R   |   PC-SAFT: m/MW, sigma, eps/k
# ---------------------------------------------------------------------------
PSEUDO_MONOMER = {
    "PSBR4.1": {"lgt": {"v_inf": 49.067, "Q": 5.0764, "Baa": 6.2829, "u0_R": -233.813},
                "pcsaft": {"m_MW": 0.045562, "sigma": 3.4622, "eps_k": 217.35}},
    "PSBR5":   {"lgt": {"v_inf": 49.232, "Q": 5.0940, "Baa": 6.2421, "u0_R": -234.223},
                "pcsaft": {"m_MW": 0.045456, "sigma": 3.4647, "eps_k": 217.86}},
    "PSBR23":  {"lgt": {"v_inf": 52.866, "Q": 5.4811, "Baa": 5.3754, "u0_R": -243.367},
                "pcsaft": {"m_MW": 0.043119, "sigma": 3.5212, "eps_k": 229.01}},
    "PSBR41":  {"lgt": {"v_inf": 57.283, "Q": 5.9518, "Baa": 4.4088, "u0_R": -254.720},
                "pcsaft": {"m_MW": 0.040280, "sigma": 3.5898, "eps_k": 242.57}},
    "PSBR45":  {"lgt": {"v_inf": 58.396, "Q": 6.0704, "Baa": 4.1805, "u0_R": -257.620},
                "pcsaft": {"m_MW": 0.039562, "sigma": 3.6072, "eps_k": 246.00}},
    "PVAE9":   {"lgt": {"v_inf": 31.368, "Q": 3.2620, "Baa": 6.7365, "u0_R": -244.147},
                "pcsaft": {"m_MW": 0.054676, "sigma": 3.4765, "eps_k": 186.36}},
    "PVAE25":  {"lgt": {"v_inf": 33.182, "Q": 3.4584, "Baa": 6.6361, "u0_R": -244.498},
                "pcsaft": {"m_MW": 0.051700, "sigma": 3.5864, "eps_k": 202.57}},
    "PVAE50":  {"lgt": {"v_inf": 37.198, "Q": 3.8934, "Baa": 6.4163, "u0_R": -245.277},
                "pcsaft": {"m_MW": 0.045107, "sigma": 3.8300, "eps_k": 238.49}},
    "PVAE70":  {"lgt": {"v_inf": 42.257, "Q": 4.4413, "Baa": 6.1447, "u0_R": -246.259},
                "pcsaft": {"m_MW": 0.036813, "sigma": 4.1364, "eps_k": 283.67}},
    "PVAVC12": {"lgt": {"v_inf": 36.576, "Q": 5.6576, "Baa": 0.4954, "u0_R": -704.204},
                "pcsaft": {"m_MW": 0.084039, "sigma": 2.6425, "eps_k": 183.57}},
}

# ---------------------------------------------------------------------------
# Tabla 3: parametros PC-SAFT de copolimeros (fracciones de segmento/enlace)
#   zib      : fraccion de segmento del monomero beta  (z_i,beta)
#   Biaia    : B_{i alpha, i alpha}
#   Biajb    : B_{i alpha, j beta}  (enlace alpha-beta)
#   Bjbjb    : B_{j beta,  j beta}
#   MW       : peso molecular del copolimero [g/mol]
#   Guiones ("-") en el paper -> None (PVT del copolimero no disponible;
#   kiajb no se pudo calcular, se trata como no-copolimero / homopolimero).
# ---------------------------------------------------------------------------
COPOLYMER_PCSAFT = {
    "PSBR4.1": {"MW": 100000, "zib": None,   "Biaia": None,   "Biajb": None,   "Bjbjb": None},
    "PSBR5":   {"MW": 100000, "zib": None,   "Biaia": None,   "Biajb": None,   "Bjbjb": None},
    "PSBR23":  {"MW": 720000, "zib": 0.1409, "Biaia": 0.7161, "Biajb": 0.2839, "Bjbjb": 0.0000},
    "PSBR41":  {"MW": 100000, "zib": None,   "Biaia": None,   "Biajb": None,   "Bjbjb": None},
    "PSBR45":  {"MW": 600000, "zib": 0.3100, "Biaia": 0.3707, "Biajb": 0.6293, "Bjbjb": 0.0000},
    "PVAE9":   {"MW": 100000, "zib": 0.0802, "Biaia": 0.8392, "Biajb": 0.1608, "Bjbjb": 0.0000},
    "PVAE25":  {"MW": 88000,  "zib": 0.2271, "Biaia": 0.5423, "Biajb": 0.4577, "Bjbjb": 0.0000},
    "PVAE50":  {"MW": 250000, "zib": 0.4686, "Biaia": 0.0487, "Biajb": 0.9513, "Bjbjb": 0.0000},
    "PVAE70":  {"MW": 280000, "zib": 0.6729, "Biaia": 0.0000, "Biajb": 0.4633, "Bjbjb": 0.5367},
    "PVAVC12": {"MW": 30000,  "zib": None,   "Biaia": None,   "Biajb": None,   "Bjbjb": None},
}


# ---------------------------------------------------------------------------
# Tabla 6: resultados / parametros de interaccion binaria (metodologia 1:
#   concepto de pseudo-monomero).
#   Cada entrada: (solvente, copolimero) -> parametros + desviacion reportada.
#   LGT:  Bba [K], uba0_R [K]   (dos parametros de interaccion)
#   PC-SAFT: kij                (un parametro, solvente-copolimero)
#   dP_lgt, dP_pcsaft: %DeltaP/P reportado por el paper (Eq. 40).
# ---------------------------------------------------------------------------
TABLE6 = {
    ("ACE",   "PSBR5"):   {"Bba": 34.080,  "uba0_R": -373.500, "kij": 0.0120,  "dP_lgt": 10.53, "dP_pcsaft": 4.78, "T": [323.15, 333.15]},
    ("ACE",   "PSBR23"):  {"Bba": -1.0193, "uba0_R": -422.469, "kij": 0.0122,  "dP_lgt": 7.18,  "dP_pcsaft": 3.04, "T": [323.15, 333.15]},
    ("ACE",   "PSBR45"):  {"Bba": -57.066, "uba0_R": -521.866, "kij": 0.0130,  "dP_lgt": 10.32, "dP_pcsaft": 5.19, "T": [323.15, 333.15]},

    ("nC4OH", "PVAE9"):   {"Bba": 81.008,  "uba0_R": -307.931, "kij": 0.0230,  "dP_lgt": 3.82,  "dP_pcsaft": 4.05, "T": [353.15]},
    ("nC4OH", "PVAVC12"): {"Bba": 53.452,  "uba0_R": -496.475, "kij": -0.0126, "dP_lgt": 4.20,  "dP_pcsaft": 5.18, "T": [353.15]},

    ("CHL",   "PVAE9"):   {"Bba": 32.860,  "uba0_R": -442.43,  "kij": -0.1208, "dP_lgt": 7.07,  "dP_pcsaft": 4.18, "T": [333.15]},
    ("CHL",   "PVAE25"):  {"Bba": 35.057,  "uba0_R": -470.868, "kij": -0.0816, "dP_lgt": 16.86, "dP_pcsaft": 5.06, "T": [333.15]},
    ("CHL",   "PVAE50"):  {"Bba": 36.124,  "uba0_R": -486.782, "kij": -0.0701, "dP_lgt": 16.23, "dP_pcsaft": 4.41, "T": [333.15]},
    ("CHL",   "PVAE70"):  {"Bba": 35.976,  "uba0_R": -489.007, "kij": -0.0516, "dP_lgt": 26.73, "dP_pcsaft": 3.28, "T": [333.15]},

    ("cC6",   "PSBR5"):   {"Bba": 51.596,  "uba0_R": -430.131, "kij": 0.0118,  "dP_lgt": 15.72, "dP_pcsaft": 3.66, "T": [333.15]},
    ("cC6",   "PSBR23"):  {"Bba": -27.743, "uba0_R": -538.326, "kij": 0.0104,  "dP_lgt": 12.72, "dP_pcsaft": 4.38, "T": [296.65, 333.15]},
    ("cC6",   "PSBR41"):  {"Bba": 50.592,  "uba0_R": -418.634, "kij": 0.0088,  "dP_lgt": 6.74,  "dP_pcsaft": 2.28, "T": [343.15]},
    ("cC6",   "PSBR45"):  {"Bba": -47.915, "uba0_R": -596.968, "kij": 0.0101,  "dP_lgt": 12.51, "dP_pcsaft": 4.31, "T": [296.65, 333.15]},

    ("cC6",   "PVAE25"):  {"Bba": 6.1817,  "uba0_R": -481.450, "kij": -0.1068, "dP_lgt": 7.51,  "dP_pcsaft": 5.36, "T": [353.15]},
    ("cC6",   "PVAE50"):  {"Bba": 6.0434,  "uba0_R": -471.136, "kij": -0.0858, "dP_lgt": 12.13, "dP_pcsaft": 6.90, "T": [353.15]},
    ("cC6",   "PVAE70"):  {"Bba": 5.9136,  "uba0_R": -469.771, "kij": -0.0636, "dP_lgt": 11.04, "dP_pcsaft": 5.73, "T": [353.15]},

    ("BZN",   "PSBR4.1"): {"Bba": 93.240,  "uba0_R": -313.260, "kij": 0.0202,  "dP_lgt": 4.92,  "dP_pcsaft": 2.44, "T": [298.15, 308.15]},

    ("TOL",   "PSBR5"):   {"Bba": 63.981,  "uba0_R": -450.566, "kij": 0.0080,  "dP_lgt": 3.37,  "dP_pcsaft": 4.07, "T": [308.15]},
    ("TOL",   "PSBR23"):  {"Bba": 69.172,  "uba0_R": -448.891, "kij": 0.0091,  "dP_lgt": 3.31,  "dP_pcsaft": 5.76, "T": [308.15]},
    ("TOL",   "PSBR41"):  {"Bba": -52.602, "uba0_R": -677.198, "kij": 0.0098,  "dP_lgt": 7.59,  "dP_pcsaft": 3.71, "T": [343.15, 373.15]},

    ("nC6",   "PSBR23"):  {"Bba": 35.724,  "uba0_R": -395.208, "kij": 0.0086,  "dP_lgt": 12.21, "dP_pcsaft": 3.92, "T": [343.15]},
    ("nC6",   "PSBR45"):  {"Bba": 42.548,  "uba0_R": -394.415, "kij": 0.0089,  "dP_lgt": 11.33, "dP_pcsaft": 4.28, "T": [343.15]},

    ("nC5",   "PSBR23"):  {"Bba": 45.408,  "uba0_R": -336.736, "kij": 0.0072,  "dP_lgt": 6.84,  "dP_pcsaft": 3.68, "T": [333.15]},
    ("nC5",   "PSBR45"):  {"Bba": 42.886,  "uba0_R": -315.799, "kij": 0.0079,  "dP_lgt": 11.49, "dP_pcsaft": 6.24, "T": [333.15]},
}

GLOBAL_DP_LGT_M1 = 10.10      # desviacion global reportada LGT (metodologia 1)
GLOBAL_DP_PCSAFT_M1 = 4.41    # desviacion global reportada PC-SAFT (metodologia 1)

# ---------------------------------------------------------------------------
# Tabla 7: resultados PC-SAFT (metodologia 2: cada co-monomero por separado).
#   kij_a : solvente-monomero1
#   kij_b : solvente-monomero2
#   kiajb : monomero1-monomero2  (por copolimero concreto)
#   dP    : %DeltaP/P reportado por familia solvente+familia copolimero.
# ---------------------------------------------------------------------------
TABLE7 = {
    ("ACE", "PSBR"):   {"kij_a": 0.0121, "kij_b": 0.0028,  "kiajb": {"PSBR5": 0.0003, "PSBR23": 0.0004, "PSBR45": 0.0001}, "dP": 2.83},
    ("nC4OH", "PVAE"): {"kij_a": 0.0160, "kij_b": 0.0131,  "kiajb": {"PVAE9": -0.0010}, "dP": 1.88},
    ("nC4OH", "PVAVC"):{"kij_a": -0.0101,"kij_b": -0.0090, "kiajb": {"PVAVC12": 0.0008}, "dP": 1.63},
    ("CHL", "PVAE"):   {"kij_a": -0.0532,"kij_b": 0.0003,  "kiajb": {"PVAE9": 0.0008, "PVAE25": 0.0005, "PVAE50": 0.0003, "PVAE70": 0.0003}, "dP": 2.59},
    ("cC6", "PSBR"):   {"kij_a": 0.0052, "kij_b": 0.0012,  "kiajb": {"PSBR5": 0.0008, "PSBR23": 0.0006, "PSBR41": 0.0005, "PSBR45": 0.0002}, "dP": 2.11},
    ("cC6", "PVAE"):   {"kij_a": -0.1085,"kij_b": 0.0101,  "kiajb": {"PVAE25": 0.0011, "PVAE50": 0.0010, "PVAE70": 0.0008}, "dP": 2.09},
    ("BZN", "PSBR"):   {"kij_a": 0.0186, "kij_b": 0.0109,  "kiajb": {"PSBR41": 0.0004}, "dP": 1.23},
    ("TOL", "PSBR"):   {"kij_a": 0.0082, "kij_b": 0.0126,  "kiajb": {"PSBR5": 0.0014, "PSBR23": 0.0011, "PSBR41": 0.0006}, "dP": 2.13},
    ("nC6", "PSBR"):   {"kij_a": 0.0080, "kij_b": 0.0091,  "kiajb": {"PSBR23": 0.0010, "PSBR45": 0.0008}, "dP": 1.84},
    ("nC5", "PSBR"):   {"kij_a": 0.0066, "kij_b": 0.0108,  "kiajb": {"PSBR23": 0.0013, "PSBR45": 0.0008}, "dP": 2.07},
}

GLOBAL_DP_PCSAFT_M2 = 2.04    # desviacion global reportada PC-SAFT (metodologia 2)
