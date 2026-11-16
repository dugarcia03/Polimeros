"""
Parameters taken EXCLUSIVELY from Arce, P.F., Mattedi, S., Aznar, M. (2006),
"Vapor-liquid equilibrium of copolymer + solvent mixtures: Thermodynamic
modeling by two theoretical equations of state", Fluid Phase Equilibria 246,
52-63. DOI: 10.1016/j.fluid.2006.05.009

Only PC-SAFT-relevant data are kept here. No values from any other publication
are mixed in (the universal a_i/b_i live in constantes_universales.py; any
external validation data live in datos_externos.py).

Transcribed tables:
  - Table 4b : PC-SAFT pure-component parameters (solvents, monomers, polymers).
  - Table 2  : systems, temperatures, w1 ranges (physical characteristics).
  - Table 3  : copolymer molecular weights and comonomer compositions (wt% of
               the FIRST monomer, i.e. the one named in the copolymer code).
  - Table 6  : binary interaction parameter kij for Methodology 1 (pseudo-monomer).
  - Table 5  : pseudo-monomer parameters (reproduced in-code from Table 4b; the
               published Table 5 values are kept in comments for verification).

Units: m/MW in mol/g ; sigma in Angstrom ; eps/k in K ; MW in g/mol ; T in K.
"""

# ---------------------------------------------------------------------------
# Table 4b - PC-SAFT pure-component parameters: SOLVENTS
# key -> (m_over_MW [mol/g], sigma [A], eps_over_k [K], MW [g/mol])
# MW of solvents taken from Table 4a (same compounds, LGT side) / standard.
# ---------------------------------------------------------------------------
SOLVENTS = {
    # code : (m/MW,     sigma,   eps/k,   MW)
    "nC5":  (0.034090, 3.8729, 240.11, 72.15),
    "nC6":  (0.034088, 3.8383, 239.57, 86.18),
    "cC6":  (0.029946, 3.8399, 279.21, 84.16),
    "BZN":  (0.032712, 3.5878, 281.55, 78.12),
    "TOL":  (0.030984, 3.6869, 282.69, 92.14),
    "CHL":  (0.021826, 3.4085, 265.94, 119.38),
    "ACE":  (0.037879, 3.4956, 272.26, 58.08),
    # nC4OH is associating; associating params (kAB=0.008752, epsAB=2489.23 K)
    # are NOT used here because Methodology 1 for PC-SAFT in the paper uses the
    # non-associating pseudo-monomer treatment for the copolymer side. The
    # associating contribution for 1-butanol is handled only if explicitly
    # enabled; see datos_externos / notes. Base 3 params below.
    "nC4OH": (0.035126, 3.7563, 263.18, 74.12),
}

# Associating parameters for nC4OH (Table 4b). Kept separate: only relevant
# when the association term (Eqs. 36-38) is activated.
SOLVENT_ASSOC = {
    "nC4OH": {"kAB": 0.008752, "epsAB_k": 2489.23, "scheme": "2B"},
}

# ---------------------------------------------------------------------------
# Table 4b - PC-SAFT pure-component parameters: MONOMERS
# key -> (m/MW [mol/g], sigma [A], eps/k [K], MW [g/mol])
# Monomer molar masses from Table 4a.
# ---------------------------------------------------------------------------
MONOMERS = {
    "STR": (0.024339, 3.9750, 318.70, 104.14),  # styrene
    "VA":  (0.011525, 5.0706, 421.42, 86.09),   # vinyl acetate
    "BTD": (0.046033, 3.4508, 215.10, 54.10),   # butadiene
    "ET":  (0.056070, 3.4250, 178.77, 28.10),   # ethylene
    "VC":  (0.091211, 2.4024, 160.05, 62.50),   # vinyl chloride
}

# ---------------------------------------------------------------------------
# Table 4b - PC-SAFT pure-component parameters: POLYMERS (homopolymers)
# key -> (m/MW [mol/g], sigma [A], eps/k [K])
# (Provided for reference / Methodology-2 style work; Methodology 1 builds the
#  copolymer from the MONOMER parameters via the pseudo-monomer rule.)
# ---------------------------------------------------------------------------
POLYMERS = {
    "PS":  (0.033238, 3.5022, 320.14),
    "PVA": (0.029905, 3.5086, 310.14),
    "PB":  (0.060552, 2.9845, 240.21),
    "PE":  (0.028193, 4.0125, 320.24),
    "PVC": (0.021369, 3.7486, 370.42),
}

# ---------------------------------------------------------------------------
# Copolymer definitions.
#   monomers   : (first_monomer, second_monomer)  matching the copolymer code
#   wt_first   : weight fraction of the FIRST monomer (Table 2/3 "number is the
#                weight percent of the first monomer-unit").
#   MW         : copolymer molecular weight [g/mol] (Table 3). 100000 assumed
#                where Table 3 lists 100000 (paper's convention for unavailable).
# ---------------------------------------------------------------------------
COPOLYMERS = {
    # poly(styrene-co-butadiene): first monomer = styrene
    "PSBR4.1": {"monomers": ("STR", "BTD"), "wt_first": 0.041, "MW": 100000.0},
    "PSBR5":   {"monomers": ("STR", "BTD"), "wt_first": 0.05,  "MW": 100000.0},
    "PSBR23":  {"monomers": ("STR", "BTD"), "wt_first": 0.23,  "MW": 720000.0},
    "PSBR41":  {"monomers": ("STR", "BTD"), "wt_first": 0.41,  "MW": 100000.0},
    "PSBR45":  {"monomers": ("STR", "BTD"), "wt_first": 0.45,  "MW": 600000.0},
    # poly(vinyl acetate-co-ethylene): first monomer = vinyl acetate
    "PVAE9":   {"monomers": ("VA", "ET"), "wt_first": 0.09, "MW": 100000.0},
    "PVAE25":  {"monomers": ("VA", "ET"), "wt_first": 0.25, "MW": 88000.0},
    "PVAE50":  {"monomers": ("VA", "ET"), "wt_first": 0.50, "MW": 250000.0},
    "PVAE70":  {"monomers": ("VA", "ET"), "wt_first": 0.70, "MW": 280000.0},
    # poly(vinyl acetate-co-vinyl chloride): first monomer = vinyl acetate
    "PVAVC12": {"monomers": ("VA", "VC"), "wt_first": 0.12, "MW": 30000.0},
}

# ---------------------------------------------------------------------------
# Table 6 - Methodology 1 (pseudo-monomer) binary interaction parameter kij.
# key -> {copolymer: kij}.  Also stores the published %dP/P for comparison.
# ---------------------------------------------------------------------------
KIJ_M1 = {
    "ACE":   {"PSBR5": 0.0120, "PSBR23": 0.0122, "PSBR45": 0.0130},
    "nC4OH": {"PVAE9": 0.0230, "PVAVC12": -0.0126},
    "CHL":   {"PVAE9": -0.1208, "PVAE25": -0.0816, "PVAE50": -0.0701, "PVAE70": -0.0516},
    # cyclohexane appears twice in Table 6 (PSBR block and PVAE block)
    "cC6":   {"PSBR5": 0.0118, "PSBR23": 0.0104, "PSBR41": 0.0088, "PSBR45": 0.0101,
              "PVAE25": -0.1068, "PVAE50": -0.0858, "PVAE70": -0.0636},
    "BZN":   {"PSBR4.1": 0.0202},
    "TOL":   {"PSBR5": 0.0080, "PSBR23": 0.0091, "PSBR41": 0.0098},
    "nC6":   {"PSBR23": 0.0086, "PSBR45": 0.0089},
    "nC5":   {"PSBR23": 0.0072, "PSBR45": 0.0079},
}

# Published Methodology-1 %dP/P (Table 6, PC-SAFT column) for comparison.
DP_M1_PUBLISHED = {
    "ACE":   {"PSBR5": 4.78, "PSBR23": 3.04, "PSBR45": 5.19},
    "nC4OH": {"PVAE9": 4.05, "PVAVC12": 5.18},
    "CHL":   {"PVAE9": 4.18, "PVAE25": 5.06, "PVAE50": 4.41, "PVAE70": 3.28},
    "cC6":   {"PSBR5": 3.66, "PSBR23": 4.38, "PSBR41": 2.28, "PSBR45": 4.31,
              "PVAE25": 5.36, "PVAE50": 6.90, "PVAE70": 5.73},
    "BZN":   {"PSBR4.1": 2.44},
    "TOL":   {"PSBR5": 4.07, "PSBR23": 5.76, "PSBR41": 3.71},
    "nC6":   {"PSBR23": 3.92, "PSBR45": 4.28},
    "nC5":   {"PSBR23": 3.68, "PSBR45": 6.24},
}
DP_M1_GLOBAL_PUBLISHED = 4.41  # Table 6 global PC-SAFT pressure deviation

# Published Methodology-2 %dP/P (Table 7, PC-SAFT) per solvent family, for the
# final comparison table. Keyed by (solvent) -> dP (one value per family).
DP_M2_PUBLISHED = {
    "ACE": 2.83, "nC4OH_PVAE9": 1.88, "nC4OH_PVAVC12": 1.63,
    "CHL": 2.59, "cC6_PSBR": 2.11, "cC6_PVAE": 2.09, "BZN": 1.23,
    "TOL": 2.13, "nC6": 1.84, "nC5": 2.07,
}
DP_M2_GLOBAL_PUBLISHED = 2.04

# ---------------------------------------------------------------------------
# Table 2 - systems with temperature(s). Each entry lists the (solvent,
# copolymer) pair and the temperature(s) [K] at which VLE data were modeled.
# This drives the per-system execution/plotting workflow.
# ---------------------------------------------------------------------------
SYSTEMS = [
    # (solvent, copolymer, [temperatures K], experimental data source tag)
    ("ACE", "PSBR5",   [323.15],          "ref1"),
    ("ACE", "PSBR23",  [323.15, 333.15],  "ref1"),
    ("ACE", "PSBR45",  [323.15, 333.15],  "ref1"),
    ("nC4OH", "PVAE9",   [353.15],        "ref1"),
    ("nC4OH", "PVAVC12", [353.15],        "ref1"),
    ("CHL", "PVAE9",   [333.15],          "ref1"),
    ("CHL", "PVAE25",  [333.15],          "ref1"),
    ("CHL", "PVAE50",  [333.15],          "ref1"),
    ("CHL", "PVAE70",  [333.15],          "ref1"),
    ("cC6", "PSBR5",   [333.15],          "ref1"),
    ("cC6", "PSBR23",  [296.65, 333.15],  "ref1"),
    ("cC6", "PSBR41",  [343.15],          "ref34"),
    ("cC6", "PSBR45",  [296.65, 333.15],  "ref1"),
    ("cC6", "PVAE25",  [353.15],          "ref1"),
    ("cC6", "PVAE50",  [353.15],          "ref1"),
    ("cC6", "PVAE70",  [353.15],          "ref1"),
    ("BZN", "PSBR4.1", [298.15, 308.15],  "ref35"),
    ("TOL", "PSBR5",   [308.15],          "ref1"),
    ("TOL", "PSBR23",  [308.15],          "ref1"),
    ("TOL", "PSBR41",  [343.15, 373.15],  "ref34"),
    ("nC6", "PSBR23",  [343.15],          "ref1"),
    ("nC6", "PSBR45",  [343.15],          "ref1"),
    ("nC5", "PSBR23",  [333.15],          "ref1"),
    ("nC5", "PSBR45",  [333.15],          "ref1"),
]
