"""
VLE con la LGT EoS (metodologia 1: pseudo-monomero).

El copolimero se representa por un unico grupo "pseudo-monomero" cuyos
parametros (v_inf, Q, Baa, u0_R) se toman de la Tabla 5 (o se calculan con las
reglas de combinacion Eqs. 13-17).  El solvente es su propio grupo.

Para la mezcla hay dos grupos:  g0 = solvente,  g1 = pseudo-monomero.
  nu[0] = [1, 0]            (una molecula de solvente = 1 grupo solvente)
  nu[1] = [0, N_psm]        (una molecula de copolimero = N pseudo-monomeros)
donde N_psm = MW_copolimero / MW_pseudomonomero.

Energias de interaccion (Eq. 12):  u_ma/R = u0_R_ma (1 + B_ma/T).
  - solvente-solvente  y  psm-psm : parametros puros (Tabla 4a / Tabla 5).
  - solvente-psm (cruzada)        : Tabla 6 (Bba, uba0_R).
"""

import math
from data.parameters import LGT_PURE, PSEUDO_MONOMER, TABLE6, COPOLYMER_PCSAFT
from data.systems import (COPOLYMERS, monomer_mole_fractions,
                          monomer_weight_fractions)
from src.lgt import LGTMixture, V_STAR


def pseudo_monomer_lgt_params(copolymer):
    """Parametros LGT del pseudo-monomero. Usa Tabla 5 si existe; de lo
    contrario aplica reglas de combinacion (Eqs. 13-16)."""
    if copolymer in PSEUDO_MONOMER:
        return dict(PSEUDO_MONOMER[copolymer]["lgt"])
    # calcular con reglas de combinacion
    cp = COPOLYMERS[copolymer]
    x1, x2 = monomer_mole_fractions(copolymer)
    m1 = LGT_PURE[cp["mon1"]]
    m2 = LGT_PURE[cp["mon2"]]
    # Eq. 14-15: aritmeticas
    v_inf = x1 * m1["v_inf"] + x2 * m2["v_inf"]
    Q = x1 * m1["Q"] + x2 * m2["Q"]
    # Eq. 16: geometrica con signo (u0_R siempre negativo)
    u0 = -((x1 * abs(m1["u0_R"]) ** 0.5 + x2 * abs(m2["u0_R"]) ** 0.5)) ** 2
    # Eq. 17: Baa con signo
    def si(v):
        return 1.0 if v >= 0 else -1.0
    sm = si(m1["Baa"] if abs(m1["Baa"]) >= abs(m2["Baa"]) else m2["Baa"])
    Baa = sm * (si(m1["Baa"]) * x1 * abs(m1["Baa"]) ** 0.5
                + si(m2["Baa"]) * x2 * abs(m2["Baa"]) ** 0.5) ** 2
    return {"v_inf": v_inf, "Q": Q, "Baa": Baa, "u0_R": u0}


def pseudo_monomer_MW(copolymer):
    """MW del grupo pseudo-monomero (Eq. 13)."""
    cp = COPOLYMERS[copolymer]
    x1, x2 = monomer_mole_fractions(copolymer)
    return x1 * LGT_PURE[cp["mon1"]]["MW"] + x2 * LGT_PURE[cp["mon2"]]["MW"]


def build_lgt_mixture(solvent, copolymer, T):
    """Construye la LGTMixture binaria para (solvente, copolimero)."""
    sp = LGT_PURE[solvent]
    pp = pseudo_monomer_lgt_params(copolymer)

    # grupos
    def mkgroup(name, p):
        va = p["v_inf"]       # A=0 => v_a = v_inf (cm3/mol)
        Ra = va / V_STAR
        return {"name": name, "Ra": Ra, "Qa": p["Q"], "va": va,
                "u0_R": p["u0_R"], "Baa": p["Baa"]}

    g_solv = mkgroup(solvent, sp)
    g_psm = mkgroup("PSM", pp)
    groups = [g_solv, g_psm]

    # nu[i][a]
    MW_cop = COPOLYMER_PCSAFT[copolymer]["MW"]
    MW_psm = pseudo_monomer_MW(copolymer)
    N_psm = MW_cop / MW_psm
    nu = [[1.0, 0.0],        # solvente
          [0.0, N_psm]]      # copolimero

    # energias de interaccion cruzadas u_ma/R [K]  (Eq. 12)
    entry = TABLE6.get((solvent, copolymer))
    # diagonales: parametros puros
    u_ss = sp["u0_R"] * (1.0 + sp["Baa"] / T)
    u_pp = pp["u0_R"] * (1.0 + pp["Baa"] / T)
    if entry is not None:
        u_cross = entry["uba0_R"] * (1.0 + entry["Bba"] / T)
    else:
        # combinacion geometrica si no hay parametro ajustado
        u0c = -math.sqrt(abs(sp["u0_R"]) * abs(pp["u0_R"]))
        Bc = math.copysign(math.sqrt(abs(sp["Baa"]) * abs(pp["Baa"])), sp["Baa"])
        u_cross = u0c * (1.0 + Bc / T)

    cross_uR = [[u_ss, u_cross],
                [u_cross, u_pp]]

    return LGTMixture(groups, nu, cross_uR)


def _solve_vred_lgt(mix, x, T, P_mpa, phase="liq"):
    """Resuelve v~ (reduced volume) tal que P(v~)=P_mpa."""
    def f(vred):
        try:
            return mix.pressure(x, vred, T) - P_mpa
        except (ValueError, ZeroDivisionError):
            return None

    if phase == "liq":
        grid = [1.02 + 0.02 * i for i in range(0, 250)]   # 1.02 .. ~6
    else:
        grid = [50.0 * (1.08 ** i) for i in range(0, 120)]  # grandes v~ (gas)

    prev_v = prev_f = None
    for v in grid:
        fv = f(v)
        if fv is None:
            prev_v = prev_f = None
            continue
        if prev_f is not None and prev_f * fv < 0:
            lo, hi, flo = prev_v, v, prev_f
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                fmid = f(mid)
                if fmid is None:
                    hi = mid
                    continue
                if abs(fmid) < 1e-9 or abs(hi - lo) < 1e-12:
                    return mid
                if flo * fmid < 0:
                    hi = mid
                else:
                    lo, flo = mid, fmid
            return 0.5 * (lo + hi)
        prev_v, prev_f = v, fv
    return None


def bubble_pressure_lgt(solvent, copolymer, T, x1):
    """Presion de burbuja [MPa] con LGT para fraccion molar de solvente x1."""
    mix = build_lgt_mixture(solvent, copolymer, T)
    x = [x1, 1.0 - x1]

    # vapor = solvente puro: construimos una mezcla de un solo grupo
    sp = LGT_PURE[solvent]
    g_solv = {"name": solvent, "Ra": sp["v_inf"] / V_STAR, "Qa": sp["Q"],
              "va": sp["v_inf"], "u0_R": sp["u0_R"], "Baa": sp["Baa"]}
    u_ss = sp["u0_R"] * (1.0 + sp["Baa"] / T)
    mix_v = LGTMixture([g_solv], [[1.0]], [[u_ss]])

    def residual(P_mpa):
        vred_l = _solve_vred_lgt(mix, x, T, P_mpa, phase="liq")
        if vred_l is None:
            return None
        lnphi_l = mix.ln_phi(0, x, vred_l, T)
        vred_v = _solve_vred_lgt(mix_v, [1.0], T, P_mpa, phase="vap")
        if vred_v is None:
            return None
        lnphi_v = mix_v.ln_phi(0, [1.0], vred_v, T)
        return (math.log(x1) + lnphi_l) - lnphi_v

    import numpy as np
    P_grid = np.logspace(math.log10(1e-5), math.log10(2.0), 80)  # MPa
    prev_P = prev_f = None
    bracket = None
    for P in P_grid:
        fv = residual(float(P))
        if fv is None:
            prev_P = prev_f = None
            continue
        if prev_f is not None and prev_f * fv < 0:
            bracket = (prev_P, float(P), prev_f)
            break
        prev_P, prev_f = float(P), fv
    if bracket is None:
        return None
    lo, hi, flo = bracket
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fmid = residual(mid)
        if fmid is None:
            hi = mid
            continue
        if abs(fmid) < 1e-9 or abs(hi - lo) < 1e-9:
            return mid
        if flo * fmid < 0:
            hi = mid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)
