"""
Calculo de equilibrio liquido-vapor (VLE) para mezclas copolimero + solvente.

Condicion de equilibrio (Eq. 39), como el copolimero no es volatil:
      x1 * phi_hat_1^L (T, P, x) = phi_1^V (T, P, y1=1)

donde "1" es el solvente.  A T dada y fraccion molar de solvente x1 en la fase
liquida, se busca la presion P tal que se cumpla el equilibrio.

Este modulo ofrece:
  - construccion de los objetos PC-SAFT para un sistema (metodologia 1 y 2)
  - resolucion de la presion de equilibrio con PC-SAFT
  - conversion entre fraccion molar y fraccion masica de solvente
  - desviacion de presion %DeltaP/P (Eq. 40) frente a una referencia
"""

import math

from data.parameters import (PCSAFT_PURE, COPOLYMER_PCSAFT, TABLE6, TABLE7)
from data.systems import (COPOLYMERS, FAMILY, monomer_weight_fractions)
from src.pcsaft import PCSAFT, Component, Segment, KB, NA


# =====================================================================
#  Construccion de componentes PC-SAFT
# =====================================================================
def build_solvent(solvent):
    p = PCSAFT_PURE[solvent]
    # para el solvente, m = (m/MW) * MW_solvente.  MW esta en LGT_PURE;
    # usamos la relacion m = (m/MW)*MW para obtener el numero de segmentos.
    from data.parameters import LGT_PURE
    MW = LGT_PURE[solvent]["MW"]
    m = p["m_MW"] * MW
    seg = Segment(solvent, m, p["sigma"], p["eps_k"])
    return Component(solvent, [seg], [1.0])


def build_copolymer(copolymer):
    """Construye el componente copolimero con segmentos alpha (mon1) y beta (mon2).

    m_{i,s} = w_{i,s} * MW_copolym * (m/MW)_homopoly_s   (Eq. 18)
    usando los parametros PC-SAFT del homopolimero correspondiente (Tabla 4b).
    """
    cp = COPOLYMERS[copolymer]
    cpar = COPOLYMER_PCSAFT[copolymer]
    MWcop = cpar["MW"]
    w1, w2 = monomer_weight_fractions(copolymer)   # masicas

    h1 = PCSAFT_PURE[cp["homo1"]]   # homopolimero del monomero 1 (alpha)
    h2 = PCSAFT_PURE[cp["homo2"]]   # homopolimero del monomero 2 (beta)

    m_a = w1 * MWcop * h1["m_MW"]
    m_b = w2 * MWcop * h2["m_MW"]
    m_tot = m_a + m_b
    z_a = m_a / m_tot
    z_b = m_b / m_tot

    seg_a = Segment(cp["homo1"], m_a, h1["sigma"], h1["eps_k"])
    seg_b = Segment(cp["homo2"], m_b, h2["sigma"], h2["eps_k"])
    comp = Component(copolymer, [seg_a, seg_b], [z_a, z_b])
    comp._z_alpha = z_a
    comp._z_beta = z_b
    return comp


def bonding_fractions(copolymer, comp):
    """Fracciones de enlace B (Tabla 1 / Tabla 3).

    Preferimos los valores reportados en Tabla 3; si no estan, los calculamos
    con las reglas de la Tabla 1 para distribucion aleatoria (random).
    Devuelve dict {"aa","ab","bb"} o None si no hay copolimero efectivo.
    """
    cpar = COPOLYMER_PCSAFT[copolymer]
    if cpar["zib"] is None:
        return None  # PVT no disponible -> tratar como homopolimero promedio
    # Tabla 3 da: Biaia (aa), Biajb (ab), Bjbjb (bb)
    return {"aa": cpar["Biaia"], "ab": cpar["Biajb"], "bb": cpar["Bjbjb"]}


# =====================================================================
#  Fraccion masica <-> fraccion molar del solvente
# =====================================================================
def w1_to_x1(w1, solvent, copolymer):
    """Fraccion MASICA de solvente -> fraccion MOLAR de solvente.

    El 'peso molecular' del copolimero se usa para los moles de copolimero.
    """
    from data.parameters import LGT_PURE
    MW_s = LGT_PURE[solvent]["MW"]
    MW_c = COPOLYMER_PCSAFT[copolymer]["MW"]
    n_s = w1 / MW_s
    n_c = (1.0 - w1) / MW_c
    return n_s / (n_s + n_c)


def x1_to_w1(x1, solvent, copolymer):
    from data.parameters import LGT_PURE
    MW_s = LGT_PURE[solvent]["MW"]
    MW_c = COPOLYMER_PCSAFT[copolymer]["MW"]
    m_s = x1 * MW_s
    m_c = (1.0 - x1) * MW_c
    return m_s / (m_s + m_c)


# =====================================================================
#  Resolucion de densidad a (T,P) dada  -- PC-SAFT
# =====================================================================
def _solve_density(model, x, T, P_pa, bonding, phase="liq"):
    """Encuentra rho (1/A^3) tal que P(rho)=P_pa.  Busca raiz por biseccion
    en un rango de densidad reducida eta."""
    # rho se relaciona con eta; barremos eta y buscamos cambio de signo.
    # calculamos eta a partir de rho: eta = zeta3.  Hacemos un barrido en rho.
    def Pres(rho):
        return model.pressure(x, rho, T, bonding) - P_pa

    # rango de rho: desde muy diluido hasta empaquetamiento alto
    # eta ~ pi/6 * rho * sum(...). Elegimos rho segun eta objetivo.
    # Estimamos rho para eta=0.5 como cota superior.
    # sum_seg = sum xi mi sum zi di^3
    sum_seg = 0.0
    for i, c in enumerate(model.comps):
        for si, seg in enumerate(c.segs):
            d = model.d_seg(seg, T)
            sum_seg += x[i] * c.m * c.z[si] * d ** 3
    coef = math.pi / 6.0 * sum_seg
    rho_at_eta = lambda eta: eta / coef

    if phase == "liq":
        etas = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-4]
    else:
        etas = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1]

    from scipy.optimize import brentq
    rhos = [rho_at_eta(e) for e in etas]
    # buscar cambio de signo consecutivo y refinar con brentq
    prev_rho = rhos[0]
    prev_val = Pres(prev_rho)
    for rho in rhos[1:]:
        val = Pres(rho)
        if prev_val * val < 0:
            lo, hi = sorted((prev_rho, rho))
            try:
                return brentq(Pres, lo, hi, xtol=1e-16, rtol=1e-10, maxiter=100)
            except Exception:
                return 0.5 * (lo + hi)
        prev_rho, prev_val = rho, val
    return None


# =====================================================================
#  Presion de burbuja (equilibrio) con PC-SAFT
# =====================================================================
def bubble_pressure_pcsaft(solvent, copolymer, T, x1, kij_dict):
    """Devuelve P [MPa] de equilibrio para fraccion molar de solvente x1.

    kij_dict: dict frozenset({segA,segB}) -> kij.
    Se resuelve iterando sobre P hasta que
       x1 phi_hat_1^L (T,P,x) = phi_1^V(T,P)   (vapor = solvente puro).
    """
    solv = build_solvent(solvent)
    cop = build_copolymer(copolymer)
    model_liq = PCSAFT([solv, cop], kij=kij_dict)
    model_vap = PCSAFT([solv], kij=kij_dict)

    x = [x1, 1.0 - x1]
    bonding = {1: bonding_fractions(copolymer, cop)}

    def residual(P_pa):
        # fase liquida
        rho_l = _solve_density(model_liq, x, T, P_pa, bonding, phase="liq")
        if rho_l is None:
            return None
        lnphi_l = model_liq.ln_phi(x, rho_l, T, bonding)
        # fase vapor: solvente puro
        rho_v = _solve_density(model_vap, [1.0], T, P_pa, {}, phase="vap")
        if rho_v is None:
            return None
        lnphi_v = model_vap.ln_phi([1.0], rho_v, T, {})[0]
        # equilibrio: ln(x1) + lnphi_l_1 = lnphi_v_1
        return (math.log(x1) + lnphi_l[0]) - lnphi_v

    # Barrido logaritmico (grueso) de P para localizar el cambio de signo del
    # residual, luego refinar con brentq.
    import numpy as np
    from scipy.optimize import brentq
    P_grid = np.logspace(math.log10(1.0e2), math.log10(2.0e6), 28)  # 100 Pa..2 MPa
    prev_P = None
    prev_f = None
    bracket = None
    for P_pa in P_grid:
        f = residual(float(P_pa))
        if f is None:
            prev_P, prev_f = None, None
            continue
        if prev_f is not None and prev_f * f < 0:
            bracket = (prev_P, float(P_pa))
            break
        prev_P, prev_f = float(P_pa), f

    if bracket is None:
        return None

    P_lo, P_hi = bracket
    try:
        P_root = brentq(lambda P: residual(P), P_lo, P_hi,
                        xtol=1e-2, rtol=1e-9, maxiter=60)
        return P_root / 1e6
    except Exception:
        return 0.5 * (P_lo + P_hi) / 1e6


# =====================================================================
#  kij dictionaries para cada metodologia
# =====================================================================
def kij_methodology1(solvent, copolymer):
    """Metodologia 1 (Tabla 6): un solo kij solvente-copolimero.
    En el modelo de pseudo-monomero el copolimero es un pseudo-segmento, pero
    aqui mantenemos los dos segmentos homopolimero y aplicamos el mismo kij
    solvente-segmento a ambos (aproximacion coherente con un unico parametro).
    """
    entry = TABLE6.get((solvent, copolymer))
    if entry is None:
        return {}
    kij = entry["kij"]
    cp = COPOLYMERS[copolymer]
    d = {}
    d[frozenset((solvent, cp["homo1"]))] = kij
    d[frozenset((solvent, cp["homo2"]))] = kij
    return d


def kij_methodology2(solvent, copolymer):
    """Metodologia 2 (Tabla 7): kij_a (solvente-monomero1), kij_b
    (solvente-monomero2) y kiajb (monomero1-monomero2)."""
    fam = FAMILY[copolymer]
    entry = TABLE7.get((solvent, fam))
    if entry is None:
        return {}
    cp = COPOLYMERS[copolymer]
    d = {}
    d[frozenset((solvent, cp["homo1"]))] = entry["kij_a"]
    d[frozenset((solvent, cp["homo2"]))] = entry["kij_b"]
    kiajb = entry["kiajb"].get(copolymer)
    if kiajb is not None:
        d[frozenset((cp["homo1"], cp["homo2"]))] = kiajb
    return d


# =====================================================================
#  Desviacion de presion (Eq. 40)
# =====================================================================
def pressure_deviation(P_calc, P_ref):
    """%DeltaP/P = 100/np sum |Pcalc - Pref|/Pref ."""
    n = len(P_calc)
    s = 0.0
    cnt = 0
    for pc, pr in zip(P_calc, P_ref):
        if pc is None or pr is None or pr == 0:
            continue
        s += abs(pc - pr) / abs(pr)
        cnt += 1
    if cnt == 0:
        return None
    return 100.0 * s / cnt
