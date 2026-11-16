#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
  VLE de mezclas COPOLIMERO + SOLVENTE con la EoS PC-SAFT
  Reproduccion FIEL del articulo -- METODOLOGIA 2 (co-monomeros separados)
===============================================================================

    Arce, Mattedi & Aznar (2006), Fluid Phase Equilibria 246, 52-63.
    "Vapor-liquid equilibrium of copolymer + solvent mixtures ..."

Este codigo implementa la SEGUNDA metodologia del paper (seccion 5, Tabla 7),
que el propio articulo identifica como la mejor (desviacion global 2.04%):

  - El copolimero se trata como DOS tipos de segmento, alpha y beta (los dos
    monomeros), NO como un unico pseudo-monomero.
  - Cada segmento usa los parametros PC-SAFT del MONOMERO correspondiente
    (Tabla 4b: STR, BTD, VA, ET, VC).
  - Fracciones de segmento z_alpha, z_beta (Ec. 20) a partir de la composicion.
  - Fracciones de enlace B_{a,a}, B_{a,b}, B_{b,b} (Tabla 3 / reglas Tabla 1),
    que entran en el termino de cadena dura (Ec. 22).
  - TRES parametros de interaccion binaria (Tabla 7):
      kij_a  = solvente - monomero alpha
      kij_b  = solvente - monomero beta
      ki_ab  = monomero alpha - monomero beta  (dentro del copolimero)
  - Condicion de equilibrio (Ec. 39): vapor = solvente puro.
  - Constantes universales aik, bik importadas de constantes_pcsaft.py
    (fuente externa: Gross & Sadowski 2001, ref. [23] del paper).

Dependencias:  numpy, scipy, matplotlib
Uso:           python pcsaft_vle.py
===============================================================================
"""

import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

import math
import numpy as np
from scipy.optimize import brentq

# Constantes universales aik, bik (fuente externa citada por el paper, ref [23])
from constantes_pcsaft import A_COEF, B_COEF


# =============================================================================
#  ////////////////////  PARAMETROS EDITABLES  //////////////////////////////
#  Valores por defecto:  cC6 (ciclohexano) + PSBR45 a 333.15 K
#  Solvente de Tabla 4b; monomeros de Tabla 4b; kij de Tabla 7 (metodologia 2).
# =============================================================================

# ---- 1) TEMPERATURA y rango de composicion --------------------------------
TEMPERATURA_K = 333.15            # Temperatura [K]
W1_MIN        = 0.011             # Fraccion masica de solvente minima
W1_MAX        = 0.301             # Fraccion masica de solvente maxima
N_PUNTOS      = 60                # Numero de puntos de la curva estimada

# ---- 2) SOLVENTE (componente 1)  --  Tabla 4b -----------------------------
#   m/MW [1/mol], sigma [Angstrom], eps/k [K], MW [g/mol].
SOLVENTE = {
    "nombre": "cC6 (ciclohexano)",
    "MW":     84.16,
    "m_MW":   0.029946,
    "sigma":  3.8399,
    "eps_k":  279.21,
}

# ---- 3) COPOLIMERO (componente 2) como DOS MONOMEROS  --  Tabla 4b + Tabla 3
#   alpha = monomero 1 (el del % del nombre),  beta = monomero 2.
#   Cada monomero con sus parametros PC-SAFT de la Tabla 4b.
#   w_alpha = fraccion MASICA del monomero alpha (= numero del nombre / 100).
#   MW = peso molecular total del copolimero (Tabla 3).
COPOLIMERO = {
    "nombre":  "PSBR45 (poly(styrene-co-butadiene), 45% estireno)",
    "MW":      600000.0,          # g/mol (peso molecular total, Tabla 3)
    "w_alpha": 0.45,              # fraccion masica del monomero alpha (estireno)

    # --- monomero ALPHA (aqui: STR, estireno) -- Tabla 4b ---
    "alpha": {
        "nombre": "STR",
        "MW":     104.14,         # g/mol del monomero
        "m_MW":   0.024339,       # (m/MW) [1/mol]
        "sigma":  3.9750,         # Angstrom
        "eps_k":  318.70,         # eps/k [K]
    },
    # --- monomero BETA (aqui: BTD, butadieno) -- Tabla 4b ---
    "beta": {
        "nombre": "BTD",
        "MW":     54.10,          # g/mol del monomero
        "m_MW":   0.046033,       # (m/MW) [1/mol]
        "sigma":  3.4508,         # Angstrom
        "eps_k":  215.10,         # eps/k [K]
    },

    # --- Fracciones de enlace B (Tabla 3 del paper).  Pon None para que el
    #     codigo las calcule con las reglas de la Tabla 1 (aleatorio).
    #     Tabla 3 para PSBR45: Biaia=0.3707, Biajb=0.6293, Bjbjb=0.0000.
    #     (el paper define alpha=estireno=el monomero del nombre)
    "bonding": {"aa": 0.3707, "ab": 0.6293, "bb": 0.0000},
}

# ---- 4) PARAMETROS DE INTERACCION BINARIA (Tabla 7, metodologia 2) ---------
#   Tres parametros:
#     KIJ_A  = solvente - monomero alpha
#     KIJ_B  = solvente - monomero beta
#     KI_AB  = monomero alpha - monomero beta (dentro del copolimero)
#   Valores por defecto: cC6+PSBR, familia PSBR (Tabla 7).
KIJ_A = 0.0052                    # cC6-STR   (Tabla 7, kij_a familia PSBR)
KIJ_B = 0.0012                    # cC6-BTD   (Tabla 7, kij_b familia PSBR)
KI_AB = 0.0002                    # STR-BTD   (Tabla 7, kiajb de PSBR45)

# ---- 5) DATOS EXPERIMENTALES (opcional) -----------------------------------
#   (w1, P[MPa]).  cC6+PSBR45 a 333.15 K (Gupta & Prausnitz [1]).
DATOS_EXPERIMENTALES = [
    (0.021, 0.0048),
    (0.050, 0.0105),
    (0.087, 0.0169),
    (0.103, 0.0199),
    (0.130, 0.0248),
    (0.163, 0.0295),
    (0.200, 0.0337),
    (0.240, 0.0376),
    (0.282, 0.0407),
]

# ---- 6) SALIDA ------------------------------------------------------------
ARCHIVO_FIGURA = "resultado_VLE.png"

# =============================================================================
#  ///////////////  FIN DE PARAMETROS EDITABLES  ////////////////////////////
#  Debajo de aqui no hay parametros del sistema; solo el motor PC-SAFT.
#  Las constantes universales aik/bik estan en constantes_pcsaft.py.
# =============================================================================

KB = 1.380649e-23        # J/K (constante de Boltzmann)


# -----------------------------------------------------------------------------
#  Segmento (un tipo): m = numero de segmentos aportado, sigma, eps/k
# -----------------------------------------------------------------------------
class Segmento:
    def __init__(self, nombre, m, sigma, eps_k):
        self.nombre = nombre
        self.m = m
        self.sigma = sigma
        self.eps_k = eps_k


# -----------------------------------------------------------------------------
#  Componente: 1 segmento (solvente) o 2 segmentos (copolimero alpha+beta).
#   z = fracciones de segmento (suman 1). m_tot = suma de m de sus segmentos.
#   bonding = dict {"aa","ab","bb"} con las fracciones de enlace (solo copol.)
# -----------------------------------------------------------------------------
class Componente:
    def __init__(self, nombre, segmentos, z, bonding=None):
        self.nombre = nombre
        self.segs = segmentos
        self.z = z
        self.m = sum(s.m for s in segmentos)
        self.bonding = bonding


# -----------------------------------------------------------------------------
#  PC-SAFT para copolimeros (Ecs. 21-35 del paper, metodologia de 2 segmentos).
# -----------------------------------------------------------------------------
class PCSAFT:
    def __init__(self, componentes, kij):
        """kij: dict con clave frozenset({nombreSeg1,nombreSeg2}) -> valor."""
        self.c = componentes
        self.nc = len(componentes)
        self.kij = kij or {}

    # --- lista plana de (componente_idx, segmento) para recorrer segmentos ---
    def _segmentos(self):
        for i, comp in enumerate(self.c):
            for si, seg in enumerate(comp.segs):
                yield i, si, seg

    def _kij(self, na, nb):
        if na == nb:
            return 0.0
        return self.kij.get(frozenset((na, nb)), 0.0)

    def d_seg(self, seg, T):
        """Diametro de segmento dependiente de T (Ec. 26)."""
        return seg.sigma * (1.0 - 0.12 * math.exp(-3.0 * seg.eps_k / T))

    def m_bar(self, x):
        """Ec. 24: m promedio = sum_i x_i m_i."""
        return sum(x[i] * self.c[i].m for i in range(self.nc))

    def zeta(self, x, rho, T):
        """Ec. 25: zeta_k = pi/6 rho sum_i x_i m_i sum_s z_is d_is^k."""
        z = [0.0, 0.0, 0.0, 0.0]
        for k in range(4):
            s = 0.0
            for i, comp in enumerate(self.c):
                for si, seg in enumerate(comp.segs):
                    s += x[i] * comp.m * comp.z[si] * self.d_seg(seg, T) ** k
            z[k] = math.pi / 6.0 * rho * s
        return z

    # ---- esfera dura ----
    def a_hs(self, z):
        z0, z1, z2, z3 = z
        return (1.0 / z0) * (3.0 * z1 * z2 / (1.0 - z3)
                             + z2 ** 3 / (z3 * (1.0 - z3) ** 2)
                             + (z2 ** 3 / z3 ** 2 - z0) * math.log(1.0 - z3))

    def g_hs(self, di, dj, z):
        """Funcion de distribucion radial (Ec. 23)."""
        z2, z3 = z[2], z[3]
        dd = di * dj / (di + dj)
        return (1.0 / (1.0 - z3)
                + dd * 3.0 * z2 / (1.0 - z3) ** 2
                + dd ** 2 * 2.0 * z2 ** 2 / (1.0 - z3) ** 3)

    # ---- cadena dura (Ec. 22) con fracciones de enlace B ----
    def a_hc(self, x, rho, T):
        z = self.zeta(x, rho, T)
        a = self.m_bar(x) * self.a_hs(z)
        for i, comp in enumerate(self.c):
            if len(comp.segs) == 1:            # solvente: cadena homonuclear
                d = self.d_seg(comp.segs[0], T)
                a -= x[i] * (comp.m - 1.0) * math.log(self.g_hs(d, d, z))
            else:                              # copolimero: dos segmentos a,b
                da = self.d_seg(comp.segs[0], T)
                db = self.d_seg(comp.segs[1], T)
                B = comp.bonding
                if B is None:
                    g = self.g_hs(da, db, z)
                    a -= x[i] * (comp.m - 1.0) * math.log(g)
                else:
                    # Ec. 22: sum sobre pares alpha-beta ponderado por B
                    term = (B["aa"] * math.log(self.g_hs(da, da, z))
                            + B["ab"] * math.log(self.g_hs(da, db, z))
                            + B["bb"] * math.log(self.g_hs(db, db, z)))
                    a -= x[i] * (comp.m - 1.0) * term
        return a

    # ---- dispersion (Ecs. 27-32) ----
    def _ab(self, mb):
        a, b = [], []
        for k in range(7):
            a.append(A_COEF[k][0] + (mb - 1.0) / mb * A_COEF[k][1]
                     + (mb - 1.0) / mb * (mb - 2.0) / mb * A_COEF[k][2])
            b.append(B_COEF[k][0] + (mb - 1.0) / mb * B_COEF[k][1]
                     + (mb - 1.0) / mb * (mb - 2.0) / mb * B_COEF[k][2])
        return a, b

    def _I1_I2(self, mb, eta):
        a, b = self._ab(mb)
        I1 = sum(a[k] * eta ** k for k in range(7))
        I2 = sum(b[k] * eta ** k for k in range(7))
        return I1, I2

    def _C1(self, mb, eta):
        return 1.0 / (1.0
                      + mb * (8.0 * eta - 2.0 * eta ** 2) / (1.0 - eta) ** 4
                      + (1.0 - mb) * (20.0 * eta - 27.0 * eta ** 2
                                      + 12.0 * eta ** 3 - 2.0 * eta ** 4)
                      / ((1.0 - eta) * (2.0 - eta)) ** 2)

    def _eps_sigma(self, sa, sb):
        """Ec. 33: reglas de combinacion con kij entre los dos segmentos."""
        sig = 0.5 * (sa.sigma + sb.sigma)
        eps = math.sqrt(sa.eps_k * sb.eps_k) * (1.0 - self._kij(sa.nombre, sb.nombre))
        return eps, sig

    def _m2es(self, x, T, power):
        """Doble suma sobre segmentos (Ecs. 28-29):
        sum_i sum_j xi xj mi mj sum_a sum_b zia zjb (eps_ab/kT)^p sigma_ab^3."""
        s = 0.0
        for i, ci in enumerate(self.c):
            for j, cj in enumerate(self.c):
                for ai, sa in enumerate(ci.segs):
                    for bj, sb in enumerate(cj.segs):
                        eps, sig = self._eps_sigma(sa, sb)
                        s += (x[i] * x[j] * ci.m * cj.m
                              * ci.z[ai] * cj.z[bj]
                              * (eps / T) ** power * sig ** 3)
        return s

    def a_disp(self, x, rho, T):
        mb = self.m_bar(x)
        eta = self.zeta(x, rho, T)[3]
        I1, I2 = self._I1_I2(mb, eta)
        C1 = self._C1(mb, eta)
        a1 = -2.0 * math.pi * rho * I1 * self._m2es(x, T, 1)       # Ec. 28
        a2 = -math.pi * rho * mb * C1 * I2 * self._m2es(x, T, 2)   # Ec. 29
        return a1 + a2

    # ---- residual, Z, presion, fugacidad ----
    def a_res(self, x, rho, T):
        return self.a_hc(x, rho, T) + self.a_disp(x, rho, T)       # Ec. 21

    def Z(self, x, rho, T):
        drho = 1e-6 * rho
        return 1.0 + rho * (self.a_res(x, rho + drho, T)
                            - self.a_res(x, rho - drho, T)) / (2.0 * drho)

    def pressure(self, x, rho, T):
        """P en Pa.  rho en 1/A^3."""
        return self.Z(x, rho, T) * (rho * 1e30) * KB * T

    def ln_phi(self, x, rho, T):
        """Coeficiente de fugacidad (Ecs. 34-35)."""
        ares = self.a_res(x, rho, T)
        Z = self.Z(x, rho, T)
        lnZ = math.log(Z)
        h = 1e-7
        dadx = []
        for i in range(self.nc):
            xp = list(x); xp[i] += h
            xm = list(x); xm[i] -= h
            dadx.append((self.a_res(xp, rho, T) - self.a_res(xm, rho, T)) / (2.0 * h))
        sdx = sum(x[j] * dadx[j] for j in range(self.nc))
        return [ares + (Z - 1.0) + dadx[i] - sdx - lnZ for i in range(self.nc)]


# -----------------------------------------------------------------------------
#  Construccion de componentes a partir del bloque de parametros
# -----------------------------------------------------------------------------
def construir_solvente():
    s = SOLVENTE
    m = s["m_MW"] * s["MW"]
    seg = Segmento("SOLV", m, s["sigma"], s["eps_k"])
    return Componente(s["nombre"], [seg], [1.0])


def _fracciones_molares_monomero():
    """Fraccion MOLAR de cada monomero (para fracciones de segmento)."""
    c = COPOLIMERO
    wA = c["w_alpha"]; wB = 1.0 - wA
    nA = wA / c["alpha"]["MW"]; nB = wB / c["beta"]["MW"]
    return nA / (nA + nB), nB / (nA + nB)


def construir_copolimero():
    c = COPOLIMERO
    MWc = c["MW"]
    wA = c["w_alpha"]; wB = 1.0 - wA
    # Ec. 18: m_{i,s} = w_{i,s} * MW_copol * (m/MW)_s   (monomeros)
    mA = wA * MWc * c["alpha"]["m_MW"]
    mB = wB * MWc * c["beta"]["m_MW"]
    m_tot = mA + mB
    zA = mA / m_tot                           # Ec. 20: fraccion de segmento
    zB = mB / m_tot
    segA = Segmento(c["alpha"]["nombre"], mA, c["alpha"]["sigma"], c["alpha"]["eps_k"])
    segB = Segmento(c["beta"]["nombre"],  mB, c["beta"]["sigma"],  c["beta"]["eps_k"])
    # Fracciones de enlace: usar las de la Tabla 3 si se dieron; si no, calcular.
    bonding = c.get("bonding")
    if bonding is None:
        bonding = _fracciones_enlace(zA, zB, m_tot)
    comp = Componente(c["nombre"], [segA, segB], [zA, zB], bonding)
    comp._zA, comp._zB = zA, zB
    return comp


def _fracciones_enlace(zA, zB, m):
    """Fracciones de enlace B para copolimero aleatorio (Tabla 1 del paper).
    alpha = monomero 1, beta = monomero 2.  (z_beta vs z_alpha deciden el caso)."""
    if m <= 1:
        return None
    if zB < zA:
        B_ab = 2.0 * (zB * m) / (m - 1.0)
        B_aa = 1.0 - B_ab
        B_bb = 0.0
    elif zB > zA:
        B_ab = 2.0 * (zA * m) / (m - 1.0)
        B_bb = 1.0 - B_ab
        B_aa = 0.0
    else:                                     # alternante
        B_ab = 1.0; B_aa = 0.0; B_bb = 0.0
    B_ab = min(max(B_ab, 0.0), 1.0)
    return {"aa": max(B_aa, 0.0), "ab": B_ab, "bb": max(B_bb, 0.0)}


def construir_kij():
    """Diccionario de kij por par de segmentos (Tabla 7):
       solvente-alpha, solvente-beta, alpha-beta."""
    solv = "SOLV"
    na = COPOLIMERO["alpha"]["nombre"]
    nb = COPOLIMERO["beta"]["nombre"]
    d = {}
    d[frozenset((solv, na))] = KIJ_A
    d[frozenset((solv, nb))] = KIJ_B
    if na != nb:
        d[frozenset((na, nb))] = KI_AB
    return d


# -----------------------------------------------------------------------------
#  Fraccion masica de solvente -> fraccion molar de solvente
# -----------------------------------------------------------------------------
def w1_a_x1(w1):
    MW_s = SOLVENTE["MW"]
    MW_c = COPOLIMERO["MW"]
    n_s = w1 / MW_s
    n_c = (1.0 - w1) / MW_c
    return n_s / (n_s + n_c)


# -----------------------------------------------------------------------------
#  Densidad a (T,P): resolver rho tal que P(rho)=P_pa
# -----------------------------------------------------------------------------
def _resolver_densidad(model, x, T, P_pa, fase):
    suma = 0.0
    for i, comp in enumerate(model.c):
        for si, seg in enumerate(comp.segs):
            suma += x[i] * comp.m * comp.z[si] * model.d_seg(seg, T) ** 3
    coef = math.pi / 6.0 * suma
    rho_de_eta = lambda eta: eta / coef
    if fase == "liq":
        etas = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1,
                0.05, 0.02, 0.01, 1e-3, 1e-4]
    else:
        etas = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1]

    def P(rho):
        return model.pressure(x, rho, T) - P_pa

    prev_r = rho_de_eta(etas[0]); prev_v = P(prev_r)
    for e in etas[1:]:
        r = rho_de_eta(e); v = P(r)
        if prev_v * v < 0:
            lo, hi = sorted((prev_r, r))
            try:
                return brentq(P, lo, hi, xtol=1e-16, rtol=1e-10, maxiter=100)
            except Exception:
                return 0.5 * (lo + hi)
        prev_r, prev_v = r, v
    return None


# -----------------------------------------------------------------------------
#  Presion de equilibrio (Ec. 39)
# -----------------------------------------------------------------------------
def presion_de_burbuja(model_liq, model_vap, T, x1):
    x = [x1, 1.0 - x1]

    def residual(P_pa):
        rl = _resolver_densidad(model_liq, x, T, P_pa, "liq")
        if rl is None:
            return None
        lnphi_l = model_liq.ln_phi(x, rl, T)[0]
        rv = _resolver_densidad(model_vap, [1.0], T, P_pa, "vap")
        if rv is None:
            return None
        lnphi_v = model_vap.ln_phi([1.0], rv, T)[0]
        return (math.log(x1) + lnphi_l) - lnphi_v          # Ec. 39

    P_grid = np.logspace(math.log10(1.0e2), math.log10(2.0e6), 28)
    prev_P = prev_f = None
    bracket = None
    for P_pa in P_grid:
        f = residual(float(P_pa))
        if f is None:
            prev_P = prev_f = None
            continue
        if prev_f is not None and prev_f * f < 0:
            bracket = (prev_P, float(P_pa))
            break
        prev_P, prev_f = float(P_pa), f
    if bracket is None:
        return None
    lo, hi = bracket
    try:
        return brentq(lambda P: residual(P), lo, hi,
                      xtol=1e-2, rtol=1e-9, maxiter=60) / 1e6
    except Exception:
        return 0.5 * (lo + hi) / 1e6


# -----------------------------------------------------------------------------
#  Programa principal
# -----------------------------------------------------------------------------
def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    solv = construir_solvente()
    cop = construir_copolimero()
    kij = construir_kij()
    model_liq = PCSAFT([solv, cop], kij)
    model_vap = PCSAFT([solv], {})

    T = TEMPERATURA_K
    ws = np.linspace(max(W1_MIN, 1e-4), W1_MAX, N_PUNTOS)

    print("=" * 66)
    print("  VLE PC-SAFT  --  METODOLOGIA 2 (co-monomeros separados)")
    print("  %s (1) + %s (2)" % (SOLVENTE["nombre"], COPOLIMERO["nombre"]))
    print("  T = %.2f K" % T)
    print("  kij_a(solv-%s)=%.4f  kij_b(solv-%s)=%.4f  ki_ab=%.4f"
          % (cop.segs[0].nombre, KIJ_A, cop.segs[1].nombre, KIJ_B, KI_AB))
    print("  z_alpha=%.4f  z_beta=%.4f  B=%s" % (cop.z[0], cop.z[1], cop.bonding))
    print("=" * 66)
    print("  %10s  %12s" % ("w1", "P (MPa)"))

    w_ok, P_ok = [], []
    for w1 in ws:
        x1 = w1_a_x1(float(w1))
        try:
            P = presion_de_burbuja(model_liq, model_vap, T, x1)
        except Exception:
            P = None
        if P is not None and P > 0:
            w_ok.append(float(w1)); P_ok.append(P)
            print("  %10.4f  %12.6f" % (w1, P))
        else:
            print("  %10.4f  %12s" % (w1, "(no converge)"))

    if not w_ok:
        print("\n[!] No se obtuvieron puntos. Revisa los parametros.")
        return

    if DATOS_EXPERIMENTALES:
        s = 0.0; n = 0
        for w1, Pexp in DATOS_EXPERIMENTALES:
            x1 = w1_a_x1(w1)
            Pc = presion_de_burbuja(model_liq, model_vap, T, x1)
            if Pc:
                s += abs(Pc - Pexp) / Pexp; n += 1
        if n:
            print("\n  %%DeltaP/P (Ec. 40) = %.2f %%   (np=%d)" % (100 * s / n, n))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(w_ok, P_ok, "-", color="#3b0a70", lw=2, label="PC-SAFT (metod. 2)")
    if DATOS_EXPERIMENTALES:
        ax.plot([d[0] for d in DATOS_EXPERIMENTALES],
                [d[1] for d in DATOS_EXPERIMENTALES],
                "o", color="#d62728", ms=6, markerfacecolor="white",
                markeredgewidth=1.4, label="Experimental")
    ax.set_xlabel("Fraccion masica de solvente  $w_1$")
    ax.set_ylabel("Presion  $P$  (MPa)")
    ax.set_title("VLE PC-SAFT (metodologia 2)\n%s (1) + %s (2)   T=%.2f K"
                 % (SOLVENTE["nombre"], COPOLIMERO["nombre"], T))
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    carpeta = os.path.dirname(os.path.abspath(__file__))
    ruta = os.path.join(carpeta, ARCHIVO_FIGURA)
    fig.savefig(ruta, dpi=130)
    print("\n  Figura guardada en: %s" % ruta)


if __name__ == "__main__":
    main()
