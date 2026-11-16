#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
  VLE de mezclas COPOLIMERO + SOLVENTE con la EoS PC-SAFT (para copolimeros)
===============================================================================

Reproduce el modelado del articulo:
    Arce, Mattedi & Aznar (2006), Fluid Phase Equilibria 246, 52-63.
    "Vapor-liquid equilibrium of copolymer + solvent mixtures ..."

Archivo UNICO y autocontenido. Calcula el equilibrio liquido-vapor
(presion P vs fraccion masica de solvente w1) para UNA sola combinacion
solvente + copolimero, la que definas en el BLOQUE DE PARAMETROS de abajo.

Para cambiar de sistema: edita a mano los valores de la seccion
"PARAMETROS EDITABLES". Ahi eliges el solvente, el copolimero (sus dos
monomeros / homopolimeros), la composicion, la temperatura, el rango de w1 y
los parametros PC-SAFT (m/MW, sigma, eps/k) y de interaccion binaria (kij).

Dependencias:  numpy, scipy, matplotlib
    pip install numpy scipy matplotlib

Uso:
    python pcsaft_vle.py
Genera la figura 'resultado_VLE.png' e imprime la tabla P(w1) en pantalla.
===============================================================================
"""

import math
import numpy as np
from scipy.optimize import brentq

# =============================================================================
#  ////////////////////  PARAMETROS EDITABLES  //////////////////////////////
#  Cambia SOLO los valores de este bloque para calcular otra combinacion.
#  Los valores por defecto corresponden a:  ACETONA (1) + PSBR5 (2) a 323.15 K
#  (poly(styrene-co-butadiene) con 5% en peso de estireno).
#  Todos los parametros estan tomados del paper (Tablas 4b, 6).
# =============================================================================

# ---- 1) TEMPERATURA y rango de composicion --------------------------------
TEMPERATURA_K   = 323.15          # Temperatura del calculo [K]
W1_MIN          = 0.001           # Fraccion masica de solvente minima
W1_MAX          = 0.12            # Fraccion masica de solvente maxima
N_PUNTOS        = 80              # Numero de puntos de la curva ESTIMADA
                                  # (mas puntos = curva mas suave y detallada)

# ---- 2) SOLVENTE (componente 1) -------------------------------------------
#   m/MW [1/mol], sigma [Angstrom], eps/k [K], y su peso molecular MW [g/mol].
#   (Tabla 4b del paper para m/MW, sigma, eps_k).
SOLVENTE = {
    "nombre":  "ACE (acetona)",
    "MW":      58.08,             # g/mol
    "m_MW":    0.037879,          # (m/MW) [1/mol]
    "sigma":   3.4956,            # Angstrom
    "eps_k":   272.26,            # eps/k [K]
}

# ---- 3) COPOLIMERO (componente 2) -----------------------------------------
#   Formado por dos monomeros (alpha y beta). Para PC-SAFT se usan los
#   parametros del HOMOPOLIMERO de cada monomero (Tabla 4b):
#     alpha -> homopolimero 1   |   beta -> homopolimero 2
#   'w_alpha' es la fraccion MASICA del monomero alpha en el copolimero
#   (= el numero que aparece en el nombre / 100; p.ej. PSBR5 -> 0.05).
#   'MW_copolimero' es el peso molecular total del copolimero [g/mol].
COPOLIMERO = {
    "nombre":        "PSBR5 (poly(styrene-co-butadiene), 5% estireno)",
    "MW_copolimero": 100000.0,    # g/mol  (peso molecular total)
    "w_alpha":       0.05,        # fraccion masica del monomero alpha (estireno)

    # --- homopolimero del monomero ALPHA (aqui: PS, poliestireno) ---
    "alpha": {
        "nombre": "PS",
        "m_MW":   0.033238,       # (m/MW) [1/mol]
        "sigma":  3.5022,         # Angstrom
        "eps_k":  320.14,         # eps/k [K]
    },
    # --- homopolimero del monomero BETA (aqui: PB, polibutadieno) ---
    "beta": {
        "nombre": "PB",
        "m_MW":   0.060552,       # (m/MW) [1/mol]
        "sigma":  2.9845,         # Angstrom
        "eps_k":  240.21,         # eps/k [K]
    },
}

# ---- 4) PARAMETROS DE INTERACCION BINARIA kij (adimensionales) ------------
#   kij corrige la energia cruzada: eps_cross = sqrt(eps_i*eps_j)*(1-kij).
#   - kij_solvente_alpha : solvente <-> monomero alpha
#   - kij_solvente_beta  : solvente <-> monomero beta
#   - kij_alpha_beta     : monomero alpha <-> monomero beta (dentro del copol.)
#
#   Metodologia 1 del paper (Tabla 6): un unico kij solvente-copolimero,
#   el mismo para alpha y beta, y kij_alpha_beta = 0.
#   Metodologia 2 del paper (Tabla 7): kij distintos para alpha y beta, mas
#   un kij_alpha_beta pequeno.
#   Valor por defecto: ACE+PSBR5 metodologia 1 -> kij = 0.0120.
KIJ_SOLVENTE_ALPHA = 0.0120
KIJ_SOLVENTE_BETA  = 0.0120
KIJ_ALPHA_BETA     = 0.0000

# ---- 5) DATOS EXPERIMENTALES (opcional) -----------------------------------
#   Pega aqui los puntos experimentales del sistema que estes calculando.
#   Se dibujan como PUNTOS sobre la curva estimada (que es solo linea).
#   Formato: una tupla (w1, P) por fila, con
#       w1 = fraccion MASICA de solvente
#       P  = presion [MPa]   (OJO: en MPa, no en kPa)
#
#   Puedes copiar los valores desde 'datos_experimentales.txt':
#     - Si ahi estan en kPa, divide entre 1000 para pasarlos a MPa.
#       (p.ej. 36.3 kPa  ->  0.0363 MPa)
#     - Si estan como presion reducida (P/Ps), multiplica por la presion de
#       saturacion Ps del solvente [MPa] a esa T y pega el resultado.
#
#   Si dejas la lista vacia ( [] ), simplemente no se dibujan puntos.
#
#   Ejemplo (ACE + PSBR5 a 323.15 K, datos de Gupta & Prausnitz, pasados a MPa):
DATOS_EXPERIMENTALES = [
    # (w1,    P [MPa])
    (0.003,  0.0137),
    (0.023,  0.0265),
    (0.040,  0.0363),
    (0.056,  0.0415),
    (0.067,  0.0453),
    (0.078,  0.0511),
    (0.080,  0.0533),
    (0.091,  0.0552),
    (0.115,  0.0596),
    (0.120,  0.0613),
]

# ---- 6) SALIDA ------------------------------------------------------------
ARCHIVO_FIGURA = "resultado_VLE.png"

# =============================================================================
#  ///////////////  FIN DE PARAMETROS EDITABLES  ////////////////////////////
#  (A partir de aqui no hace falta editar nada para cambiar de sistema.)
# =============================================================================


# -----------------------------------------------------------------------------
#  Constantes universales de PC-SAFT (Gross & Sadowski, 2001)
# -----------------------------------------------------------------------------
A_COEF = [
    [0.9105631789, -0.3084016918, -0.0906148351],
    [0.6361281449,  0.1860531159,  0.4527842806],
    [2.6861347891, -2.5030047259,  0.5962700728],
    [-26.547362491, 21.419793629, -1.7241829131],
    [97.759208784, -65.255885330, -4.1302112531],
    [-159.59154087, 83.318680481, 13.776631870],
    [91.297774084, -33.746922930, -8.6728470368],
]
B_COEF = [
    [0.7240946941, -0.5755498075,  0.0976883116],
    [2.2382791861,  0.6995095521, -0.2557574982],
    [-4.0025849485, 3.8925673390, -9.1558561530],
    [-21.003576815, -17.215471648, 20.642075974],
    [26.855641363,  192.67226447, -38.804430052],
    [206.55133841, -161.82646165,  93.626774077],
    [-355.60235612, -165.20769346, -29.666905585],
]

KB = 1.380649e-23        # J/K  (constante de Boltzmann)


# -----------------------------------------------------------------------------
#  Representacion interna de segmentos y componentes
# -----------------------------------------------------------------------------
class Segmento:
    """Un tipo de segmento PC-SAFT."""
    def __init__(self, nombre, m, sigma, eps_k):
        self.nombre = nombre
        self.m = m            # numero de segmentos aportado
        self.sigma = sigma    # Angstrom
        self.eps_k = eps_k    # K


class Componente:
    """Componente con 1 (solvente) o 2 (copolimero) tipos de segmento."""
    def __init__(self, nombre, segmentos, z):
        self.nombre = nombre
        self.segs = segmentos         # lista de Segmento
        self.z = z                    # fracciones de segmento (suman 1)
        self.m = sum(s.m for s in segmentos)


# -----------------------------------------------------------------------------
#  Motor PC-SAFT para una mezcla
# -----------------------------------------------------------------------------
class PCSAFT:
    def __init__(self, componentes, kij=None):
        self.comps = componentes
        self.kij = kij or {}     # dict frozenset({nombreA,nombreB}) -> kij

    def _kij(self, a, b):
        if a == b:
            return 0.0
        return self.kij.get(frozenset((a, b)), 0.0)

    def _sigma_cross(self, sa, sb):
        return 0.5 * (sa.sigma + sb.sigma)

    def _eps_cross(self, sa, sb):
        k = self._kij(sa.nombre, sb.nombre)
        return math.sqrt(sa.eps_k * sb.eps_k) * (1.0 - k)

    def d_seg(self, seg, T):
        """Diametro de segmento dependiente de T (Ec. 26)."""
        return seg.sigma * (1.0 - 0.12 * math.exp(-3.0 * seg.eps_k / T))

    def m_bar(self, x):
        """Numero de segmentos promedio (Ec. 24)."""
        return sum(x[i] * self.comps[i].m for i in range(len(self.comps)))

    def _zeta(self, x, rho, T):
        """zeta_k, k=0..3 (Ec. 25). rho en 1/A^3."""
        z = [0.0, 0.0, 0.0, 0.0]
        for k in range(4):
            s = 0.0
            for i, c in enumerate(self.comps):
                for si, seg in enumerate(c.segs):
                    d = self.d_seg(seg, T)
                    s += x[i] * c.m * c.z[si] * d ** k
            z[k] = math.pi / 6.0 * rho * s
        return z

    def eta_of(self, x, rho, T):
        return self._zeta(x, rho, T)[3]

    # ---- contribucion hard-chain (Ecs. 22-23) ----
    def _ghs(self, di, dj, zeta):
        z2, z3 = zeta[2], zeta[3]
        dd = (di * dj) / (di + dj)
        return (1.0 / (1.0 - z3)
                + dd * 3.0 * z2 / (1.0 - z3) ** 2
                + dd ** 2 * 2.0 * z2 ** 2 / (1.0 - z3) ** 3)

    def _a_hs(self, zeta):
        z0, z1, z2, z3 = zeta
        return (1.0 / z0) * (3.0 * z1 * z2 / (1.0 - z3)
                             + z2 ** 3 / (z3 * (1.0 - z3) ** 2)
                             + (z2 ** 3 / z3 ** 2 - z0) * math.log(1.0 - z3))

    def a_hc(self, x, rho, T, bonding):
        """Termino de cadena dura (Ec. 22)."""
        zeta = self._zeta(x, rho, T)
        a = self.m_bar(x) * self._a_hs(zeta)
        for i, c in enumerate(self.comps):
            if len(c.segs) == 1:
                d = self.d_seg(c.segs[0], T)
                a -= x[i] * (c.m - 1.0) * math.log(self._ghs(d, d, zeta))
            else:
                B = bonding.get(i)
                da = self.d_seg(c.segs[0], T)
                db = self.d_seg(c.segs[1], T)
                if B is None:
                    g = self._ghs(da, db, zeta)
                    a -= x[i] * (c.m - 1.0) * math.log(g)
                else:
                    term = (B["aa"] * math.log(self._ghs(da, da, zeta))
                            + B["ab"] * math.log(self._ghs(da, db, zeta))
                            + B["bb"] * math.log(self._ghs(db, db, zeta)))
                    a -= x[i] * (c.m - 1.0) * term
        return a

    # ---- contribucion de dispersion (Ecs. 27-32) ----
    def _ab_k(self, mbar):
        a, b = [], []
        for k in range(7):
            a.append(A_COEF[k][0] + (mbar - 1.0) / mbar * A_COEF[k][1]
                     + (mbar - 1.0) / mbar * (mbar - 2.0) / mbar * A_COEF[k][2])
            b.append(B_COEF[k][0] + (mbar - 1.0) / mbar * B_COEF[k][1]
                     + (mbar - 1.0) / mbar * (mbar - 2.0) / mbar * B_COEF[k][2])
        return a, b

    def _I1_I2(self, mbar, eta):
        a, b = self._ab_k(mbar)
        I1 = sum(a[k] * eta ** k for k in range(7))
        I2 = sum(b[k] * eta ** k for k in range(7))
        return I1, I2

    def _C1(self, mbar, eta):
        return 1.0 / (1.0
                      + mbar * (8.0 * eta - 2.0 * eta ** 2) / (1.0 - eta) ** 4
                      + (1.0 - mbar) * (20.0 * eta - 27.0 * eta ** 2
                                        + 12.0 * eta ** 3 - 2.0 * eta ** 4)
                      / ((1.0 - eta) * (2.0 - eta)) ** 2)

    def _m2es3(self, x, T, power):
        s = 0.0
        for i, ci in enumerate(self.comps):
            for j, cj in enumerate(self.comps):
                for ai, sa in enumerate(ci.segs):
                    for bj, sb in enumerate(cj.segs):
                        eps = self._eps_cross(sa, sb)
                        sig = self._sigma_cross(sa, sb)
                        s += (x[i] * x[j] * ci.m * cj.m
                              * ci.z[ai] * cj.z[bj]
                              * (eps / T) ** power * sig ** 3)
        return s

    def a_disp(self, x, rho, T):
        mbar = self.m_bar(x)
        eta = self.eta_of(x, rho, T)
        I1, I2 = self._I1_I2(mbar, eta)
        C1 = self._C1(mbar, eta)
        a1 = -2.0 * math.pi * rho * I1 * self._m2es3(x, T, 1)
        a2 = -math.pi * rho * mbar * C1 * I2 * self._m2es3(x, T, 2)
        return a1 + a2

    # ---- residual, Z, presion y fugacidad ----
    def a_res(self, x, rho, T, bonding):
        return self.a_hc(x, rho, T, bonding) + self.a_disp(x, rho, T)

    def Z(self, x, rho, T, bonding):
        drho = 1e-6 * rho
        ap = self.a_res(x, rho + drho, T, bonding)
        am = self.a_res(x, rho - drho, T, bonding)
        return 1.0 + rho * (ap - am) / (2.0 * drho)

    def pressure(self, x, rho, T, bonding):
        """Presion en Pa. rho en 1/A^3."""
        Z = self.Z(x, rho, T, bonding)
        return Z * (rho * 1e30) * KB * T

    def ln_phi(self, x, rho, T, bonding):
        nc = len(self.comps)
        ares = self.a_res(x, rho, T, bonding)
        Z = self.Z(x, rho, T, bonding)
        lnZ = math.log(Z)
        h = 1e-7
        dadx = []
        for i in range(nc):
            xp = list(x); xp[i] += h
            xm = list(x); xm[i] -= h
            dadx.append((self.a_res(xp, rho, T, bonding)
                         - self.a_res(xm, rho, T, bonding)) / (2.0 * h))
        sdx = sum(x[j] * dadx[j] for j in range(nc))
        return [ares + (Z - 1.0) + dadx[i] - sdx - lnZ for i in range(nc)]


# -----------------------------------------------------------------------------
#  Construccion de los componentes a partir del bloque de parametros
# -----------------------------------------------------------------------------
def construir_solvente():
    s = SOLVENTE
    m = s["m_MW"] * s["MW"]
    seg = Segmento(s["nombre_seg"] if "nombre_seg" in s else "SOLV",
                   m, s["sigma"], s["eps_k"])
    return Componente(s["nombre"], [seg], [1.0])


def construir_copolimero():
    c = COPOLIMERO
    MWc = c["MW_copolimero"]
    w_a = c["w_alpha"]
    w_b = 1.0 - w_a
    # Ec. 18: m_{i,s} = w_{i,s} * MW_copol * (m/MW)_s
    m_a = w_a * MWc * c["alpha"]["m_MW"]
    m_b = w_b * MWc * c["beta"]["m_MW"]
    m_tot = m_a + m_b
    z_a = m_a / m_tot
    z_b = m_b / m_tot
    seg_a = Segmento(c["alpha"]["nombre"], m_a, c["alpha"]["sigma"], c["alpha"]["eps_k"])
    seg_b = Segmento(c["beta"]["nombre"],  m_b, c["beta"]["sigma"],  c["beta"]["eps_k"])
    return Componente(c["nombre"], [seg_a, seg_b], [z_a, z_b])


def construir_kij():
    """Diccionario de kij por par de segmentos."""
    sa = SOLVENTE.get("nombre_seg", "SOLV")
    na = COPOLIMERO["alpha"]["nombre"]
    nb = COPOLIMERO["beta"]["nombre"]
    d = {}
    d[frozenset((sa, na))] = KIJ_SOLVENTE_ALPHA
    d[frozenset((sa, nb))] = KIJ_SOLVENTE_BETA
    if na != nb:
        d[frozenset((na, nb))] = KIJ_ALPHA_BETA
    return d


def bonding_fractions(comp):
    """Fracciones de enlace B para distribucion aleatoria (Tabla 1 del paper).

    Se calculan a partir de la fraccion de segmento beta y el numero de
    segmentos del copolimero.  Si no es copolimero efectivo (un solo tipo),
    devuelve None.
    """
    if len(comp.segs) < 2:
        return None
    z_a, z_b = comp.z[0], comp.z[1]
    m = comp.m
    if m <= 1:
        return None
    # Tabla 1: distribucion aleatoria (random)
    if z_b < z_a:
        B_ab = 2.0 * (z_b * m) / (m - 1.0)
        B_aa = 1.0 - B_ab
        B_bb = 0.0
    elif z_b > z_a:
        B_ab = 2.0 * (z_a * m) / (m - 1.0)
        B_bb = 1.0 - B_ab
        B_aa = 0.0
    else:  # alternante
        B_ab = 1.0
        B_aa = 0.0
        B_bb = 0.0
    # normalizacion defensiva
    B_ab = min(max(B_ab, 0.0), 1.0)
    return {"aa": max(B_aa, 0.0), "ab": B_ab, "bb": max(B_bb, 0.0)}


# -----------------------------------------------------------------------------
#  Conversiones y resolucion del equilibrio
# -----------------------------------------------------------------------------
def w1_a_x1(w1):
    """Fraccion masica de solvente -> fraccion molar de solvente."""
    MW_s = SOLVENTE["MW"]
    MW_c = COPOLIMERO["MW_copolimero"]
    n_s = w1 / MW_s
    n_c = (1.0 - w1) / MW_c
    return n_s / (n_s + n_c)


def _resolver_densidad(model, x, T, P_pa, bonding, fase):
    """Encuentra rho (1/A^3) tal que P(rho) = P_pa."""
    sum_seg = 0.0
    for i, c in enumerate(model.comps):
        for si, seg in enumerate(c.segs):
            d = model.d_seg(seg, T)
            sum_seg += x[i] * c.m * c.z[si] * d ** 3
    coef = math.pi / 6.0 * sum_seg
    rho_de_eta = lambda eta: eta / coef

    if fase == "liq":
        etas = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1,
                0.05, 0.02, 0.01, 1e-3, 1e-4]
    else:
        etas = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1]

    def Pres(rho):
        return model.pressure(x, rho, T, bonding) - P_pa

    prev_rho = rho_de_eta(etas[0])
    prev_val = Pres(prev_rho)
    for e in etas[1:]:
        rho = rho_de_eta(e)
        val = Pres(rho)
        if prev_val * val < 0:
            lo, hi = sorted((prev_rho, rho))
            try:
                return brentq(Pres, lo, hi, xtol=1e-16, rtol=1e-10, maxiter=100)
            except Exception:
                return 0.5 * (lo + hi)
        prev_rho, prev_val = rho, val
    return None


def presion_de_burbuja(solv, cop, model_liq, model_vap, T, x1, bonding):
    """Presion de equilibrio [MPa] para fraccion molar de solvente x1.

    Condicion (Ec. 39):  x1 * phi_hat_1^L = phi_1^V  (vapor = solvente puro).
    """
    x = [x1, 1.0 - x1]

    def residual(P_pa):
        rho_l = _resolver_densidad(model_liq, x, T, P_pa, bonding, "liq")
        if rho_l is None:
            return None
        lnphi_l = model_liq.ln_phi(x, rho_l, T, bonding)
        rho_v = _resolver_densidad(model_vap, [1.0], T, P_pa, {}, "vap")
        if rho_v is None:
            return None
        lnphi_v = model_vap.ln_phi([1.0], rho_v, T, {})[0]
        return (math.log(x1) + lnphi_l[0]) - lnphi_v

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

    # construir componentes y modelos
    solv = construir_solvente()
    cop = construir_copolimero()
    kij = construir_kij()
    model_liq = PCSAFT([solv, cop], kij=kij)
    model_vap = PCSAFT([solv], kij=kij)
    bonding = {1: bonding_fractions(cop)}

    T = TEMPERATURA_K
    ws = np.linspace(max(W1_MIN, 1e-4), W1_MAX, N_PUNTOS)

    print("=" * 66)
    print("  VLE PC-SAFT  |  %s (1) + %s (2)" % (SOLVENTE["nombre"],
                                                 COPOLIMERO["nombre"]))
    print("  T = %.2f K" % T)
    print("  copolimero:  m_total=%.1f  z_alpha=%.4f  z_beta=%.4f"
          % (cop.m, cop.z[0], cop.z[1]))
    print("  kij: solv-alpha=%.4f  solv-beta=%.4f  alpha-beta=%.4f"
          % (KIJ_SOLVENTE_ALPHA, KIJ_SOLVENTE_BETA, KIJ_ALPHA_BETA))
    print("=" * 66)
    print("  %10s  %12s" % ("w1", "P (MPa)"))

    w_ok, P_ok = [], []
    for w1 in ws:
        x1 = w1_a_x1(float(w1))
        try:
            P = presion_de_burbuja(solv, cop, model_liq, model_vap, T, x1, bonding)
        except Exception:
            P = None
        if P is not None and P > 0:
            w_ok.append(float(w1))
            P_ok.append(P)
            print("  %10.4f  %12.6f" % (w1, P))
        else:
            print("  %10.4f  %12s" % (w1, "(no converge)"))

    if not w_ok:
        print("\n[!] No se obtuvieron puntos. Revisa los parametros.")
        return

    # grafica
    fig, ax = plt.subplots(figsize=(7, 5))

    # curva ESTIMADA: solo linea (sin marcadores)
    ax.plot(w_ok, P_ok, "-", color="#3b0a70", lw=2,
            label="PC-SAFT (estimado)")

    # datos EXPERIMENTALES: solo puntos (si hay)
    if DATOS_EXPERIMENTALES:
        w_exp = [d[0] for d in DATOS_EXPERIMENTALES]
        P_exp = [d[1] for d in DATOS_EXPERIMENTALES]
        ax.plot(w_exp, P_exp, "o", color="#d62728", ms=6,
                markerfacecolor="white", markeredgewidth=1.4,
                label="Experimental")

    ax.set_xlabel("Fraccion masica de solvente  $w_1$")
    ax.set_ylabel("Presion  $P$  (MPa)")
    ax.set_title("VLE PC-SAFT\n%s (1) + %s (2)   T=%.2f K"
                 % (SOLVENTE["nombre"], COPOLIMERO["nombre"], T))
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    # Guardar SIEMPRE en la misma carpeta donde esta este archivo .py
    # (no en el directorio desde el que se lanza el comando).
    import os
    carpeta_script = os.path.dirname(os.path.abspath(__file__))
    ruta_figura = os.path.join(carpeta_script, ARCHIVO_FIGURA)
    fig.savefig(ruta_figura, dpi=130)
    print("\n  Figura guardada en: %s" % ruta_figura)


if __name__ == "__main__":
    main()
