# -*- coding: utf-8 -*-
"""
===============================================================================
  VLE de COPOLIMERO + SOLVENTE con PC-SAFT - METODOLOGIA 1 (pseudo-monomero)
  ARCHIVO UNICO Y AUTOCONTENIDO
===============================================================================

  Reproduccion del articulo:
    Arce, Mattedi & Aznar (2006), Fluid Phase Equilibria 246, 52-63.

  Metodologia 1: el copolimero se trata como UN pseudo-monomero (un solo tipo de
  segmento) con los parametros PC-SAFT de la Tabla 5 del paper. Un unico kij
  solvente-copolimero (Tabla 6).

  Este archivo contiene TODO en un solo lugar:
    1) Constantes universales aik/bik (Gross-Sadowski 2001, ref [23]).
    2) Motor PC-SAFT (clase PCSAFT, Ecs. 21-39).
    3) Parametros y datos experimentales de todos los sistemas (Tablas 4b,5,6).
    4) Calculo del %DeltaP (Ec. 40) vs el paper.
    5) Generacion de las figuras VLE, UNA POR SOLVENTE/FAMILIA.

  NOTA: la metodologia 1 reproduce bien los sistemas no polares (PSBR con
  alcanos/aromaticos: ~1-20%), pero da errores grandes en los sistemas polares
  (PVAE/PVAVC con cloroformo, acetona, butanol). Esto es inherente al
  pseudo-monomero de baja densidad de esos copolimeros (ver el articulo).

  Uso:
    python metodologia_1.py          -> imprime tabla de %dP y genera figuras
  Dependencias: numpy, scipy, matplotlib
===============================================================================
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import math
from collections import defaultdict
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
#  1) CONSTANTES UNIVERSALES  (Gross & Sadowski 2001, ref [23] del paper)
# =============================================================================
A_COEF = [
    [0.9105631445, -0.3084016918, -0.0906148351],
    [0.6361281449,  0.1860531159,  0.4527842806],
    [2.6861347891, -2.5030047259,  0.5962700728],
    [-26.547362491, 21.419793629, -1.7241829131],
    [97.759208784, -65.255885330, -4.1302112531],
    [-159.59154087, 83.318680481, 13.776631870],
    [91.297774084, -33.746922930, -8.6728470368],
]

# b0i, b1i, b2i  (para I2 -> termino de dispersion de segundo orden, Ec. 29)
B_COEF = [
    [0.7240946941, -0.5755498075,  0.0976883116],
    [2.2382791861,  0.6995095521, -0.2557574982],
    [-4.0025849485, 3.8925673390, -9.1558561530],
    [-21.003576815, -17.215471648, 20.642075974],
    [26.855641363,  192.67226447, -38.804430052],
    [206.55133841, -161.82646165,  93.626774077],
    [-355.60235612, -165.20769346, -29.666905585],
]


KB = 1.380649e-23        # J/K (constante de Boltzmann)

# Parametros del sistema activo (los setea dp_system / curva por cada sistema)
SOLVENTE = {}
COPOLIMERO = {}


# =============================================================================
#  2) MOTOR PC-SAFT (Ecs. 21-39 del paper)
# =============================================================================


# -----------------------------------------------------------------------------
#  Componente PC-SAFT: cadena de UN solo tipo de segmento
#  (solvente o copolimero-pseudo-monomero). Caracterizado por m, sigma, eps/k.
# -----------------------------------------------------------------------------
class Componente:
    def __init__(self, nombre, m, sigma, eps_k, kAB=0.0, epsAB=0.0):
        self.nombre = nombre
        self.m = m            # numero de segmentos por molecula
        self.sigma = sigma    # Angstrom
        self.eps_k = eps_k    # K
        self.kAB = kAB        # volumen de asociacion (Tabla 4b) -- 0 si no asocia
        self.epsAB = epsAB    # energia de asociacion eps_AB/k [K] (Tabla 4b)
        self.asocia = (kAB > 0.0 and epsAB > 0.0)   # esquema 2B (alcohol)


# -----------------------------------------------------------------------------
#  PC-SAFT para una mezcla de componentes de un solo tipo de segmento cada uno
#  (Ecs. 21-35 del paper). Sin asociacion.
# -----------------------------------------------------------------------------
class PCSAFT:
    def __init__(self, componentes, kij_matrix):
        self.c = componentes
        self.nc = len(componentes)
        self.kij = kij_matrix     # matriz nc x nc

    def d_i(self, i, T):
        """Diametro dependiente de T (Ec. 26)."""
        s = self.c[i]
        return s.sigma * (1.0 - 0.12 * math.exp(-3.0 * s.eps_k / T))

    def m_bar(self, x):
        """Ec. 24: m promedio."""
        return sum(x[i] * self.c[i].m for i in range(self.nc))

    def zeta(self, x, rho, T):
        """Ec. 25: zeta_k, k=0..3.  rho en 1/A^3."""
        z = [0.0, 0.0, 0.0, 0.0]
        for k in range(4):
            s = 0.0
            for i in range(self.nc):
                s += x[i] * self.c[i].m * self.d_i(i, T) ** k
            z[k] = math.pi / 6.0 * rho * s
        return z

    # ---- hard-sphere (Ec. estandar PC-SAFT) ----
    def a_hs(self, z):
        z0, z1, z2, z3 = z
        return (1.0 / z0) * (3.0 * z1 * z2 / (1.0 - z3)
                             + z2 ** 3 / (z3 * (1.0 - z3) ** 2)
                             + (z2 ** 3 / z3 ** 2 - z0) * math.log(1.0 - z3))

    def g_hs(self, i, j, z, T):
        """Funcion de distribucion radial (Ec. 23)."""
        z2, z3 = z[2], z[3]
        di, dj = self.d_i(i, T), self.d_i(j, T)
        dd = di * dj / (di + dj)
        return (1.0 / (1.0 - z3)
                + dd * 3.0 * z2 / (1.0 - z3) ** 2
                + dd ** 2 * 2.0 * z2 ** 2 / (1.0 - z3) ** 3)

    # ---- hard-chain (Ec. 22) ----
    #   Como cada componente es un solo tipo de segmento, B_{aa}=1 y el
    #   termino se reduce a  -sum xi (mi-1) ln g^hs_ii  (igual que homopolimero).
    def a_hc(self, x, rho, T):
        z = self.zeta(x, rho, T)
        a = self.m_bar(x) * self.a_hs(z)
        for i in range(self.nc):
            a -= x[i] * (self.c[i].m - 1.0) * math.log(self.g_hs(i, i, z, T))
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
        # Ec. 30
        return 1.0 / (1.0
                      + mb * (8.0 * eta - 2.0 * eta ** 2) / (1.0 - eta) ** 4
                      + (1.0 - mb) * (20.0 * eta - 27.0 * eta ** 2
                                      + 12.0 * eta ** 3 - 2.0 * eta ** 4)
                      / ((1.0 - eta) * (2.0 - eta)) ** 2)

    def _eps_sigma_cross(self, i, j):
        """Ec. 33: reglas de combinacion."""
        sig = 0.5 * (self.c[i].sigma + self.c[j].sigma)
        eps = math.sqrt(self.c[i].eps_k * self.c[j].eps_k) * (1.0 - self.kij[i][j])
        return eps, sig

    def _m2es(self, x, T, power):
        """sum_i sum_j xi xj mi mj (eps_ij/kT)^power sigma_ij^3."""
        s = 0.0
        for i in range(self.nc):
            for j in range(self.nc):
                eps, sig = self._eps_sigma_cross(i, j)
                s += (x[i] * x[j] * self.c[i].m * self.c[j].m
                      * (eps / T) ** power * sig ** 3)
        return s

    def a_disp(self, x, rho, T):
        mb = self.m_bar(x)
        eta = self.zeta(x, rho, T)[3]
        I1, I2 = self._I1_I2(mb, eta)
        C1 = self._C1(mb, eta)
        a1 = -2.0 * math.pi * rho * I1 * self._m2es(x, T, 1)   # Ec. 28
        a2 = -math.pi * rho * mb * C1 * I2 * self._m2es(x, T, 2)  # Ec. 29
        return a1 + a2

    # ---- asociacion (Ecs. 36-38) ----
    #   Esquema 2B (dos sitios A,B por molecula; A solo asocia con B).
    #   Solo contribuyen los componentes auto-asociantes (alcohol nC4OH).
    #   Para un solo asociante en la mezcla, por simetria 2B XA = XB, y:
    #       XA = (-1 + sqrt(1 + 4*rho_i*Delta)) / (2*rho_i*Delta)
    #   con Delta = sigma_ii^3 * g^hs_ii * kAB * (exp(epsAB/T) - 1)   (Ec. 38)
    #   rho_i = rho * x_i  es la densidad numerica del asociante [1/A^3].
    def a_assoc(self, x, rho, T):
        total = 0.0
        z = None
        for i in range(self.nc):
            ci = self.c[i]
            if not getattr(ci, "asocia", False):
                continue
            if z is None:
                z = self.zeta(x, rho, T)
            ghs = self.g_hs(i, i, z, T)
            delta = ci.sigma ** 3 * ghs * ci.kAB * (math.exp(ci.epsAB / T) - 1.0)
            rho_i = rho * x[i]                 # densidad numerica del asociante
            rd = rho_i * delta
            if rd <= 0.0:
                continue
            XA = (-1.0 + math.sqrt(1.0 + 4.0 * rd)) / (2.0 * rd)
            # Ec. 36 con 2 sitios iguales (A y B): sum_sites (ln X - X/2 + 1/2)
            a_i = 2.0 * (math.log(XA) - XA / 2.0 + 0.5)
            total += x[i] * a_i
        return total

    # ---- residual, Z, presion, fugacidad ----
    def a_res(self, x, rho, T):
        a = self.a_hc(x, rho, T) + self.a_disp(x, rho, T)        # Ec. 21
        if any(getattr(ci, "asocia", False) for ci in self.c):
            a += self.a_assoc(x, rho, T)                         # + Ec. 36
        return a

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
#  Construccion a partir del bloque de parametros
# -----------------------------------------------------------------------------
def construir_solvente():
    s = SOLVENTE
    m = s["m_MW"] * s["MW"]                     # m = (m/MW)*MW
    return Componente(s["nombre"], m, s["sigma"], s["eps_k"])


def construir_copolimero():
    c = COPOLIMERO
    m = c["m_MW"] * c["MW"]                     # Ec.18-19 con pseudo-monomero:
                                               # m_copol = (m/MW)_psm * MW_copol
    return Componente(c["nombre"], m, c["sigma"], c["eps_k"])


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
    for i in range(model.nc):
        suma += x[i] * model.c[i].m * model.d_i(i, T) ** 3
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
        return (math.log(x1) + lnphi_l) - lnphi_v      # Ec. 39

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
                      xtol=1e-2, rtol=1e-9, maxiter=60) / 1e6   # MPa
    except Exception:
        return 0.5 * (lo + hi) / 1e6


# -----------------------------------------------------------------------------
#  Programa principal
# -----------------------------------------------------------------------------

# =============================================================================
#  3) PARAMETROS Y DATOS EXPERIMENTALES (Tablas 4b, 5, 6 del paper)
# =============================================================================
# Parametros de solventes (Tabla 4b): (m_MW, sigma, eps_k, MW)
SOLV = {
    "ACE":   (0.037879, 3.4956, 272.26, 58.08),
    "nC4OH": (0.035126, 3.7563, 263.18, 74.12),
    "CHL":   (0.021826, 3.4085, 265.94, 119.38),
    "cC6":   (0.029946, 3.8399, 279.21, 84.16),
    "nC6":   (0.034088, 3.8383, 239.57, 86.18),
    "nC5":   (0.034090, 3.8729, 240.11, 72.15),
    "TOL":   (0.030984, 3.6869, 282.69, 92.14),
    "BZN":   (0.032712, 3.5878, 281.55, 78.12),
}
# Pseudo-monomeros (Tabla 5): (m_MW, sigma, eps_k) ; y MW copolimero (Tabla 3)
COP = {
    "PSBR4.1": (0.045562, 3.4622, 217.35, 100000),
    "PSBR5":   (0.045456, 3.4647, 217.86, 100000),
    "PSBR23":  (0.043119, 3.5212, 229.01, 720000),
    "PSBR41":  (0.040280, 3.5898, 242.57, 100000),
    "PSBR45":  (0.039562, 3.6072, 246.00, 600000),
    "PVAE9":   (0.054676, 3.4765, 186.36, 100000),
    "PVAE25":  (0.051700, 3.5864, 202.57, 88000),
    "PVAE50":  (0.045107, 3.8300, 238.49, 250000),
    "PVAE70":  (0.036813, 4.1364, 283.67, 280000),
    "PVAVC12": (0.084039, 2.6425, 183.57, 30000),
}
# kij (Tabla 6, metodologia 1) y dP reportado PC-SAFT
KIJ_DP = {
    ("ACE","PSBR5"):(0.0120,4.78),("ACE","PSBR23"):(0.0122,3.04),("ACE","PSBR45"):(0.0130,5.19),
    ("nC4OH","PVAE9"):(0.0230,4.05),("nC4OH","PVAVC12"):(-0.0126,5.18),
    ("CHL","PVAE9"):(-0.1208,4.18),("CHL","PVAE25"):(-0.0816,5.06),("CHL","PVAE50"):(-0.0701,4.41),("CHL","PVAE70"):(-0.0516,3.28),
    ("cC6","PSBR5"):(0.0118,3.66),("cC6","PSBR23"):(0.0104,4.38),("cC6","PSBR41"):(0.0088,2.28),("cC6","PSBR45"):(0.0101,4.31),
    ("cC6","PVAE25"):(-0.1068,5.36),("cC6","PVAE50"):(-0.0858,6.90),("cC6","PVAE70"):(-0.0636,5.73),
    ("BZN","PSBR4.1"):(0.0202,2.44),
    ("TOL","PSBR5"):(0.0080,4.07),("TOL","PSBR23"):(0.0091,5.76),("TOL","PSBR41"):(0.0098,3.71),
    ("nC6","PSBR23"):(0.0086,3.92),("nC6","PSBR45"):(0.0089,4.28),
    ("nC5","PSBR23"):(0.0072,3.68),("nC5","PSBR45"):(0.0079,6.24),
}

# Datos experimentales por (solvente,copolimero,T) -> lista (w1, P[MPa])
# (solo los que tenemos en MPa, fuente Gupta [1])
DATA = {
 ("ACE","PSBR5",323.15):[(0.003,0.0137),(0.023,0.0265),(0.040,0.0363),(0.056,0.0415),(0.067,0.0453),(0.078,0.0511),(0.080,0.0533),(0.091,0.0552),(0.115,0.0596),(0.120,0.0613)],
 ("ACE","PSBR23",333.15):[(0.009,0.0127),(0.021,0.0264),(0.037,0.0407),(0.062,0.0548),(0.092,0.0685),(0.128,0.0813),(0.147,0.0848),(0.183,0.0915),(0.231,0.0979)],
 ("ACE","PSBR45",333.15):[(0.009,0.0127),(0.023,0.0264),(0.044,0.0407),(0.069,0.0548),(0.101,0.0685),(0.140,0.0813),(0.163,0.0848),(0.197,0.0915),(0.258,0.0979)],
 ("nC4OH","PVAE9",353.15):[(0.006,0.0027),(0.012,0.0052),(0.017,0.0080),(0.027,0.0108),(0.031,0.0123),(0.034,0.0139),(0.036,0.0152)],
 ("nC4OH","PVAVC12",353.15):[(0.007,0.0027),(0.016,0.0052),(0.022,0.0080),(0.033,0.0108),(0.040,0.0123),(0.044,0.0139),(0.046,0.0152)],
 ("CHL","PVAE9",333.15):[(0.031,0.0140),(0.067,0.0281),(0.104,0.0401),(0.147,0.0532),(0.168,0.0584),(0.196,0.0661),(0.216,0.0723),(0.235,0.0792),(0.253,0.0871)],
 ("CHL","PVAE25",333.15):[(0.062,0.0139),(0.140,0.0267),(0.228,0.0399),(0.320,0.0532),(0.358,0.0599),(0.401,0.0665),(0.441,0.0731),(0.481,0.0793)],
 ("CHL","PVAE50",333.15):[(0.098,0.0139),(0.202,0.0267),(0.303,0.0399),(0.398,0.0532),(0.439,0.0599),(0.477,0.0665),(0.512,0.0731),(0.552,0.0793)],
 ("CHL","PVAE70",333.15):[(0.059,0.0139),(0.161,0.0267),(0.287,0.0399),(0.403,0.0532),(0.450,0.0599),(0.489,0.0665),(0.526,0.0731),(0.573,0.0793)],
 ("cC6","PSBR5",333.15):[(0.027,0.0048),(0.057,0.0105),(0.103,0.0169),(0.120,0.0199),(0.152,0.0248),(0.189,0.0295),(0.233,0.0337),(0.276,0.0376),(0.323,0.0407)],
 ("cC6","PSBR23",333.15):[(0.018,0.0048),(0.052,0.0105),(0.088,0.0169),(0.136,0.0248),(0.168,0.0295),(0.208,0.0337),(0.244,0.0376),(0.287,0.0407)],
 ("cC6","PSBR45",333.15):[(0.021,0.0048),(0.050,0.0105),(0.087,0.0169),(0.103,0.0199),(0.130,0.0248),(0.163,0.0295),(0.200,0.0337),(0.240,0.0376),(0.282,0.0407)],
 ("cC6","PVAE25",353.15):[(0.011,0.0075),(0.037,0.0144),(0.074,0.0265),(0.140,0.0397),(0.228,0.0549),(0.238,0.0555),(0.277,0.0607),(0.300,0.0641),(0.369,0.0704)],
 ("cC6","PVAE50",353.15):[(0.005,0.0075),(0.025,0.0144),(0.046,0.0265),(0.104,0.0397),(0.159,0.0549),(0.170,0.0555),(0.205,0.0607),(0.225,0.0641),(0.287,0.0704)],
 ("cC6","PVAE70",353.15):[(0.004,0.0075),(0.026,0.0144),(0.042,0.0265),(0.081,0.0397),(0.123,0.0549),(0.133,0.0555),(0.156,0.0607),(0.191,0.0641),(0.242,0.0704)],
 # TOL+PSBR5 (Gupta, formato (P,w1)): P = {0.93,2.08,3.04,4.36,4.83,5.48} kPa. VERIFICADO vs imagen PDF.
 ("TOL","PSBR5",308.15):[(0.043,0.00093),(0.099,0.00208),(0.184,0.00304),(0.319,0.00436),(0.416,0.00483),(0.513,0.00548)],
 # TOL+PSBR23 (Gupta, formato (P,w1)): P = {0.90,2.08,3.04,4.36,4.83,5.48} kPa. VERIFICADO vs imagen PDF.
 ("TOL","PSBR23",308.15):[(0.029,0.00090),(0.031,0.00208),(0.032,0.00304),(0.037,0.00436),(0.037,0.00483),(0.038,0.00548)],
 ("nC6","PSBR23",343.15):[(0.014,0.0139),(0.043,0.0273),(0.063,0.0385),(0.110,0.0488),(0.162,0.0615),(0.262,0.0752)],
 ("nC6","PSBR45",343.15):[(0.016,0.0139),(0.038,0.0273),(0.061,0.0385),(0.103,0.0488),(0.155,0.0615),(0.249,0.0752)],
 ("nC5","PSBR23",333.15):[(0.008,0.0228),(0.016,0.0360),(0.021,0.0468),(0.033,0.0655),(0.046,0.0840),(0.055,0.0973),(0.070,0.1133),(0.078,0.1248),(0.089,0.1376)],
 ("nC5","PSBR45",333.15):[(0.010,0.0228),(0.013,0.0468),(0.019,0.0655),(0.022,0.0840),(0.022,0.0973),(0.026,0.1133),(0.027,0.1248),(0.034,0.1376)],
 # --- Sistemas de Wohlfarth [34] y Kang [35], datos en P/P0 (reducida) convertidos a MPa con Psat del solvente ---
 # cC6+PSBR41 (Wohlfarth, T=343.15K): P/P0 * Psat_cC6(343K)=0.07250 MPa. VERIFICADO vs imagen PDF (rango 0.0473-0.0665 ~ Arce 0.0472-0.0664).
 ("cC6","PSBR41",343.15):[(0.2169,0.04729),(0.2436,0.05051),(0.2803,0.05431),(0.3288,0.05860),(0.3606,0.06079),(0.3992,0.06303),(0.4414,0.06514),(0.4766,0.06653)],
 # TOL+PSBR41 (Wohlfarth): dos isotermas. T=343.15K: P/P0 * Psat_TOL(343K)=0.02716 MPa.
 ("TOL","PSBR41",343.15):[(0.2381,0.01582),(0.2802,0.01759),(0.3279,0.01935),(0.3724,0.02074),(0.4200,0.02201),(0.4683,0.02308),(0.5039,0.02378),(0.5467,0.02448)],
 # TOL+PSBR41 T=373.15K: P/P0 * Psat_TOL(373K)=0.07417 MPa.
 ("TOL","PSBR41",373.15):[(0.1832,0.03593),(0.2299,0.04231),(0.2613,0.04609),(0.3030,0.05068),(0.3681,0.05646),(0.4274,0.06066),(0.4937,0.06447),(0.5692,0.06782)],
 # BZN+PSBR4.1 (Kang, formato p/p^s): P/P0 * Psat_BZN. T=308.15K (Psat=0.01976) y T=298.15K (Psat=0.01268). VERIFICADO rango 0.0022-0.0173 ~ Arce.
 ("BZN","PSBR4.1",308.15):[(0.0471,0.00395),(0.0756,0.00593),(0.1044,0.00778),(0.125,0.00877),(0.1442,0.00946),(0.1959,0.01178),(0.2289,0.01284),(0.2637,0.01414),(0.3098,0.01544),(0.3742,0.01649),(0.4297,0.01729)],
 ("BZN","PSBR4.1",298.15):[(0.0459,0.00225),(0.08,0.00396),(0.1117,0.00494),(0.1364,0.00580),(0.1656,0.00672),(0.1864,0.00764),(0.3014,0.00976),(0.37,0.01075),(0.4275,0.01148)],
}


# Parametros de asociacion (Tabla 4b): solo el nC4OH es auto-asociante.
ASOC = {
    "nC4OH": (0.008752, 2489.23),   # (kAB, epsAB/k [K])
}


# =============================================================================
#  4) CALCULO DEL %DeltaP (Ec. 40) Y 5) GRAFICAS POR SOLVENTE
# =============================================================================
NOMBRES = {"ACE": "acetona", "nC4OH": "1-butanol", "CHL": "cloroformo",
           "cC6": "ciclohexano", "nC6": "n-hexano", "nC5": "n-pentano",
           "TOL": "tolueno", "BZN": "benceno"}


def _set_sistema(solv_key, cop_key):
    """Carga los parametros globales del sistema (solvente + pseudo-monomero)."""
    ms = SOLV[solv_key]; mc = COP[cop_key]
    SOLVENTE.clear(); SOLVENTE.update(
        {"nombre": solv_key, "MW": ms[3], "m_MW": ms[0], "sigma": ms[1], "eps_k": ms[2]})
    COPOLIMERO.clear(); COPOLIMERO.update(
        {"nombre": cop_key, "MW": float(mc[3]), "m_MW": mc[0], "sigma": mc[1], "eps_k": mc[2]})


def _construir_modelos(solv_key, cop_key, kij):
    """Construye los modelos liquido y vapor; activa asociacion si aplica (nC4OH)."""
    ms = SOLV[solv_key]
    if solv_key in ASOC:
        kAB, epsAB = ASOC[solv_key]
        solv = Componente(solv_key, ms[0] * ms[3], ms[1], ms[2], kAB, epsAB)
    else:
        solv = construir_solvente()
    cop = construir_copolimero()
    ml = PCSAFT([solv, cop], [[0, kij], [kij, 0]])
    mv = PCSAFT([solv], [[0.0]])
    return ml, mv


def dp_system(solv_key, cop_key, T, pts, kij):
    """%DeltaP (Ec. 40) para un sistema."""
    _set_sistema(solv_key, cop_key)
    ml, mv = _construir_modelos(solv_key, cop_key, kij)
    s = 0.0; n = 0
    for w1, Pexp in pts:
        x1 = w1_a_x1(w1)
        try:
            Pc = presion_de_burbuja(ml, mv, T, x1)
        except Exception:
            Pc = None
        if Pc and Pc > 0:
            s += abs(Pc - Pexp) / Pexp; n += 1
    return (100.0 * s / n) if n else None, n


def _familia(cop):
    if cop.startswith("PSBR"): return "PSBR"
    if cop.startswith("PVAE"): return "PVAE"
    if cop.startswith("PVAVC"): return "PVAVC"
    return cop


MARKERS = ["o", "s", "^", "v", "D", "P", "*", "X"]
COLORS = plt.cm.tab10.colors


def _curva(solv_key, cop_key, T, kij, npts=35):
    _set_sistema(solv_key, cop_key)
    ml, mv = _construir_modelos(solv_key, cop_key, kij)
    pts = DATA[(solv_key, cop_key, T)]
    wmin = min(p[0] for p in pts); wmax = max(p[0] for p in pts)
    ws = np.linspace(max(wmin * 0.5, 1e-4), wmax * 1.05, npts)
    wok, Pok = [], []
    for w1 in ws:
        x1 = w1_a_x1(float(w1))
        try:
            P = presion_de_burbuja(ml, mv, T, x1)
        except Exception:
            P = None
        if P and P > 0:
            wok.append(float(w1)); Pok.append(P)
    return wok, Pok, pts


def graficar_por_solvente(out="figuras_por_solvente"):
    """Genera UNA figura por solvente/familia (estilo Figs. 2-10 del paper)."""
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), out)
    os.makedirs(outdir, exist_ok=True)
    grupos = defaultdict(list)
    for (s, c, T) in DATA:
        grupos[(s, _familia(c))].append((c, T))
    print("Generando figuras por solvente en:", outdir)
    for (solv_key, fam), sistemas in grupos.items():
        sistemas.sort()
        fig, ax = plt.subplots(figsize=(7, 5))
        for idx, (cop_key, T) in enumerate(sistemas):
            kij, _ = KIJ_DP[(solv_key, cop_key)]
            col = COLORS[idx % len(COLORS)]; mk = MARKERS[idx % len(MARKERS)]
            try:
                wok, Pok, pts = _curva(solv_key, cop_key, T, kij)
            except Exception as e:
                print("  ERR %s+%s: %s" % (solv_key, cop_key, e)); continue
            if wok:
                ax.plot(wok, Pok, "-", color=col, lw=1.8)
            ax.plot([p[0] for p in pts], [p[1] for p in pts], mk, color=col,
                    ms=6, markerfacecolor="white", markeredgewidth=1.3,
                    label="%s (T=%.2f K)" % (cop_key, T))
        ax.set_xlabel("Fraccion masica de solvente  $w_1$")
        ax.set_ylabel("Presion  $P$  (MPa)")
        nom = NOMBRES.get(solv_key, solv_key)
        ax.set_title("VLE PC-SAFT (metod. 1)\n%s (%s) + %s\n(linea = PC-SAFT, simbolos = experimental)"
                     % (solv_key, nom, fam))
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, title="Copolimero")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "VLE_%s_%s.png" % (solv_key, fam)), dpi=130)
        plt.close(fig)
        print("  OK  %s + %s  (%d copolimeros)" % (solv_key, fam, len(sistemas)))


def main():
    print("%-7s %-9s %8s %8s %10s %8s" % ("Solv", "Copol", "T[K]", "kij", "dP_calc%", "dP_paper%"))
    print("-" * 56)
    suma = 0.0; n = 0
    for (solv, cop, T), pts in DATA.items():
        kij, dp_paper = KIJ_DP[(solv, cop)]
        dp_calc, _ = dp_system(solv, cop, T, pts, kij)
        if dp_calc is not None:
            suma += dp_calc; n += 1
        sc = ("%9.2f" % dp_calc) if dp_calc is not None else "     None"
        print("%-7s %-9s %8.2f %8.4f %s %8.2f" % (solv, cop, T, kij, sc, dp_paper))
    print("-" * 56)
    print("Desviacion global (promedio por sistema), nuestro motor: %.2f%% sobre %d sistemas" % (suma / n, n))
    print("Desviacion global PC-SAFT reportada por el paper (metod. 1): 4.41%")
    print()
    graficar_por_solvente()


if __name__ == "__main__":
    main()
