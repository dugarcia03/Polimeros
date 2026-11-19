# -*- coding: utf-8 -*-
"""
===============================================================================
  PC-SAFT para VLE copolimero + solvente  -  METODOLOGIAS 1 y 2 UNIFICADAS
===============================================================================

  Arce, Mattedi & Aznar (2006), Fluid Phase Equilibria 246, 52-63.
  Extension a copolimeros: Gross, Spuhl, Tumakaka & Sadowski (2003), ref [24].
  Constantes universales aik/bik: Gross & Sadowski (2001), ref [23].

  Este archivo UNIFICA en un solo codigo las dos metodologias del paper:

    METODOLOGIA 1  (pseudo-monomero)
      El copolimero se representa como UN solo tipo de segmento "promedio",
      con parametros de la Tabla 5 y un unico kij solvente-copolimero (Tabla 6).

    METODOLOGIA 2  (copolimero de 2 segmentos)
      El copolimero conserva sus DOS monomeros como segmentos distintos
      (alpha/beta), con parametros de la Tabla 4b, fracciones de enlace B de
      la Tabla 3 y tres kij por par de segmento (Tabla 7).

  Un unico motor termodinamico (clase PCSAFT) sirve para AMBAS: cuando el
  copolimero tiene 1 segmento se reduce a la metodologia 1; cuando tiene 2 se
  comporta como la metodologia 2.  La diferencia esta solo en como se CONSTRUYE
  el copolimero (construir_cop_m1 vs construir_cop_m2) y en los kij.

  Salidas:
    1) Una TABLA de errores comparando %dP de M1 y M2 lado a lado (stdout).
    2) Una sola IMAGEN por solvente/familia con las curvas de M1 y M2
       superpuestas sobre los puntos experimentales (figuras_m1_vs_m2/).
    3) Una sola IMAGEN RESUMEN con la tabla de errores M1 vs M2 (tabla_errores_m1_vs_m2.png).
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


# =============================================================================
#  2) MOTOR PC-SAFT UNIFICADO  (Ecs. 18-39 del paper)
#     Vale para componentes de 1 segmento (M1) o de varios segmentos (M2).
# =============================================================================
class Segmento:
    """Un tipo de segmento: parametros PC-SAFT de un monomero o de un solvente."""
    def __init__(self, nombre, m, sigma, eps_k, kAB=0.0, epsAB=0.0):
        self.nombre = nombre
        self.m = m            # numero de segmentos de este tipo en la molecula
        self.sigma = sigma    # Angstrom
        self.eps_k = eps_k    # K
        self.kAB = kAB        # volumen de asociacion (Tabla 4b)
        self.epsAB = epsAB    # energia de asociacion eps_AB/k [K] (Tabla 4b)
        self.asocia = (kAB > 0.0 and epsAB > 0.0)   # esquema 2B (alcohol)


class Componente:
    """
    Molecula = cadena de uno o varios tipos de segmento.
      segmentos: lista de Segmento.
      B: matriz (ns x ns) de fracciones de enlace B_{a,b} (Tabla 3).
         Para 1 segmento (solvente o pseudo-monomero M1): B=[[1.0]].
    """
    def __init__(self, nombre, segmentos, B=None):
        self.nombre = nombre
        self.seg = segmentos
        self.ns = len(segmentos)
        self.m = sum(s.m for s in segmentos)              # Ec. 19: m_i total
        self.z = [s.m / self.m for s in segmentos]        # Ec. 20: z_{i,a}
        self.B = B if B is not None else [[1.0]]


class PCSAFT:
    """
    Mezcla de 'Componente'. Cada componente puede tener uno o varios segmentos.
      kseg: dict {(i,a,j,b): kij} con el parametro de interaccion entre el
            segmento a del componente i y el segmento b del componente j.
            Simetrico; si falta un par, se asume 0.
    """
    def __init__(self, componentes, kseg=None):
        self.c = componentes
        self.nc = len(componentes)
        self.kseg = kseg or {}

    # ---- helpers ----
    def _k(self, i, a, j, b):
        return self.kseg.get((i, a, j, b), self.kseg.get((j, b, i, a), 0.0))

    def d_seg(self, i, a, T):
        """Ec. 26: diametro dependiente de T del segmento a del componente i."""
        s = self.c[i].seg[a]
        return s.sigma * (1.0 - 0.12 * math.exp(-3.0 * s.eps_k / T))

    def m_bar(self, x):
        """Ec. 24: m promedio = sum_i x_i m_i."""
        return sum(x[i] * self.c[i].m for i in range(self.nc))

    def zeta(self, x, rho, T):
        """Ec. 25: zeta_k = pi/6 rho sum_i x_i m_i sum_a z_{i,a} d_{i,a}^k."""
        z = [0.0, 0.0, 0.0, 0.0]
        for k in range(4):
            s = 0.0
            for i in range(self.nc):
                ci = self.c[i]
                sa = 0.0
                for a in range(ci.ns):
                    sa += ci.z[a] * self.d_seg(i, a, T) ** k
                s += x[i] * ci.m * sa
            z[k] = math.pi / 6.0 * rho * s
        return z

    # ---- hard-sphere ----
    def a_hs(self, z):
        z0, z1, z2, z3 = z
        return (1.0 / z0) * (3.0 * z1 * z2 / (1.0 - z3)
                             + z2 ** 3 / (z3 * (1.0 - z3) ** 2)
                             + (z2 ** 3 / z3 ** 2 - z0) * math.log(1.0 - z3))

    def g_seg(self, di, dj, z):
        """Ec. 23: g^hs entre dos segmentos de diametros di, dj."""
        z2, z3 = z[2], z[3]
        dd = di * dj / (di + dj)
        return (1.0 / (1.0 - z3)
                + dd * 3.0 * z2 / (1.0 - z3) ** 2
                + dd ** 2 * 2.0 * z2 ** 2 / (1.0 - z3) ** 3)

    # ---- hard-chain (Ec. 22, con fracciones de enlace B) ----
    def a_hc(self, x, rho, T):
        z = self.zeta(x, rho, T)
        a = self.m_bar(x) * self.a_hs(z)
        for i in range(self.nc):
            ci = self.c[i]
            if ci.ns == 1:
                # homopolimero/solvente/pseudo-monomero: -x_i (m_i-1) ln g_ii
                di = self.d_seg(i, 0, T)
                a -= x[i] * (ci.m - 1.0) * math.log(self.g_seg(di, di, z))
            else:
                # copolimero: -x_i (m_i-1) sum_a sum_b B_{ab} ln g_{ab}
                suma = 0.0
                for aa in range(ci.ns):
                    da = self.d_seg(i, aa, T)
                    for bb in range(ci.ns):
                        db = self.d_seg(i, bb, T)
                        Bab = ci.B[aa][bb]
                        if Bab != 0.0:
                            suma += Bab * math.log(self.g_seg(da, db, z))
                a -= x[i] * (ci.m - 1.0) * suma
        return a

    # ---- dispersion (Ecs. 27-32, doble suma sobre segmentos) ----
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

    def _m2es(self, x, T, power):
        """
        Ecs. 28-29 extendidas a copolimeros (ref [24]):
          sum_i sum_j x_i x_j m_i m_j sum_a sum_b z_ia z_jb (eps_ab/kT)^p sigma_ab^3
        """
        s = 0.0
        for i in range(self.nc):
            ci = self.c[i]
            for j in range(self.nc):
                cj = self.c[j]
                sab = 0.0
                for a in range(ci.ns):
                    sa = ci.seg[a]
                    for b in range(cj.ns):
                        sb = cj.seg[b]
                        kij = self._k(i, a, j, b)
                        eps = math.sqrt(sa.eps_k * sb.eps_k) * (1.0 - kij)
                        sig = 0.5 * (sa.sigma + sb.sigma)
                        sab += ci.z[a] * cj.z[b] * (eps / T) ** power * sig ** 3
                s += x[i] * x[j] * ci.m * cj.m * sab
        return s

    def a_disp(self, x, rho, T):
        mb = self.m_bar(x)
        eta = self.zeta(x, rho, T)[3]
        I1, I2 = self._I1_I2(mb, eta)
        C1 = self._C1(mb, eta)
        a1 = -2.0 * math.pi * rho * I1 * self._m2es(x, T, 1)       # Ec. 28
        a2 = -math.pi * rho * mb * C1 * I2 * self._m2es(x, T, 2)   # Ec. 29
        return a1 + a2

    # ---- asociacion (Ecs. 36-38), esquema 2B ----
    def a_assoc(self, x, rho, T):
        total = 0.0
        z = None
        for i in range(self.nc):
            ci = self.c[i]
            for a in range(ci.ns):
                sa = ci.seg[a]
                if not sa.asocia:
                    continue
                if z is None:
                    z = self.zeta(x, rho, T)
                da = self.d_seg(i, a, T)
                ghs = self.g_seg(da, da, z)
                delta = sa.sigma ** 3 * ghs * sa.kAB * (math.exp(sa.epsAB / T) - 1.0)
                rho_i = rho * x[i] * ci.z[a]
                rd = rho_i * delta
                if rd <= 0.0:
                    continue
                XA = (-1.0 + math.sqrt(1.0 + 4.0 * rd)) / (2.0 * rd)
                total += x[i] * ci.z[a] * 2.0 * (math.log(XA) - XA / 2.0 + 0.5)
        return total

    # ---- residual, Z, presion, fugacidad ----
    def a_res(self, x, rho, T):
        a = self.a_hc(x, rho, T) + self.a_disp(x, rho, T)         # Ec. 21
        if any(s.asocia for ci in self.c for s in ci.seg):
            a += self.a_assoc(x, rho, T)                          # + Ec. 36
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
#  Densidad a (T,P): resolver rho tal que P(rho)=P_pa
# -----------------------------------------------------------------------------
def resolver_densidad(model, x, T, P_pa, fase):
    suma = 0.0
    for i in range(model.nc):
        ci = model.c[i]
        for a in range(ci.ns):
            suma += x[i] * ci.m * ci.z[a] * model.d_seg(i, a, T) ** 3
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
                return brentq(P, lo, hi, rtol=1e-12, maxiter=100)
            except Exception:
                return None   # si brentq falla, NO devolver un punto medio falso
        prev_r, prev_v = r, v
    return None


# -----------------------------------------------------------------------------
#  Presion de burbuja (Ec. 39): x1 phi_L1 = phi_V1(y1=1)
# -----------------------------------------------------------------------------
def presion_de_burbuja(model_liq, model_vap, T, x1):
    x = [x1, 1.0 - x1]

    def residual(P_pa):
        rl = resolver_densidad(model_liq, x, T, P_pa, "liq")
        if rl is None:
            return None
        lnphi_l = model_liq.ln_phi(x, rl, T)[0]
        rv = resolver_densidad(model_vap, [1.0], T, P_pa, "vap")
        if rv is None:
            return None
        lnphi_v = model_vap.ln_phi([1.0], rv, T)[0]
        return (math.log(x1) + lnphi_l) - lnphi_v      # Ec. 39

    # Rango de presion amplio (10 Pa a 2 MPa), 36 puntos log-espaciados.
    P_grid = np.logspace(math.log10(1.0e1), math.log10(2.0e6), 36)
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
        return None   # si brentq falla, NO devolver un punto medio falso


# =============================================================================
#  3) PARAMETROS Y DATOS EXPERIMENTALES (Tablas 3, 4b, 5, 6, 7 del paper)
# =============================================================================
# --- Solventes (Tabla 4b): (m_MW, sigma, eps_k, MW) ---
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

# --- Asociacion (Tabla 4b): solo el nC4OH es auto-asociante (esquema 2B) ---
ASOC = {"nC4OH": (0.008752, 2489.23)}   # (kAB, epsAB/k [K])

# =====================  METODOLOGIA 1  =======================================
# Pseudo-monomeros (Tabla 5): (m_MW, sigma, eps_k) ; MW copolimero (Tabla 3)
COP_M1 = {
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
# kij metodologia 1 (Tabla 6) y %dP reportado por el paper
KIJ_M1 = {
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

# =====================  METODOLOGIA 2  =======================================
# NOTA DE AUDITORIA (ambiguedad del paper en M2 -- leer antes de tocar):
#   El texto del paper (sec. 3, pag. 55) dice que los parametros de un sistema M2
#   son los de "solvent, poly(alpha) and poly(beta)" (POLIMEROS: PS, PB, PVA, PE),
#   y con esos parametros la fraccion de segmento ziβ reproduce EXACTAMENTE la
#   Tabla 3 (PSBR23: 0.1409; PSBR45: 0.3100).  SIN EMBARGO, al usar los parametros
#   de POLIMERO con los kij publicados (Tabla 7) el error global sube a ~70%,
#   mientras que usando los parametros de MONOMERO (STR, BTD, VA, ET de la Tabla
#   4b) el error global es ~33%.  Esto indica que los kij de la Tabla 7 fueron
#   ajustados por los autores junto con los parametros de MONOMERO (o con una
#   construccion no documentada), de modo que NO es posible reproducir el 2.04%
#   del paper sin reajustar los kij.  Se mantienen los MONOMEROS porque dan el
#   mejor acuerdo con los kij tal como estan publicados.  (Las fracciones de
#   segmento resultantes no igualan la Tabla 3; es una inconsistencia del paper.)
#
# Monomeros (Tabla 4b, seccion "Monomers"): (m_MW, sigma, eps_k, MW)
MON = {
    "STR": (0.024339, 3.9750, 318.70, 104.14), "VA": (0.011525, 5.0706, 421.42, 86.09),
    "BTD": (0.046033, 3.4508, 215.10, 54.10),  "ET": (0.056070, 3.4250, 178.77, 28.10),
    "VC":  (0.091211, 2.4024, 160.05, 62.50),
}
# Copolimero -> (mon1, mon2, wt%_mon1, MW_copol, (B_aa_t3, B_ab_t3, B_bb_t3))
#   mon1 = constituyente del nombre (STR en PSBR, VA en PVAE); mon2 = el otro.
#   Mapeo Tabla 3: beta=mon1 (del wt%), alpha=mon2.  Motor: seg 0=mon1, seg 1=mon2.
COP_M2 = {
    "PSBR23": ("STR","BTD", 23, 720000, (0.7161, 0.2839, 0.0000)),
    "PSBR45": ("STR","BTD", 45, 600000, (0.3707, 0.6293, 0.0000)),
    "PVAE9":  ("VA","ET",    9, 100000, (0.8392, 0.1608, 0.0000)),
    "PVAE25": ("VA","ET",   25,  88000, (0.5423, 0.4577, 0.0000)),
    "PVAE50": ("VA","ET",   50, 250000, (0.0487, 0.9513, 0.0000)),
    "PVAE70": ("VA","ET",   70, 280000, (0.0000, 0.4633, 0.5367)),
}
# Tabla 7: kappa_a (solv-mon1), kappa_b (solv-mon2) por FAMILIA
KAB_FAM = {
    ("ACE","PSBR"): (0.0121, 0.0028), ("nC4OH","PVAE"): (0.0160, 0.0131),
    ("CHL","PVAE"): (-0.0532, 0.0003), ("cC6","PSBR"): (0.0052, 0.0012),
    ("cC6","PVAE"): (-0.1085, 0.0101), ("TOL","PSBR"): (0.0082, 0.0126),
    ("nC6","PSBR"): (0.0080, 0.0091), ("nC5","PSBR"): (0.0066, 0.0108),
}
# Tabla 7: kappa_ab (mon1-mon2) por SISTEMA individual
KAB_AB = {
    ("ACE","PSBR23"):0.0004,("ACE","PSBR45"):0.0001,
    ("nC4OH","PVAE9"):-0.0010,
    ("CHL","PVAE9"):0.0008,("CHL","PVAE25"):0.0005,("CHL","PVAE50"):0.0003,("CHL","PVAE70"):0.0003,
    ("cC6","PSBR23"):0.0006,("cC6","PSBR45"):0.0002,
    ("cC6","PVAE25"):0.0011,("cC6","PVAE50"):0.0010,("cC6","PVAE70"):0.0008,
    ("TOL","PSBR23"):0.0011,
    ("nC6","PSBR23"):0.0010,("nC6","PSBR45"):0.0008,
    ("nC5","PSBR23"):0.0013,("nC5","PSBR45"):0.0008,
}
# %dP reportado por el paper (Tabla 7), es por FAMILIA
DP_FAM_M2 = {("ACE","PSBR"):2.83,("nC4OH","PVAE"):1.88,("CHL","PVAE"):2.59,
             ("cC6","PSBR"):2.11,("cC6","PVAE"):2.09,("TOL","PSBR"):2.13,
             ("nC6","PSBR"):1.84,("nC5","PSBR"):2.07}

# --- Datos experimentales por (solvente, copolimero, T) -> [(w1, P[MPa])] ---
#     (fuentes: Gupta [1], Wohlfarth [34], Kang [35])
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
 ("TOL","PSBR5",308.15):[(0.043,0.00093),(0.099,0.00208),(0.184,0.00304),(0.319,0.00436),(0.416,0.00483),(0.513,0.00548)],
 ("TOL","PSBR23",308.15):[(0.029,0.00090),(0.031,0.00208),(0.032,0.00304),(0.037,0.00436),(0.037,0.00483),(0.038,0.00548)],
 ("nC6","PSBR23",343.15):[(0.014,0.0139),(0.043,0.0273),(0.063,0.0385),(0.110,0.0488),(0.162,0.0615),(0.262,0.0752)],
 ("nC6","PSBR45",343.15):[(0.016,0.0139),(0.038,0.0273),(0.061,0.0385),(0.103,0.0488),(0.155,0.0615),(0.249,0.0752)],
 ("nC5","PSBR23",333.15):[(0.008,0.0228),(0.016,0.0360),(0.021,0.0468),(0.033,0.0655),(0.046,0.0840),(0.055,0.0973),(0.070,0.1133),(0.078,0.1248),(0.089,0.1376)],
 ("nC5","PSBR45",333.15):[(0.010,0.0228),(0.013,0.0468),(0.019,0.0655),(0.022,0.0840),(0.022,0.0973),(0.026,0.1133),(0.027,0.1248),(0.034,0.1376)],
 ("cC6","PSBR41",343.15):[(0.2169,0.04729),(0.2436,0.05051),(0.2803,0.05431),(0.3288,0.05860),(0.3606,0.06079),(0.3992,0.06303),(0.4414,0.06514),(0.4766,0.06653)],
 ("TOL","PSBR41",343.15):[(0.2381,0.01582),(0.2802,0.01759),(0.3279,0.01935),(0.3724,0.02074),(0.4200,0.02201),(0.4683,0.02308),(0.5039,0.02378),(0.5467,0.02448)],
 ("TOL","PSBR41",373.15):[(0.1832,0.03593),(0.2299,0.04231),(0.2613,0.04609),(0.3030,0.05068),(0.3681,0.05646),(0.4274,0.06066),(0.4937,0.06447),(0.5692,0.06782)],
 ("BZN","PSBR4.1",308.15):[(0.0471,0.00395),(0.0756,0.00593),(0.1044,0.00778),(0.125,0.00877),(0.1442,0.00946),(0.1959,0.01178),(0.2289,0.01284),(0.2637,0.01414),(0.3098,0.01544),(0.3742,0.01649),(0.4297,0.01729)],
 ("BZN","PSBR4.1",298.15):[(0.0459,0.00225),(0.08,0.00396),(0.1117,0.00494),(0.1364,0.00580),(0.1656,0.00672),(0.1864,0.00764),(0.3014,0.00976),(0.37,0.01075),(0.4275,0.01148)],
}

NOMBRES = {"ACE": "acetona", "nC4OH": "1-butanol", "CHL": "cloroformo",
           "cC6": "ciclohexano", "nC6": "n-hexano", "nC5": "n-pentano",
           "TOL": "tolueno", "BZN": "benceno"}


def familia(cop):
    if cop.startswith("PSBR"):  return "PSBR"
    if cop.startswith("PVAE"):  return "PVAE"
    if cop.startswith("PVAVC"): return "PVAVC"
    return cop


# =============================================================================
#  4) CONSTRUCCION DE MODELOS PARA CADA METODOLOGIA
# =============================================================================
def _solvente(solv_key):
    """Solvente como componente de 1 segmento; activa asociacion si aplica."""
    ms = SOLV[solv_key]
    kAB, epsAB = ASOC.get(solv_key, (0.0, 0.0))
    seg = Segmento(solv_key, ms[0] * ms[3], ms[1], ms[2], kAB, epsAB)
    return Componente(solv_key, [seg])


def _w1_a_x1(w1, MW_s, MW_c):
    ns = w1 / MW_s
    nc = (1.0 - w1) / MW_c
    return ns / (ns + nc)


# ---------- METODOLOGIA 1: pseudo-monomero (1 segmento) ----------
def modelos_m1(solv_key, cop_key):
    """Devuelve (modelo_liquido, modelo_vapor, MW_copol) o None si falta kij."""
    if (solv_key, cop_key) not in KIJ_M1:
        return None
    kij, _ = KIJ_M1[(solv_key, cop_key)]
    mc = COP_M1[cop_key]
    MW_c = float(mc[3])
    seg_c = Segmento(cop_key, mc[0] * MW_c, mc[1], mc[2])   # m = (m/MW)*MW
    cop = Componente(cop_key, [seg_c])
    solv = _solvente(solv_key)
    # kseg: solvente(comp 0, seg 0) - copolimero(comp 1, seg 0)
    kseg = {(0, 0, 1, 0): kij}
    ml = PCSAFT([solv, cop], kseg)
    mv = PCSAFT([solv], {})
    return ml, mv, MW_c


# ---------- METODOLOGIA 2: copolimero de 2 segmentos ----------
def modelos_m2(solv_key, cop_key):
    """Devuelve (modelo_liquido, modelo_vapor, MW_copol) o None si no aplica."""
    if cop_key not in COP_M2:
        return None
    fam = (solv_key, familia(cop_key))
    if fam not in KAB_FAM:
        return None
    A, B, wtA, MW_c, Bfr = COP_M2[cop_key]
    wA = wtA / 100.0; wB = 1.0 - wA
    mA = MON[A]; mB = MON[B]                                # parametros de MONOMERO (ver nota)
    segA = Segmento(A, wA * MW_c * mA[0], mA[1], mA[2])     # Ec. 18
    segB = Segmento(B, wB * MW_c * mB[0], mB[1], mB[2])
    B_aa, B_ab, B_bb = Bfr                                  # Tabla 3
    Bmat = [[B_bb,       B_ab / 2.0],                       # seg0=mon1(beta)
            [B_ab / 2.0, B_aa]]                             # seg1=mon2(alpha)
    cop = Componente(cop_key, [segA, segB], Bmat)
    solv = _solvente(solv_key)
    ka, kb = KAB_FAM[fam]
    kab = KAB_AB.get((solv_key, cop_key), 0.0)
    kseg = {(0, 0, 1, 0): ka,    # solvente - mon1
            (0, 0, 1, 1): kb,    # solvente - mon2
            (1, 0, 1, 1): kab}   # mon1 - mon2
    ml = PCSAFT([solv, cop], kseg)
    mv = PCSAFT([solv], {})
    return ml, mv, MW_c


# =============================================================================
#  5) %dP (Ec. 40) Y CURVAS, PARA CUALQUIER METODOLOGIA
# =============================================================================
def dp_system(modelos, solv_key, T, pts):
    """%dP (Ec. 40) de un sistema dado un constructor de modelos (m1 o m2)."""
    res = modelos
    if res is None:
        return None, 0
    ml, mv, MW_c = res
    MW_s = SOLV[solv_key][3]
    s = 0.0; n = 0
    for w1, Pexp in pts:
        x1 = _w1_a_x1(w1, MW_s, MW_c)
        try:
            Pc = presion_de_burbuja(ml, mv, T, x1)
        except Exception:
            Pc = None
        if Pc and Pc > 0:
            s += abs(Pc - Pexp) / Pexp; n += 1
    return (100.0 * s / n if n else None), n


def curva(modelos, solv_key, T, pts, npts=24):
    """Curva P(w1) calculada en una malla fina, para graficar."""
    res = modelos
    if res is None:
        return [], []
    ml, mv, MW_c = res
    MW_s = SOLV[solv_key][3]
    wmin = min(p[0] for p in pts); wmax = max(p[0] for p in pts)
    ws = np.linspace(max(wmin * 0.5, 1e-4), wmax * 1.05, npts)
    wok, Pok = [], []
    for w1 in ws:
        x1 = _w1_a_x1(float(w1), MW_s, MW_c)
        try:
            P = presion_de_burbuja(ml, mv, T, x1)
        except Exception:
            P = None
        if P and P > 0:
            wok.append(float(w1)); Pok.append(P)
    return wok, Pok


# Cache de curvas: evita recalcular la misma curva en la imagen resumen y en
# las figuras individuales.  clave = (metodologia, solv, cop, T).
_CURVAS = {}


def curva_cacheada(metod, solv_key, cop_key, T, pts):
    key = (metod, solv_key, cop_key, T)
    if key not in _CURVAS:
        mods = modelos_m1(solv_key, cop_key) if metod == "M1" \
            else modelos_m2(solv_key, cop_key)
        _CURVAS[key] = curva(mods, solv_key, T, pts) if mods is not None else ([], [])
    return _CURVAS[key]


# =============================================================================
#  6) GRAFICAS: UNA IMAGEN POR SOLVENTE/FAMILIA CON M1 Y M2 SUPERPUESTAS
# =============================================================================
MARKERS = ["o", "s", "^", "v", "D", "P", "*", "X"]
COLORS = plt.cm.tab10.colors


def grafica_resumen_unica(out="resultados_m1_vs_m2.png"):
    """UNA sola imagen (grid de subplots) con TODOS los sistemas: M1 vs M2."""
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), out)
    # un subplot por (solvente, familia), ordenados
    grupos = defaultdict(list)
    for (s, c, T) in DATA:
        grupos[(s, familia(c))].append((c, T))
    claves = sorted(grupos.keys())
    ncols = 3
    nrows = math.ceil(len(claves) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.2 * nrows))
    axes = np.atleast_1d(axes).ravel()
    print("\nGenerando imagen resumen unica:", outpath)
    for ax_i, (solv_key, fam) in enumerate(claves):
        ax = axes[ax_i]
        sistemas = sorted(grupos[(solv_key, fam)])
        for idx, (cop_key, T) in enumerate(sistemas):
            col = COLORS[idx % len(COLORS)]; mk = MARKERS[idx % len(MARKERS)]
            pts = DATA[(solv_key, cop_key, T)]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], mk, color=col,
                    ms=5, markerfacecolor="white", markeredgewidth=1.1,
                    label="%s (%.0f K)" % (cop_key, T))
            w1, P1 = curva_cacheada("M1", solv_key, cop_key, T, pts)
            if w1:
                ax.plot(w1, P1, "-", color=col, lw=1.5)
            w2, P2 = curva_cacheada("M2", solv_key, cop_key, T, pts)
            if w2:
                ax.plot(w2, P2, "--", color=col, lw=1.5)
        ax.set_title("%s (%s) + %s" % (solv_key, NOMBRES.get(solv_key, solv_key), fam),
                     fontsize=10)
        ax.set_xlabel("$w_1$", fontsize=9)
        ax.set_ylabel("$P$ (MPa)", fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=6.5)
        ax.tick_params(labelsize=8)
    # apaga ejes sobrantes
    for k in range(len(claves), len(axes)):
        axes[k].axis("off")
    fig.suptitle("VLE PC-SAFT por sistema: metodologia 1 (linea continua) vs "
                 "metodologia 2 (discontinua)\nsimbolos = datos experimentales | "
                 "Arce, Mattedi & Aznar (2006)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outpath, dpi=140)
    plt.close(fig)
    print("  OK  imagen resumen guardada")


def graficar_m1_vs_m2(out="figuras_m1_vs_m2"):
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), out)
    os.makedirs(outdir, exist_ok=True)
    grupos = defaultdict(list)
    for (s, c, T) in DATA:
        grupos[(s, familia(c))].append((c, T))
    print("\nGenerando figuras M1 vs M2 en:", outdir)
    for (solv_key, fam), sistemas in grupos.items():
        sistemas.sort()
        fig, ax = plt.subplots(figsize=(8, 6))
        for idx, (cop_key, T) in enumerate(sistemas):
            col = COLORS[idx % len(COLORS)]; mk = MARKERS[idx % len(MARKERS)]
            pts = DATA[(solv_key, cop_key, T)]
            # puntos experimentales (marcadores huecos)
            ax.plot([p[0] for p in pts], [p[1] for p in pts], mk, color=col,
                    ms=6, markerfacecolor="white", markeredgewidth=1.3,
                    label="%s exp (T=%.1f K)" % (cop_key, T))
            # curva M1 (linea continua)
            w1, P1 = curva_cacheada("M1", solv_key, cop_key, T, pts)
            if w1:
                ax.plot(w1, P1, "-", color=col, lw=1.8)
            # curva M2 (linea discontinua), solo si el sistema esta en M2
            w2, P2 = curva_cacheada("M2", solv_key, cop_key, T, pts)
            if w2:
                ax.plot(w2, P2, "--", color=col, lw=1.8)
        ax.set_xlabel("Fraccion masica de solvente  $w_1$")
        ax.set_ylabel("Presion  $P$  (MPa)")
        nom = NOMBRES.get(solv_key, solv_key)
        ax.set_title("VLE PC-SAFT: metodologia 1 vs 2\n%s (%s) + %s\n"
                     "simbolos = experimental | linea continua = M1 | discontinua = M2"
                     % (solv_key, nom, fam))
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, title="Copolimero")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "VLE_%s_%s_M1vsM2.png" % (solv_key, fam)), dpi=130)
        plt.close(fig)
        print("  OK  %s + %s" % (solv_key, fam))


# =============================================================================
#  7) TABLA DE ERRORES M1 vs M2  (stdout + imagen)
# =============================================================================
def calcular_tabla():
    """Devuelve lista de filas: (solv, cop, T, dp1, dp1_paper, dp2, dp2_paper)."""
    filas = []
    for (solv, cop, T), pts in DATA.items():
        # M1
        dp1, _ = dp_system(modelos_m1(solv, cop), solv, T, pts)
        dp1_paper = KIJ_M1[(solv, cop)][1] if (solv, cop) in KIJ_M1 else None
        # M2
        dp2, _ = dp_system(modelos_m2(solv, cop), solv, T, pts)
        dp2_paper = DP_FAM_M2.get((solv, familia(cop)), None)
        if cop not in COP_M2:      # sistema sin metodologia 2 en el paper
            dp2 = None; dp2_paper = None
        filas.append((solv, cop, T, dp1, dp1_paper, dp2, dp2_paper))
    return filas


def imprimir_tabla(filas):
    print("\n" + "=" * 78)
    print("  TABLA DE ERRORES  %dP  -  METODOLOGIA 1 vs METODOLOGIA 2")
    print("=" * 78)
    print("%-7s %-9s %7s | %8s %8s | %8s %8s"
          % ("Solv", "Copol", "T[K]", "M1 calc", "M1 pap", "M2 calc", "M2 pap"))
    print("-" * 78)
    s1 = n1 = 0.0; s2 = n2 = 0.0
    for (solv, cop, T, dp1, dp1p, dp2, dp2p) in filas:
        f1 = ("%7.2f" % dp1) if dp1 is not None else "   None"
        f1p = ("%7.2f" % dp1p) if dp1p is not None else "      -"
        f2 = ("%7.2f" % dp2) if dp2 is not None else "   None"
        f2p = ("%7.2f" % dp2p) if dp2p is not None else "      -"
        print("%-7s %-9s %7.2f | %8s %8s | %8s %8s"
              % (solv, cop, T, f1, f1p, f2, f2p))
        if dp1 is not None: s1 += dp1; n1 += 1
        if dp2 is not None: s2 += dp2; n2 += 1
    print("-" * 78)
    g1 = s1 / n1 if n1 else float("nan")
    g2 = s2 / n2 if n2 else float("nan")
    # M2: promedio por familia (criterio con el que el paper reporta 2.04%)
    g2_fam = _global_m2_por_familia(filas)
    print("GLOBAL:")
    print("  Metodologia 1 (prom. por sistema): %.2f%%  sobre %d sistemas   (paper 4.41%%)"
          % (g1, int(n1)))
    print("  Metodologia 2 (prom. por sistema): %.2f%%  sobre %d sistemas   (paper 2.04%%)"
          % (g2, int(n2)))
    print("  Metodologia 2 (prom. por familia): %.2f%%   (criterio del paper, 2.04%%)"
          % g2_fam)
    print("=" * 78)
    return g1, g2, g2_fam, int(n1), int(n2)


def _global_m2_por_familia(filas):
    """Promedio de los promedios por familia (solvente, familia), como el paper."""
    fam_acc = defaultdict(lambda: [0.0, 0])
    for (solv, cop, T, dp1, dp1p, dp2, dp2p) in filas:
        if dp2 is None:
            continue
        key = (solv, familia(cop))
        fam_acc[key][0] += dp2; fam_acc[key][1] += 1
    medias = [s / n for (s, n) in fam_acc.values() if n]
    return sum(medias) / len(medias) if medias else float("nan")


def imagen_tabla(filas, g1, g2, g2_fam, n1, n2, out="tabla_errores_m1_vs_m2.png"):
    """Guarda la tabla de errores como una imagen PNG."""
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), out)
    encabezado = ["Solvente", "Copolimero", "T [K]",
                  "M1 calc %", "M1 paper %", "M2 calc %", "M2 paper %"]
    celdas = []
    for (solv, cop, T, dp1, dp1p, dp2, dp2p) in filas:
        celdas.append([
            solv, cop, "%.2f" % T,
            "%.2f" % dp1 if dp1 is not None else "None",
            "%.2f" % dp1p if dp1p is not None else "-",
            "%.2f" % dp2 if dp2 is not None else "None",
            "%.2f" % dp2p if dp2p is not None else "-",
        ])
    # fila resumen (M2 se reporta con promedio por familia, criterio del paper)
    celdas.append(["GLOBAL", "(prom.)", "",
                   "%.2f" % g1, "4.41", "%.2f" % g2_fam, "2.04"])

    fig_h = 0.42 * (len(celdas) + 2)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.axis("off")
    ax.set_title("Desviacion media en presion  %dP  (Ec. 40)\n"
                 "PC-SAFT metodologia 1 vs 2 - Arce et al. (2006)",
                 fontsize=13, fontweight="bold", pad=14)
    tabla = ax.table(cellText=celdas, colLabels=encabezado,
                     cellLoc="center", loc="center")
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(9)
    tabla.scale(1.0, 1.35)
    # estilo encabezado
    for j in range(len(encabezado)):
        c = tabla[0, j]
        c.set_facecolor("#2f5d62"); c.set_text_props(color="white", fontweight="bold")
    # estilo fila global (ultima)
    r = len(celdas)
    for j in range(len(encabezado)):
        c = tabla[r, j]
        c.set_facecolor("#dbe9e9"); c.set_text_props(fontweight="bold")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Imagen de la tabla de errores guardada en:", outpath)


# =============================================================================
#  8) PROGRAMA PRINCIPAL
# =============================================================================
def main():
    filas = calcular_tabla()
    g1, g2, g2_fam, n1, n2 = imprimir_tabla(filas)
    imagen_tabla(filas, g1, g2, g2_fam, n1, n2)   # tabla de errores -> PNG
    grafica_resumen_unica()                        # UNA imagen con todos los sistemas
    graficar_m1_vs_m2()                            # ademas, una figura por solvente


if __name__ == "__main__":
    main()
