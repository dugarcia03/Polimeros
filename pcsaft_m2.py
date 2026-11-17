# -*- coding: utf-8 -*-
"""
===============================================================================
  PC-SAFT para COPOLIMEROS - METODOLOGIA 2 (segmentos alpha/beta)
===============================================================================

  Arce, Mattedi & Aznar (2006), Fluid Phase Equilibria 246, 52-63.
  Extension a copolimeros de Gross, Spuhl, Tumakaka & Sadowski (2003), ref [24].

  A diferencia de la metodologia 1 (pseudo-monomero = un solo segmento), aqui el
  copolimero se modela como una cadena con DOS tipos de segmento (alpha y beta,
  los dos monomeros), con:

    - Ec. 18:  m_{i,a} = w_{i,a} * MW_copol * (m/MW)_{i,a}   (segmentos tipo a)
    - Ec. 19:  m_i = sum_a m_{i,a}                            (segmentos totales)
    - Ec. 20:  z_{i,a} = m_{i,a} / m_i                        (fraccion de segmento)
    - Ec. 22:  a_hc con doble suma sum_a sum_b B_{ia,ib} ln g_{ia,ib}
               (fracciones de enlace B de la Tabla 3)
    - Ecs. 28-29: dispersion con doble suma sobre segmentos  z_ia z_jb
    - Ec. 33:  reglas cruzadas con kij POR PAR DE SEGMENTOS (Tabla 7):
                 k_solv-mon1 (kappa_a), k_solv-mon2 (kappa_b), k_mon1-mon2 (kappa_ab)

  Las constantes universales aik/bik y la constante de Boltzmann se reutilizan de
  constantes_pcsaft.py (fuente externa: Gross & Sadowski 2001, ref [23]).

  Estructura de datos:
    Un "Componente" es una lista de segmentos. Cada segmento = (m, sigma, eps_k).
    Para un solvente: 1 segmento.  Para un copolimero: 2 segmentos (alpha, beta).
    z[i] = lista de fracciones de segmento de cada segmento del componente i
           (para solvente z=[1.0]).
    B[i] = matriz de fracciones de enlace B_{a,b} del componente i
           (para solvente B=[[1.0]] ... en realidad no hay enlaces, m puede ser 1).
===============================================================================
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import math
import numpy as np
from scipy.optimize import brentq

from constantes_pcsaft import A_COEF, B_COEF

KB = 1.380649e-23        # J/K


class Segmento:
    """Un tipo de segmento: parametros PC-SAFT de un monomero o de un solvente."""
    def __init__(self, nombre, m, sigma, eps_k, kAB=0.0, epsAB=0.0):
        self.nombre = nombre
        self.m = m            # numero de segmentos de este tipo en la molecula
        self.sigma = sigma    # Angstrom
        self.eps_k = eps_k    # K
        self.kAB = kAB
        self.epsAB = epsAB
        self.asocia = (kAB > 0.0 and epsAB > 0.0)


class Componente:
    """
    Molecula = cadena de uno o varios tipos de segmento.
    segmentos: lista de Segmento
    B: matriz (ns x ns) de fracciones de enlace B_{a,b} (Tabla 3).
       Para un homopolimero/solvente de un segmento: B=[[1.0]] no se usa (m-1 term).
    """
    def __init__(self, nombre, segmentos, B=None):
        self.nombre = nombre
        self.seg = segmentos
        self.ns = len(segmentos)
        self.m = sum(s.m for s in segmentos)              # Ec. 19: m_i total
        self.z = [s.m / self.m for s in segmentos]        # Ec. 20: z_{i,a}
        if B is None:
            # un solo segmento: no hay enlaces alpha-beta; B_aa=1
            B = [[1.0]]
        self.B = B


class PCSAFT_M2:
    """
    Mezcla de 'Componente'. Cada componente puede tener varios tipos de segmento.
    kseg: diccionario {(i, a, j, b): kij} con el parametro de interaccion entre
          el segmento a del componente i y el segmento b del componente j.
          Simetrico. Si falta un par, se asume 0.
    """
    def __init__(self, componentes, kseg=None):
        self.c = componentes
        self.nc = len(componentes)
        self.kseg = kseg or {}

    # ---- helpers de segmento ----
    def _k(self, i, a, j, b):
        return self.kseg.get((i, a, j, b),
               self.kseg.get((j, b, i, a), 0.0))

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
                # homopolimero/solvente: -x_i (m_i-1) ln g_ii
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
        return 1.0 / (1.0
                      + mb * (8.0 * eta - 2.0 * eta ** 2) / (1.0 - eta) ** 4
                      + (1.0 - mb) * (20.0 * eta - 27.0 * eta ** 2
                                      + 12.0 * eta ** 3 - 2.0 * eta ** 4)
                      / ((1.0 - eta) * (2.0 - eta)) ** 2)

    def _m2es(self, x, T, power):
        """
        Ec. 28/29 (Arce) y A.12/A.13 (GS) extendidas a copolimeros (ref [24]):
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
        a1 = -2.0 * math.pi * rho * I1 * self._m2es(x, T, 1)
        a2 = -math.pi * rho * mb * C1 * I2 * self._m2es(x, T, 2)
        return a1 + a2

    # ---- asociacion (Ecs. 36-38), esquema 2B, segmento asociante ----
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
        a = self.a_hc(x, rho, T) + self.a_disp(x, rho, T)
        if any(s.asocia for ci in self.c for s in ci.seg):
            a += self.a_assoc(x, rho, T)
        return a

    def Z(self, x, rho, T):
        drho = 1e-6 * rho
        return 1.0 + rho * (self.a_res(x, rho + drho, T)
                            - self.a_res(x, rho - drho, T)) / (2.0 * drho)

    def pressure(self, x, rho, T):
        return self.Z(x, rho, T) * (rho * 1e30) * KB * T

    def ln_phi(self, x, rho, T):
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
#  Densidad a (T,P)
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
                return brentq(P, lo, hi, xtol=1e-16, rtol=1e-10, maxiter=100)
            except Exception:
                return 0.5 * (lo + hi)
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
        return (math.log(x1) + lnphi_l) - lnphi_v

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
