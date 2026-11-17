#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
  VLE de mezclas COPOLIMERO + SOLVENTE con la EoS PC-SAFT
  Reproduccion FIEL del articulo (metodologia 1: pseudo-monomero)
===============================================================================

    Arce, Mattedi & Aznar (2006), Fluid Phase Equilibria 246, 52-63.
    "Vapor-liquid equilibrium of copolymer + solvent mixtures ..."

Este codigo sigue PASO POR PASO lo que dice el paper, sin anadir nada extra:

  - PC-SAFT estandar de Gross & Sadowski (Ecs. 21-35 del paper).
  - El copolimero se trata como UN pseudo-monomero (concepto usado en la
    seccion 5 y Tabla 5 del paper para obtener la Tabla 6, metodologia 1).
    Es decir, el copolimero es un UNICO segmento con los parametros PC-SAFT
    (m/MW, sigma, eps/k) que el paper ya da en la Tabla 5.
  - Como el pseudo-monomero es un solo tipo de segmento, NO hay fracciones de
    enlace B ni termino alpha-beta: el copolimero es una cadena de un solo
    tipo (igual que un homopolimero en PC-SAFT).
  - Un unico parametro de interaccion binaria solvente-copolimero (kij),
    tal como en la Tabla 6 (metodologia 1).
  - Condicion de equilibrio (Ec. 39): el copolimero no es volatil, el vapor
    es solvente puro:   x1 * phi_hat_1^L = phi_1^V(y1=1).
  - Las constantes universales aik, bik son las de Gross & Sadowski (2001),
    ref. [23] del paper (el paper remite a esa fuente, Ec. 31-32).

NOTA: el 1-butanol (nC4OH) es autoasociante; el paper le anade un termino de
asociacion (Ecs. 36-38). Ese termino NO se incluye aqui (por eso nC4OH quedaria
aproximado). Para los demas solventes el modelo es el completo del paper.

Dependencias:  numpy, scipy, matplotlib
Uso:           python pcsaft_vle.py
===============================================================================
"""

import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

import math
import numpy as np
from scipy.optimize import brentq


# =============================================================================
#  ////////////////////  PARAMETROS EDITABLES  //////////////////////////////
#  Valores por defecto:  ACETONA (1) + PSBR5 (2) a 323.15 K   (Fig. 2)
#  Todos tomados del paper (Tabla 4b para el solvente, Tabla 5 para el
#  pseudo-monomero del copolimero, Tabla 6 para kij).
# =============================================================================

# ---- 1) TEMPERATURA y rango de composicion --------------------------------
TEMPERATURA_K = 333.15            # Temperatura [K]
W1_MIN        = 0.011             # Fraccion masica de solvente minima
W1_MAX        = 0.301             # Fraccion masica de solvente maxima
N_PUNTOS      = 80                # Numero de puntos de la curva estimada

# ---- 2) SOLVENTE (componente 1)  --  Tabla 4b -----------------------------
#   m/MW [1/mol], sigma [Angstrom], eps/k [K], MW [g/mol].
SOLVENTE = {
    "nombre": "cC6 (ciclohexano)",
    "MW":     84.16,              # g/mol
    "m_MW":   0.029946,           # (m/MW) [1/mol]
    "sigma":  3.8399,             # Angstrom
    "eps_k":  279.21,             # eps/k [K]
}

# ---- 3) COPOLIMERO como PSEUDO-MONOMERO (componente 2)  --  Tabla 5 --------
#   El paper representa el copolimero por UN pseudo-monomero con estos tres
#   parametros PC-SAFT (Tabla 5). MW es el peso molecular total del copolimero
#   (Tabla 3). El pseudo-monomero se trata como un solo tipo de segmento.
COPOLIMERO = {
    "nombre": "PSBR45 (poly(styrene-co-butadiene), 45% estireno)",
    "MW":     600000.0,           # g/mol  (peso molecular total del copolimero)
    "m_MW":   0.039562,           # (m/MW) del pseudo-monomero [1/mol]  (Tabla 5)
    "sigma":  3.6072,             # sigma del pseudo-monomero [Angstrom] (Tabla 5)
    "eps_k":  246.00,             # eps/k del pseudo-monomero [K]        (Tabla 5)
}

# ---- 4) PARAMETRO DE INTERACCION BINARIA kij  --  Tabla 6 (metodologia 1) --
#   Un unico kij solvente-copolimero (Ec. 33: eps_12 = sqrt(e1 e2)(1-kij)).
KIJ = 0.0101                      # cC6+PSBR45, Tabla 6

# ---- 5) DATOS EXPERIMENTALES (opcional) -----------------------------------
#   Pega aqui los puntos del sistema (de datos_experimentales.txt), en MPa.
#   Formato (w1, P[MPa]). Si kPa, divide entre 1000. Lista vacia = sin puntos.
#   Ejemplo: cC6+PSBR45 a 333.15 K (Gupta & Prausnitz [1]), kPa -> MPa.
DATOS_EXPERIMENTALES = [
    (0.021,  0.0048),
    (0.050,  0.0105),
    (0.087,  0.0169),
    (0.103,  0.0199),
    (0.130,  0.0248),
    (0.163,  0.0295),
    (0.200,  0.0337),
    (0.240,  0.0376),
    (0.282,  0.0407),
]

# ---- 6) SALIDA ------------------------------------------------------------
ARCHIVO_FIGURA = "resultado_VLE.png"

# =============================================================================
#  ///////////////  FIN DE PARAMETROS EDITABLES  ////////////////////////////
#
#  A PARTIR DE AQUI NO HAY NINGUN PARAMETRO DEL SISTEMA.
#  Todo lo que depende del solvente/copolimero/temperatura esta ARRIBA (1-6).
#  Lo que sigue son CONSTANTES FIJAS (universales de PC-SAFT y de la fisica)
#  y el codigo de calculo. NO hace falta editar nada de aqui para abajo.
# =============================================================================


# -----------------------------------------------------------------------------
#  CONSTANTES FIJAS -- NO EDITAR
#  Constantes universales aik, bik importadas de constantes_pcsaft.py
#  (fuente externa: Gross & Sadowski 2001, ref. [23] del paper). Universales.
# -----------------------------------------------------------------------------
from constantes_pcsaft import A_COEF, B_COEF

KB = 1.380649e-23        # J/K (constante de Boltzmann) -- CONSTANTE FIJA, NO EDITAR


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
def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    solv = construir_solvente()
    cop = construir_copolimero()
    kij_m = [[0.0, KIJ], [KIJ, 0.0]]
    model_liq = PCSAFT([solv, cop], kij_m)
    model_vap = PCSAFT([solv], [[0.0]])

    T = TEMPERATURA_K
    ws = np.linspace(max(W1_MIN, 1e-4), W1_MAX, N_PUNTOS)

    print("=" * 66)
    print("  VLE PC-SAFT (pseudo-monomero, metodologia 1 del paper)")
    print("  %s (1) + %s (2)" % (SOLVENTE["nombre"], COPOLIMERO["nombre"]))
    print("  T = %.2f K   kij = %.4f" % (T, KIJ))
    print("  m_solvente = %.3f   m_copolimero(pseudo) = %.1f"
          % (solv.m, cop.m))
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

    # ---- desviacion de presion %DeltaP/P (Ec. 40) si hay experimentales ----
    if DATOS_EXPERIMENTALES:
        import numpy as _np
        wE = _np.array([d[0] for d in DATOS_EXPERIMENTALES])
        PE = _np.array([d[1] for d in DATOS_EXPERIMENTALES])
        s = 0.0; n = 0
        for w1, Pexp in DATOS_EXPERIMENTALES:
            x1 = w1_a_x1(w1)
            Pc = presion_de_burbuja(model_liq, model_vap, T, x1)
            if Pc:
                s += abs(Pc - Pexp) / Pexp; n += 1
        if n:
            print("\n  %%DeltaP/P (Ec. 40) = %.2f %%   (np=%d)" % (100 * s / n, n))

    # ---- grafica ----
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(w_ok, P_ok, "-", color="#3b0a70", lw=2, label="PC-SAFT (estimado)")
    if DATOS_EXPERIMENTALES:
        ax.plot([d[0] for d in DATOS_EXPERIMENTALES],
                [d[1] for d in DATOS_EXPERIMENTALES],
                "o", color="#d62728", ms=6, markerfacecolor="white",
                markeredgewidth=1.4, label="Experimental")
    ax.set_xlabel("Fraccion masica de solvente  $w_1$")
    ax.set_ylabel("Presion  $P$  (MPa)")
    ax.set_title("VLE PC-SAFT\n%s (1) + %s (2)   T=%.2f K"
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
