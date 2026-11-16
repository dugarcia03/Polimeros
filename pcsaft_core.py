"""
PC-SAFT core engine.

Implements the perturbed-chain SAFT equation of state exactly as summarized in
Appendix A of:
    Gross, J.; Sadowski, G. Ind. Eng. Chem. Res. 2001, 40, 1244-1260
    (reference [23] of Arce et al. 2006),
and equivalently written in Eqs. 18-40 of Arce et al. (2006).

This module works for a MIXTURE of homonuclear chains, each component i having
pure parameters (m_i, sigma_i, eps_i/k). That is exactly what Methodology 1
(pseudo-monomer) requires: the solvent is one such chain and the copolymer is a
single pseudo-monomer chain. One binary interaction parameter kij corrects the
unlike dispersion energy (Eq. 24 / A.15 / Eq. 33 of Arce).

EQUATION MAP (Arce 2006  <->  Gross&Sadowski 2001 Appendix A  <->  code):
  a_res = a_hc + a_disp                         Eq.21 / A.3           -> a_res()
  a_hc  = mbar*a_hs - sum xi(mi-1) ln g_ii       Eq.22 / A.4           -> a_hc()
  mbar  = sum xi mi                              Eq.24 / A.5           -> mbar
  a_hs  (per segment)                            A.6                   -> a_hs()
  g_ij^hs                                        Eq.23 / A.7           -> g_hs()
  zeta_n = (pi/6) rho sum xi mi di^n             Eq.25 / A.8           -> zeta()
  d_i = sigma_i[1 - 0.12 exp(-3 eps_i/kT)]       Eq.26 / A.9           -> d_i()
  a_disp = a1 + a2                               Eq.27 / A.10          -> a_disp()
  a1 = -2 pi rho I1 (m2eps sig3)                 Eq.28 / A.10,A.12     -> a_disp()
  a2 = -pi rho mbar C1 I2 (m2eps2 sig3)          Eq.29 / A.10,A.13     -> a_disp()
  C1 (compressibility)                           A.11 / Eq.30          -> C1()
  m2es3 , m2e2s3 (VDW one-fluid sums)            A.12,A.13 / Eq.28,29  -> m2es3()
  sigma_ij=(si+sj)/2 ; eps_ij=sqrt(ei ej)(1-kij) Eq.23,24 / A.14,A.15  -> combine
  I1,I2 power series in eta                      A.16,A.17            -> I1I2()
  ai(mbar), bi(mbar)                             Eq.31,32 / A.18,A.19 -> ab_coeff()
  Z = 1 + Zhc + Zdisp                            A.24                  -> Z_from_eta
  Zhc, Zhs, rho dg/drho                          A.25,A.26,A.27        -> Z_hc()
  Zdisp, d(eta I1)/deta, d(eta I2)/deta, C2      A.28-A.31             -> Z_disp()
  P = Z k T rho (1e10)^3                          A.23                  -> pressure()
  density iteration: solve Pcalc(eta)=Psys       Appendix A "Density"  -> solve_eta
  ln phi_k  (fugacity coefficient)               A.32-A.45             -> ln_phi()

All reduced-density / packing-fraction relations use eta = zeta_3.
"""

import numpy as np
from constantes_universales import A_COEFFS, B_COEFFS, N_AV, K_B


# ----------------------------------------------------------------------------
# Low-level building blocks.  Component arrays are 1-D numpy arrays.
# x   : mole fractions            (n,)
# m   : segment numbers           (n,)
# sig : segment diameters [A]     (n,)
# eps : eps/k [K]                 (n,)
# ----------------------------------------------------------------------------

def d_i(m, sig, eps, T):
    """Temperature-dependent segment diameter d_i(T). Eq.26 / A.9."""
    return sig * (1.0 - 0.12 * np.exp(-3.0 * eps / T))


def mbar(x, m):
    """Mean segment number mbar = sum xi mi. Eq.24 / A.5."""
    return np.sum(x * m)


def zeta_n(x, m, d, rho):
    """zeta_0..zeta_3 = (pi/6) rho sum xi mi di^n. Eq.25 / A.8.
    rho in 1/A^3. Returns array [z0,z1,z2,z3]."""
    out = np.empty(4)
    for n in range(4):
        out[n] = (np.pi / 6.0) * rho * np.sum(x * m * d**n)
    return out


def ab_coeff(mb):
    """ai(mbar), bi(mbar) power-series coefficients. Eq.31,32 / A.18,A.19."""
    f1 = (mb - 1.0) / mb
    f2 = (mb - 1.0) / mb * (mb - 2.0) / mb
    a = A_COEFFS[:, 0] + f1 * A_COEFFS[:, 1] + f2 * A_COEFFS[:, 2]
    b = B_COEFFS[:, 0] + f1 * B_COEFFS[:, 1] + f2 * B_COEFFS[:, 2]
    return a, b


def I1I2(eta, a, b):
    """I1,I2 power series in eta. A.16,A.17."""
    powers = eta ** np.arange(7)
    return np.sum(a * powers), np.sum(b * powers)


def dI1I2_deta(eta, a, b):
    """d(eta I1)/deta and d(eta I2)/deta. A.29,A.30."""
    j = np.arange(7)
    powers = eta ** j
    d1 = np.sum(a * (j + 1) * powers)
    d2 = np.sum(b * (j + 1) * powers)
    return d1, d2


def C1_C2(mb, eta):
    """C1 (A.11) and C2 = dC1/deta (A.31)."""
    omm = 1.0 - mb
    C1_inv = (1.0
              + mb * (8.0 * eta - 2.0 * eta**2) / (1.0 - eta)**4
              + omm * (20.0 * eta - 27.0 * eta**2 + 12.0 * eta**3 - 2.0 * eta**4)
              / ((1.0 - eta) * (2.0 - eta))**2)
    C1 = 1.0 / C1_inv
    C2 = -C1**2 * (mb * (-4.0 * eta**2 + 20.0 * eta + 8.0) / (1.0 - eta)**5
                   + omm * (2.0 * eta**3 + 12.0 * eta**2 - 48.0 * eta + 40.0)
                   / ((1.0 - eta) * (2.0 - eta))**3)
    return C1, C2


def m2es3_terms(x, m, sig, eps, T, kij):
    """VDW one-fluid dispersion sums:
       m2es3  = sum_i sum_j xi xj mi mj (eps_ij/kT) sigma_ij^3     (A.12)
       m2e2s3 = sum_i sum_j xi xj mi mj (eps_ij/kT)^2 sigma_ij^3   (A.13)
    with combining rules sigma_ij=(si+sj)/2, eps_ij=sqrt(ei ej)(1-kij) (A.14,A.15).
    Returns (m2es3, m2e2s3)."""
    n = len(x)
    s1 = 0.0
    s2 = 0.0
    for i in range(n):
        for j in range(n):
            sig_ij = 0.5 * (sig[i] + sig[j])
            eps_ij = np.sqrt(eps[i] * eps[j]) * (1.0 - kij[i, j])
            base = x[i] * x[j] * m[i] * m[j] * sig_ij**3
            s1 += base * (eps_ij / T)
            s2 += base * (eps_ij / T)**2
    return s1, s2


def g_hs_ii(d, zeta):
    """Hard-sphere radial distribution function at contact for i-i pairs.
    Eq.23 / A.7. Returns array g_ii."""
    z2, z3 = zeta[2], zeta[3]
    om = 1.0 - z3
    # for i-i: di*dj/(di+dj) = di/2
    r = d / 2.0
    return 1.0 / om + r * 3.0 * z2 / om**2 + r**2 * 2.0 * z2**2 / om**3


# ----------------------------------------------------------------------------
# Helmholtz energy contributions (reduced, per molecule: a = A/NkT)
# ----------------------------------------------------------------------------

def a_hs_per_segment(zeta):
    """Hard-sphere Helmholtz energy per segment. A.6."""
    z0, z1, z2, z3 = zeta
    om = 1.0 - z3
    return (1.0 / z0) * (3.0 * z1 * z2 / om
                         + z2**3 / (z3 * om**2)
                         + (z2**3 / z3**2 - z0) * np.log(om))


def a_hc(x, m, d, zeta):
    """Hard-chain contribution. Eq.22 / A.4."""
    mb = np.sum(x * m)
    ahs = a_hs_per_segment(zeta)
    g = g_hs_ii(d, zeta)
    return mb * ahs - np.sum(x * (m - 1.0) * np.log(g))


def a_disp(x, m, sig, eps, T, kij, rho, eta):
    """Dispersion contribution. Eq.27-29 / A.10."""
    mb = np.sum(x * m)
    a, b = ab_coeff(mb)
    I1, I2 = I1I2(eta, a, b)
    C1, _ = C1_C2(mb, eta)
    m2es3, m2e2s3 = m2es3_terms(x, m, sig, eps, T, kij)
    a1 = -2.0 * np.pi * rho * I1 * m2es3
    a2 = -np.pi * rho * mb * C1 * I2 * m2e2s3
    return a1 + a2


# ----------------------------------------------------------------------------
# Compressibility factor Z and pressure
# ----------------------------------------------------------------------------

def Z_hc(x, m, d, zeta):
    """Hard-chain compressibility factor contribution. A.25,A.26,A.27."""
    mb = np.sum(x * m)
    z0, z1, z2, z3 = zeta
    om = 1.0 - z3
    Zhs = (z3 / om
           + 3.0 * z1 * z2 / (z0 * om**2)
           + (3.0 * z2**3 - z3 * z2**3) / (z0 * om**3))
    g = g_hs_ii(d, zeta)
    r = d / 2.0  # di*dj/(di+dj) for i-i
    # rho * d g_ii / d rho  (A.27); note zeta_n propto rho so rho d/drho scales each zeta term
    rho_dg = (z3 / om**2
              + r * (3.0 * z2 / om**2 + 6.0 * z2 * z3 / om**3)
              + r**2 * (4.0 * z2**2 / om**3 + 6.0 * z2**2 * z3 / om**4))
    return mb * Zhs - np.sum(x * (m - 1.0) * (1.0 / g) * rho_dg)


def Z_disp(x, m, sig, eps, T, kij, rho, eta):
    """Dispersion compressibility factor contribution. A.28-A.31."""
    mb = np.sum(x * m)
    a, b = ab_coeff(mb)
    I2, = (I1I2(eta, a, b)[1],)
    detaI1, detaI2 = dI1I2_deta(eta, a, b)
    C1, C2 = C1_C2(mb, eta)
    m2es3, m2e2s3 = m2es3_terms(x, m, sig, eps, T, kij)
    Zdisp1 = -2.0 * np.pi * rho * detaI1 * m2es3
    Zdisp2 = -np.pi * rho * mb * (C1 * detaI2 + C2 * eta * I2) * m2e2s3
    return Zdisp1 + Zdisp2


def Z_total(x, m, sig, eps, T, kij, eta):
    """Z = 1 + Zhc + Zdisp. A.24. Computes rho from eta internally."""
    d = d_i(m, sig, eps, T)
    # rho from eta = zeta_3 (A.20): rho = (6/pi) eta / sum xi mi di^3
    rho = (6.0 / np.pi) * eta / np.sum(x * m * d**3)
    zeta = zeta_n(x, m, d, rho)
    zhc = Z_hc(x, m, d, zeta)
    zdisp = Z_disp(x, m, sig, eps, T, kij, rho, eta)
    return 1.0 + zhc + zdisp, rho


def pressure_Pa(x, m, sig, eps, T, kij, eta):
    """Pressure [Pa]. P = Z k T rho (1e10)^3. A.23 (rho in 1/A^3)."""
    Z, rho = Z_total(x, m, sig, eps, T, kij, eta)
    return Z * K_B * T * rho * (1e10)**3


def solve_eta(x, m, sig, eps, T, kij, P_target, phase="liquid"):
    """Find eta such that pressure(eta) = P_target [Pa].
    Starting guesses per Appendix A: liquid eta~0.5, vapor eta~1e-10.
    The liquid root is the DENSEST eta with Pcalc=P_target; the vapor root is
    the most dilute one."""
    from scipy.optimize import brentq

    def f(eta):
        return pressure_Pa(x, m, sig, eps, T, kij, eta) - P_target

    if phase == "liquid":
        etas = np.linspace(0.05, 0.7405, 80)
        vals = np.array([f(e) for e in etas])
        sign = np.sign(vals)
        idx = np.where(np.diff(sign) != 0)[0]
        if len(idx) == 0:
            return etas[np.argmin(np.abs(vals))]
        i = idx[-1]  # densest sign change = liquid
        return brentq(f, etas[i], etas[i + 1], xtol=1e-12, rtol=1e-12)
    else:
        etas = np.linspace(1e-10, 0.3, 80)
        vals = np.array([f(e) for e in etas])
        sign = np.sign(vals)
        idx = np.where(np.diff(sign) != 0)[0]
        if len(idx) == 0:
            return etas[np.argmin(np.abs(vals))]
        i = idx[0]  # most dilute sign change = vapor
        return brentq(f, etas[i], etas[i + 1], xtol=1e-14, rtol=1e-12)


# ----------------------------------------------------------------------------
# Fugacity coefficient  (A.32 - A.45)
# ----------------------------------------------------------------------------

def a_res_reduced(x, m, sig, eps, T, kij, eta):
    """Residual Helmholtz energy a_res = A_res/(N k T) (per molecule).
    a_res = a_hc + a_disp  (Eq.21 / A.3). eta = zeta_3 is the packing fraction."""
    d = d_i(m, sig, eps, T)
    rho = (6.0 / np.pi) * eta / np.sum(x * m * d**3)
    zeta = zeta_n(x, m, d, rho)
    return a_hc(x, m, d, zeta) + a_disp(x, m, sig, eps, T, kij, rho, eta)


def ln_phi(x, m, sig, eps, T, kij, eta):
    """ln phi_k for all components at converged eta. Returns array (n,).

    Uses the standard thermodynamic identity (A.32-A.33 of Gross & Sadowski):
        ln phi_k = mu_k^res/kT - ln Z,
        mu_k^res/kT = d(A_res/kT)/d n_k  at constant T, V, n_{j!=k}.
    The derivative wrt mole number n_k is evaluated by central finite
    differences on A_res/kT = N * a_res, with the system volume V held fixed
    (so the molecular number density rho = N/V changes as n_k changes, and
    eta tracks rho accordingly). This is mathematically equivalent to the
    closed-form A.33-A.45 derivatives but avoids any transcription risk; the
    Helmholtz energy itself (Eqs. 21-32) is the faithfully transcribed core.
    """
    n = len(x)
    d = d_i(m, sig, eps, T)
    rho = (6.0 / np.pi) * eta / np.sum(x * m * d**3)  # 1/A^3
    V = 1.0 / rho                                     # volume for N=1 molecule
    nvec = x.copy()                                   # mole-number vector (N=1)

    Z, _ = Z_total(x, m, sig, eps, T, kij, eta)

    def A_over_kT(nv):
        N = np.sum(nv)
        xx = nv / N
        r = N / V                                     # number density, 1/A^3
        eta_l = (np.pi / 6.0) * r * np.sum(xx * m * d**3)
        return N * a_res_reduced(xx, m, sig, eps, T, kij, eta_l)

    mu_res = np.empty(n)
    for k in range(n):
        h = 1e-7 * max(nvec[k], 1e-6)
        np_ = nvec.copy(); np_[k] += h
        nm_ = nvec.copy(); nm_[k] -= h
        mu_res[k] = (A_over_kT(np_) - A_over_kT(nm_)) / (2.0 * h)

    return mu_res - np.log(Z)


def ln_phi_pure(m_i, sig_i, eps_i, T, P_target, phase):
    """ln phi for a pure component (single chain). Convenience wrapper."""
    x = np.array([1.0])
    m = np.array([m_i]); sig = np.array([sig_i]); eps = np.array([eps_i])
    kij = np.zeros((1, 1))
    eta = solve_eta(x, m, sig, eps, T, kij, P_target, phase)
    return ln_phi(x, m, sig, eps, T, kij, eta)[0], eta


# ----------------------------------------------------------------------------
# Binary VLE: solvent(1) + nonvolatile copolymer(2).  Vapor = pure solvent.
# Equilibrium condition (Arce Eq. 39):  x1 phiL_1(T,P,x) = phiV_1(T,P,y1=1).
# Taking logs:  ln(x1) + ln phiL_1(T,P) = ln phiV_1(T,P).
# ----------------------------------------------------------------------------

def bubble_pressure(x1, T, solvent_params, copolymer_params, kij,
                    P_hint=None):
    """Bubble pressure [Pa] of a binary solvent(1)+copolymer(2) mixture.

    solvent_params   = (m1, sigma1, eps1, MW1)
    copolymer_params = (m2, sigma2, eps2, MW2)   (pseudo-monomer chain)
    x1               = solvent MOLE fraction in the liquid.
    kij              = binary interaction parameter (Eq.24/A.15/Eq.33).
    P_hint           = optional pressure bracket hint [Pa] (e.g. pure Psat).
    """
    from scipy.optimize import brentq

    m1, s1, e1, _ = solvent_params
    m2, s2, e2, _ = copolymer_params
    m = np.array([m1, m2]); sig = np.array([s1, s2]); eps = np.array([e1, e2])
    kijm = np.array([[0.0, kij], [kij, 0.0]])
    xL = np.array([x1, 1.0 - x1])

    def resid(P):
        etaL = solve_eta(xL, m, sig, eps, T, kijm, P, "liquid")
        lnphiL = ln_phi(xL, m, sig, eps, T, kijm, etaL)[0]
        lnphiV = ln_phi_pure(m1, s1, e1, T, P, "vapor")[0]
        return np.log(x1) + lnphiL - lnphiV

    # Bracket: the bubble pressure lies between ~0 and the pure-solvent Psat.
    # Use P_hint to size the bracket; otherwise scan logarithmically.
    if P_hint is None:
        P_hint = 1.0e5
    lo = max(P_hint * 1e-4, 1.0)
    hi = P_hint * 2.0
    # ensure a sign change by expanding if needed
    flo, fhi = resid(lo), resid(hi)
    tries = 0
    while flo * fhi > 0 and tries < 40:
        hi *= 1.5
        fhi = resid(hi)
        tries += 1
    if flo * fhi > 0:
        # fall back to coarse scan
        Ps = np.logspace(np.log10(lo), np.log10(hi), 60)
        vals = np.array([resid(P) for P in Ps])
        idx = np.where(np.diff(np.sign(vals)) != 0)[0]
        if len(idx) == 0:
            raise RuntimeError("bubble_pressure: no bracket found")
        i = idx[0]
        return brentq(resid, Ps[i], Ps[i + 1], xtol=1.0, rtol=1e-10)
    return brentq(resid, lo, hi, xtol=1.0, rtol=1e-10)
