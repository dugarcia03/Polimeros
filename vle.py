"""
VLE driver for Methodology 1 (pseudo-monomer) PC-SAFT, Arce et al. (2006).

Builds the copolymer pseudo-monomer chain from the monomer PC-SAFT parameters
(Table 4b) using mole-fraction-weighted combining (Eqs. 13-15 of Arce; these
reproduce the published Table 5 pseudo-monomer values exactly), then computes
the bubble pressure P vs solvent mass fraction w1 and the pressure deviation
%dP/P (Eq. 40).
"""

import numpy as np
import pcsaft_core as pc
from parametros_articulo import SOLVENTS, MONOMERS, COPOLYMERS


def pseudomonomer_params(copolymer):
    """Return (m_total, sigma, eps/k, MW) for the copolymer pseudo-monomer.

    Methodology 1:
      x_i  = mole fraction of monomer i (from its weight fraction, Eq. 13 area)
      m/MW, sigma, eps/k : mole-fraction arithmetic averages of the monomers
                           (reproduces Table 5 of Arce 2006 exactly)
      m_total = (m/MW)_psm * MW_copolymer   (Eq. 18-19 collapsed to one chain)
    """
    info = COPOLYMERS[copolymer]
    mon1, mon2 = info["monomers"]
    w1 = info["wt_first"]
    w2 = 1.0 - w1
    mMW1, sig1, eps1, MW1 = MONOMERS[mon1]
    mMW2, sig2, eps2, MW2 = MONOMERS[mon2]
    n1 = w1 / MW1
    n2 = w2 / MW2
    x1 = n1 / (n1 + n2)
    x2 = n2 / (n1 + n2)
    mMW = x1 * mMW1 + x2 * mMW2
    sig = x1 * sig1 + x2 * sig2
    eps = x1 * eps1 + x2 * eps2
    MWc = info["MW"]
    m_total = mMW * MWc
    return m_total, sig, eps, MWc


def solvent_params(solvent):
    """Return (m, sigma, eps/k, MW) for a solvent (Table 4b)."""
    mMW, sig, eps, MW = SOLVENTS[solvent]
    return mMW * MW, sig, eps, MW


def w1_to_x1(w1, MW_solvent, MW_copolymer):
    """Convert solvent MASS fraction w1 to solvent MOLE fraction x1."""
    n1 = w1 / MW_solvent
    n2 = (1.0 - w1) / MW_copolymer
    return n1 / (n1 + n2)


def bubble_P_MPa(w1, T, solvent, copolymer, kij, P_hint_Pa=None):
    """Bubble pressure [MPa] at solvent mass fraction w1 and temperature T."""
    sp = solvent_params(solvent)
    cp = pseudomonomer_params(copolymer)
    x1 = w1_to_x1(w1, sp[3], cp[3])
    P = pc.bubble_pressure(x1, T, sp, cp, kij, P_hint=P_hint_Pa)
    return P / 1.0e6


def model_curve(T, solvent, copolymer, kij, w1_min, w1_max, npts=60,
                P_hint_Pa=None):
    """Smooth P(w1) model curve [MPa] over [w1_min, w1_max]."""
    w1s = np.linspace(w1_min, w1_max, npts)
    Ps = np.empty(npts)
    hint = P_hint_Pa
    for i, w1 in enumerate(w1s):
        Ps[i] = bubble_P_MPa(w1, T, solvent, copolymer, kij, hint)
        hint = Ps[i] * 1e6 * 1.5 + 1.0  # warm-start next point
    return w1s, Ps


def deviation_percent(w1_arr, Pexp_MPa, T, solvent, copolymer, kij):
    """%dP/P per Eq. 40 over the experimental points (same T)."""
    sp = solvent_params(solvent)
    cp = pseudomonomer_params(copolymer)
    tot = 0.0
    n = len(w1_arr)
    Pcalc = np.empty(n)
    for i in range(n):
        x1 = w1_to_x1(w1_arr[i], sp[3], cp[3])
        P = pc.bubble_pressure(x1, T, sp, cp, kij,
                               P_hint=Pexp_MPa[i] * 1e6)
        Pcalc[i] = P / 1.0e6
        tot += abs(Pcalc[i] - Pexp_MPa[i]) / Pexp_MPa[i]
    return 100.0 * tot / n, Pcalc
