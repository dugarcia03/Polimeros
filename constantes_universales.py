"""
Universal PC-SAFT model constants.

SOURCE (the ONLY permitted external numerical source besides the paper):
    Gross, J.; Sadowski, G. "Perturbed-Chain SAFT: An Equation of State Based
    on a Perturbation Theory for Chain Molecules."
    Ind. Eng. Chem. Res. 2001, 40, 1244-1260. DOI: 10.1021/ie0003887
    -> Table 1 ("Universal Model Constants for Equations 18 and 19").

These are the a0i,a1i,a2i and b0i,b1i,b2i coefficients that define the
chain-length dependence of the dispersion power-series coefficients
ai(mbar) and bi(mbar) (Eqs. 18-19 / A.18-A.19 of Gross & Sadowski 2001,
equivalent to Eqs. 31-32 of Arce et al. 2006).

Index i runs 0..6. Values are transcribed EXACTLY as published.
"""

import numpy as np

# Columns: a0i, a1i, a2i   (Eq. 18 / A.18)
A_COEFFS = np.array([
    [ 0.9105631445, -0.3084016918, -0.0906148351],
    [ 0.6361281449,  0.1860531159,  0.4527842806],
    [ 2.6861347891, -2.5030047259,  0.5962700728],
    [-26.547362491, 21.419793629,  -1.7241829131],
    [ 97.759208784, -65.255885330, -4.1302112531],
    [-159.59154087,  83.318680481, 13.776631870],
    [ 91.297774084, -33.746922930, -8.6728470368],
])

# Columns: b0i, b1i, b2i   (Eq. 19 / A.19)
B_COEFFS = np.array([
    [ 0.7240946941, -0.5755498075,  0.0976883116],
    [ 2.2382791861,  0.6995095521, -0.2557574982],
    [-4.0025849485,  3.8925673390, -9.1558561530],
    [-21.003576815, -17.215471648, 20.642075974],
    [ 26.855641363, 192.67226447, -38.804430052],
    [206.55133841, -161.82646165, 93.626774077],
    [-355.60235612, -165.20769346, -29.666905585],
])

# Avogadro constant (standard physical constant; used in P = Z k T rho).
N_AV = 6.022e23  # 1/mol  (value used by Gross & Sadowski 2001, Appendix A)
# Boltzmann constant (J/K).
K_B = 1.380649e-23
