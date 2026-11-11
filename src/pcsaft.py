"""
PC-SAFT equation of state extended to copolymers.

Implementacion de las Ecs. (18)-(38) del paper de Arce, Mattedi y Aznar (2006),
que a su vez sigue a Gross & Sadowski (2001) y Gross et al. (2003).

La mezcla binaria es:  solvente (componente 1) + copolimero (componente 2).
El copolimero contiene dos tipos de segmento: alpha (monomero 1) y beta (monomero 2).
El solvente es un componente de un solo tipo de segmento.

Convenciones internas
---------------------
Trabajamos con una lista de "especies de segmento". Cada componente i aporta
sus segmentos con una fraccion de segmento z_{i,s}.  Para el solvente hay un
unico segmento con z = 1.  Para el copolimero hay dos segmentos (alpha, beta)
con z_alpha + z_beta = 1.

Parametros por segmento:  m (num. de segmentos), sigma [A], eps/k [K].
Para el copolimero, m_{i,alpha} = w_{i,alpha} * MW_copolym * (m/MW)_alpha  (Eq.18)
y m_i = m_{i,alpha} + m_{i,beta}  (Eq.19),  z_{i,s} = m_{i,s}/m_i  (Eq.20).

Las constantes universales a0k..a2k, b0k..b2k son las de Gross & Sadowski (2001).
"""

import math

# constantes universales de Gross & Sadowski (2001), Tablas 1 y 2
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

KB = 1.380649e-23        # J/K
NA = 6.02214076e23       # 1/mol


class Segment:
    def __init__(self, name, m, sigma, eps_k, kAB=None, epsAB_k=None):
        self.name = name
        self.m = m            # numero de segmentos aportado por este segmento
        self.sigma = sigma    # Angstrom
        self.eps_k = eps_k    # K


class Component:
    """Componente PC-SAFT. Puede tener 1 (solvente) o 2 (copolimero) tipos
    de segmento.  m_total es la suma de m de cada tipo."""

    def __init__(self, name, segs, z):
        self.name = name
        self.segs = segs           # lista de Segment
        self.z = z                 # fracciones de segmento por tipo (suman 1)
        self.m = sum(s.m for s in segs)


class PCSAFT:
    """Mezcla PC-SAFT. kij: diccionario de parametros de interaccion binaria
    por par de tipos de segmento (clave = frozenset de nombres o tupla)."""

    def __init__(self, components, kij=None):
        self.comps = components
        self.kij = kij or {}
        # catalogo plano de segmentos (por componente)
        # usamos indices (i, s) para recorrer

    def _kij(self, name_a, name_b):
        if name_a == name_b:
            return 0.0
        key = frozenset((name_a, name_b))
        return self.kij.get(key, 0.0)

    def _sigma_cross(self, sa, sb):
        return 0.5 * (sa.sigma + sb.sigma)

    def _eps_cross(self, sa, sb):
        k = self._kij(sa.name, sb.name)
        return math.sqrt(sa.eps_k * sb.eps_k) * (1.0 - k)

    def d_seg(self, seg, T):
        """Diametro de segmento dependiente de T (Eq. 26)."""
        return seg.sigma * (1.0 - 0.12 * math.exp(-3.0 * seg.eps_k / T))

    def m_bar(self, x):
        """Eq. 24: m promedio."""
        return sum(x[i] * self.comps[i].m for i in range(len(self.comps)))

    def _zeta(self, x, rho, T):
        """Eq. 25: zeta_k, k=0..3.  rho en 1/A^3 (densidad numerica de moleculas)."""
        zeta = [0.0, 0.0, 0.0, 0.0]
        for k in range(4):
            s = 0.0
            for i, c in enumerate(self.comps):
                for si, seg in enumerate(c.segs):
                    d = self.d_seg(seg, T)
                    s += x[i] * c.m * c.z[si] * (d ** k)
            zeta[k] = math.pi / 6.0 * rho * s
        return zeta

    def eta_of(self, x, rho, T):
        return self._zeta(x, rho, T)[3]

    # ---------------- contribucion hard-chain (Eqs. 22-23) ----------------
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
        """Eq. 22. bonding: dict por componente copolimero con fracciones de
        enlace B_{alpha,beta}; para solventes se ignora (m-1 cadena normal).

        Para un componente de 1 segmento (solvente) el termino de cadena es
        -(m-1) ln g^hs(d).  Para el copolimero se usa la suma ponderada por
        las fracciones de enlace B (Eq. 22)."""
        zeta = self._zeta(x, rho, T)
        mbar = self.m_bar(x)
        a = mbar * self._a_hs(zeta)
        for i, c in enumerate(self.comps):
            if len(c.segs) == 1:
                d = self.d_seg(c.segs[0], T)
                g = self._ghs(d, d, zeta)
                a -= x[i] * (c.m - 1.0) * math.log(g)
            else:
                # copolimero: suma sobre pares de segmentos con B
                B = bonding.get(i, None)
                da = self.d_seg(c.segs[0], T)
                db = self.d_seg(c.segs[1], T)
                g_aa = self._ghs(da, da, zeta)
                g_ab = self._ghs(da, db, zeta)
                g_bb = self._ghs(db, db, zeta)
                if B is None:
                    # sin datos de enlace -> tratar como homopolimero promedio
                    g = self._ghs(da, db, zeta)
                    a -= x[i] * (c.m - 1.0) * math.log(g)
                else:
                    term = (B["aa"] * math.log(g_aa)
                            + B["ab"] * math.log(g_ab)
                            + B["bb"] * math.log(g_bb))
                    a -= x[i] * (c.m - 1.0) * term
        return a

    # ---------------- contribucion de perturbacion (Eqs. 27-32) ----------------
    def _ab_k(self, mbar):
        a = []
        b = []
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
        """sum sum xi xj mi mj sum_a sum_b zia zjb (eps_ab/kT)^power sigma_ab^3."""
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
        """Eqs. 28-29 (dispersion)."""
        mbar = self.m_bar(x)
        eta = self.eta_of(x, rho, T)
        I1, I2 = self._I1_I2(mbar, eta)
        C1 = self._C1(mbar, eta)
        m2es3_1 = self._m2es3(x, T, 1)
        m2es3_2 = self._m2es3(x, T, 2)
        a1 = -2.0 * math.pi * rho * I1 * m2es3_1
        a2 = -math.pi * rho * mbar * C1 * I2 * m2es3_2
        return a1 + a2

    # ---------------- residual y presion ----------------
    def a_res(self, x, rho, T, bonding):
        return self.a_hc(x, rho, T, bonding) + self.a_disp(x, rho, T)

    def Z(self, x, rho, T, bonding):
        """Factor de compresibilidad via Z = 1 + eta d(a_res)/d(eta) (derivada numerica)."""
        eta = self.eta_of(x, rho, T)
        deta = 1e-6 * max(eta, 1e-6)
        # rho y eta proporcionales -> d a/d eta = d a/d rho * d rho/d eta
        # usamos derivada respecto a rho directamente: Z = 1 + rho d a/d rho
        drho = 1e-6 * rho
        a_plus = self.a_res(x, rho + drho, T, bonding)
        a_minus = self.a_res(x, rho - drho, T, bonding)
        da_drho = (a_plus - a_minus) / (2.0 * drho)
        return 1.0 + rho * da_drho

    def pressure(self, x, rho, T, bonding):
        """Presion en Pa.  rho en 1/A^3.
        P = Z rho kB T  con rho en 1/m^3.
        """
        Z = self.Z(x, rho, T, bonding)
        rho_m3 = rho * 1e30          # 1/A^3 -> 1/m^3
        return Z * rho_m3 * KB * T   # Pa

    def ln_phi(self, x, rho, T, bonding):
        """Coeficientes de fugacidad de todos los componentes (Eqs. 34-35).

        Devuelve lista ln_phi_i.  Derivadas de a_res respecto a x_i por
        diferencias finitas (manteniendo T, V=cte via rho).
        """
        nc = len(self.comps)
        ares = self.a_res(x, rho, T, bonding)
        Z = self.Z(x, rho, T, bonding)
        lnZ = math.log(Z)

        # d a_res / d x_i a T, V (rho) constantes
        dadx = []
        h = 1e-7
        for i in range(nc):
            xp = list(x)
            xp[i] += h
            ap = self.a_res(xp, rho, T, bonding)
            xm = list(x)
            xm[i] -= h
            am = self.a_res(xm, rho, T, bonding)
            dadx.append((ap - am) / (2.0 * h))

        sum_xj_dadxj = sum(x[j] * dadx[j] for j in range(nc))
        lnphi = []
        for i in range(nc):
            mu_res = ares + (Z - 1.0) + dadx[i] - sum_xj_dadxj
            lnphi.append(mu_res - lnZ)
        return lnphi
