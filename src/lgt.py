"""
LGT (Lattice Gas Theory) equation of state -- group-contribution form.

Implementacion de las Ecs. (1)-(12) del paper de Arce, Mattedi y Aznar (2006).

Modelo de contribucion de grupos para una mezcla binaria:
   componente 1 = solvente  (1 grupo: el propio solvente)
   componente 2 = copolimero (N grupos pseudo-monomero)

Constantes del paper:
   Zc    = 10.0   (numero de coordinacion)
   Psi   = 1.0    (constante universal)
   v*    = 5.0 cm3/mol  (volumen molar de celda)
   R     = 8.314 J/mol/K   (en las unidades del paper u0/R y B estan en K,
                            por lo que las energias entran como u/R [K])

Parametros de grupo (Tabla 4a / pseudo-monomero Tabla 5):
   v_inf [1e-6 m3/mol], Q, Baa [K], u0_R [K]
Reglas (Chen-Kreglewski, Eqs. 11-12):
   v_a   = v_inf * (1 - A*exp(-3 u0_R/ (R T)))^3     con A = 0  -> v_a = v_inf
   u_aa/R = (u0_R) * (1 + Baa/T)

El parametro de segmentos R_a (numero de segmentos del grupo) se toma igual a
   R_a = v_inf_a / v*     (coherente con v_a = R_a * v* a A=0),
de modo que r = sum x_i sum nu_i^a R_a reproduce rv* = sum x_i sum nu_i^a v_a.

Para grupos puros (solvente) R_a se define de la misma forma.
"""

import math

ZC = 10.0
PSI = 1.0
V_STAR = 5.0            # cm3/mol  -> usar v_inf en cm3/mol = 1e-6 m3/mol * 1e6 = valor tal cual
# Nota: v_inf esta en 1e-6 m3/mol = cm3/mol, asi que v* y v_inf quedan en las
# mismas unidades (cm3/mol).


def _group_params(group, params, T):
    """Devuelve (Ra, Qa, va, uaa_R) para un grupo a temperatura T.

    params: dict con keys v_inf, Q, Baa, u0_R  (A fijado a 0 en el paper).
    """
    v_inf = params["v_inf"]          # cm3/mol
    Q = params["Q"]
    Baa = params["Baa"]
    u0_R = params["u0_R"]
    # A = 0  =>  v_a = v_inf   (Eq. 11)
    va = v_inf
    Ra = va / V_STAR                 # numero de segmentos del grupo
    # Eq. 12: u_aa/R = u0_R * (1 + Baa/T)
    uaa_R = u0_R * (1.0 + Baa / T)
    return Ra, Qa_of(Q), va, uaa_R


def Qa_of(Q):
    return Q


class LGTMixture:
    """Mezcla binaria solvente(1) + copolimero(2) para la LGT EoS.

    groups: lista de dicts, uno por grupo distinto presente en la mezcla:
        {"name", "Ra", "Qa", "va", "uaa_R"}
    nu: matriz nu[i][a] = numero de grupos 'a' en una molecula del componente i.
    cross: matriz de energias de interaccion u_ma/R [K] entre grupos m,a.
    """

    def __init__(self, groups, nu, cross_uR):
        self.groups = groups
        self.ng = len(groups)
        self.nu = nu                 # nu[i][a]
        self.nc = len(nu)
        self.cross_uR = cross_uR     # cross_uR[m][a] en K
        self.Ra = [g["Ra"] for g in groups]
        self.Qa = [g["Qa"] for g in groups]
        self.va = [g["va"] for g in groups]

    # ---- cantidades dependientes de composicion (Eqs. 2-6) ----
    def r_i(self, i):
        return sum(self.nu[i][a] * self.Ra[a] for a in range(self.ng))

    def q_i(self, i):
        return sum(self.nu[i][a] * self.Qa[a] for a in range(self.ng))

    def r_mix(self, x):
        return sum(x[i] * self.r_i(i) for i in range(self.nc))

    def q_mix(self, x):
        return sum(x[i] * self.q_i(i) for i in range(self.nc))

    def rvstar(self, x):
        # Eq. 5/6 (con A=0, V_a = v_a): rv* = sum_i x_i sum_a nu_i^a v_a
        return sum(x[i] * sum(self.nu[i][a] * self.va[a] for a in range(self.ng))
                   for i in range(self.nc))

    # ---- factores de interaccion (Eqs. 7-9) ----
    def gamma(self, m, a, T):
        """gamma_ma = exp(-u_ma/(R T)).  cross_uR en K."""
        return math.exp(-self.cross_uR[m][a] / T)

    def S(self, x, m):
        """Eq. 8: S^m = (sum_i nu_i^m x_i Q^m)/q."""
        q = self.q_mix(x)
        num = sum(self.nu[i][m] * x[i] * self.Qa[m] for i in range(self.nc))
        return num / q

    def Gamma(self, x, a, T):
        """Eq. 7: Gamma^a = sum_m S^m gamma_ma."""
        return sum(self.S(x, m) * self.gamma(m, a, T) for m in range(self.ng))

    # ---- factor de compresibilidad (Eq. 1) ----
    def Z(self, x, v_red, T):
        """v_red = v~ (reduced volume)."""
        r = self.r_mix(x)
        q = self.q_mix(x)
        qr = q / r
        # l en la Eq.1: l = (Zc/2)*sum... termino de referencia; el paper define
        # l como el factor de Staverman   l_i = (Zc/2)(r_i - q_i) - (r_i - 1).
        # En la forma de la Ec.1 el termino 'l' acompana a la combinatoria;
        # usamos la definicion estandar de la red (l de mezcla):
        l = (ZC / 2.0) * (r - q) - (r - 1.0)

        term1 = v_red * r * math.log(v_red / (v_red - 1.0))
        term2 = (ZC / 2.0) * v_red * r * math.log((v_red - 1.0 + qr) / v_red)
        term3 = l
        # termino energetico (doble suma sobre grupos)
        s = 0.0
        for a in range(self.ng):
            Ga = self.Gamma(x, a, T)
            for i in range(self.nc):
                s += x[i] * self.nu[i][a] * self.Qa[a] * (Ga - 1.0) / (v_red - 1.0 + qr * Ga)
        term4 = - v_red * PSI * qr * s
        return term1 + term2 + term3 + term4

    def pressure(self, x, v_red, T):
        """P a partir de Z = P v / (R T),  v = v_red * rv* (cm3/mol).

        Devuelve P en MPa.
        R = 8.314 J/mol/K = 8.314 cm3*MPa/mol/K.
        """
        Rgas = 8.314            # cm3*MPa/(mol*K)
        rvstar = self.rvstar(x)       # cm3/mol
        v = v_red * rvstar            # cm3/mol
        Z = self.Z(x, v_red, T)
        return Z * Rgas * T / v

    def ln_phi(self, i, x, v_red, T):
        """Coeficiente de fugacidad del componente i (Eq. 10)."""
        r = self.r_mix(x)
        q = self.q_mix(x)
        qr = q / r
        ri = self.r_i(i)
        qi = self.q_i(i)
        li = (ZC / 2.0) * (ri - qi) - (ri - 1.0)

        t1 = -ri * math.log((v_red - 1.0) / (v_red - 1.0 + qr))
        t2 = (1.0 - li) * math.log(v_red / (v_red - 1.0 + qr))
        t3 = PSI * qr * (qi - ri) / (v_red - 1.0 + qr)

        t4 = 0.0
        for a in range(self.ng):
            Ga = self.Gamma(x, a, T)
            t4 += self.nu[i][a] * self.Qa[a] * math.log((v_red - 1.0 + qr) / (v_red - 1.0 + qr * Ga))
        t4 *= PSI

        t5 = 0.0
        for a in range(self.ng):
            Ga = self.Gamma(x, a, T)
            inner = 0.0
            for k in range(self.nc):
                num = sum(self.nu[i][e] * self.Qa[e] * self.gamma(e, a, T) for e in range(self.ng))
                inner += x[k] * self.nu[k][a] * self.Qa[a] * (num - ri) / (v_red - 1.0 + qr * Ga)
            t5 += inner
        t5 *= (PSI / r)

        Z = self.Z(x, v_red, T)
        return t1 + t2 + t3 + t4 - t5 - math.log(Z)
