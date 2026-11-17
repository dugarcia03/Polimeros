#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRUEBA NUEVA: anade el termino de asociacion (Ecs. 36-38 del paper) al nC4OH
usando SOLO los parametros del paper (Tabla 4b): kAB=0.008752, epsAB=2489.23 K.
Modelo 2B (dos sitios, esquema estandar para alcoholes).
Objetivo: ver si el Psat del 1-butanol puro baja de 565% a ~pocos %.
Si funciona, confirma que el termino faltante es LA causa para nC4OH.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import math, importlib.util
spec = importlib.util.spec_from_file_location('pv', 'pcsaft_vle.py')
pv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)
from scipy.optimize import brentq
import numpy as np

# nC4OH Tabla 4b
mMW,sig,eps_k,MW = 0.035126, 3.7563, 263.18, 74.12
kAB, epsAB = 0.008752, 2489.23   # Tabla 4b
m = mMW*MW
NAV = 6.02214076e23

class C: pass
c = C(); c.m=m; c.sigma=sig; c.eps_k=eps_k
ml = pv.PCSAFT([pv.Componente("nC4OH",m,sig,eps_k)],[[0.0]])

def d1(T): return sig*(1-0.12*math.exp(-3*eps_k/T))

def a_assoc(rho, T):
    # esquema 2B: un sitio A y un sitio B, asociacion A-B cruzada.
    # XA = XB por simetria 2B. Delta = sigma^3 * ghs(d) * kAB * (exp(epsAB/T)-1)
    z = ml.zeta([1.0], rho, T)
    ghs = ml.g_hs(0,0,z,T)
    d = d1(T)
    delta = (sig**3) * ghs * kAB * (math.exp(epsAB/T)-1.0)   # [A^3]
    # rho aqui en 1/A^3 (densidad numerica de moleculas). ec 37: rho_j molar...
    # usamos densidad numerica: XA = (-1+sqrt(1+4*rho*delta))/(2*rho*delta)
    rd = rho*delta
    XA = (-1.0 + math.sqrt(1.0+4.0*rd))/(2.0*rd)
    # 2 sitios iguales -> suma sobre sitios A y B: 2*(lnXA - XA/2 + 1/2)
    return 2.0*(math.log(XA) - XA/2.0 + 0.5), XA

# parche: a_res con asociacion
_orig_ares = ml.a_res
def a_res_assoc(x, rho, T):
    base = _orig_ares(x, rho, T)
    aa,_ = a_assoc(rho, T)
    return base + aa
ml.a_res = a_res_assoc

def psat(T):
    def eq(P):
        rl = pv._resolver_densidad(ml,[1.0],T,P,"liq")
        rv = pv._resolver_densidad(ml,[1.0],T,P,"vap")
        if rl is None or rv is None: return None
        return ml.ln_phi([1.0],rl,T)[0]-ml.ln_phi([1.0],rv,T)[0]
    grid=np.logspace(2,6.5,50); pP=pf=None
    for P in grid:
        f=eq(float(P))
        if f is None: pP=pf=None; continue
        if pf is not None and pf*f<0:
            try: return brentq(lambda p:eq(p),pP,float(P),xtol=1,maxiter=60)/1e6
            except: return None
        pP,pf=float(P),f
    return None

# Antoine 1-butanol
def ant(T): return 10**(4.55100-1351.555/(T-93.34))*0.1
for T in [353.15, 330.0, 370.0]:
    pm=psat(T); pa=ant(T)
    print(f"T={T}: Psat_mod={pm} MPa  Antoine={pa:.5f} MPa  err={abs(pm-pa)/pa*100:.1f}%" if pm else f"T={T}: sin raiz")
