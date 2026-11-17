#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRUEBA DE ORO: valida nuestro motor PC-SAFT contra los parametros ORIGINALES
de Gross-Sadowski (2001, Tabla 2) para n-alcanos puros. GS reporta AAD(Psat)<1%.
Si nuestro motor reproduce la Psat experimental con error pequeno, el motor esta
100% verificado (independiente del paper de Arce).
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import math, importlib.util
import numpy as np
from scipy.optimize import brentq
spec = importlib.util.spec_from_file_location('pv','pcsaft_vle.py'); pv=importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)

# Parametros GS 2001 (m, sigma, eps/k) -- m es numero de segmentos directo (no m/MW)
GS = {
 'n-pentane':(2.6896,3.7729,231.20,72.146),
 'n-hexane' :(3.0576,3.7983,236.77,86.177),
 'n-heptane':(3.4831,3.8049,238.40,100.203),
 'n-octane' :(3.8176,3.8373,242.78,114.231),
}
# Antoine (bar), NIST
ANT = {
 'n-pentane':(3.97786,1064.840,-41.136),
 'n-hexane' :(4.00266,1171.530,-48.784),
 'n-heptane':(4.02832,1268.636,-56.199),
 'n-octane' :(4.04867,1355.126,-63.633),
}
def antoine(s,T): A,B,C=ANT[s]; return 10**(A-B/(T+C))*0.1  # MPa

def psat(m,sig,eps,T):
    comp = pv.Componente('x', m, sig, eps)
    ml = pv.PCSAFT([comp],[[0.0]])
    def eq(P):
        rl=pv._resolver_densidad(ml,[1.0],T,P,'liq'); rv=pv._resolver_densidad(ml,[1.0],T,P,'vap')
        if rl is None or rv is None: return None
        return ml.ln_phi([1.0],rl,T)[0]-ml.ln_phi([1.0],rv,T)[0]
    g=np.logspace(2,6.5,60); pP=pf=None
    for P in g:
        f=eq(float(P))
        if f is None: pP=pf=None; continue
        if pf is not None and pf*f<0:
            try: return brentq(lambda p:eq(p),pP,float(P),xtol=1,maxiter=80)/1e6
            except: return None
        pP,pf=float(P),f
    return None

print(f"{'substancia':>11} {'T[K]':>7} {'Psat_motor':>11} {'Psat_exp':>10} {'%err':>7}")
for s,(m,sig,eps,MW) in GS.items():
    for T in [300.0, 330.0, 360.0]:
        pm=psat(m,sig,eps,T); pa=antoine(s,T)
        e=abs(pm-pa)/pa*100 if pm else float('nan')
        pms=f"{pm:.5f}" if pm else "None"
        print(f"{s:>11} {T:>7.1f} {pms:>11} {pa:>10.5f} {e:>7.2f}")
