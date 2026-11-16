#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tabla comparativa METODOLOGIA 1 (pseudo-monomero) vs METODOLOGIA 2 (co-monomeros)
para cada sistema con datos experimentales, frente a lo reportado por el paper.

- Met 1: copolimero = 1 pseudo-monomero (Tabla 5) + 1 kij (Tabla 6).
- Met 2: copolimero = 2 monomeros alpha/beta (Tabla 4b) + B (Tabla 3) + 3 kij (Tabla 7).
Usa el motor de pcsaft_vle.py (metodologia 2) y, para met 1, un motor de 1 segmento.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import math, importlib.util
import numpy as np
from scipy.optimize import brentq
from constantes_pcsaft import A_COEF, B_COEF

KB = 1.380649e-23

# ---- motor PC-SAFT generico (1 o 2 segmentos), copia minima ----
spec = importlib.util.spec_from_file_location('pv', 'pcsaft_vle.py')
pv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)
Segmento, Componente, PCSAFT = pv.Segmento, pv.Componente, pv.PCSAFT

# ====== DATOS DE TABLAS DEL PAPER ======
# Solventes (Tabla 4b): m_MW, sigma, eps_k, MW
SOLV = {
 'ACE':(0.037879,3.4956,272.26,58.08),'nC4OH':(0.035126,3.7563,263.18,74.12),
 'CHL':(0.021826,3.4085,265.94,119.38),'cC6':(0.029946,3.8399,279.21,84.16),
 'nC6':(0.034088,3.8383,239.57,86.18),'nC5':(0.034090,3.8729,240.11,72.15),
 'TOL':(0.030984,3.6869,282.69,92.14),'BZN':(0.032712,3.5878,281.55,78.12)}
# Pseudo-monomeros (Tabla 5): m_MW, sigma, eps_k  (para met 1)
PSM = {
 'PSBR4.1':(0.045562,3.4622,217.35),'PSBR5':(0.045456,3.4647,217.86),
 'PSBR23':(0.043119,3.5212,229.01),'PSBR41':(0.040280,3.5898,242.57),'PSBR45':(0.039562,3.6072,246.00),
 'PVAE9':(0.054676,3.4765,186.36),'PVAE25':(0.051700,3.5864,202.57),
 'PVAE50':(0.045107,3.8300,238.49),'PVAE70':(0.036813,4.1364,283.67),'PVAVC12':(0.084039,2.6425,183.57)}
# Monomeros (Tabla 4b): m_MW, sigma, eps_k, MW  (para met 2)
MON = {
 'STR':(0.024339,3.9750,318.70,104.14),'BTD':(0.046033,3.4508,215.10,54.10),
 'VA':(0.011525,5.0706,421.42,86.09),'ET':(0.056070,3.4250,178.77,28.10),
 'VC':(0.091211,2.4024,160.05,62.50)}
# copolimero -> (monomero_alpha, monomero_beta, w_alpha, MW_copol)
COPDEF = {
 'PSBR4.1':('STR','BTD',0.041,100000),'PSBR5':('STR','BTD',0.05,100000),
 'PSBR23':('STR','BTD',0.23,720000),'PSBR41':('STR','BTD',0.41,100000),'PSBR45':('STR','BTD',0.45,600000),
 'PVAE9':('VA','ET',0.09,100000),'PVAE25':('VA','ET',0.25,88000),
 'PVAE50':('VA','ET',0.50,250000),'PVAE70':('VA','ET',0.70,280000),'PVAVC12':('VA','VC',0.12,30000)}
# Fracciones de enlace B (Tabla 3): aa, ab, bb ; None si no estan en Tabla 3
BOND = {
 'PSBR23':(0.7161,0.2839,0.0),'PSBR45':(0.3707,0.6293,0.0),
 'PVAE9':(0.8392,0.1608,0.0),'PVAE25':(0.5423,0.4577,0.0),
 'PVAE50':(0.0487,0.9513,0.0),'PVAE70':(0.0,0.4633,0.5367)}
# kij Tabla 6 (met1) y dP_pcsaft reportado
T6 = {
 ('ACE','PSBR5'):(0.0120,4.78),('ACE','PSBR23'):(0.0122,3.04),('ACE','PSBR45'):(0.0130,5.19),
 ('nC4OH','PVAE9'):(0.0230,4.05),('nC4OH','PVAVC12'):(-0.0126,5.18),
 ('CHL','PVAE9'):(-0.1208,4.18),('CHL','PVAE25'):(-0.0816,5.06),('CHL','PVAE50'):(-0.0701,4.41),('CHL','PVAE70'):(-0.0516,3.28),
 ('cC6','PSBR5'):(0.0118,3.66),('cC6','PSBR23'):(0.0104,4.38),('cC6','PSBR45'):(0.0101,4.31),
 ('cC6','PVAE25'):(-0.1068,5.36),('cC6','PVAE50'):(-0.0858,6.90),('cC6','PVAE70'):(-0.0636,5.73),
 ('TOL','PSBR5'):(0.0080,4.07),('TOL','PSBR23'):(0.0091,5.76),
 ('nC6','PSBR23'):(0.0086,3.92),('nC6','PSBR45'):(0.0089,4.28),
 ('nC5','PSBR23'):(0.0072,3.68),('nC5','PSBR45'):(0.0079,6.24)}
# kij Tabla 7 (met2): (solvente,familia)->(kij_a, kij_b) ; kiajb por copolimero
T7 = {
 ('ACE','PSBR'):(0.0121,0.0028),('nC4OH','PVAE'):(0.0160,0.0131),('nC4OH','PVAVC'):(-0.0101,-0.0090),
 ('CHL','PVAE'):(-0.0532,0.0003),('cC6','PSBR'):(0.0052,0.0012),('cC6','PVAE'):(-0.1085,0.0101),
 ('TOL','PSBR'):(0.0082,0.0126),('nC6','PSBR'):(0.0080,0.0091),('nC5','PSBR'):(0.0066,0.0108)}
KIAB = {'PSBR5':0.0003,'PSBR23':0.0004,'PSBR45':0.0001,'PVAE9':-0.0010,'PVAVC12':0.0008,
 'PVAE25':{'CHL':0.0005,'cC6':0.0011},'PVAE50':{'CHL':0.0003,'cC6':0.0010},
 'PVAE70':{'CHL':0.0003,'cC6':0.0008},'PSBR41':0.0006}
# dP met2 reportado (Tabla 7, por familia)
DP7 = {('ACE','PSBR'):2.83,('nC4OH','PVAE'):1.88,('nC4OH','PVAVC'):1.63,('CHL','PVAE'):2.59,
 ('cC6','PSBR'):2.11,('cC6','PVAE'):2.09,('TOL','PSBR'):2.13,('nC6','PSBR'):1.84,('nC5','PSBR'):2.07}
FAM = {'PSBR4.1':'PSBR','PSBR5':'PSBR','PSBR23':'PSBR','PSBR41':'PSBR','PSBR45':'PSBR',
 'PVAE9':'PVAE','PVAE25':'PVAE','PVAE50':'PVAE','PVAE70':'PVAE','PVAVC12':'PVAVC'}

# datos experimentales (w1, P[MPa]) -- de Gupta [1]
DATA = {
 ('ACE','PSBR5',323.15):[(0.003,0.0137),(0.023,0.0265),(0.040,0.0363),(0.056,0.0415),(0.067,0.0453),(0.078,0.0511),(0.080,0.0533),(0.091,0.0552),(0.115,0.0596),(0.120,0.0613)],
 ('ACE','PSBR23',333.15):[(0.009,0.0127),(0.021,0.0264),(0.037,0.0407),(0.062,0.0548),(0.092,0.0685),(0.128,0.0813),(0.147,0.0848),(0.183,0.0915),(0.231,0.0979)],
 ('ACE','PSBR45',333.15):[(0.009,0.0127),(0.023,0.0264),(0.044,0.0407),(0.069,0.0548),(0.101,0.0685),(0.140,0.0813),(0.163,0.0848),(0.197,0.0915),(0.258,0.0979)],
 ('nC4OH','PVAE9',353.15):[(0.006,0.0027),(0.012,0.0052),(0.017,0.0080),(0.027,0.0108),(0.031,0.0123),(0.034,0.0139),(0.036,0.0152)],
 ('nC4OH','PVAVC12',353.15):[(0.007,0.0027),(0.016,0.0052),(0.022,0.0080),(0.033,0.0108),(0.040,0.0123),(0.044,0.0139),(0.046,0.0152)],
 ('CHL','PVAE9',333.15):[(0.031,0.0140),(0.067,0.0281),(0.104,0.0401),(0.147,0.0532),(0.168,0.0584),(0.196,0.0661),(0.216,0.0723),(0.235,0.0792),(0.253,0.0871)],
 ('CHL','PVAE25',333.15):[(0.062,0.0139),(0.140,0.0267),(0.228,0.0399),(0.320,0.0532),(0.358,0.0599),(0.401,0.0665),(0.441,0.0731),(0.481,0.0793)],
 ('CHL','PVAE50',333.15):[(0.098,0.0139),(0.202,0.0267),(0.303,0.0399),(0.398,0.0532),(0.439,0.0599),(0.477,0.0665),(0.512,0.0731),(0.552,0.0793)],
 ('CHL','PVAE70',333.15):[(0.059,0.0139),(0.161,0.0267),(0.287,0.0399),(0.403,0.0532),(0.450,0.0599),(0.489,0.0665),(0.526,0.0731),(0.573,0.0793)],
 ('cC6','PSBR5',333.15):[(0.027,0.0048),(0.057,0.0105),(0.103,0.0169),(0.120,0.0199),(0.152,0.0248),(0.189,0.0295),(0.233,0.0337),(0.276,0.0376),(0.323,0.0407)],
 ('cC6','PSBR23',333.15):[(0.018,0.0048),(0.052,0.0105),(0.088,0.0169),(0.136,0.0248),(0.168,0.0295),(0.208,0.0337),(0.244,0.0376),(0.287,0.0407)],
 ('cC6','PSBR45',333.15):[(0.021,0.0048),(0.050,0.0105),(0.087,0.0169),(0.103,0.0199),(0.130,0.0248),(0.163,0.0295),(0.200,0.0337),(0.240,0.0376),(0.282,0.0407)],
 ('cC6','PVAE25',353.15):[(0.011,0.0075),(0.037,0.0144),(0.074,0.0265),(0.140,0.0397),(0.228,0.0549),(0.238,0.0555),(0.277,0.0607),(0.300,0.0641),(0.369,0.0704)],
 ('cC6','PVAE50',353.15):[(0.005,0.0075),(0.025,0.0144),(0.046,0.0265),(0.104,0.0397),(0.159,0.0549),(0.170,0.0555),(0.205,0.0607),(0.225,0.0641),(0.287,0.0704)],
 ('cC6','PVAE70',353.15):[(0.004,0.0075),(0.026,0.0144),(0.042,0.0265),(0.081,0.0397),(0.123,0.0549),(0.133,0.0555),(0.156,0.0607),(0.191,0.0641),(0.242,0.0704)],
 ('TOL','PSBR5',308.15):[(0.043,0.00093),(0.099,0.00208),(0.184,0.00304),(0.319,0.00436),(0.416,0.00483),(0.513,0.00548)],
 ('TOL','PSBR23',308.15):[(0.029,0.00090),(0.031,0.00208),(0.032,0.00304),(0.037,0.00436),(0.037,0.00483),(0.038,0.00548)],
 ('nC6','PSBR23',343.15):[(0.014,0.0139),(0.043,0.0273),(0.063,0.0385),(0.110,0.0488),(0.162,0.0615),(0.262,0.0752)],
 ('nC6','PSBR45',343.15):[(0.016,0.0139),(0.038,0.0273),(0.061,0.0385),(0.103,0.0488),(0.155,0.0615),(0.249,0.0752)],
 ('nC5','PSBR23',333.15):[(0.008,0.0228),(0.016,0.0360),(0.021,0.0468),(0.033,0.0655),(0.046,0.0840),(0.055,0.0973),(0.070,0.1133),(0.078,0.1248),(0.089,0.1376)],
 ('nC5','PSBR45',333.15):[(0.010,0.0228),(0.013,0.0468),(0.019,0.0655),(0.022,0.0840),(0.022,0.0973),(0.026,0.1133),(0.027,0.1248),(0.034,0.1376)],
}


def w1_to_x1(w1, MWs, MWc):
    ns=w1/MWs; nc=(1-w1)/MWc
    return ns/(ns+nc)


def dp_met1(solv_k, cop_k, T, pts):
    ms = SOLV[solv_k]; ps = PSM[cop_k]; MWc = COPDEF[cop_k][3]
    kij,_ = T6[(solv_k,cop_k)]
    segS = Segmento('S', ms[0]*ms[3], ms[1], ms[2])
    cS = Componente('S',[segS],[1.0])
    segP = Segmento('P', ps[0]*MWc, ps[1], ps[2])
    cP = Componente('P',[segP],[1.0])
    kd = {frozenset(('S','P')):kij}
    ml = PCSAFT([cS,cP],kd); mv = PCSAFT([cS],{})
    return _dp(ml,mv,T,pts,ms[3],MWc)


def dp_met2(solv_k, cop_k, T, pts):
    ms = SOLV[solv_k]; mA_k,mB_k,wA,MWc = COPDEF[cop_k]
    fam = FAM[cop_k]
    if (solv_k,fam) not in T7: return None,0
    ka,kb = T7[(solv_k,fam)]
    kiab_raw = KIAB.get(cop_k,0.0)
    kiab = kiab_raw[solv_k] if isinstance(kiab_raw,dict) else kiab_raw
    mA = MON[mA_k]; mB = MON[mB_k]
    segS = Segmento('S', ms[0]*ms[3], ms[1], ms[2]); cS = Componente('S',[segS],[1.0])
    maA = wA*MWc*mA[0]; maB = (1-wA)*MWc*mB[0]; mt=maA+maB
    zA=maA/mt; zB=maB/mt
    segA=Segmento(mA_k,maA,mA[1],mA[2]); segB=Segmento(mB_k,maB,mB[1],mB[2])
    B = BOND.get(cop_k)
    bonding = {'aa':B[0],'ab':B[1],'bb':B[2]} if B else pv._fracciones_enlace(zA,zB,mt)
    cP = Componente('P',[segA,segB],[zA,zB],bonding)
    kd = {frozenset(('S',mA_k)):ka, frozenset(('S',mB_k)):kb, frozenset((mA_k,mB_k)):kiab}
    ml = PCSAFT([cS,cP],kd); mv = PCSAFT([cS],{})
    return _dp(ml,mv,T,pts,ms[3],MWc)


def _dp(ml,mv,T,pts,MWs,MWc):
    s=0.0; n=0
    for w1,Pe in pts:
        x1=w1_to_x1(w1,MWs,MWc)
        try:
            Pc = pv.presion_de_burbuja.__wrapped__ if False else bubble(ml,mv,T,x1)
        except Exception:
            Pc=None
        if Pc and Pc>0:
            s+=abs(Pc-Pe)/Pe; n+=1
    return (100*s/n if n else None), n


def bubble(ml,mv,T,x1):
    x=[x1,1-x1]
    def res(P):
        rl=pv._resolver_densidad(ml,x,T,P,'liq')
        if rl is None: return None
        ll=ml.ln_phi(x,rl,T)[0]
        rv=pv._resolver_densidad(mv,[1.0],T,P,'vap')
        if rv is None: return None
        lv=mv.ln_phi([1.0],rv,T)[0]
        return (math.log(x1)+ll)-lv
    grid=np.logspace(2,math.log10(2e6),28); pP=pf=None; br=None
    for P in grid:
        f=res(float(P))
        if f is None: pP=pf=None; continue
        if pf is not None and pf*f<0: br=(pP,float(P)); break
        pP,pf=float(P),f
    if br is None: return None
    try: return brentq(lambda P:res(P),br[0],br[1],xtol=1e-2,rtol=1e-9,maxiter=60)/1e6
    except Exception: return 0.5*(br[0]+br[1])/1e6


if __name__ == "__main__":
    print("%-6s %-9s %7s | %8s %8s | %8s %8s"%(
        "Solv","Copol","T[K]","M1_calc","M1_pap","M2_calc","M2_pap"))
    print("-"*62)
    for (solv,cop,T),pts in DATA.items():
        fam=FAM[cop]
        _,dpm1_pap = T6[(solv,cop)]
        dpm2_pap = DP7.get((solv,fam),None)
        d1,_ = dp_met1(solv,cop,T,pts)
        d2,_ = dp_met2(solv,cop,T,pts)
        f1 = "%.1f"%d1 if d1 else "-"
        f2 = "%.1f"%d2 if d2 else "-"
        f2p = "%.2f"%dpm2_pap if dpm2_pap else "-"
        print("%-6s %-9s %7.2f | %8s %8.2f | %8s %8s"%(
            solv,cop,T,f1,dpm1_pap,f2,f2p))
    print("-"*62)
    print("Global paper:  M1 PC-SAFT = 4.41%   |   M2 PC-SAFT = 2.04%")
