#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PRUEBA A FONDO de la teoria: el m del copolimero por Ec.18 (fraccion MASICA)
difiere del (m/MW) de Tabla5 (fraccion MOLAR). Se prueban TODAS las variantes
de construir el pseudo-monomero y se compara %dP vs paper (Tabla 6, metod.1).

Variantes:
  M1_tabla5  : (m/MW,sigma,eps) directo de Tabla5 (lo que hace el codigo actual).
  M2_ec18    : m por Ec.18 (masa); sigma/eps por fraccion de SEGMENTO.
  M3_molar   : m, sigma, eps por fraccion MOLAR (deberia == Tabla5).
  M4_masa_all: m, sigma, eps TODOS por fraccion MASICA.
  M5_ec18_molarSE: m por Ec.18 (masa); sigma/eps por fraccion MOLAR.

Solo parametros del paper (Tabla 4b monomeros).
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import math, importlib.util
spec = importlib.util.spec_from_file_location('pv','pcsaft_vle.py'); pv=importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)
spec2 = importlib.util.spec_from_file_location('ct','correr_todos.py'); ct=importlib.util.module_from_spec(spec2); spec2.loader.exec_module(ct)

# Tabla 4b monomeros: (m/MW, sigma, eps/k, MW)
MON = {
 'STR':(0.024339,3.9750,318.70,104.14),'VA':(0.011525,5.0706,421.42,86.09),
 'BTD':(0.046033,3.4508,215.10,54.10),'ET':(0.056070,3.4250,178.77,28.10),
 'VC':(0.091211,2.4024,160.05,62.50),
}
# copolimero -> (mon1, mon2, wt%_mon1, MW_copol)
COPDEF = {
 'PSBR5':('STR','BTD',5,100000),'PSBR23':('STR','BTD',23,720000),'PSBR45':('STR','BTD',45,600000),
 'PVAE9':('VA','ET',9,100000),'PVAE25':('VA','ET',25,88000),'PVAE50':('VA','ET',50,250000),
 'PVAE70':('VA','ET',70,280000),'PVAVC12':('VA','VC',12,30000),
}

def frac(cop):
    A,B,wtA,MW = COPDEF[cop]; wA=wtA/100.0; wB=1-wA
    mA,mB = MON[A],MON[B]
    nA=wA/mA[3]; nB=wB/mB[3]; xA=nA/(nA+nB); xB=1-xA      # fraccion molar
    return A,B,wA,wB,xA,xB,MW,mA,mB

def build(cop, modo):
    A,B,wA,wB,xA,xB,MW,mA,mB = frac(cop)
    if modo=='M1_tabla5':
        c=ct.COP[cop]; return c[0]*c[3], c[1], c[2], float(c[3])
    if modo=='M3_molar':
        mMW=xA*mA[0]+xB*mB[0]; sig=xA*mA[1]+xB*mB[1]; eps=xA*mA[2]+xB*mB[2]
        return mMW*MW, sig, eps, MW
    if modo=='M4_masa_all':
        mMW=wA*mA[0]+wB*mB[0]; sig=wA*mA[1]+wB*mB[1]; eps=wA*mA[2]+wB*mB[2]
        return mMW*MW, sig, eps, MW
    # Ec.18: m por segmentos (masa)
    mA_seg=wA*MW*mA[0]; mB_seg=wB*MW*mB[0]; m=mA_seg+mB_seg
    zA=mA_seg/m; zB=mB_seg/m     # fraccion de SEGMENTO
    if modo=='M2_ec18':
        sig=zA*mA[1]+zB*mB[1]; eps=zA*mA[2]+zB*mB[2]
    elif modo=='M5_ec18_molarSE':
        sig=xA*mA[1]+xB*mB[1]; eps=xA*mA[2]+xB*mB[2]
    return m, sig, eps, MW

def dp(solv,cop,T,kij,modo):
    ms=ct.SOLV[solv]
    m,sig,eps,MW = build(cop,modo)
    pv.SOLVENTE={'nombre':solv,'MW':ms[3],'m_MW':ms[0],'sigma':ms[1],'eps_k':ms[2]}
    pv.COPOLIMERO={'MW':MW}
    if solv in ('nC4OH',):
        solvo=pv.Componente(solv,ms[0]*ms[3],ms[1],ms[2],0.008752,2489.23)
    else:
        solvo=pv.construir_solvente()
    copo=pv.Componente(cop,m,sig,eps)
    ml=pv.PCSAFT([solvo,copo],[[0,kij],[kij,0]]); mv=pv.PCSAFT([solvo],[[0.0]])
    s=0.0;n=0
    for w1,Pe in ct.DATA[(solv,cop,T)]:
        x1=pv.w1_a_x1(w1)
        try: Pc=pv.presion_de_burbuja(ml,mv,T,x1)
        except: Pc=None
        if Pc and Pc>0: s+=abs(Pc-Pe)/Pe;n+=1
    return (100*s/n if n else None)

modos=['M1_tabla5','M3_molar','M2_ec18','M5_ec18_molarSE','M4_masa_all']
print(f"{'sistema':>14} {'paper':>6} " + " ".join(f"{m:>12}" for m in modos))
tot={m:[0.0,0] for m in modos}
for (solv,cop,T),pts in ct.DATA.items():
    kij,dpp=ct.KIJ_DP[(solv,cop)]
    row=[]
    for m in modos:
        d=dp(solv,cop,T,kij,m); row.append(d)
        if d is not None: tot[m][0]+=d; tot[m][1]+=1
    print(f"{solv+'+'+cop:>14} {dpp:>6.2f} " + " ".join(f"{(d if d is not None else -1):>12.2f}" for d in row))
print("-"*90)
print(f"{'GLOBAL(prom)':>14} {4.41:>6.2f} " + " ".join(f"{tot[m][0]/tot[m][1]:>12.2f}" for m in modos))
