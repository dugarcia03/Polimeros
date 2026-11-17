#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
METODOLOGIA 2: corre PC-SAFT con copolimero de 2 segmentos (alpha/beta) para los
sistemas donde el paper da fracciones de enlace (Tabla 3) y tres kij (Tabla 7).
Compara %dP calculado vs reportado (Tabla 7, global paper = 2.04%).

Todos los parametros del paper:
  - Tabla 4b: parametros PC-SAFT de solventes y MONOMEROS.
  - Tabla 3 : fracciones de enlace B_{aa}, B_{ab}, B_{bb} de cada copolimero.
  - Tabla 7 : kappa_a (solv-mon1), kappa_b (solv-mon2), kappa_ab (mon1-mon2).
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import importlib.util
spec = importlib.util.spec_from_file_location('m2', 'pcsaft_m2.py')
m2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m2)
spec2 = importlib.util.spec_from_file_location('ct', 'correr_todos.py')
ct = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(ct)

# Tabla 4b: solventes (m/MW, sigma, eps_k, MW [, kAB, epsAB])
SOLV = dict(ct.SOLV)
ASOC = {"nC4OH": (0.008752, 2489.23)}

# Tabla 4b: MONOMEROS (m/MW, sigma, eps_k, MW)
MON = {
 'STR':(0.024339,3.9750,318.70,104.14),'VA':(0.011525,5.0706,421.42,86.09),
 'BTD':(0.046033,3.4508,215.10,54.10),'ET':(0.056070,3.4250,178.77,28.10),
 'VC':(0.091211,2.4024,160.05,62.50),
}

# Copolimero -> (mon1, mon2, wt%_mon1, MW_copol, (B_aa_t3, B_ab_t3, B_bb_t3))
# mon1 = monomero del nombre (STR en PSBR, VA en PVAE).  mon2 = el otro.
#
# MAPEO Tabla 3 -> segmentos:
#   La Tabla 3 lista: ziβ | B_{iα,iα} | B_{iα,jβ} | B_{jβ,jβ}
#   donde β = mon1 (el del wt%, minoritario en PVAE9) y α = mon2.
#   Verificado via ziβ: PVAE9 ziβ=0.0802 ~ fraccion del segmento VA(mon1).
#   En el motor: segmento indice 0 = mon1 (=β_tabla3), indice 1 = mon2 (=α_tabla3).
#   Por tanto la matriz B del motor es:
#       B[0][0] = B_{ββ} = B_bb_t3   (mon1-mon1)
#       B[1][1] = B_{αα} = B_aa_t3   (mon2-mon2)
#       B[0][1] = B[1][0] = B_{αβ} = B_ab_t3
# Tupla guardada = (B_aa_t3, B_ab_t3, B_bb_t3) tal como aparece en Tabla 3.
COP = {
 'PSBR23':  ('STR','BTD',23,720000,(0.7161,0.2839,0.0000)),
 'PSBR45':  ('STR','BTD',45,600000,(0.3707,0.6293,0.0000)),
 'PVAE9':   ('VA','ET', 9,100000,(0.8392,0.1608,0.0000)),
 'PVAE25':  ('VA','ET',25, 88000,(0.5423,0.4577,0.0000)),
 'PVAE50':  ('VA','ET',50,250000,(0.0487,0.9513,0.0000)),
 'PVAE70':  ('VA','ET',70,280000,(0.0000,0.4633,0.5367)),
}

# Tabla 7: kappa_a (solv-mon1), kappa_b (solv-mon2) por FAMILIA (solv, tipo).
KAB_FAM = {
 ('ACE','PSBR'): (0.0121, 0.0028),
 ('nC4OH','PVAE'):(0.0160, 0.0131),
 ('CHL','PVAE'): (-0.0532, 0.0003),
 ('cC6','PSBR'): (0.0052, 0.0012),
 ('cC6','PVAE'): (-0.1085, 0.0101),
 ('TOL','PSBR'): (0.0082, 0.0126),
 ('nC6','PSBR'): (0.0080, 0.0091),
 ('nC5','PSBR'): (0.0066, 0.0108),
}
# Tabla 7: kappa_ab (mon1-mon2) por SISTEMA individual.
KAB_AB = {
 ('ACE','PSBR23'):0.0004,('ACE','PSBR45'):0.0001,
 ('nC4OH','PVAE9'):-0.0010,
 ('CHL','PVAE9'):0.0008,('CHL','PVAE25'):0.0005,('CHL','PVAE50'):0.0003,('CHL','PVAE70'):0.0003,
 ('cC6','PSBR23'):0.0006,('cC6','PSBR45'):0.0002,
 ('cC6','PVAE25'):0.0011,('cC6','PVAE50'):0.0010,('cC6','PVAE70'):0.0008,
 ('TOL','PSBR23'):0.0011,
 ('nC6','PSBR23'):0.0010,('nC6','PSBR45'):0.0008,
 ('nC5','PSBR23'):0.0013,('nC5','PSBR45'):0.0008,
}
# %dP reportado por el paper (Tabla 7) - es por FAMILIA (en fila del 1er copol)
DP_FAM = {('ACE','PSBR'):2.83,('nC4OH','PVAE'):1.88,('CHL','PVAE'):2.59,
          ('cC6','PSBR'):2.11,('cC6','PVAE'):2.09,('TOL','PSBR'):2.13,
          ('nC6','PSBR'):1.84,('nC5','PSBR'):2.07}

def familia(cop):
    if cop.startswith('PSBR'): return 'PSBR'
    if cop.startswith('PVAE'): return 'PVAE'
    if cop.startswith('PVAVC'): return 'PVAVC'

def construir_solvente(solv):
    ms = SOLV[solv]
    kAB,epsAB = ASOC.get(solv,(0.0,0.0))
    return m2.Componente(solv, [m2.Segmento(solv, ms[0]*ms[3], ms[1], ms[2], kAB, epsAB)])

def construir_copolimero(cop):
    A,B,wtA,MW,Bfr = COP[cop]
    wA = wtA/100.0; wB = 1-wA
    mA = MON[A]; mB = MON[B]
    # Ec.18: m_a = w_a * MW * (m/MW)_a
    segA = m2.Segmento(A, wA*MW*mA[0], mA[1], mA[2])
    segB = m2.Segmento(B, wB*MW*mB[0], mB[1], mB[2])
    # Bfr = (B_aa_t3, B_ab_t3, B_bb_t3) de Tabla 3, con α=mon2, β=mon1.
    # Los tres valores de Tabla 3 SUMAN 1 (son las fracciones de los 3 tipos de enlace).
    # En la doble suma Σ_a Σ_b B[a][b] ln g, el enlace cruzado αβ aparece DOS veces
    # (posiciones [0][1] y [1][0]); por eso se reparte B_αβ/2 en cada una, para que
    # la doble suma total sea exactamente B_αα + B_αβ + B_ββ = 1.
    # Motor: seg 0 = mon1(A) = β_tabla3 ; seg 1 = mon2(B) = α_tabla3.
    B_aa_t3, B_ab_t3, B_bb_t3 = Bfr
    Bmat = [[B_bb_t3,       B_ab_t3 / 2.0],   # [0][0]=B_ββ ; [0][1]=B_αβ/2
            [B_ab_t3 / 2.0, B_aa_t3]]         # [1][1]=B_αα ; [1][0]=B_αβ/2
    return m2.Componente(cop, [segA, segB], Bmat), MW

def dp_sistema(solv, cop, T, pts):
    fam = (solv, familia(cop))
    if fam not in KAB_FAM: return None, 0
    ka, kb = KAB_FAM[fam]
    kab = KAB_AB.get((solv,cop), 0.0)
    solvo = construir_solvente(solv)
    copo, MWc = construir_copolimero(cop)
    # kseg: componente 0 = solvente (1 seg), componente 1 = copolimero (2 seg: 0=mon1,1=mon2)
    # solv-mon1 = ka ; solv-mon2 = kb ; mon1-mon2 = kab
    kseg = {
        (0,0,1,0): ka,   # solvente - mon1(alpha)
        (0,0,1,1): kb,   # solvente - mon2(beta)
        (1,0,1,1): kab,  # mon1 - mon2
    }
    ml = m2.PCSAFT_M2([solvo, copo], kseg)
    mv = m2.PCSAFT_M2([solvo], {})
    # conversion w1->x1 con MW del copolimero
    MWs = SOLV[solv][3]
    s=0.0;n=0
    for w1,Pexp in pts:
        ns = w1/MWs; ncop = (1-w1)/MWc; x1 = ns/(ns+ncop)
        try: Pc = m2.presion_de_burbuja(ml,mv,T,x1)
        except Exception: Pc=None
        if Pc and Pc>0: s+=abs(Pc-Pexp)/Pexp; n+=1
    return (100*s/n if n else None), n

if __name__ == "__main__":
    print("%-7s %-9s %8s %8s %8s %9s %8s" % ("Solv","Copol","ka","kb","kab","dP_calc%","dP_paper%"))
    print("-"*62)
    # agrupar por familia para promedio
    fam_acc = {}
    for (solv,cop,T),pts in ct.DATA.items():
        if cop not in COP:   # solo los que tienen Tabla 3
            continue
        fam = (solv, familia(cop))
        if fam not in KAB_FAM: continue
        d,n = dp_sistema(solv,cop,T,pts)
        ka,kb = KAB_FAM[fam]; kab = KAB_AB.get((solv,cop),0.0)
        dpp = DP_FAM.get(fam, float('nan'))
        print("%-7s %-9s %8.4f %8.4f %8.4f %9s %8.2f" %
              (solv,cop,ka,kb,kab, (f"{d:.2f}" if d else "None"), dpp))
        if d is not None:
            fam_acc.setdefault(fam,[0.0,0])
            fam_acc[fam][0]+=d; fam_acc[fam][1]+=1
    print("-"*62)
    print("Promedio por familia (metod.2):")
    allv=[]
    for fam,(s,n) in fam_acc.items():
        print("  %-14s calc=%.2f  paper=%.2f" % (str(fam), s/n, DP_FAM.get(fam,float('nan'))))
        allv.append(s/n)
    if allv:
        print("GLOBAL metod.2 (prom familias): %.2f%%   (paper 2.04%%)" % (sum(allv)/len(allv)))
