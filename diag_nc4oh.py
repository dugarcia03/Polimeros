#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Corre los sistemas nC4OH (PVAE9, PVAVC12) CON el termino de asociacion
activado para el solvente (kAB=0.008752, epsAB=2489.23, Tabla 4b).
Compara %dP con y sin asociacion, vs paper (metod.1).
Usa SOLO parametros del paper (Tablas 4b, 5, 6).
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import math, importlib.util
spec = importlib.util.spec_from_file_location('pv', 'pcsaft_vle.py')
pv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)
spec2 = importlib.util.spec_from_file_location('ct', 'correr_todos.py')
ct = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(ct)

# params asociacion del nC4OH (Tabla 4b)
KAB_NC4OH, EPSAB_NC4OH = 0.008752, 2489.23

def dp(solv, cop, T, kij, con_asoc):
    ms = ct.SOLV[solv]; mc = ct.COP[cop]
    if con_asoc and solv == "nC4OH":
        solvo = pv.Componente(solv, ms[0]*ms[3], ms[1], ms[2], KAB_NC4OH, EPSAB_NC4OH)
    else:
        solvo = pv.Componente(solv, ms[0]*ms[3], ms[1], ms[2])
    copo = pv.Componente(cop, mc[0]*float(mc[3]), mc[1], mc[2])
    ml = pv.PCSAFT([solvo,copo],[[0,kij],[kij,0]]); mv = pv.PCSAFT([solvo],[[0.0]])
    # necesitamos w1_a_x1 con MW correctos
    pv.SOLVENTE = {"MW":ms[3]}; pv.COPOLIMERO = {"MW":float(mc[3])}
    pts = ct.DATA[(solv,cop,T)]
    s=0.0;n=0; detalle=[]
    for w1,Pexp in pts:
        x1 = pv.w1_a_x1(w1)
        try: Pc = pv.presion_de_burbuja(ml,mv,T,x1)
        except Exception as e: Pc=None
        if Pc and Pc>0:
            s+=abs(Pc-Pexp)/Pexp; n+=1; detalle.append((w1,Pexp,Pc))
    return (100*s/n if n else None), n, detalle

for solv,cop,T,kij,dpp in [("nC4OH","PVAE9",353.15,0.0230,4.05),
                           ("nC4OH","PVAVC12",353.15,-0.0126,5.18)]:
    d0,n0,_ = dp(solv,cop,T,kij,False)
    d1,n1,det = dp(solv,cop,T,kij,True)
    print(f"\n{solv}+{cop}  kij={kij}  (paper {dpp}%)")
    print(f"  SIN asociacion: %dP={d0:.2f}  (n={n0})")
    print(f"  CON asociacion: %dP={d1:.2f}  (n={n1})")
    print(f"  {'w1':>7} {'Pexp':>9} {'Pcalc':>10} {'ratio':>7}")
    for w1,pe,pc in det:
        print(f"  {w1:>7.3f} {pe:>9.4f} {pc:>10.5f} {pc/pe:>7.2f}")
