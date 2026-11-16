#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnostico punto-a-punto: Pcalc vs Pexp para sistemas reveladores."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import importlib.util
spec = importlib.util.spec_from_file_location('pv', 'pcsaft_vle.py')
pv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pv)
import importlib.util as iu
spec2 = iu.spec_from_file_location('ct', 'correr_todos.py')
ct = iu.module_from_spec(spec2); spec2.loader.exec_module(ct)

def diag(solv, cop, T, kij):
    ms = ct.SOLV[solv]; mc = ct.COP[cop]
    pv.SOLVENTE = {"nombre":solv,"MW":ms[3],"m_MW":ms[0],"sigma":ms[1],"eps_k":ms[2]}
    pv.COPOLIMERO = {"nombre":cop,"MW":float(mc[3]),"m_MW":mc[0],"sigma":mc[1],"eps_k":mc[2]}
    solvo = pv.construir_solvente(); copo = pv.construir_copolimero()
    ml = pv.PCSAFT([solvo,copo],[[0,kij],[kij,0]]); mv = pv.PCSAFT([solvo],[[0.0]])
    pts = ct.DATA[(solv,cop,T)]
    print(f"\n=== {solv}+{cop} T={T} kij={kij} ===")
    print(f"{'w1':>7} {'x1':>9} {'Pexp[MPa]':>10} {'Pcalc[MPa]':>11} {'ratio':>7}")
    for w1,Pexp in pts:
        x1 = pv.w1_a_x1(w1)
        try:
            Pc = pv.presion_de_burbuja(ml,mv,T,x1)
        except Exception as e:
            Pc = None
        r = (Pc/Pexp) if (Pc and Pexp) else float('nan')
        pcs = f"{Pc:.5f}" if Pc else "None"
        print(f"{w1:>7.3f} {x1:>9.5f} {Pexp:>10.4f} {pcs:>11} {r:>7.2f}")

diag("cC6","PSBR45",333.15, 0.0101)   # control bueno (2%)
diag("cC6","PVAE25",353.15,-0.1068)   # malo (84%), mismo solvente apolar
diag("cC6","PVAE70",353.15,-0.0636)   # malo (81%)
