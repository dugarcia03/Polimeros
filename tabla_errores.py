#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tabla consolidada de errores %dP: M1 y M2 (calculado) vs Paper.
  M1 = metodologia 1 (pseudo-monomero), paper = Tabla 6.
  M2 = metodologia 2 (segmentos alpha/beta), paper = Tabla 7 (por familia).
M2 solo aplica a los 6 copolimeros con fracciones de enlace (Tabla 3):
  PSBR23, PSBR45, PVAE9, PVAE25, PVAE50, PVAE70.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import importlib.util as iu

# --- M1 ---
spec=iu.spec_from_file_location('ct','correr_todos.py');ct=iu.module_from_spec(spec);spec.loader.exec_module(ct)
# --- M2 ---
specm=iu.spec_from_file_location('m2','pcsaft_m2.py');m2=iu.module_from_spec(specm);specm.loader.exec_module(m2)
specc=iu.spec_from_file_location('cm','correr_m2.py');cm=iu.module_from_spec(specc);specc.loader.exec_module(cm)

def m1_dp(solv,cop,T,pts,kij):
    return ct.dp_system(solv,cop,T,pts,kij)[0]

def m2_dp(solv,cop,T):
    if cop not in cm.COP: return None
    fam=(solv,cm.familia(cop))
    if fam not in cm.KAB_FAM: return None
    d,_=cm.dp_sistema(solv,cop,T,ct.DATA[(solv,cop,T)])
    return d

# paper M2 por familia (Tabla 7)
DP2_PAPER=cm.DP_FAM

print("%-6s %-9s | %8s %8s | %8s %8s" % ("Solv","Cop","M1_calc","M1_pap","M2_calc","M2_pap"))
print("-"*60)
m1c=[];m1p=[];m2c=[];m2p=[]
for (solv,cop,T),pts in ct.DATA.items():
    kij,dp1_paper = ct.KIJ_DP[(solv,cop)]
    d1 = m1_dp(solv,cop,T,pts,kij)
    d2 = m2_dp(solv,cop,T)
    fam=(solv,cm.familia(cop))
    dp2_paper = DP2_PAPER.get(fam)
    s1c = f"{d1:.2f}" if d1 is not None else "-"
    s2c = f"{d2:.2f}" if d2 is not None else "-"
    s2p = f"{dp2_paper:.2f}" if dp2_paper is not None else "-"
    print("%-6s %-9s | %8s %8.2f | %8s %8s" % (solv,cop,s1c,dp1_paper,s2c,s2p))
    if d1 is not None: m1c.append(d1); m1p.append(dp1_paper)
    if d2 is not None: m2c.append(d2); m2p.append(dp2_paper)
print("-"*60)
print("PROMEDIOS (sobre sistemas calculados):")
print("  M1: calc=%.2f%%  paper(Tabla6)=%.2f%%  (n=%d)" % (sum(m1c)/len(m1c), sum(m1p)/len(m1p), len(m1c)))
print("  M2: calc=%.2f%%  paper(Tabla7)=%.2f%%  (n=%d)" % (sum(m2c)/len(m2c), sum(m2p)/len(m2p), len(m2c)))
print()
print("Global paper reportado: M1=4.41%, M2=2.04% (sobre los 24 sistemas del paper)")
