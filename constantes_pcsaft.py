# -*- coding: utf-8 -*-
"""
===============================================================================
  CONSTANTES UNIVERSALES DE PC-SAFT  (fuente EXTERNA al paper de Arce)
===============================================================================

  Estas constantes NO aparecen impresas en una tabla del paper de Arce, Mattedi
  & Aznar (2006). El paper las usa y las referencia como [23] en sus Ecs. 31-32
  ("the constants aik and bik were fitted with thermophysical properties of pure
  n-alkanes [23]"), pero los valores numericos provienen de:

      J. Gross, G. Sadowski, "Perturbed-Chain SAFT: An Equation of State Based
      on a Perturbation Theory for Chain Molecules",
      Ind. Eng. Chem. Res. 40 (2001) 1244-1260.   DOI: 10.1021/ie0003887
      (Tabla 1 de ese articulo).

  Se mantienen en este archivo aparte para dejar claro que son de una fuente
  externa (aunque citada por el paper de Arce). Son universales: iguales para
  todos los sistemas, no ajustables.

  Verificadas fila por fila contra el PDF original de Gross & Sadowski (2001).
===============================================================================
"""

# a0i, a1i, a2i  (para I1 -> termino de dispersion de primer orden, Ec. 28)
A_COEF = [
    [0.9105631445, -0.3084016918, -0.0906148351],
    [0.6361281449,  0.1860531159,  0.4527842806],
    [2.6861347891, -2.5030047259,  0.5962700728],
    [-26.547362491, 21.419793629, -1.7241829131],
    [97.759208784, -65.255885330, -4.1302112531],
    [-159.59154087, 83.318680481, 13.776631870],
    [91.297774084, -33.746922930, -8.6728470368],
]

# b0i, b1i, b2i  (para I2 -> termino de dispersion de segundo orden, Ec. 29)
B_COEF = [
    [0.7240946941, -0.5755498075,  0.0976883116],
    [2.2382791861,  0.6995095521, -0.2557574982],
    [-4.0025849485, 3.8925673390, -9.1558561530],
    [-21.003576815, -17.215471648, 20.642075974],
    [26.855641363,  192.67226447, -38.804430052],
    [206.55133841, -161.82646165,  93.626774077],
    [-355.60235612, -165.20769346, -29.666905585],
]
