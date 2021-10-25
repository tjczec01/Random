# -*- coding: utf-8 -*-
"""

Created on Thu Oct  7 00:10:15 2021

@author: Travis J Czechorski

Github: https://github.com/tjczec01

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com

Advisor: thomas.schwartz@maine.edu

Github:  https://github.com/tjczec01
         https://github.com/TravisCzechorskiUMaine

"""

import math as mt

ekv = [119.8, 164.0]
R = 8.31446261815324 # kg * m**2 / s**2 * K mol
T_star = 2.0
T = [T_star*i for i in ekv]
print(T)
sig = [3.41, 3.38]
sigm = [i*1E-10 for i in sig]
mm = [39.95, 83.798]
mmkg= [i/1000 for i in mm]  #  kg/mol
abg = 6.02E23
M = [i/(1000 * abg) for i in mm]
KB=1.38064852E-23
en = [i*KB for i in ekv]
t0 = [sigm[i]*mt.sqrt(M[i]/en[i]) for i in range(0, len(ekv), 1)]
print(t0)
print(en)
print(M, sigm)
print(1.09E-14/5)
print("Ar: {}K\nKypton: {}K".format(*T))
densr = 2000.0*sigm[0]**3
den = densr/sigm[1]**3
print(den)
TT = 298/ekv[0]
print(TT)
P = 2000.0*R*mmkg[0]*298
Ps = (P * (sigm[0]**3))/ en[0] 
print(P, Ps, 2.4874791318864777 * ekv[1], (Ps*en[1])/(sigm[1]**3))
