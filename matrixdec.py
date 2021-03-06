# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 19:17:48 2020

@author: Travis Czechorski tjczec01@gmail.com
"""

import numpy as np
import pprint as pp
import sympy as sp
import scipy.linalg as la
import scipy as sc
import decimal
from decimal import Decimal, getcontext, DefaultContext
De = decimal.Decimal



def transposed(M, pr):
       columns = len(M[0][:])
       rows = len(M)
       vf = 0.0
       a = '{}:.{}f{}'.format('{', pr,'}')
       a1 = '{}'.format(a)
       a2 = str(a1)
       a3 = str(a2.format(vf))
       a4 = De(a3)
       tmat = [[a4 for i in range(rows)] for j in range(columns)]
       for i in range(rows):
              for j in range(columns):
                     vvf = M[i][j]
                     ad = '{}:.{}f{}'.format('{', pr,'}')
                     a1d = '{}'.format(ad)
                     a2d = str(a1d)
                     a3d = str(a2d.format(vvf))
                     a4d = De(a3d)
                     tmat[j][i] = a4d
       return tmat
                     
def dotd(v1, v2, pr):
       vv = sum([x*y for x, y in zip(v1, v2)])
       aa = '{}:.{}f{}'.format('{', pr,'}')
       aa1 = '{}'.format(aa)
       aa2 = str(aa1)
       aa3 = str(aa2.format(vv))
       aa4 = De(aa3)
       return aa4

def zerod(n, pr):
       mm = []
       for ni in range(n):
               vf = 0.0
               a = '{}:.{}f{}'.format('{', pr,'}')
               a1 = '{}'.format(a)
               a2 = str(a1)
               a3 = str(a2.format(vf))
               a4 = De(a3)
               mm.append([a4 for mi in range(n)])
       
       for nn in range(n):
               vfb = 0.0
               ab = '{}:.{}f{}'.format('{', pr,'}')
               a1b = '{}'.format(ab)
               a2b = str(a1b)
               a3b = str(a2b.format(vfb))
               a4b = De(a3b)
               mm[nn][nn] = a4b
       return mm 


def LU_decompositiond(A, pr):
    """Perform LU decomposition using the Doolittle factorisation."""

    L = zerod(len(A), pr)
    U = zerod(len(A), pr)
    N = len(A)
    
    def uvals(Um, k, n):
            ulist = []
            for i in range(k):
                  uu = Um[i][n]
                  ulist.append(uu)
            return ulist
     
    def lvals(Lm, k, n):
            llist = []
            lu = Lm[k]
            lul = lu[0:k]
            return lul
    
    def uvalsd(Um, k, n, pr):
            ulist = []
            for i in range(k):
                  uu = Um[i][n]
                  aa = '{}:.{}f{}'.format('{', pr,'}')
                  aa1 = '{}'.format(aa)
                  aa2 = str(aa1)
                  aa3 = str(aa2.format(uu))
                  aa4 = De(aa3)
                  ulist.append(aa4)
            return ulist
     
    def lvalsd(Lm, k, n, pr):
            llist = []
            lu = Lm[k]
            lul = lu[0:k]
            for i in range(len(lul)):
                     val_ij = lul[i]
                     aa = '{}:.{}f{}'.format('{' , pr,'}')
                     aa1 = '{}'.format(aa)
                     aa2 = str(aa1)
                     aa3 = str(aa2.format(val_ij))
                     aa4 = De(aa3)
                     llist.append(aa4)
            return lul
       
    for k in range(N):
        v1 = 1.0
        a = '{}:.{}f{}'.format('{', pr,'}')
        a1 = '{}'.format(a)
        a2 = str(a1)
        a3 = str(a2.format(v1))
        a4 = De(a3)
        L[k][k] = a4
        v2 = (A[k][k] - dotd(lvalsd(L, k, k, pr), uvalsd(U, k, k, pr), pr)) / L[k][k]
        ab = '{}:.{}f{}'.format('{', pr,'}')
        ab1 = '{}'.format(ab)
        ab2 = str(ab1)
        ab3 = str(ab2.format(v2))
        ab4 = De(ab3)
        U[k][k] = ab4
        for j in range(k+1, N):
              val_i = float((A[k][j] - dotd(lvalsd(L, k, k, pr), uvalsd(U, j, j, pr), pr)) / L[k][k])
              aa = '{}:.{}f{}'.format('{', pr,'}')
              aa1 = '{}'.format(aa)
              aa2 = str(aa1)
              aa3 = str(aa2.format(val_i))
              aa4 = De(aa3)
              U[k][j] = aa4
        for i in range(k+1, N):
              val_ib = float((A[i][k] - dotd(lvalsd(L, i, i, pr), uvalsd(U, k, k, pr), pr)) / U[k][k])
              aab = '{}:.{}f{}'.format('{', pr,'}')
              aa1b = '{}'.format(aab)
              aa2b = str(aa1b)
              aa3b = str(aa2b.format(val_ib))
              aa4b = De(aa3b)
              L[i][k] = aa4b

    return L, U

def forward_subd(L, b, pr):
    """Given a lower triangular matrix L and right-side vector b,
    compute the solution vector y solving Ly = b."""

    y = []
    for i in range(len(b)):
        y.append(b[i])
        for j in range(i):
              val_i = y[i]-(L[i][j]*y[j])
              aa = '{}:.{}f{}'.format('{', pr,'}')
              aa1 = '{}'.format(aa)
              aa2 = str(aa1)
              aa3 = str(aa2.format(val_i))
              aa4 = De(aa3)
              y[i]= aa4
        y[i] = De(y[i])/De(L[i][i])

    return y

def backward_subd(U, y, pr):
    """Given a lower triangular matrix U and right-side vector y,
    compute the solution vector x solving Ux = y."""

    x = [De(0.0) for ix in y]

    for i in range(len(x), 0, -1):
       val_i = (y[i-1] - dotd(U[i-1][i:], x[i:], pr)) / U[i-1][i-1]
       aa = '{}:.{}f{}'.format('{', pr,'}')
       aa1 = '{}'.format(aa)
       aa2 = str(aa1)
       aa3 = str(aa2.format(val_i))
       aa4 = De(aa3)
       x[i-1] = aa4
       
    return x

def lu_solved(L, U, b, pr):
    # Step 1: Solve Uy = b using forward substitution
    # Step 2: Solve Lx = y using backward substitution
    y = forward_subd(L, b, pr)
    yv = []
    for i in range(len(y)):
       val_yi = float(y[i])
       aa = '{}:.{}f{}'.format('{', pr,'}')
       aa1 = '{}'.format(aa)
       aa2 = str(aa1)
       aa3 = str(aa2.format(val_yi))
       aa4 = De(aa3)
       yv.append(aa4)
    x = backward_subd(U, yv, pr)
    return x

def linear_solved(Ad, bd, pr):
    Ld, Ud = LU_decompositiond(Ad, pr)
    x = lu_solved(Ld, Ud, bd, pr)
    return x

def decfunc(Ax, p):
       rows = len(Ax)
       cols = len(Ax[0])
       AD = [[De(ij) for ij in Ax[ii]] for ii in range(len(Ax))]
       for i in range(rows):
              for j in range(cols):
                     decimal.getcontext().prec = 25
                     val_ij = Ax[i][j]
                     aa = '{}:.{}f{}'.format('{',p,'}')
                     aa1 = '{}'.format(aa)
                     aa2 = str(aa1)
                     aa3 = str(aa2.format(val_ij))
                     aa4 = De(aa3)
                     AD[i][j] = aa4
                     Ax[i][j] = aa4
       return AD, Ax
       
def decfuncl(Axl, p):
       vals = len(Axl)
       ADl = [De(il) for il in Axl]
       for i in range(vals):
              val_i = Axl[i]
              aa = '{}:.{}f{}'.format('{',p,'}')
              aa1 = '{}'.format(aa)
              aa2 = str(aa1)
              aa3 = str(aa2.format(val_i))
              aa4 = De(aa3)
              ADl[i] = aa4
       return ADl


pre = 10
AAd = np.array([[3.0, 2.0, 3.0],     [4.0, 6.0, 6.0],     [7.0, 4.0, 9.0]])
AAD, AAX = decfunc(AAd, pre)
ATX = transposed(AAD, pre)
bd = decfuncl([6.0, -4.0, 27.0], pre)
LX, UX = LU_decompositiond(AAD, pre)
xd = linear_solved(AAD, bd, pre)
print(np.array(LX))
print(np.array(UX))
print(xd)
