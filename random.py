# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 01:02:01 2020

@author: Travis Czechorski tjczec01@gmail.com
"""

import sympy as sp
from sympy import *
import math
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, CubicHermiteSpline, Akima1DInterpolator, BPoly
from scipy.linalg import solve_triangular, solve_banded, solve, toeplitz
from scipy.optimize import root
from scipy.sparse import spdiags, kron
from scipy.sparse.linalg import spilu, LinearOperator
from numpy import cosh, zeros_like, mgrid, zeros, eye 

# parameters
nx, ny = 75, 75
hx, hy = 1./(nx-1), 1./(ny-1)

P_left, P_right = 0, 0
P_top, P_bottom = 1, 0

def get_preconditioner():
    """Compute the preconditioner M"""
    diags_x = zeros((3, nx))
    diags_x[0,:] = 1/hx/hx
    diags_x[1,:] = -2/hx/hx
    diags_x[2,:] = 1/hx/hx
    Lx = spdiags(diags_x, [-1,0,1], nx, nx)

    diags_y = zeros((3, ny))
    diags_y[0,:] = 1/hy/hy
    diags_y[1,:] = -2/hy/hy
    diags_y[2,:] = 1/hy/hy
    Ly = spdiags(diags_y, [-1,0,1], ny, ny)

    J1 = kron(Lx, eye(ny)) + kron(eye(nx), Ly)

    # Now we have the matrix `J_1`. We need to find its inverse `M` --
    # however, since an approximate inverse is enough, we can use
    # the *incomplete LU* decomposition

    J1_ilu = spilu(J1)

    # This returns an object with a method .solve() that evaluates
    # the corresponding matrix-vector product. We need to wrap it into
    # a LinearOperator before it can be passed to the Krylov methods:

    M = LinearOperator(shape=(nx*ny, nx*ny), matvec=J1_ilu.solve)
    return M


def diagonal_form(a, upper=1, lower=1):

    """
    a is a numpy square matrix
    this function converts a square matrix to diagonal ordered form
    returned matrix in ab shape which can be used directly for scipy.linalg.solve_banded
    """
    nf = a.shape[1]
    ALb = [] # [[0*kk for kk in a[ii]] for ii in range(len(a))]
    assert(np.all(a.shape ==(nf,nf)))
    
    abd = np.zeros((2*nf-1, nf))
    
    for i in range(nf):
        xxl = list(np.diagonal(a,(nf-1)-i))
        xln = len(xxl)
        while len(xxl) < nf:
               xxl.append(0.0)
        nvn = list(np.diagonal(a,(nf-1)-i))
        abd[i,(nf-1)-i:] = np.diagonal(a,(nf-1)-i)
        if all(v == 0 for v in xxl):
               pass
        else:
               ALb.append(xxl)
        
    for t in range(nf-1): 
        abd[(2*nf-2)-t,:t+1] = np.diagonal(a,t-(nf-1))
        xxlb =  list(np.diagonal(a,t-(n-1)))
        xlnb = len(xxlb)
        while len(xxlb) < nf:
               xxlb.append(0.0)
        if all(vb == 0 for vb in xxlb):
               pass
        else:
               ALb.append(xxlb[::-1])
        # ALb[i] = xxlb[:]

    mid_row_inx = int(abd.shape[0]/2)
    upper_rows = [mid_row_inx - ijj for ijj in range(1, upper+1)]
    upper_rows.reverse()
    upper_rows.append(mid_row_inx)
    lower_rows = [mid_row_inx + ij for ij in range(1, lower+1)]
    keep_rows = upper_rows+lower_rows
    abd = abd[keep_rows,:]
    print(ALb)


    return ALb #abd

def divd(x, y):
       def d2(x, y):
              x1, x2 = x
              y1, y2 = y
              try:
                     d2v = (y2 - y1)/(x2 - x1)
              except:
                     d2v = 0
              return d2v
       def d3(x, y):
              x1, x2, x3 = x 
              y1, y2, y3 = y
              # v1 = (y3 - y2)/((x3 - x2)*(x2 - x1)) - (y2 - y1)/((x2 - x1)*(x3 - x1))
              v1 = d2([x2, x3], [y2, y3])
              v2 = d2([x1, x2], [y1, y2])
              v3b = v1 - v2
              v3 = x3 - x1
              vf = v3b/v3
              return vf
       def d4(x, y):
              x1, x2, x3, x4 = x
              y1, y2, y3, y4 = y
              v1 = d3([x2, x3, x4], [y2, y3, y4])
              v2 = d3([x1, x2, x3], [y1, y2, y3])
              v3 = x4 - x1
              vf = (v1 - v2)/v3
              return vf

       xl = len(x)
       if xl == 2:
              return d2(x, y)
       elif xl == 3:
              return d3(x, y)
       elif xl == 4:
              return d4(x, y)


def cfuncf1(x, y, Sl):
       xl = len(x)
       dvs = []
       lams = [1]
       mus = []
       hs = list(np.diff(x))
       x0 = x[0]
       y0 = y[0]
       xx1 = x[1]
       yy1 = y[1]
       def mf(h1, h2):
              hf = h1/(h1 + h2)
              return hf
       def lf(h1, h2):
              hf = h2/(h1 + h2)
              return hf
       for ii in range(0, xl-2):
              hm = mf(hs[ii], hs[ii+1])
              mus.append(hm)
              hl = lf(hs[ii], hs[ii+1])
              lams.append(hl)
       for i in range(xl):
              if i == 0:
                     dva = divd([x[0], x[0], x[1]], [y[0], y[0], y[1]]) 
                     dvaa = dva - Sl[0]
                     dvs.append(6*dvaa)
              elif i > 0 and i < xl - 1:
                     # print(i)
                     dvb = divd([x[i-1], x[i], x[i+1]], [y[i-1], y[i], y[i+1]])
                     dvs.append(6*dvb)
              else:
                     dvc = divd([x[-2], x[-1], x[-1]], [y[-2], y[-1], y[-1]]) + S[-1]
                     dvs.append(6*dvc)
       fc = [2, mus[0]]
       fr = [2, lams[0]]
       while len(fr) < xl:
              fc.append(0)
              fr.append(0)
       mus.append(lams[0])
       dia = toeplitz(fc, fr)
       dia2 = np.zeros((xl, xl))  # This is a banded matrix representation.
       dia2b = [list(dia2[ii]) for ii in range(len(dia2))]
       
       for i in range(xl):
              for j in range(xl):
                     ii = i + 1
                     if i >= 1:
                            ij = i - 1
                     else:
                            ij = 0
                     if i == j:
                            dia2b[i][j] = 2
                            
                     elif ii == j:
                            dia2b[i][j] = lams[i]
                     elif ij == j:
                            dia2b[i][i - 1] = mus[i - 1]
       dia[-1][-2] = fr[-1]
       Ms = solve(dia2b, np.array(dvs))
       
       def cfunc(xi, yi, m, h):
              xs = Symbol('x')
              xl = len(xi)
              CS = []
              for i in range(1, xl-1):
                     t1 = (m[i - 1]*(xi[i] - xs)**3)/(6*h[i])
                     t2 = (m[i]*(xs - xi[i - 1])**3)/(6*h[i])
                     t3 = (yi[i - 1] - ((m[i - 1]*(h[i]**2))/6))*((xi[i] - xs)/h[i])
                     t4 = (yi[i] - ((m[i]*(h[i]**2))/6))*((xs - xi[i - 1])/h[i])
                     tf = t1 + t2 + t3 + t4
                     print(tf)
                     eq = tf.expand()
                     # tff = sp.simplify(tf)
                     # print(tff)
                     # eq = tff.expand()
                     print(eq)
                     cl = [eq.coeff(xs, i) for i in reversed(range(4))]
                     print(cl)
                     CS.append(cl)
              return CS
       
       ansf = cfunc(x, y, Ms, hs)
       return ansf

x3i = [0, 1, 2, 3]
y3i = [0, 0.5, 2, 1.5]
S = [0.2, 0.0, 0.0, -1]
ff = cfuncf1([0, 1, 2], [0.15505102572168228, 0.6449489742783174, 1.0], [0.0, 0])
print(ff)
n = 3
A = np.zeros((3, n))
AL = [list(A[ii]) for ii in range(len(A))]
# print(AL)
dx = np.diff([1, 2, 3])
# print(list(dx))
x2i = [-1, 0, 3]
y2i = [0.5, 0, 3]
def func1(xi, yi):
       nn = len(xi)
       u = nn - 2
       l = nn - 2
       dx = np.diff(xi)
       # dx = np.diff([0, 1, 2])
       dxr = dx.reshape([dx.shape[0]] + [1] * (yi.ndim - 1))
       slope = np.diff(yi, axis=0) / dxr
       A = np.zeros((nn, nn))  # This is a banded matrix representation.
       AL = [list(A[ii]) for ii in range(len(A))]
       # # b = np.empty((nn,) + yi.shape[1:], dtype=yi.dtype)
       # b = [0, 0, 0]
       # abb = [list(A[ii][:]) for ii in range(len(A))] # list(ab)
       # BL = list(b)
       # A[1, 1:-1] = 2 * (dx[:-1] + dx[1:])  # The diagonal
       # A[0, 2:] = dx[:-1]                   # The upper diagonal
       # A[-1, :-2] = dx[1:]                  # The lower diagonal
       # # b[1:-1] = 3 * (dxr[1:] * slope[:-1] + dxr[:-1] * slope[1:])
       A = np.zeros((3, 3))  # This is a standard matrix.
       b = np.empty((3,) + yi.shape[1:], dtype=yi.dtype)

       A[0, 0] = 2/(xi[1] - xi[0])
       A[0, 1] = 1/(xi[1] - xi[0])
       A[1, 0] = 1/(xi[1] - xi[0])
       A[1, 1] = 2 * (1/(xi[1] - xi[0]) + 1/(xi[2] - xi[1]))
       A[1, 2] = 1/(xi[2] - xi[1])
       A[2, 1] = 1/(xi[2] - xi[1])
       A[2, 2] = 2/(xi[2] - xi[1])

       b[0] = 3 * (yi[1] - yi[0])/((xi[1] - xi[0])**2)
       b[1] = 3 * ((yi[1] - yi[0])/((xi[1] - xi[0])**2) + (yi[2] - yi[1])/((xi[2] - xi[1])**2))
       b[2] = 3 * (yi[2] - yi[1])/((xi[2] - xi[1])**2)
       # abss = diagonal_form(A, u, l)
       print(A)
       # print(np.array(b))
       # print(abss)
       # for i in range(nn):
       #        for j in range(nn):
       #               val = A[i-1][j-1]
       #               abb[u + i - j - 1][j-1] = val
       #               # abss[u + i - j - 1, j - 1] = A[i - 0,j - 1]
       # print(np.array(abb))
       km = solve(A, b)
       # km = solve_banded((l, u), abss, b)
       print(km)
       return km
       
# func1(np.array([1, 2, 3]), np.array([0.15505102572168228, 0.6449489742783174, 1.0]))
# ks = func1(np.array(x2), np.array(y2))

def acfunc(x, y, k):
       xl = len(list(x))
       n = len(x)
       A = np.zeros((3, n))
       AL = list(A)
       yl = len(list(y))
       al = []
       bl = []
       slope = (y[-1] - y[0])/(x[-1] - x[0])
       for i in range(1, xl):
                     av = k[i-1]*(x[i] - x[i-1]) - (y[i] - y[i-1])
                     al.append(av)
                     bv = -k[i]*(x[i] - x[i-1]) + (y[i] - y[i-1])
                     bl.append(bv)
                     # aal = [0*k for k in range(yl)]
                     # al.append(aal)
                     # al.append(aal)
                     
       
       # if xl == yl:
       #        for i in range(1, xl):
       #               aa = []
       #               kv = k[i-1]
       #               for j in range(1, yl):
       #                      # kv = k[j-1]
       #                      aav = kv*(x[i] - x[i - 1]) - (y[i] - y[i - 1])
       #                      aa.append(aav)
       #                      if i != j:
       #                             al[i][j] = aav
       #                      elif i == j:
       #                             al[i][j] = 0.0
                                   
                     # al.append(aa)
                     # aa.clear()
       return al, bl
  # np.array(x2, np.array(y2)                   
ks = func1(np.array([1, 2, 3]), np.array([0.15505102572168228, 0.6449489742783174, 1.0]))
print(acfunc(np.array([1, 2, 3]), np.array([0.15505102572168228, 0.6449489742783174, 1.0]), ks))
# ks = func1(np.array(x2), np.array(y2))
# print(acfunc(np.array(x2), np.array(y2), ks))
S6 = 6**0.5
x = Symbol('x')
Eq = (x**3*sp.cos(x/2) + 1/2)*sp.sqrt(4 - x**2)
ans = sp.Integral(Eq, (x, -2, 2)).evalf(10)
# print(ans)
Pc = [[13/3 + 7*S6/3, -23/3 - 22*S6/3, 10/3 + 5 * S6],
     [13/3 - 7*S6/3, -23/3 + 22*S6/3, 10/3 - 5 * S6],
     [1/3, -8/3, 10/3]]
X = [0, 1, 2]
XX = [0, 1, 2, 3]
# print(Pc)
# X = [1, 2, 3]

Y = [0.15505102572168228, 0.6449489742783174, 1.0]

cs = CubicSpline(X, Y, axis=1)
css = PchipInterpolator(X, Y)
csss = CubicHermiteSpline(X, Y, [0, 0, 0])
cs4 = Akima1DInterpolator(X, Y)
cs5 = BPoly(Pc, XX)
# print(cs.c)
# print(css.c)
# print(csss.c)
# print(cs4.c)
# print(cs5)
# print(-0.9797959   + 1.46969385)