# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:00:34 2020

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
DefaultContext.prec = 25

pp = pp.pprint
EPS = np.finfo(float).eps
# print("{:.52f}".format(EPS))
As = sp.Matrix([[3, 2, 3],
     [4, 6, 6],
     [7, 4, 9]])
Bs = sp.Matrix([[5, 5], [6, 7], [9, 9]])
AM = [[3, 2, 3],      [4, 6, 6],      [7, 4, 9]]
B = [[5, 5], [6, 7], [9, 9]]
AN = np.array(AM)
# print(np.array(A).T)
# print(len("0000000000000002220446049250313080847263336181640625"))
flatten = lambda l: [item for sublist in l for item in sublist] 

def shape(Ax):
       rows = len(Ax)
       cols = len(Ax[0])
       shape = list((rows, cols))
       ts = tuple((rows, cols))
       print("{} Rows x {} Columns".format(ts[0], ts[1]))
       print(ts)
       return shape
   
def zeros_matrix(rows, cols):
    """
    Creates a matrix filled with zeros.
        :param rows: the number of rows the matrix should have
        :param cols: the number of columns the matrix should have
        :return: list of lists that form the matrix
    """
    M = []
    while len(M) < rows:
        M.append([])
        while len(M[-1]) < cols:
            M[-1].append(0.0)

    return M


def identity_matrix(n):
    """
    Creates and returns an identity matrix.
        :param n: the square size of the matrix
        :return: a square identity matrix
    """
    IdM = zeros_matrix(n, n)
    for i in range(n):
        IdM[i][i] = 1.0

    return IdM    
   
def copy_matrix(M):
    """
    Creates and returns a copy of a matrix.
        :param M: The matrix to be copied
        :return: A copy of the given matrix
    """
    # Section 1: Get matrix dimensions
    rows = len(M)
    cols = len(M[0])

    # Section 2: Create a new matrix of zeros
    MC = zeros_matrix(rows, cols)

    # Section 3: Copy values of M into the copy
    for i in range(rows):
        for j in range(cols):
            MC[i][j] = M[i][j]

    return MC

def check_matrix_equality(A, B, tol=None):
    """
    Checks the equality of two matrices.
        :param A: The first matrix
        :param B: The second matrix
        :param tol: The decimal place tolerance of the check
        :return: The boolean result of the equality check
    """
    # Section 1: First ensure matrices have same dimensions
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        return False

    # Section 2: Check element by element equality
    #            use tolerance if given
    for i in range(len(A)):
        for j in range(len(A[0])):
            if tol is None:
                if A[i][j] != B[i][j]:
                    return False
            else:
                if round(A[i][j], tol) != round(B[i][j], tol):
                    return False

    return True

def check_squareness(A):
    """
    Makes sure that a matrix is square
        :param A: The matrix to be checked.
    """
    if len(A) != len(A[0]):
        raise ArithmeticError("Matrix must be square to inverse.")
        
def check_non_singular(A):
    """
    Ensure matrix is NOT singular
        :param A: The matrix under consideration
        :return: determinant of A - nonzero is positive boolean
                  otherwise, raise ArithmeticError
    """
    det = determinant_fast(A)
    if det != 0:
        return det
    else:
        raise ArithmeticError("Singular Matrix!")
  
def matrix_multiply(A, B):
    """
    Returns the product of the matrix A * B
        :param A: The first matrix - ORDER MATTERS!
        :param B: The second matrix
        :return: The product of the two matrices
    """
    # Section 1: Ensure A & B dimensions are correct for multiplication
    rowsA = len(A)
    colsA = len(A[0])
    rowsB = len(B)
    colsB = len(B[0])
    if colsA != rowsB:
        raise ArithmeticError(
            'Number of A columns must equal number of B rows.')

    # Section 2: Store matrix multiplication in a new matrix
    C = zeros_matrix(rowsA, colsB)
    for i in range(rowsA):
        for j in range(colsB):
            total = 0
            for ii in range(colsA):
                total += A[i][ii] * B[ii][j]
            C[i][j] = total

    return C     
  
def determinant_recursive(A, total=0):
    """
    Find determinant of a square matrix using full recursion
        :param A: the matrix to find the determinant for
        :param total=0: safely establish a total at each recursion level
        :returns: the running total for the levels of recursion
    """
    # Section 1: store indices in list for flexible row referencing
    indices = list(range(len(A)))

    # Section 2: when at 2x2 submatrices recursive calls end
    if len(A) == 2 and len(A[0]) == 2:
        val = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        return val

    # Section 3: define submatrix for focus column and call this function
    for fc in indices:  # for each focus column, find the submatrix ...
        As = copy_matrix(A)  # make a copy, and ...
        As = As[1:]  # ... remove the first row
        height = len(As)

        for i in range(height):  # for each remaining row of submatrix ...
            As[i] = As[i][0:fc] + As[i][fc+1:]  # zero focus column elements

        sign = (-1) ** (fc % 2)  # alternate signs for submatrix multiplier
        sub_det = determinant_recursive(As)  # pass submatrix recursively
        total += sign * A[0][fc] * sub_det  # total all returns from recursion

    return total

def determinant_fast(A):
    # Section 1: Establish n parameter and copy A
    n = len(A)
    AM = copy_matrix(A)
 
    # Section 2: Row ops on A to get in upper triangle form
    for fd in range(n): # A) fd stands for focus diagonal
        for i in range(fd+1,n): # B) only use rows below fd row
            if AM[fd][fd] == 0: # C) if diagonal is zero ...
                AM[fd][fd] == 1.0e-18 # change to ~zero
            # D) cr stands for "current row"
            crScaler = AM[i][fd] / AM[fd][fd] 
            # E) cr - crScaler * fdRow, one element at a time
            for j in range(n): 
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
     
    # Section 3: Once AM is in upper triangle form ...
    product = 1.0
    for i in range(n):
        # ... product of diagonals is determinant
        product *= AM[i][i] 
 
    return product

def print_matrix(M, decimals=3):
    """
    Print a matrix one row at a time
        :param M: The matrix to be printed
    """
    for row in M:
        print([round(x, decimals)+0 for x in row])

def invert_matrix(A, tol=None):
    """
    Returns the inverse of the passed in matrix.
        :param A: The matrix to be inversed
 
        :return: The inverse of the matrix A
    """
    # Section 1: Make sure A can be inverted.
    check_squareness(A)
    check_non_singular(A)
 
    # Section 2: Make copies of A & I, AM & IM, to use for row ops
    n = len(A)
    AM = copy_matrix(A)
    I = identity_matrix(n)
    IM = copy_matrix(I)
 
    # Section 3: Perform row operations
    indices = list(range(n)) # to allow flexible row referencing ***
    for fd in range(n): # fd stands for focus diagonal
        fdScaler = 1.0 / AM[fd][fd]
        # FIRST: scale fd row with fd inverse. 
        for j in range(n): # Use j to indicate column looping.
            AM[fd][j] *= fdScaler
            IM[fd][j] *= fdScaler
        # SECOND: operate on all rows except fd row as follows:
        for i in indices[0:fd] + indices[fd+1:]: 
            # *** skip row with fd in it.
            crScaler = AM[i][fd] # cr stands for "current row".
            for j in range(n): 
                # cr - crScaler * fdRow, but one element at a time.
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
                IM[i][j] = IM[i][j] - crScaler * IM[fd][j]
 
    # Section 4: Make sure IM is an inverse of A with specified tolerance
    if check_matrix_equality(I,matrix_multiply(A,IM),tol):
        return IM
    else:
           # return IM
            raise ArithmeticError("Matrix inverse out of tolerance.")

def permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = list(range(n))
    cycles = list(range(n, n-r, -1))
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return
     
def rotate_left(x, y):
    """
    Left rotates a list x by the number of steps specified
    in y.
    Examples
    ========
    >>> from sympy.utilities.iterables import rotate_left
    >>> a = [0, 1, 2]
    >>> rotate_left(a, 1)
    [1, 2, 0]
    """
    if len(x) == 0:
        return []
    y = y % len(x)
    return x[y:] + x[:y]


def rotate_right(x, y):
    """
    Right rotates a list x by the number of steps specified
    in y.
    Examples
    ========
    >>> from sympy.utilities.iterables import rotate_right
    >>> a = [0, 1, 2]
    >>> rotate_right(a, 1)
    [2, 0, 1]
    """
    if len(x) == 0:
        return []
    y = len(x) - y % len(x)
    return x[y:] + x[:y]  

def minlex(seq, directed=True, is_set=False, small=None):
    """
    Return a tuple where the smallest element appears first; if
    ``directed`` is True (default) then the order is preserved, otherwise
    the sequence will be reversed if that gives a smaller ordering.
    If every element appears only once then is_set can be set to True
    for more efficient processing.
    If the smallest element is known at the time of calling, it can be
    passed and the calculation of the smallest element will be omitted.
    Examples
    ========
    >>> from sympy.combinatorics.polyhedron import minlex
    >>> minlex((1, 2, 0))
    (0, 1, 2)
    >>> minlex((1, 0, 2))
    (0, 2, 1)
    >>> minlex((1, 0, 2), directed=False)
    (0, 1, 2)
    >>> minlex('11010011000', directed=True)
    '00011010011'
    >>> minlex('11010011000', directed=False)
    '00011001011'
    """
    is_str = type(seq)
    seq = list(seq)
    if small is None:
       small = min(seq)
    if is_set:
        i = seq.index(small)
        if not directed:
            n = len(seq)
            p = (i + 1) % n
            m = (i - 1) % n
            if seq[p] > seq[m]:
                seq = list(reversed(seq))
                i = n - i - 1
        if i:
            seq = rotate_left(seq, i)
        best = seq
    else:
        count = seq.count(small)
        if count == 1 and directed:
            best = rotate_left(seq, seq.index(small))
        else:
            # if not directed, and not a set, we can't just
            # pass this off to minlex with is_set True since
            # peeking at the neighbor may not be sufficient to
            # make the decision so we continue...
            best = seq
            for i in range(count):
                seq = rotate_left(seq, seq.index(small, count != 1))
                if seq < best:
                    best = seq
                # it's cheaper to rotate now rather than search
                # again for these in reversed order so we test
                # the reverse now
                if not directed:
                    seq = rotate_left(seq, 1)
                    seq = list(reversed(seq))
                    if seq < best:
                        best = seq
                    seq = list(reversed(seq))
                    seq = rotate_right(seq, 1)
    # common return
    if is_str == str:
        return ''.join(best)
    return tuple(best)       

def cyclic_form(Pl):
        """
        This is used to convert to the cyclic notation
        from the canonical notation. Singletons are omitted.
        Examples
        ========
        >>> from sympy.combinatorics.permutations import Permutation
        >>> Permutation.print_cyclic = False
        >>> p = Permutation([0, 3, 1, 2])
        >>> p.cyclic_form
        [[1, 3, 2]]
        >>> Permutation([1, 0, 2, 4, 3, 5]).cyclic_form
        [[0, 1], [3, 4]]
        See Also
        ========
        array_form, full_cyclic_form
        """
        pt = type(Pl)
        if pt == tuple:
               Pl = list(Pl)
               # return list(Pl)
        elif pt == str:
               raise Exception('Given Value must be either tuple or list')
        elif pt == int:
               raise Exception('Given Value must be either tuple or list')
        elif pt == float:
               raise Exception('Given Value must be either tuple or list')
        elif pt == list:
               Pl = Pl
                # return Pl
        array_form = Pl
        unchecked = [True] * len(Pl)
        cyclic_form = []
        for i in range(len(Pl)):
            if unchecked[i]:
                cycle = []
                cycle.append(i)
                unchecked[i] = False
                j = i
                while unchecked[array_form[j]]:
                    j = array_form[j]
                    cycle.append(j)
                    unchecked[j] = False
                if len(cycle) > 1:
                    cyclic_form.append(cycle)
                    assert cycle == list(minlex(cycle, is_set=True))
        cyclic_form.sort()
        cyclic_form = cyclic_form[:]
        return cyclic_form
 
def transpositions(Pl):
        """
        Return the permutation decomposed into a list of transpositions.
        It is always possible to express a permutation as the product of
        transpositions, see [1]
        Examples
        ========
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[1, 2, 3], [0, 4, 5, 6, 7]])
        >>> t = p.transpositions()
        >>> t
        [(0, 7), (0, 6), (0, 5), (0, 4), (1, 3), (1, 2)]
        >>> print(''.join(str(c) for c in t))
        (0, 7)(0, 6)(0, 5)(0, 4)(1, 3)(1, 2)
        >>> Permutation.rmul(*[Permutation([ti], size=p.size) for ti in t]) == p
        True
        References
        ==========
        .. [1] https://en.wikipedia.org/wiki/Transposition_%28mathematics%29#Properties
        """
        al = [i for i in range(len(Pl))]
        if al == Pl is True:
               return [0, 0]    
        else:
               a = cyclic_form(Pl)
               res = []
               for x in a:
                   nx = len(x)
                   if nx == 2:
                       res.append(tuple(x))
                   elif nx > 2:
                       first = x[0]
                       for y in x[nx - 1:0:-1]:
                           res.append((first, y))
               return res      


def transpose(M):
       columns = len(M[0][:])
       rows = len(M)
       tmat = [[i*0.0 for i in range(rows)] for j in range(columns)]
       for i in range(rows):
              for j in range(columns):
                     tmat[j][i] = M[i][j]
       return tmat

# def transposed(M):
#        columns = len(M[0][:])
#        rows = len(M)
#        vf = 0.0
#        a = '{}:.{}f{}'.format('{', pr,'}')
#        a1 = '{}'.format(a)
#        a2 = str(a1)
#        a3 = str(a2.format(vf))
#        a4 = De(a3)
#        tmat = [[a4 for i in range(rows)] for j in range(columns)]
#        for i in range(rows):
#               for j in range(columns):
#                      vvf = M[i][j]
#                      ad = '{}:.{}f{}'.format('{', pr,'}')
#                      a1d = '{}'.format(ad)
#                      a2d = str(a1d)
#                      a3d = str(a2d.format(vvf))
#                      a4d = De(a3d)
#                      tmat[j][i] = afd
       # return tmat

def factorial(n):
     if n == 0:
            return 1
     else:
            nn = n
            ne = n - 1
            while ne >= 1:
                   nn = nn*ne
                   ne -= 1
            return nn

def perm(n, r):
       pp = factorial(n)/factorial(n-r)
       return int(pp)

def comb(n, r):
       cc = factorial(n)/(factorial(n-r)*factorial(r))
       return cc

def Gamma(n):
       if n <= 0:
              return None
       else:
              return factorial(n-1)

def dot(v1, v2):
     return sum([x*y for x,y in zip(v1,v2)])

def dotd(v1, v2, pr):
       vv = sum([x*y for x, y in zip(v1, v2)])
       aa = '{}:.{}f{}'.format('{', pr,'}')
       aa1 = '{}'.format(aa)
       aa2 = str(aa1)
       aa3 = str(aa2.format(vv))
       aa4 = De(aa3)
       return aa4

def matmul(A, B):
       acolumns = len(A[0][:])
       arows = len(A)
       bcolumns = len(B[0][:])
       brows = len(B)
       if acolumns == brows:
              nmat = [[i*0.0 for i in range(bcolumns)] for j in range(arows)]
              for i in range(arows):
                     Ar = A[i][:]
                     for j in range(bcolumns):
                            Bc = [B[i][j] for i in range(brows)]
                            Cij = dot(Ar, Bc) #sum([i*j for i, j in zip(Ar,Bc)])
                            nmat[i][j] = Cij
              return nmat
       
       elif acolumns != brows:
              raise Exception('Columns of matrix A ({}) needs to equal the Rows of Matrix B ({}) {} != {}'.format(acolumns, brows, acolumns, brows))
       
def kcycle(S, k):
       n = len(S)
       kn = factorial(n + 1)/(factorial(n - k + 1)*k)
       return kn
       

def multiply(n):
    total = 1
    for i in n:
        total *= i
    return total       

def sgnf(m):
       sgn = (-1)**m
       return sgn

def det(A):
       det = []
       colsl = len(A[0][:])
       p = [list(i) for i in list(permutations([pi for pi in range(colsl)]))]
       ts = []
       tns = []
       ss = []
       ais = []
       for pi in range(len(p)):
              tl = transpositions(p[pi])
              ts.append(tl)
              tn = len(transpositions(p[pi]))
              tns.append(tn)
       for i in tns:
              ss.append(sgnf(i))
       for fi in range(len(p)):
                            σ = [i + 1 for i in p[fi]]
                            sig = ss[fi]
                            for i, j in enumerate(σ):
                                   ai = A[i-1][j-1]
                                   ais.append(ai)
                            fin = sig*multiply(ais)
                            det.append(fin)
                            ais.clear()
       return sum(det)


def iden(n):
       mm = []
       for ni in range(n):
              mm.append([mi*0.0 for mi in range(n)])
       
       for nn in range(n):
              mm[nn][nn] = 1
       return mm

def zero(n):
       mm = []
       for ni in range(n):
              mm.append([mi*0.0 for mi in range(n)])
       
       for nn in range(n):
              mm[nn][nn] = 0.0
       return mm              

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

def rowred(A):
       rows = len(A)
       cols = len(A[0])
       itern = 0
       II = iden(rows)
       for i in range(cols):
              start = A[i][i]
              vals = [A[j][i]/start for j in range(1 + itern, rows)]
              nr = [[vals[v]*rv for rv in A[i][:]] for v in range(len(vals))]
              ni = [[vals[v]*rv for rv in II[i][:]] for v in range(len(vals))]
              rrows = [A[iv] for iv in range(1 + itern, cols)]
              rrowsi = [II[iv] for iv in range(1 + itern, cols)]
              # print(np.array(A))
              for g in range(len(rrows)):
                     vv = nr[g]
                     nn = rrows[g]
                     nnr = [i - j for i, j in zip(nn, vv)]
                     vvi = ni[g]
                     nni = rrowsi[g]
                     nnii = [i - j for i, j in zip(nni, vvi)]
                     A[g+1+itern] = nnr
                     II[g+1+itern] = nnii
              itern += 1
       return A, II

def solsys(A, II):
       rows = len(A)
       cols = len(A[0])
       itern = 0
       for i in range(cols):
              start = A[-1-i][-1-i]
              vals = [A[j][-1 - itern]/start for j in range(0, rows-itern)]
              nr = [[vals[v]*rv for rv in A[-1 - itern][:]] for v in range(len(vals))]
              rrows = [A[iv] for iv in range(0, cols-itern-1)]
              ni = [[vals[v]*rv for rv in II[-1 - itern][:]] for v in range(len(vals))]
              rrowsi = [II[iv] for iv in range(0, cols-itern-1)]
              for g in range(len(rrows)):
                     vv = nr[g]
                     nn = rrows[g]
                     nnr = [round(i - j, 5) for i, j in zip(nn, vv)]
                     vvi = ni[g]
                     nni = rrowsi[g]
                     nnii = [round(i - j, 5) for i, j in zip(nni, vvi)]
                     A[g] = nnr
                     II[g] = nnii
              itern += 1
       for i in range(rows):
              start = A[i][i]
              IIv = [iv/start for iv in II[i]]
              AAv = [av/start for av in A[i]]
              A[i] = AAv
              II[i] = IIv
       return A, II

# def LUf(A):
#        rows = len(A)
#        cols = len(A[0])
       
#        def U_mat(row, col):
              
#               BMat = [[0.0*j for j in range(col)] for i in range(row)]
              
#               for i in range(row):
                     
#                      for j in range(col):
                            
#                             if i == j:
#                                    strv = "U_{}{}".format(i, j)
#                                    BMat[i][j] = strv
#                             elif i > j:
#                                    BMat[i][j] = 0.0
#                             elif i < j:
#                                    strv = "U_{}{}".format(i, j)
#                                    BMat[i][j] = strv
#               return BMat
                            
                            
#        def L_mat(row, col):
              
#               BMatl = [[0.0*j for j in range(col)] for i in range(row)]
              
#               for i in range(row):
                     
#                      for j in range(col):
                            
#                             if i == j:
#                                    strv = "L_{}{}".format(i, j)
#                                    BMatl[i][j] = 1.0
#                             elif i < j:
#                                    BMatl[i][j] = 0.0
#                             elif i > j:
#                                    strv = "L_{}{}".format(i, j)
#                                    BMatl[i][j] = strv
#               return BMatl

def matmullu(A, B):
              acolumns = len(A[0][:])
              arows = len(A)
              bcolumns = len(B[0][:])
              brows = len(B)
              if acolumns == brows:
                     nmat = [[i*0.0 for i in range(bcolumns)] for j in range(arows)]
                     for i in range(arows):
                            Ar = A[i][:]
                            for j in range(bcolumns):
                                   Bc = [B[i][j] for i in range(brows)]
                                   print(Bc)
                                   Cc = [A[ia][j] for ia in range(brows)]
                                   print(Cc)
                                   Cij = sum(["{}*{}".format(Bc[i], Cc[i])  for i in range(acolumns)])
                                   print(Cij)
                                   nmat[i][j] = Cij
                     return nmat
              
              elif acolumns != brows:
                     raise Exception('Columns of matrix A ({}) needs to equal the Rows of Matrix B ({}) {} != {}'.format(acolumns, brows, acolumns, brows))
              

def LU_decomposition(A):
    """Perform LU decomposition using the Doolittle factorisation."""

    L = zero(len(A))
    U = zero(len(A))
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
       
    for k in range(N):
        L[k][k] = 1
        U[k][k] = (A[k][k] - dot(lvals(L, k, k), uvals(U, k, k))) / L[k][k]
        for j in range(k+1, N):
            U[k][j] = (A[k][j] - dot(lvals(L, k, k), uvals(U, j, j))) / L[k][k]
        for i in range(k+1, N):
            L[i][k] = (A[i][k] - dot(lvals(L, i, i), uvals(U, k, k))) / U[k][k]

    return L, U

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
            # print(lul)
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
            # print(lul)
            for i in range(len(lul)):
                     val_ij = lul[i]
                     aa = '{}:.{}f{}'.format('{' , pr,'}')
                     aa1 = '{}'.format(aa)
                     aa2 = str(aa1)
                     aa3 = str(aa2.format(val_ij))
                     aa4 = De(aa3)
                     llist.append(aa4)
            # print(lul)
            # print(llist)
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
        # print(ab4)
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

def backward_sub(U, y):
    """Given a lower triangular matrix U and right-side vector y,
    compute the solution vector x solving Ux = y."""

    # x = zero(len(y))
    x = [0.0 for ix in y]

    for i in range(len(x), 0, -1):
      x[i-1] = De((y[i-1] - dot(U[i-1][i:], x[i:])) / U[i-1][i-1])

    return x

def forward_sub(L, b):
    """Given a lower triangular matrix L and right-side vector b,
    compute the solution vector y solving Ly = b."""

    y = []
    for i in range(len(b)):
        y.append(b[i])
        for j in range(i):
            y[i]=y[i]-(L[i][j]*y[j])
        y[i] = y[i]/L[i][i]

    return y



    return x

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

    # x = zerod(len(y))
    x = [De(0.0) for ix in y]

    for i in range(len(x), 0, -1):
       val_i = (y[i-1] - dot(U[i-1][i:], x[i:])) / U[i-1][i-1]
       aa = '{}:.{}f{}'.format('{', pr,'}')
       aa1 = '{}'.format(aa)
       aa2 = str(aa1)
       aa3 = str(aa2.format(val_i))
       aa4 = De(aa3)
       x[i-1] = aa4
       
    return x



def lu_solve(L, U, b):
    # Step 1: Solve Uy = b using forward substitution
    # Step 2: Solve Lx = y using backward substitution
    y = forward_sub(L, b)
    x = backward_sub(U, y)
    return x

def linear_solve(A, b):
    L, U = LU_decomposition(A)
    x = lu_solve(L, U, b)
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

def LUdecomp(a):
    n = len(a)
    for k in range(0,n-1):
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a[i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                a[i,k] = lam
    return a

def sparsity(A):
       rows = len(A)
       cols = len(A[0])
       Zs = 0
       Nzs = 0
       tot = 0
       for i in range(rows):
              for j in range(cols):
                     val_ij = A[i][j]
                     if val_ij == 0 or val_ij == 0.0:
                            Zs += 1
                            tot += 1
                     elif val_ij != 0 or val_ij != 0.0:
                            Nzs += 1
                            tot += 1
       return Zs, Nzs, tot, Zs/tot

def csc_groups(A):
       rows = len(A)
       cols = len(A[0])
       rows_l = []
       cols_l = []
       data_l = []
       Zs = []
       Nzs = []
       tot = 0
       for i in range(rows):
              for j in range(cols):
                     val_ij = A[i][j]
                     if val_ij == 0 or val_ij == 0.0:
                            pass
                     elif val_ij != 0 or val_ij != 0.0:
                            rows_l.append(int(i))
                            cols_l.append(int(j))
                            data_l.append(A[i][j])
       return rows_l, cols_l, data_l

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

Av, IIa = rowred(AM)
B, IF = solsys(Av, IIa)
# pp((np.array(IF)).tolist())
# print(det(A))
# pp(np.array(Av))
# pp(np.array(IIa))
# pp(np.array(B))
# pp(np.array(IF))
AA = np.array([[3, 2, 3],     [4, 6, 6],     [7, 4, 9]])
AAd = np.array([[3.0, 2.0, 3.0],     [4.0, 6.0, 6.0],     [7.0, 4.0, 9.0]])
Ab = [[3.0, 2.0, 3.0, 4.0],     
      [4.0, 6.0, 6.0, 8.0],     
      [7.0, 4.0, 9.0, 11.0],
      [5.0, 4.0, 9.0, 16.0]]

As = shape(Ab)
Adet = determinant_fast(Ab)
Abinv = invert_matrix(Ab, 1)
Abinv2 = la.inv(np.array(Ab))
print(Adet, np.linalg.det(np.array(Ab)), la.det(np.array(Ab)))
print_matrix(Abinv)
print(Abinv2)
AAD, AAX = decfunc(AAd, 25)
# print(AAX, AAD)
ATX = transpose(AAD)
aa = transpose([[3, 4, 7], [2, 6, 4], [3, 6, 9]])
bd = decfuncl([6.0, -4.0, 27.0], 25)
# print(bd)
L, U = LU_decomposition(AA)
Lx, UX = LU_decompositiond(AAD, 25)
# print(Lx, UX)
Px, Llx, UuX = la.lu(ATX)
LUX = la.lu_factor(AAD)
# print(U)
xd = linear_solved(AAD, bd, 25)
# print(xd)
# LUf(AM)
# print(P)
# print(np.array(L))
# print(np.array(U))

P, LL, UU = la.lu(AA.T)
# print(np.array(P))
# print(np.array(LL))
# print(np.array(UU))

b = [6, -4, 27]
As = [[3, 2, 3],     [4, 6, 6],     [7, 4, 9]]
LU = la.lu_factor(As)
x = la.lu_solve(LU, b)
# print(linear_solve(As,b))
# print(x)

ZS, NZS, TOT, sps = sparsity([[3, 0, 3],     [4, 0.0, 0],     [0, 4, 9]])
# print(ZS, NZS, TOT, sps)

S6 = 6**0.5
S6d = De(S6)
# print(S6d)
P = [[13/3 + 7*S6/3, -23/3 - 22*S6/3, 10/3 + 5 * S6],
    [13/3 - 7*S6/3, -23/3 + 22*S6/3, 10/3 - 5 * S6],
    [1/3, -8/3, 10/3]]
# print(P)
row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
data = np.array([3, 2, 3, 4, 6, 6, 7, 4, 9])
# print(sc.sparse.csc_matrix((data,(row,col)), shape=(3,3)).todense() )
rowsl, colsl, datal = csc_groups(As)
# print(rowsl, colsl, datal)
# print(sc.sparse.csc_matrix((datal,(rowsl,colsl)), shape=(3,3)).todense() )
