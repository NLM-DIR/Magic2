"""
sl21_b.py -- Kac modules of sl(2|1) with b kept as a FORMAL VARIABLE.

Companion to sl21.py.  Same algebra, same generators, same checks; the one
difference is that the Kac-Dynkin weight ``b`` is never given a number.  It is
the SymPy symbol ``B`` throughout, so every matrix, every bracket and every
invariant comes out as an exact function of b.  Run it as

        python sl21_b.py -a <int>                 (e.g. -a 1: R(1, b))
        python sl21_b.py -a <int> -N <k>          (k-jet Matryoshka, 1 <= k <= 4)
        python sl21_b.py -h                       (help)

There is no -b and no -y: that is the whole point of this file.

Why a formal b
--------------
Two things become possible that a numeric b cannot give:

  * **Derivatives of every order.**  d^p R / db^p is taken by SymPy's ``diff``
    on each entry, for any p.  sl21.py could only ever reach p = 1, and only
    by the finite difference R(a,1) - R(a,0), which is exact solely because
    the present generators happen to be affine in b.  Here nothing is assumed:
    the derivative is the real thing, and it stays correct the moment the
    generators stop being affine.

  * **The zeta-Hermitian real form**, which mixes the odd pairs (u,v) and
    (w,x) with coefficients in sqrt(b).  sqrt(b) needs b as a symbol; and once
    sqrt(b) is present the b-dependence is no longer affine, so the higher
    derivatives above stop being a formality and start carrying content.
    (Not built in this file yet -- this file is the ground it stands on.)

``B`` is declared ``positive``, which is what makes ``sqrt(B)**2 -> B``
simplify by itself.  It also means the generic module is typical: SymPy knows
b != 0, and b = a+1 stays an open condition rather than an assumed one.

The N-layer Matryoshka is now a genuine TAYLOR JET
--------------------------------------------------
        block (I, J) of the N x N layer grid  =  (1/p!) d^p R/db^p ,  p = I - J

i.e. R(a, b + n) expanded in the nilpotent Jordan shift n (n^N = 0), which is
exactly the Taylor series truncated at order N-1.  sl21.py fills only the
p = 1 sub-diagonal; that is the special case of this, correct there because
every p >= 2 jet of an affine matrix vanishes.  The run reports which jet
orders are actually nonzero, so the day the generators cease to be affine the
report says so instead of quietly dropping terms.

``a`` and ``b`` are the Kac-Dynkin weights, in Kac's own convention

        b = (a + y) / 2            equivalently        y = 2b - a

where y is the module's central-charge eigenvalue on the top layer.  The Kac
label is the one in which typicality reads off directly: R(a,b) is atypical
exactly at

        b = 0        and        b = a + 1

which are the vanishing points of the raising scales in ``_scales``.  In terms
of the rho-shifted weights mu = a+1, lambda = y-1,

        2b = mu + lambda           2(a + 1 - b) = mu - lambda

so the atypicality polynomial factorises as

        C_2 = mu^2 - lambda^2 = 4 b (a + 1 - b) = -4 b (b - a - 1)

The adjoint is R(1,1); the two atypical fundamentals are R(1,0) and R(0,1).

The module splits into four sl(2) layers, block-diagonal for e,f,h and scalar
for Y (eigenvalues y, y-1, y-1, y-2 with y = 2b - a):

        L1 = R[a]      (top,    Y = y)          dim a+1
        L2 = R[a+1]    (middle, Y = y-1)        dim a+2
        L3 = R[a-1]    (middle, Y = y-1)        dim a     (empty when a=0)
        L4 = R[a]      (bottom, Y = y-2)        dim a+1

The odd generators map between layers. The two lowering operators v, x
(lower-triangular) are b-independent; all b-dependence sits in the two raising
operators u, w (upper-triangular), through the four reduced matrix elements
below. Each odd block is an sl(2) doublet intertwiner in the integer basis used
by sl2.py, so all entries stay rational.

Structure constants:

        [h,e]=2e  [h,f]=-2f  [e,f]=h
        [Y,u]=u   [Y,v]=-v   [Y,w]=w   [Y,x]=-x
        [h,u]=-u  [h,v]=v     [h,w]=w   [h,x]=-x
        [e,u]=w  [e,x]=-v     [f,v]=-x   [f,w]=-u
        {u,v}=(Y+h)/2   {w,x}=(Y-h)/2   {v,w}=-e   {u,x}=-f

The three invariants
--------------------
C_2, C_3 and T are UNIVERSAL: each is one fixed element of U(g), written down
once and then evaluated in whatever module is asked for.  No module under
study ever supplies the tensors used to build them.

        C_2 = h^2 - Y^2 + 2(ef + fe) + 2 N_2
        C_3 = 8 . Sum d^ABC(ref) M_A M_B M_C
        T   = (2/3) (N_4 - N_2)

with N_2 = uv - vu + wx - xw, N_4 the fully antisymmetric quartic in the odd
generators, and d^ABC(ref) the fixed cubic tensor of the anchor R(0,0).  In the
rho-shifted weights mu = a+1, lambda = y-1 their eigenvalues are

        C_2 = mu^2 - lambda^2            [the atypicality polynomial, monic]
        C_3 = lambda (mu^2 - lambda^2)   = lambda . C_2
        T   = C_2 . chi                  [operator identity]
        T^2 = C_2^2                      [central: Gorelik's theorem]

C_2 and T carry no metric at all.  C_3 needs one tensor, and the master
equation d_ABC(R) = Tr(Y)_R . c_ABC makes that tensor fixed: it is read off
once, on an anchor with Tr(Y) != 0, and never again.

The anchors R(0,0) and R(1,1/2) stay NUMERIC even here, and should: the fixed
tensor is a tensor of numbers, computed once and cached.  The symbol lives only
in the module under study.  --casimirs is therefore still exact, but it is much
slower than in sl21.py, because the 512 d-components and the metric inverse are
now symbolic in b.  Without it a run does structure constants + chi + the
anticenter T and stays quick.
"""

import sys
import argparse
from math import factorial
from itertools import permutations
from collections import namedtuple

from sympy import (Rational, Integer, Symbol, sympify, simplify, diff,
                   sqrt, cancel)

from matrix import Matrix
from algebra import Algebra, EVEN, ODD, is_zero_scalar

# The formal weight.  ``positive`` is not decoration: it is what lets
# sqrt(B)**2 collapse back to B in the zeta-Hermitian real form, and it tells
# every zero-test that the generic module is typical (b != 0).  b = a+1 is left
# open -- SymPy carries it as an unresolved condition, which is correct.
B = Symbol('b', positive=True)

# canonical numbers, fixed once for the whole family
Y, E, F, H, U, V, W, X = 0, 1, 2, 3, 4, 5, 6, 7
PARITY = {Y: EVEN, E: EVEN, F: EVEN, H: EVEN, U: ODD, V: ODD, W: ODD, X: ODD}

_half = Rational(1, 2)
# independent nonzero (super-)brackets, in canonical order i<=j
INDEPENDENT = {
    (E, F): {H: 1}, (E, H): {E: -2}, (F, H): {F: 2},
    (Y, U): {U: 1}, (Y, V): {V: -1}, (Y, W): {W: 1}, (Y, X): {X: -1},
    (H, U): {U: -1}, (H, V): {V: 1}, (H, W): {W: 1}, (H, X): {X: -1},
    (E, U): {W: 1}, (E, X): {V: -1}, (F, V): {X: -1}, (F, W): {U: 1},
    (U, V): {Y: _half, H: _half}, (W, X): {Y: _half, H: -_half},
    (V, W): {E: -1}, (U, X): {F: -1},
}


# -- 1. structure constants f(i, j, k), initialized explicitly ---------------

def structure_constants():
    """Build f(i, j, k) from INDEPENDENT, adding graded partners:
    [B,A} = -(-1)^(|A||B|) [A,B}  (sign +1 for two odds, -1 otherwise)."""
    table = {}
    for (i, j), terms in INDEPENDENT.items():
        table[(i, j)] = dict(terms)
        sign = 1 if (PARITY[i] == ODD and PARITY[j] == ODD) else -1
        table[(j, i)] = {k: sign * c for k, c in terms.items()}

    def f(i, j, k):
        return table.get((i, j), {}).get(k, 0)

    return f


# Module-level f so callers can do  from sl21 import f
f = structure_constants()


# -- 2. construction of R(a, b) ----------------------------------------------

def _sl2(j):
    """e, f, h of the (j+1)-dimensional sl(2) irrep R[j] (integer basis)."""
    n = j + 1
    e = Matrix(n, n, lambda i, k: (i + 1) * (j - i) if k == i + 1 else 0)
    fm = Matrix(n, n, lambda i, k: 1 if i == k + 1 else 0)
    h = Matrix(n, n, lambda i, k: (j - 2 * i) if i == k else 0)
    return e, fm, h


def _shape(kind, js, letter):
    """Odd doublet-intertwiner block on source R[js]: list of (row, col, val).

    kind 'up' maps R[js]->R[js+1], 'down' maps R[js]->R[js-1]. v and w share the
    weight-raising component; x, u are the weight-lowering ones (with u carrying
    the sign that makes the pair a genuine sl(2) doublet)."""
    out = []
    if kind == 'up':                          # target R[js+1]
        if letter in ('v', 'w'):
            for m in range(js + 1):
                out.append((m, m, js + 1 - m))
        if letter == 'x':
            for m in range(js + 1):
                out.append((m + 1, m, 1))
        if letter == 'u':
            for m in range(js + 1):
                out.append((m + 1, m, -1))
    else:                                     # target R[js-1]
        if letter in ('v', 'w'):
            for m in range(1, js + 1):
                out.append((m - 1, m, m))
        if letter == 'x':
            for m in range(js):
                out.append((m, m, -1))
        if letter == 'u':
            for m in range(js):
                out.append((m, m, 1))
    return out


def _scales(a, b):
    """The eight block scales: lowering (v,x) b-independent, raising (u,w)
    carrying all the b-dependence.

    In the Kac label the two atypicality conditions are manifest:
    m21, m43 vanish at b = 0;  m31, m42 vanish at b = a+1."""
    D = 4 * (a + 1)
    s12, s13, s24, s34 = 1, 1, -1, 1
    m21 = Rational(4, D) * b
    m43 = Rational(4, D) * b
    m31 = Rational(4, D) * (a + 1 - b)
    m42 = -Rational(4, D) * (a + 1 - b)
    return (s12, s13, s24, s34, m21, m31, m42, m43)


def _build(a, b, scales=None):
    """Assemble the eight 4(a+1)-square generator matrices for R(a, b).

    ``scales`` overrides ``_scales(a, b)``.  The zeta ground uses that to hand
    in the same eight numbers written in alpha and gamma instead of in b, so
    that the gauge below cancels exactly instead of leaving ``(a+1-b)/gamma``
    debris behind."""
    if scales is None:
        scales = _scales(a, b)
    s12, s13, s24, s34, m21, m31, m42, m43 = scales
    y = 2 * b - a

    d = [a + 1, a + 2, a, a + 1]                     # layer dims (L3 empty if a=0)
    off = [0, d[0], d[0] + d[1], d[0] + d[1] + d[2]]
    N = sum(d)
    reps = [a, a + 1, a - 1, a]                      # sl(2) label per layer
    yof = [y, y - 1, y - 1, y - 2]                   # Y eigenvalue per layer

    Emat = [[0] * N for _ in range(N)]
    Fmat = [[0] * N for _ in range(N)]
    Hmat = [[0] * N for _ in range(N)]
    Ymat = [[0] * N for _ in range(N)]
    for L in range(4):
        if d[L] == 0:
            continue
        e, fm, h = _sl2(reps[L])
        o = off[L]
        for i in range(d[L]):
            for k in range(d[L]):
                Emat[o + i][o + k] = e[i, k]
                Fmat[o + i][o + k] = fm[i, k]
                Hmat[o + i][o + k] = h[i, k]
            Ymat[o + i][o + i] = yof[L]

    Umat = [[0] * N for _ in range(N)]
    Vmat = [[0] * N for _ in range(N)]
    Wmat = [[0] * N for _ in range(N)]
    Xmat = [[0] * N for _ in range(N)]

    def put(M, Lt, Ls, triples, scale):
        for (r, c, v) in triples:
            M[off[Lt] + r][off[Ls] + c] += scale * v

    # lowering blocks carry v, x  (b-independent)
    for (Lt, Ls, kind, js, sc) in [(1, 0, 'up',   a,     s12),
                                   (2, 0, 'down', a,     s13),
                                   (3, 1, 'down', a + 1, s24),
                                   (3, 2, 'up',   a - 1, s34)]:
        if d[Lt] == 0 or d[Ls] == 0:
            continue
        put(Vmat, Lt, Ls, _shape(kind, js, 'v'), sc)
        put(Xmat, Lt, Ls, _shape(kind, js, 'x'), -sc)   # x sign flipped

    # raising blocks carry u, w  (all the b-dependence)
    for (Lt, Ls, kind, js, sc) in [(0, 1, 'down', a + 1, m21),
                                   (0, 2, 'up',   a - 1, m31),
                                   (1, 3, 'up',   a,     m42),
                                   (2, 3, 'down', a,     m43)]:
        if d[Lt] == 0 or d[Ls] == 0:
            continue
        put(Umat, Lt, Ls, _shape(kind, js, 'u'), sc)
        put(Wmat, Lt, Ls, _shape(kind, js, 'w'), -sc)   # w sign flipped

    def toM(grid):
        return Matrix(N, N, lambda i, k: grid[i][k])

    return (N, toM(Ymat), toM(Emat), toM(Fmat), toM(Hmat),
            toM(Umat), toM(Vmat), toM(Wmat), toM(Xmat))


def Rsl21(a, b=B):
    """Return R(a, b): the 4(a+1)-dimensional Kac module of sl(2|1).

    ``a`` is a non-negative integer.  ``b`` defaults to the formal symbol ``B``,
    which is how this file is meant to be used; passing a number still works and
    is what the numeric anchors of the cubic tensor do."""
    if a < 0:
        raise ValueError("a must be a non-negative integer")
    N, Ym, Em, Fm, Hm, Um, Vm, Wm, Xm = _build(a, b)
    g = Algebra(name=f"R({a},{b})")
    g.add(Ym, "Y", EVEN, number=Y)
    g.add(Em, "e", EVEN, number=E)
    g.add(Fm, "f", EVEN, number=F)
    g.add(Hm, "h", EVEN, number=H)
    g.add(Um, "u", ODD, number=U)
    g.add(Vm, "v", ODD, number=V)
    g.add(Wm, "w", ODD, number=W)
    g.add(Xm, "x", ODD, number=X)
    return g


# -- 2b. the two grounds -----------------------------------------------------
#
# A "ground" is the answer to four questions about the scalar entries: which
# symbols carry the weight, how to differentiate with respect to y, how to
# conjugate, and how to decide zero.  There are exactly two.
#
#   Chevalley ground : entries are rational in the symbol B.  d/dy = (1/2)d/db.
#                      Conjugation is the identity (everything is real in b).
#   zeta ground      : entries live in Q(sqrt ints)(AL, GA), the two radicals
#                      being carried as INDEPENDENT positive symbols
#
#                          AL = sqrt(b) ,      GA = sqrt(b - a - 1)
#
#                      with the relation GA^2 = AL^2 - a - 1 imposed only where
#                      it is needed, namely in the zero test.  Independence is
#                      what makes conjugation a one-line substitution instead of
#                      a fight with SymPy's radicals; the relation is what makes
#                      the brackets close.  Since b = AL^2 and b-a-1 = GA^2,
#
#                          dAL/dy = 1/(4 AL) ,   dGA/dy = 1/(4 GA)
#
#                      so d/dy stays inside the field and the derivatives are
#                      negative powers of the radicals, as they must be.

AL = Symbol('alpha', positive=True)      # alpha = sqrt(b)          = sqrt(K_1)
GA = Symbol('gamma', positive=True)      # gamma = sqrt(b - a - 1)  = sqrt(K_2)


def dy_chevalley(x):
    """d/dy on the Chevalley ground:  y = 2b - a, so d/dy = (1/2) d/db."""
    return diff(sympify(x), B) / 2


def dy_zeta(x):
    """d/dy on the zeta ground, by the chain rule through AL and GA."""
    x = sympify(x)
    return (diff(x, AL) / AL + diff(x, GA) / GA) / 4


# -- 2c. jets and the Matryoshka ---------------------------------------------
#
# JETS ARE TAKEN IN y = 2b - a, not in b (paper eq. 5.1: M_Y then carries a
# plain identity on the first sub-diagonal, not 2I).  The two conventions are
# similar via P = diag(1, 1/2, 1/4, ...) (x) Id, which scales band p by 2^-p.

def jets(rep, N, dy):
    """[J_0, ..., J_{N-1}] with J_p = {k: (1/p!) (d/dy)^p M_k}.

    The p-th Taylor coefficient of every generator, in y.  This is a REAL
    derivative -- entrywise ``diff`` through ``Matrix.applyfunc``, which is
    T-blind, so matrix.py and algebra.py know nothing about it.  It replaces
    the finite difference R(a,1) - R(a,0) of sl21.py, which was exact only
    because the Chevalley generators happen to be affine in b and which could
    never reach p >= 2.  After the zeta gauge the entries are square roots and
    every order is alive."""
    nums = rep.numbers()
    raw = [{k: rep[k] for k in nums}]               # raw[p] = (d/dy)^p M
    for _ in range(1, N):
        prev = raw[-1]
        raw.append({k: prev[k].applyfunc(dy) for k in nums})
    return [{k: Rational(1, factorial(p)) * M for k, M in J.items()}
            for p, J in enumerate(raw)]


def jet_orders(rep, N, dy, is_zero=is_zero_scalar):
    """Which of the orders p = 0 .. N-1 are actually nonzero.

    Reported on every run.  In the Chevalley basis the answer is [0, 1] and the
    Matryoshka is a two-band matrix; in the zeta basis it is every order, which
    is the whole reason the jet construction had to be made real."""
    return [p for p, J in enumerate(jets(rep, N, dy))
            if any(not is_zero(M[i, j]) for M in J.values()
                   for i in range(M.rows) for j in range(M.cols))]


def matryoshka(rep, N, dy):
    """The N-layer Matryoshka of ``rep``: block Toeplitz, lower triangular.

        [ mu       0      0    0 ]
        [ mu'      mu     0    0 ]        band p = I - J carries
        [ mu''/2   mu'    mu   0 ]        (1/p!) (d/dy)^p mu
        [ mu'''/6  mu''/2 mu'  mu ]

    This is R at y*Id_N + n with n the strictly-lower Jordan shift.  Since
    n^N = 0 the Taylor series terminates, so the substitution is exact for ANY
    y-dependence -- no affine assumption anywhere -- and the Matryoshka
    satisfies the same y-independent structure constants, checked by the same f.

    The corner blocks survive at every N because band p of {M_i, M_j} is
    (1/p!)(d/dy)^p (d^a_ij mu_a) and every mu_a is at most affine in y, so it
    vanishes identically for p >= 2.  Paper eq. 5.3 is the p = 2 instance of a
    statement with no upper bound on N.  N = 1 returns ``rep`` unchanged."""
    if N < 1:
        raise ValueError("N must be a positive integer")
    J = jets(rep, N, dy)
    d = rep.dim
    D = N * d
    g = Algebra(name=f"M[{rep.name}, N={N}]")
    for k in rep.numbers():
        band = [J[p][k] for p in range(N)]

        def entry(i, j, band=band, d=d):
            I, Jc = i // d, j // d
            if I < Jc:                              # strictly upper: empty
                return 0
            return band[I - Jc][i % d, j % d]       # band p = I - J

        g.add(Matrix(D, D, entry), rep.title(k), rep.parity(k), number=k)
    return g


def Rmatryoshka(a, N, b=B):
    """The Chevalley Matryoshka MR(a,b,N), jets in y."""
    return matryoshka(Rsl21(a, b), N, dy_chevalley)


# -- 2d. the zeta-Hermitian gauge --------------------------------------------
#
# The eight Chevalley generators are conjugated by a diagonal D so that the odd
# blocks come out symmetric in the two radicals
#
#       alpha = sqrt(K_1) = sqrt(b) ,     gamma = sqrt(K_2) = sqrt(b - a - 1)
#
# K_1 = <Lambda+rho|beta> and K_2 = <Lambda+rho|beta_2> being the two Kac
# atypicality conditions.  This is the observation the whole construction rests
# on: THE MATRIX ENTRIES ARE THE SQUARE ROOTS OF THE KAC CONDITIONS, so their
# reality is governed by the sign of those conditions.
#
# The atypical points b(b-a-1) = 0 are out of scope: D is singular there, and
# the Chevalley and zeta-Hermitian pictures are genuinely different at those
# points (paper sec. 7).
#
# Which layer pair carries which condition is NOT the layer distance:
#
#       L1-L2 and L3-L4 carry K_1 = b            (scales m21, m43)
#       L1-L3 and L2-L4 carry K_2 = b - a - 1    (scales m31, m42)
#
# so a zeta constant on layers exists, and is forced up to overall sign:
#
#       zeta = diag( 1, sign K_1, sign K_2, sign K_1 . sign K_2 )   per layer
#
# The two consistency conditions zeta_L1 zeta_L4 read through L2 and through L3
# agree, which is why this closes at all.  At a = 0 the layer L3 is empty and
# this collapses to the paper's diag(1, sign b, sign b, sign b sign(b-1)).

CASES = {                       # case -> (sign K_1, sign K_2)
    1: (+1, +1),                # b > a+1     : both conditions positive
    2: (+1, -1),                # 0 < b < a+1 : the quark window
    3: (-1, -1),                # b < 0
}

CASE_REGION = {
    1: "b > a+1        (K_1 > 0, K_2 > 0)",
    2: "0 < b < a+1    (K_1 > 0, K_2 < 0)   the quark window",
    3: "b < 0          (K_1 < 0, K_2 < 0)",
}

# zeta M_k^dagger zeta = M_{PAIR[k]}.  Not "each generator is zeta-Hermitian":
# in the Chevalley basis the raising and lowering partners are exchanged.  It is
# the Hermitian combinations lambda_1..lambda_7 of the paper that are fixed
# individually; e = (l1+i l2)/2 and f = (l1-i l2)/2 give zeta e^dag zeta = f,
# and likewise (u,v) <-> l6,l7 and (w,x) <-> l4,l5.
PAIR = {Y: Y, H: H, E: F, F: E, U: V, V: U, W: X, X: W}


def _scales_zeta(a):
    """``_scales`` written in AL and GA:  b -> AL^2,  a+1-b -> -GA^2."""
    s12, s13, s24, s34 = 1, 1, -1, 1
    m21 = AL**2 / (a + 1)
    m43 = AL**2 / (a + 1)
    m31 = -GA**2 / (a + 1)
    m42 = GA**2 / (a + 1)
    return (s12, s13, s24, s34, m21, m31, m42, m43)


def _wt(j, m):
    """sqrt(m! . j!/(j-m)!) -- the intra-layer sl(2) orthonormalisation.

    Within one sl(2) layer of label j the integer basis of ``_sl2`` has
    e[m,m+1] = (m+1)(j-m) against f[m+1,m] = 1, so e is not the adjoint of f.
    The ratio that repairs it is d_{m+1}/d_m = sqrt((m+1)(j-m)), whose product
    is this.  It is b-INDEPENDENT, so it never touches a derivative or a jet.

    This is why a > 0 is not the paper's a = 0 case with letters changed, and
    it first bites at a = 1, not a = 2: the middle layer is R[a+1] = R[2], where
    e[0,1] = e[1,2] = 2.  At a = 0 the only nontrivial layer is R[1], whose
    integer basis is already orthonormal, which is why the paper never meets
    this."""
    return sqrt(Integer(factorial(m)) * Integer(factorial(j))
                / Integer(factorial(j - m)))


def gauge_diagonal(a):
    """The diagonal D of the zeta gauge, as a list of 4(a+1) entries.

    HOW THIS WAS OBTAINED IS NOT PART OF THE ARGUMENT (MATRIX_CONTEXT sec. 14).
    It is written down here as a formula and then checked a posteriori by
    ``zeta_relations`` on the finished matrices, entry by entry, in exact
    arithmetic.  For the record, the route was: zeta-Hermiticity forces

        (d_i/d_j)^2 = zeta_i zeta_j c_i c_j M_k[j,i] / M_{PAIR(k)}[i,j]

    on every transposed pair of nonzero entries; propagating that from d_0 = 1
    fixes D, and the surviving equations are then constraints rather than
    definitions.  All four bands close, and the sl(2) index m drops out of each
    one -- that cancellation is the content, and it is what ``zeta_relations``
    re-establishes without reference to any of this.

        layer   sl(2) label   dim     constant
        L1      a             a+1     (a+1)
        L2      a+1           a+2     alpha
        L3      a-1           a       gamma sqrt(a(a+1))
        L4      a             a+1     alpha gamma

    times the intra-layer ``_wt``.  At a = 0, L3 is empty and this reduces to
    diag(1, alpha, alpha, alpha gamma), which is the paper's S of eq. 7.3 up to
    the sign of the last layer.  That sign is free: flipping D on a whole layer
    flips both directions of every band touching it, and zeta-Hermiticity only
    ever sees (d_i/d_j)^2."""
    labels = [a, a + 1, a - 1, a]
    dims = [a + 1, a + 2, a, a + 1]
    consts = [Integer(a + 1), AL, GA * sqrt(Integer(a * (a + 1))), AL * GA]
    D = []
    for L in range(4):
        for m in range(dims[L]):
            D.append(consts[L] * _wt(labels[L], m))
    return D


def Rzeta(a):
    """R(a,b) in the zeta-Hermitian basis: D R D^-1, entries in AL and GA.

    Same eight matrices for all three cases -- no factor of i appears anywhere.
    The case enters only through the two sign rules, ``zeta_matrix`` and
    ``bar``.  Since D is a similarity, the structure constants are untouched;
    that is checked rather than assumed."""
    if a < 0:
        raise ValueError("a must be a non-negative integer")
    n, Ym, Em, Fm, Hm, Um, Vm, Wm, Xm = _build(a, AL**2, _scales_zeta(a))
    D = gauge_diagonal(a)

    def gauge(M):
        return Matrix(n, n, lambda i, j: cancel(D[i] * M[i, j] / D[j]))

    g = Algebra(name=f"Rzeta({a})")
    for M, title, parity, num in [(Ym, "Y", EVEN, Y), (Em, "e", EVEN, E),
                                  (Fm, "f", EVEN, F), (Hm, "h", EVEN, H),
                                  (Um, "u", ODD, U), (Vm, "v", ODD, V),
                                  (Wm, "w", ODD, W), (Xm, "x", ODD, X)]:
        g.add(gauge(M), title, parity, number=num)
    return g


def zeta_is_zero(a):
    """Zero test on the zeta ground: impose GA^2 = AL^2 - a - 1, then decide.

    AL and GA are carried as independent symbols so that conjugation is a
    substitution; the relation between them is real and has to be put back
    before anything can be called zero.  This is the only place it enters, and
    it is the reason the brackets close."""
    rel = sqrt(AL**2 - a - 1)

    def test(x):
        return is_zero_scalar(sympify(x).subs(GA, rel))

    return test


def zeta_reduce(a):
    """Eliminate GA in favour of AL, for code that predates the zeta ground.

    ``_casimir_report`` and its helpers use bare ``simplify`` and structural
    ``!= 0`` in a dozen places.  Rather than thread a zero test through all of
    them, the Casimir work is done on a copy of the representation with GA
    substituted out, where the ordinary ``is_zero_scalar`` is already correct.
    C_2, C_3 and T are polynomials in the generators, so this changes nothing
    about what is computed."""
    rel = sqrt(AL**2 - a - 1)

    def reduce_rep(rep):
        g = Algebra(name=rep.name + " [GA eliminated]")
        for k in rep.numbers():
            g.add(rep[k].applyfunc(lambda x: sympify(x).subs(GA, rel)),
                  rep.title(k), rep.parity(k), number=k)
        return g

    return reduce_rep


def bar(x, case):
    """Complex conjugation on the zeta ground:  AL -> c1 AL,  GA -> c2 GA.

    Every entry is a real rational times square roots of positive integers
    times a monomial in AL and GA, so flipping the sign of whichever radical is
    imaginary IS complex conjugation.  SymPy's own ``conjugate`` cannot do this
    -- it has no way to know the signs of b and b-a-1 -- and must never be
    called on these entries.

    The conjugation signs are the zeta layer signs: alpha is real exactly when
    K_1 > 0, gamma exactly when K_2 > 0.  That coincidence is what makes one
    set of matrices serve all three cases, since zeta-Hermiticity on the band
    carrying K_j needs epsilon_j c_j = 1, and epsilon_j = c_j = sign K_j."""
    c1, c2 = CASES[case]
    return sympify(x).subs({AL: c1 * AL, GA: c2 * GA}, simultaneous=True)


def dagger(M, case):
    """M^dagger on the zeta ground: transpose, then ``bar`` each entry."""
    return Matrix(M.cols, M.rows, lambda i, j: bar(M[j, i], case))


def zeta_matrix(a, case, N=1):
    """zeta, or zeta_N = J (x) zeta for the Matryoshka.

    zeta itself is diagonal, constant on each sl(2) layer, with the signs of
    ``CASES``.  For N > 1 a block-diagonal zeta CANNOT work: the Matryoshka is
    block lower triangular and its adjoint is block upper triangular, and no
    block-diagonal conjugation changes that.  J is the reversal permutation of
    the N generations, J_pq = delta_{p,N-1-q}, J^2 = 1.  Since
    J (n^dag)^p J = n^p for the generation shift n, zeta_N = J (x) zeta sends
    band p back to band p (paper sec. 6), provided every Taylor coefficient is
    separately zeta-Hermitian -- which ``zeta_jet_relations`` checks.

    Consequence: the invariant form <p|zeta_N|q> pairs generation p with
    generation N+1-p, so it is anti-diagonal in the generation index."""
    s1, s2 = CASES[case]
    per_layer = [1, s1, s2, s1 * s2]
    dims = [a + 1, a + 2, a, a + 1]
    diag = []
    for s, dl in zip(per_layer, dims):
        diag.extend([s] * dl)
    d = len(diag)

    def entry(i, j):
        if i % d != j % d:                       # zeta is diagonal internally
            return 0
        if i // d + j // d != N - 1:             # J reverses the generation
            return 0
        return diag[i % d]

    return Matrix(N * d, N * d, entry)


# One reported zeta relation: its label, whether it holds, and the leftover.
ZetaCheck = namedtuple("ZetaCheck", "label ok residual")


def _is_zero_matrix(M, is_zero):
    return all(is_zero(M[i, j])
               for i in range(M.rows) for j in range(M.cols))


def zeta_relations(rep, Z, case, is_zero):
    """The reality condition, checked on the finished matrices:

        zeta M_k^dagger zeta = M_{PAIR(k)}        for all eight
        zeta^2 = 1

    This is the a-posteriori proof.  Nothing about how the gauge diagonal was
    derived enters it: it reads the eight matrices as they stand and multiplies
    them out.  Returns data only, never printed."""
    out = []
    for k in rep.numbers():
        partner = PAIR[k]
        r = Z @ dagger(rep[k], case) @ Z - rep[partner]
        lhs = f"zeta {rep.title(k)}^dag zeta"
        out.append(ZetaCheck(f"{lhs} = {rep.title(partner)}",
                             _is_zero_matrix(r, is_zero), r))
    sq = Z @ Z - Matrix.one(Z.rows)
    out.append(ZetaCheck("zeta^2 = 1", _is_zero_matrix(sq, is_zero), sq))
    return out


def zeta_jet_relations(rep, a, case, N, is_zero):
    """Each Taylor coefficient must be separately zeta-Hermitian.

    This is the hypothesis zeta_N = J (x) zeta rests on, so it is checked
    rather than assumed.  It holds because ``bar`` commutes with d/dy: the
    derivative of a real multiple of a radical is a real multiple of the same
    radical to a negative power, and differentiation therefore never moves an
    entry across the real/imaginary divide."""
    Z1 = zeta_matrix(a, case, 1)
    out = []
    for p, J in enumerate(jets(rep, N, dy_zeta)):
        ok = True
        worst = None
        for k in rep.numbers():
            r = Z1 @ dagger(J[k], case) @ Z1 - J[PAIR[k]]
            if not _is_zero_matrix(r, is_zero):
                ok = False
                worst = r
        out.append(ZetaCheck(f"jet p={p}: zeta mu_p^dag zeta = mu_p(pair)",
                             ok, worst))
    return out


def Rchi(a, N=1):
    """The grading operator chi (super identity / Klein operator) for R(a,b)
    (N=1) or the Matryoshka MR(a,b,N).

    chi is diagonal and constant on each sl(2) layer of R(a,b):

        L1 -> -1,   L2 -> +1,   L3 -> +1,   L4 -> -1
        dims  a+1        a+2         a          a+1

    That sign pattern is forced (up to overall sign) by which layers the odd
    generators connect: chi carries opposite signs on any two layers linked by
    an odd block, so it commutes with the even generators (block-diagonal) and
    anticommutes with the odd ones (off-diagonal between opposite-sign layers).
    It is +1 on one Z2-parity class and -1 on the other, and chi^2 = 1: a
    Clifford involution. Treated as ODD, its super-bracket with a generator g is
    a commutator when g is even and an anticommutator when g is odd -- exactly
    the two relations we test.

    chi is independent of b (b only rescales blocks within their layers, never
    moving weight between layers), so in the Matryoshka it is simply the same
    diagonal repeated block-diagonally on each of the N layers (no sub-diagonal).
    """
    signs_per_layer = [-1, +1, +1, -1]
    dims_per_layer = [a + 1, a + 2, a, a + 1]
    diag_one = []
    for s, dl in zip(signs_per_layer, dims_per_layer):
        diag_one.extend([s] * dl)
    diag = diag_one * N                     # repeat on each Matryoshka layer
    D = len(diag)
    return Matrix(D, D, lambda i, j: diag[i] if i == j else 0)


# One reported chi relation: its label, whether it holds, and the leftover.
ChiCheck = namedtuple("ChiCheck", "label ok residual")


def verify_metric_inverse(g_lower, g_upper):
    """Verify that g_upper is the left inverse of g_lower: g^ij g_jk = delta^i_k.
    
    Returns a list of (i, k) pairs where the inverse property fails."""
    n = 8
    failures = []
    for i in range(n):
        for k in range(n):
            delta = 1 if i == k else 0
            entry = simplify(sum(g_upper[i, j] * g_lower[j, k] for j in range(n)))
            if entry != delta:
                failures.append((i, k, entry, delta))
    return failures


def chi_relations(rep, chi, is_zero=is_zero_scalar):
    """The chi relations as plain data (never printed).

    For every generator g: the super-bracket of chi (odd) with g, i.e.
    [chi, g] = chi g - g chi when g is even and {chi, g} = chi g + g chi when g
    is odd; each must be the zero matrix. Plus the involution chi^2 = 1. Zero-
    testing goes through is_zero_scalar (invariant 4), same policy as check()."""
    D = chi.rows

    def is_zero_matrix(M):
        return all(is_zero(M[i, j])
                   for i in range(M.rows) for j in range(M.cols))

    out = []
    for k in rep.numbers():
        M = rep[k]
        if rep.parity(k) == ODD:
            br = chi @ M + M @ chi                  # {chi, odd}
            op, cl = "{", "}"
        else:
            br = chi @ M - M @ chi                  # [chi, even]
            op, cl = "[", "]"
        out.append(ChiCheck(f"{op}chi, {rep.title(k)}{cl} = 0",
                            is_zero_matrix(br), br))

    sq = chi @ chi - Matrix.one(D)
    out.append(ChiCheck("chi^2 = 1", is_zero_matrix(sq), sq))
    return out


# -- the lower-index Killing metric ------------------------------------------

def supertrace(chi, M):
    """STr(M) = Trace(chi . M), the chi-graded trace (exact, sympified)."""
    P = chi @ M
    return sympify(sum((P[i, i] for i in range(P.rows)), 0))


def killing_metric(rep, chi, N):
    """Lower-index metric on the generators, as a Matrix over the canonical
    numbers, computed directly from the matrices in the current basis:

        even a,b :  g_ab = (1/2N) STr(a b + b a)
        odd  i,j :  g_ij = (1/2N) STr(i j - j i)

    with STr(M) = Trace(chi . M) and the ordinary matrix product throughout.
    N is the Matryoshka number."""
    nums = rep.numbers()
    two_N = 2 * N

    def entry(p, q):
        i, j = nums[p], nums[q]
        A, B = rep[i], rep[j]
        if rep.parity(i) == ODD and rep.parity(j) == ODD:
            P = A @ B - B @ A
        else:
            P = A @ B + B @ A
        return simplify(supertrace(chi, P) / two_N)

    return Matrix(len(nums), len(nums), entry)


def cubic_d(rep, chi, N):
    """The 512 constants d_ABC = (1/2) STr(A [[B,C]]), keyed (A,B,C) over the
    canonical numbers, with [[B,C]] = B C - C B when B,C both odd else B C + C B.

    The 1/2 is the proper normalisation: [[B,C]] carries two terms.  It fixes
    the scale of the constant tensor c_ABC = -d_ABC(R(0,0))/4, so it is the
    convention to quote wherever c_ABC appears.  It also gives

        d_YYY = STr(Y^3) = -6 (a+1) (y-1)

    which the verification section checks directly."""
    nums = rep.numbers()
    d = {}
    for A in nums:
        MA = rep[A]
        for B in nums:
            MB = rep[B]
            for C in nums:
                MC = rep[C]
                if rep.parity(B) == ODD and rep.parity(C) == ODD:
                    bc = MB @ MC - MC @ MB
                else:
                    bc = MB @ MC + MC @ MB
                value = simplify(supertrace(chi, MA @ bc) / 2)
                
                d[(A, B, C)] = value / N
    return d


def trace(M):
    """Ordinary matrix trace -- NOT the supertrace.

    Tr(Y) is the object the anomaly is proportional to; STr(Y) vanishes
    identically on every Kac module (Y is supertraceless), so the two must not
    be confused here."""
    return simplify(sum(M[i, i] for i in range(M.rows)))


# -- the master equation:  d_ABC is proportional to Tr(Y) --------------------
#
# The cubic tensor of a module is one FIXED tensor times a single number, and
# that number is Tr(Y):
#
#       d_ABC(R) = Tr(Y)_R . c_ABC ,      c_ABC = -d_ABC(R(0,0))/4
#
# with c_ABC carrying no dependence on (a, y) whatsoever.  All the module
# dependence of the cubic tensor sits in the trace of the hypercharge, exactly
# as a gauge anomaly does.  Two consequences used throughout:
#
#   * d_ABC is strictly ADDITIVE over a composite module, since Tr(Y) is;
#   * d_ABC vanishes identically on any composite with total Tr(Y) = 0.
#
# On the Kac module R(a,b), summing y, y-1, y-1, y-2 over the four sl(2)
# layers of dimensions a+1, a+2, a, a+1:
#
#       Tr(Y) = 4 (a+1) (y-1) = 4 mu lambda
#
# NOTE the (a+1).  Tr(Y) = 4(y-1) is the a = 0 case only.
#
# The anchor here is R(0,0), dimension 4, where Tr(Y) = -4 != 0.  Any module
# with Tr(Y) != 0 would serve; the adjoint R(1,1) would not, since Tr(Y) = 0
# there and its whole d-tensor vanishes with it.

_D_REF_CACHE = {}


def cubic_d_reference():
    """The d-tensor and Tr(Y) of the anchor module R(0,0). Computed once."""
    if not _D_REF_CACHE:
        rep = Rsl21(0, 0)
        _D_REF_CACHE['d'] = cubic_d(rep, Rchi(0, 1), 1)
        _D_REF_CACHE['trY'] = trace(rep[Y])          # = -4
    return _D_REF_CACHE['d'], _D_REF_CACHE['trY']


def cubic_d_proportionality(d_lower, tr_Y):
    """Check d_ABC(R) = (constant tensor) . Tr(Y), componentwise over all 512.

    Cross-multiplied rather than divided, so the adjoint (Tr(Y) = 0, whole
    tensor zero) is tested on the same footing as everything else instead of
    being a special case:

        d_ABC(R) . Tr(Y)_ref  ==  d_ABC(ref) . Tr(Y)_R

    ``tr_Y`` is Tr(Y) per layer, i.e. the raw trace divided by N, matching the
    1/N already carried by ``cubic_d``.

    Returns (ratio, failures): the proportionality constant Tr(Y)_R/Tr(Y)_ref
    -- which is the anomaly coefficient A(R) = -mu*lambda -- and the list of
    components that fail, as (key, got, expected).  Data only, never printed.
    An empty failure list is the whole content of the claim."""
    d_ref, trY_ref = cubic_d_reference()

    failures = []
    for k, v in d_lower.items():
        lhs = simplify(v * trY_ref)
        rhs = simplify(d_ref[k] * tr_Y)
        if not is_zero_scalar(lhs - rhs):
            failures.append((k, v, simplify(rhs / trY_ref)))

    return simplify(tr_Y / trY_ref), failures


def cubic_d_upper(d_lower, g_upper):
    """The 512 upper-index constants d^ABC = g^AA' g^BB' g^CC' d_A'B'C',
    keyed (A,B,C) over the canonical numbers."""
    nums = range(8)
    d_up = {}
    for A in nums:
        for B in nums:
            for C in nums:
                total = 0
                for Ap in nums:
                    gA = g_upper[A, Ap]
                    if gA == 0:
                        continue
                    for Bp in nums:
                        gB = g_upper[B, Bp]
                        if gB == 0:
                            continue
                        for Cp in nums:
                            gC = g_upper[C, Cp]
                            if gC == 0:
                                continue
                            total += gA * gB * gC * d_lower[(Ap, Bp, Cp)]
                d_up[(A, B, C)] = simplify(total)
    return d_up


# -- the FIXED cubic tensor, read off once from a y = 0 anchor ---------------
#
# The master equation makes C_3 universal.  Since d_ABC(R) = Tr(Y)_R . c_ABC,
# fixing the cubic tensor takes nothing more than one module with Tr(Y) != 0.
# Every y = 0 Kac module R(a, a/2) qualifies: there lambda = -1 and
# Tr(Y) = -4(a+1).  The tensor is read off once, cached, and evaluated in
# whatever module is asked for.
#
# Which y = 0 anchor is used does not matter, once one factor is divided out.
# Going from anchor a to anchor 0:
#
#     d_ABC   picks up  Tr(Y) ratio          = (a+1)
#     raising picks up  three inverse forms  = (a+1)^-3   [g_AB = I(R) g_AB(ref), I = a+1]
#     -----------------------------------------------------------------------
#     d^ABC   picks up                         (a+1)^-2
#
# so multiplying by (a+1)^2 makes every y = 0 anchor deliver the IDENTICAL
# tensor.  ``reference_anchor_check`` verifies this componentwise rather than
# taking the argument on trust.

REF_ANCHOR_A = 0            # anchor R(0,0): y = 0, dimension 4, Tr(Y) = -4

_REF_D_CACHE = {}


def reference_cubic_tensor(anchor_a=REF_ANCHOR_A):
    """d^ABC of the y = 0 anchor R(anchor_a, anchor_a/2), times (anchor_a+1)^2.

    The normalisation is what makes the choice of anchor immaterial.  Cached:
    computed once per anchor per process."""
    if anchor_a not in _REF_D_CACHE:
        rep = Rsl21(anchor_a, Rational(anchor_a, 2))       # y = 2b - a = 0
        chi = Rchi(anchor_a, 1)
        g_up = upper_killing_metric(killing_metric(rep, chi, 1), rep)
        d_up = cubic_d_upper(cubic_d(rep, chi, 1), g_up)
        norm = (anchor_a + 1)**2
        _REF_D_CACHE[anchor_a] = {k: simplify(norm * v) for k, v in d_up.items()}
    return _REF_D_CACHE[anchor_a]


def reference_anchor_check(anchors=(0, 1)):
    """Different y = 0 anchors must give the same fixed tensor, componentwise.

    Returns the list of components that disagree -- empty is the claim."""
    base = reference_cubic_tensor(anchors[0])
    failures = []
    for anchor_a in anchors[1:]:
        other = reference_cubic_tensor(anchor_a)
        for k, v in base.items():
            if not is_zero_scalar(v - other[k]):
                failures.append((anchor_a, k, other[k], v))
    return failures


def casimir_quadratic(rep, g_upper, coef=1):
    """C_2 by contraction: coef * (1/2) g^{AB} M_A M_B.

    A diagnostic on the metric section, not the reported C_2 -- that is
    ``casimir_quadratic_direct``.  Contracting with a module's own inverse form
    carries that module's Dynkin index I(R) = mu, so

        g^{AB} from the module's own supertrace form  ->  C_2 / (4 mu)
        g^{AB} from R(0,0)                            ->  C_2 / 4

    where C_2 is the universal operator.  ``coef`` rescales on the way out.

    Returns a dict with:
        'even_sector': contribution from even generators only
        'odd_sector':  contribution from odd generators only
        'total':       even + odd
        'is_scalar_multiple': True if result is a multiple of the identity
        'eigenvalue':  that multiple when it is, else None
    """
    nums = rep.numbers()
    d = rep.dim
    
    # Start with zero matrix
    C2_even = Matrix.zero(d, d)
    C2_odd = Matrix.zero(d, d)
    
    # Contributions from even-even pairs (canonical numbers 0, 1, 2, 3)
    even_indices = [0, 1, 2, 3]
    for a in even_indices:
        for b in even_indices:
            c_ab = g_upper[a, b]  # g^{ab} from the 8x8 metric matrix
            if c_ab != 0:
                term = (c_ab / 2) * (rep[a] @ rep[b])
                C2_even = C2_even + term
    
    # Contributions from odd-odd pairs (canonical numbers 4, 5, 6, 7)
    odd_indices = [4, 5, 6, 7]
    for a in odd_indices:
        for b in odd_indices:
            # Access the metric at the correct position for odd generators
            c_ab = g_upper[a, b]
            if c_ab != 0:
                term = (c_ab / 2) * (rep[a] @ rep[b])
                C2_odd = C2_odd + term
    
    # Note: cross terms (even-odd) should be zero due to metric structure
    
    return _c2_result(coef * C2_even, coef * C2_odd)


def _c2_result(C2_even, C2_odd):
    """Package a C_2 built in two sectors: total, scalarity, eigenvalue.

    Shared by ``casimir_quadratic`` (contracted) and
    ``casimir_quadratic_direct`` (explicit polynomial) so the two routes are
    reported identically and the scalarity test lives in exactly one place.

    Scalarity is a property of the module, not of the Casimir: it holds on the
    irreps by Schur and fails on the Matryoshka N > 1, where C_2 picks up a
    nilpotent off-diagonal part.  Commutation with every generator, checked by
    ``casimir_commutes``, is the property that survives there."""
    C2_total = (C2_even + C2_odd).applyfunc(simplify)
    d = C2_total.rows

    off_diag_nonzero = any(not is_zero_scalar(C2_total[i, j])
                           for i in range(d) for j in range(d) if i != j)

    diag = [simplify(C2_total[i, i]) for i in range(d)]
    diag_constant = all(is_zero_scalar(x - diag[0]) for x in diag) if diag else True

    is_scalar_multiple = (not off_diag_nonzero) and diag_constant

    return {
        'even_sector': C2_even,
        'odd_sector': C2_odd,
        'total': C2_total,
        'is_scalar_multiple': is_scalar_multiple,
        'eigenvalue': diag[0] if (is_scalar_multiple and d > 0) else None,
    }


# One reported Casimir-commutation relation: generator label, whether the
# bracket vanished, and the leftover matrix.
CasimirCheck = namedtuple("CasimirCheck", "label ok residual")


def casimir_cubic(rep, d_upper, coef=1, debug=False):
    """C_3 = coef * Sum d^ABC M_A M_B M_C.

    Called with ``d_upper = reference_cubic_tensor()`` and ``coef = 8`` this is
    the universal cubic Casimir, one fixed element of U(g), with eigenvalue

        C_3 = lambda (mu^2 - lambda^2) = lambda . C_2      [monic]

    linear in lambda as an element of U(g) must be.  The 8 sets the monic
    scale; it absorbs the internal 1/6 below, so changing that convention moves
    the 8 and leaves C_3 = lambda . C_2 untouched.

    Passing a module's own raised d-tensor instead gives a per-module operator
    carrying A(R)/I(R)^3 = -lambda/mu^2, i.e. eigenvalue
    -lambda^2 (mu^2 - lambda^2) / (8 mu^2) -- quadratic in lambda, and not
    comparable across modules.  That route is a diagnostic only.

    If debug=True, prints each term's contribution to C_3[0,0] and C_3[3,3].
    """
    nums = rep.numbers()
    d = rep.dim
    titles = ['Y', 'e', 'f', 'h', 'u', 'v', 'w', 'x']
    
    # Initialize C_3 to zero matrix
    C3 = Matrix.zero(d, d)
    
    if debug:
        print("CASIMIR C_3 CUMULATIVE BUILD (debug=True):")
        print("FORMULA: C_3 = Σ d^ABC M_A M_B M_C (with sign-corrected d^ABC)")
        print("=" * 160)
        print(f"{'Triple':<12} {'d^ABC':<15} {'Prod[0,0]':<18} {'Prod[3,3]':<18} {'C3[0,0] after':<20} {'C3[3,3] after':<20} {'Diff [0,0]-[3,3]':<20}")
        print("=" * 160)
    
    # Sum over all d^ABC terms
    for A in nums:
        for B in nums:
            for C in nums:
                coeff = d_upper.get((A, B, C), 0)/6
                
                # Skip zero terms
                if coeff == 0:
                    continue
                
                # Compute M_A @ M_B @ M_C
                M_A = rep[A]
                M_B = rep[B]
                M_C = rep[C]
                product = M_C @ M_B @ M_A
                
                # Add this term to C_3
                term = coeff * product 
                C3 = C3 + term
                
                if debug:
                    triple_label = f"{titles[A]}{titles[B]}{titles[C]}"
                    prod_00 = product[0, 0]
                    prod_33 = product[3, 3]
                    c3_00 = C3[0, 0]
                    c3_33 = C3[3, 3]
                    diff = c3_00 - c3_33
                    print(f"{triple_label:<12} {str(coeff):<15} {str(prod_00):<18} {str(prod_33):<18} {str(c3_00):<20} {str(c3_33):<20} {str(diff):<20}")
    
    if debug:
        print("=" * 160)
        print()
    
    return (coef * C3).applyfunc(simplify)


def casimir_commutes(rep, C2, nam, is_zero=is_zero_scalar):
    """Check that the quadratic Casimir C2 commutes with every generator.

    Being a Casimir means exactly this: [C2, M_k] = 0 for all k. C2 is even,
    so its super-bracket with any generator -- even OR odd -- is the ordinary
    commutator C2 @ M_k - M_k @ C2 (the graded sign (-1)^(|C2||M_k|) is +1
    because |C2| = 0). So no anticommutator ever arises here.

    For N = 1, C2 is a scalar multiple of the identity (Schur on the irrep),
    so commutation is automatic; the check is trivially satisfied. For the
    Matryoshka N > 1 the module is indecomposable, C2 acquires an off-diagonal
    (nilpotent) part and is no longer scalar -- yet a true Casimir must still
    commute with every generator. That is the substantive content here, and the
    proper N > 1 replacement for the ``is_scalar_multiple`` test.

    Returns a list of CasimirCheck records -- data only, never printed.
    Zero-testing goes through is_zero_scalar (invariant 4), same policy as
    check() and chi_relations()."""
    def is_zero_matrix(M):
        return all(is_zero(M[i, j])
                   for i in range(M.rows) for j in range(M.cols))

    out = []
    for k in rep.numbers():
        M = rep[k]
        br = C2 @ M - M @ C2                    # [C2, M_k]  (C2 even -> commutator)
        out.append(CasimirCheck(f"[{nam}, {rep.title(k)}] = 0",
                                is_zero_matrix(br), br))
    return out


def export_metrics(rep, chi, N, g_lower=None, g_upper=None):
    """Export the Killing metrics (lower and upper index) as a structured dict.
    
    If g_lower/g_upper not provided, they are computed. Returns a dict with
    the sector structure for easy access:
    {
        'lower': g_lower_matrix,
        'upper': g_upper_matrix,
        'sectors': {
            'Y-Y': (g_lower[0,0], g_upper[0,0]),
            'E-F': (2x2 blocks for lower and upper),
            'H-H': (g_lower[3,3], g_upper[3,3]),
            'U-V': (2x2 blocks for lower and upper),
            'W-X': (2x2 blocks for lower and upper),
        }
    }
    """
    if g_lower is None:
        g_lower = killing_metric(rep, chi, N)
    if g_upper is None:
        g_upper = upper_killing_metric(g_lower, rep)
    
    sectors = {
        'Y-Y': (g_lower[0,0], g_upper[0,0]),
        'E-F': {
            'lower': [[g_lower[1,1], g_lower[1,2]], [g_lower[2,1], g_lower[2,2]]],
            'upper': [[g_upper[1,1], g_upper[1,2]], [g_upper[2,1], g_upper[2,2]]]
        },
        'H-H': (g_lower[3,3], g_upper[3,3]),
        'U-V': {
            'lower': [[g_lower[4,4], g_lower[4,5]], [g_lower[5,4], g_lower[5,5]]],
            'upper': [[g_upper[4,4], g_upper[4,5]], [g_upper[5,4], g_upper[5,5]]]
        },
        'W-X': {
            'lower': [[g_lower[6,6], g_lower[6,7]], [g_lower[7,6], g_lower[7,7]]],
            'upper': [[g_upper[6,6], g_upper[6,7]], [g_upper[7,6], g_upper[7,7]]]
        }
    }
    
    return {
        'lower': g_lower,
        'upper': g_upper,
        'sectors': sectors
    }


def upper_killing_metric(g_lower, rep):
    """Upper-index Killing metric g^ij (inverse of g_ij), computed sector-by-sector.

    The canonical numbers are Y(0), E(1), F(2), H(3), U(4), V(5), W(6), X(7).
    The structure is:
        - Y-Y block:   1x1 diagonal  (even)
        - E-F block:   2x2 symmetric (even)
        - H-H block:   1x1 diagonal  (even)
        - U-V block:   2x2 antisymmetric (odd)
        - W-X block:   2x2 antisymmetric (odd)

    For antisymmetric 2x2 blocks (u,v) with g_uv = -g_vu:
        g^uv = -1/g_uv,   g^vu = 1/g_uv
    satisfies g^ui * g_iv = delta^u_i."""
    n = 8  # exactly 8 generators in sl(2|1)

    # Initialize the upper-index metric
    def g_upper_entry(p, q):
        if p == 0 and q == 0:  # Y-Y (1D diagonal)
            g_yy = g_lower[0, 0]
            return 1 / g_yy if g_yy != 0 else 0

        elif (p, q) in [(1, 1), (1, 2), (2, 1), (2, 2)]:  # E-F block (2x2 symmetric)
            g_ee = g_lower[1, 1]
            g_ef = g_lower[1, 2]
            g_ff = g_lower[2, 2]
            det = simplify(g_ee * g_ff - g_ef * g_ef)
            if det == 0:
                return 0
            inv_det = 1 / det
            if p == 1 and q == 1:
                return simplify(inv_det * g_ff)
            elif (p, q) in [(1, 2), (2, 1)]:
                return simplify(-inv_det * g_ef)
            else:  # (2, 2)
                return simplify(inv_det * g_ee)

        elif p == 3 and q == 3:  # H-H (1D diagonal)
            g_hh = g_lower[3, 3]
            return 1 / g_hh if g_hh != 0 else 0

        elif (p, q) in [(4, 4), (4, 5), (5, 4), (5, 5)]:  # U-V block (2x2 antisymmetric)
            g_uv = g_lower[4, 5]  # the only nonzero off-diagonal
            if g_uv == 0:
                return 0
            inv = 1 / g_uv
            if p == 4 and q == 5:
                return simplify(-inv)  # g^uv = -1/g_uv
            elif p == 5 and q == 4:
                return simplify(inv)   # g^vu = 1/g_uv
            else:
                return 0  # diagonals are zero

        elif (p, q) in [(6, 6), (6, 7), (7, 6), (7, 7)]:  # W-X block (2x2 antisymmetric)
            g_wx = g_lower[6, 7]  # the only nonzero off-diagonal
            if g_wx == 0:
                return 0
            inv = 1 / g_wx
            if p == 6 and q == 7:
                return simplify(-inv)  # g^wx = -1/g_wx
            elif p == 7 and q == 6:
                return simplify(inv)   # g^xw = 1/g_wx
            else:
                return 0  # diagonals are zero

        else:  # off-diagonal cross-sector blocks are zero
            return 0

    return Matrix(n, n, g_upper_entry)


# -- the Gorelik anticenter element T_4 (built from N_2 and N_4) -------------

# the four odd generators, in canonical order
ODD_GENS = [U, V, W, X]


def _perm_sign(perm):
    """Sign (+1/-1) of a permutation given as a sequence, by inversion count
    against ascending order (so (U,V,W,X) itself is +1)."""
    perm = list(perm)
    n = len(perm)
    s = 1
    for i in range(n):
        for j in range(i + 1, n):
            if perm[i] > perm[j]:
                s = -s
    return s


def anticenter_N2(rep):
    """N_2 = UV - VU + WX - XW  (ordinary matrix products of the odd gens).

    U,V,W,X are u,v,w,x (canonical 4,5,6,7). Returns a Matrix (simplified)."""
    u, v, w, x = rep[U], rep[V], rep[W], rep[X]
    N2 = u @ v - v @ u + w @ x - x @ w
    return N2.applyfunc(simplify)


def anticenter_N4(rep):
    """N_4 = fully antisymmetric sum over the four odd generators:

        N_4 = Σ_{σ ∈ S_4} sign(σ) · M_{σ(1)} M_{σ(2)} M_{σ(3)} M_{σ(4)}

    where σ permutes (u, v, w, x): 24 signed quartic matrix products. Returns a
    Matrix (simplified)."""
    d = rep.dim
    N4 = Matrix.zero(d, d)
    for perm in permutations(ODD_GENS):
        s = _perm_sign(perm)
        prod = rep[perm[0]] @ rep[perm[1]] @ rep[perm[2]] @ rep[perm[3]]
        N4 = N4 + s * prod
    return N4.applyfunc(simplify)


def anticenter_T(rep):
    """T = (2/3) (N_4 - N_2)  --  the MONIC normalisation.

    Gorelik's theorem identifies the Harish-Chandra image of the anticenter
    with the atypicality polynomial

        prod_{alpha in Delta_1^+} <Lambda + rho | alpha>  =  mu^2 - lambda^2

    taken monic.  N_4 - N_2 has eigenvalue (3/2)(mu^2 - lambda^2), so the 2/3
    buys the leading coefficient 1.  That scale is fixed by the root data
    alone: no module and no choice of invariant form enter it.

    Paired with the C_2 of ``casimir_quadratic_direct`` this gives the operator
    identity

        T = C_2 . chi

    verified on the full matrices by ``anticenter_casimir_relations``."""
    N2 = anticenter_N2(rep)
    N4 = anticenter_N4(rep)
    return (Rational(2, 3) * (N4 - N2)).applyfunc(simplify)


# -- the quadratic Casimir, built directly (no metric at all) ----------------


def casimir_quadratic_direct(rep):
    """C_2 = h^2 - Y^2 + 2(ef + fe) + 2 N_2,  the universal quadratic Casimir.

    An explicit polynomial in the generators, using no metric and no inverse:
    a Casimir is an element of U(g), and this writes it down as one.  Built the
    same way T is.

    Eigenvalue on the Kac module R(a,b):

        C_2 = mu^2 - lambda^2 = -(y+a)(y-a-2),    mu = a+1,  lambda = y-1

    monic in the atypicality polynomial, matching the normalisation of T.
    For comparison with the metric section, (1/2) g^{AB}(R(0,0)) M_A M_B is
    C_2 / 4, and a module's own form adds a further 1/I(R) = 1/mu.

    The split reported below is by sector, as in ``casimir_quadratic``:

        even sector:  h^2 - Y^2 + 2(ef + fe)
        odd  sector:  2 N_2 = 2(uv - vu + wx - xw)

    Returns the same dict shape as ``casimir_quadratic``."""
    Ym, e, f, h = rep[Y], rep[E], rep[F], rep[H]

    C2_even = h @ h - Ym @ Ym + 2 * (e @ f + f @ e)
    C2_odd = 2 * anticenter_N2(rep)

    return _c2_result(C2_even, C2_odd)


# One reported T relation: its label, whether it holds, and the leftover.
TCheck = namedtuple("TCheck", "label ok residual")


def anticenter_relations(rep, T, is_zero=is_zero_scalar):
    """T must commute with every even generator and anticommute with every odd
    one -- checked NUMERICALLY on the full matrices, never theorized.

    For each generator g:
        even g:  [T, g] = T g - g T   must be the zero matrix
        odd  g:  {T, g} = T g + g T   must be the zero matrix

    This is the odd-element behaviour (same shape as the chi relations): T's
    super-bracket with g is a commutator when g is even and an anticommutator
    when g is odd. Zero-testing goes through is_zero_scalar (invariant 4).
    Returns a list of TCheck records -- data only, never printed."""
    def is_zero_matrix(M):
        return all(is_zero(M[i, j])
                   for i in range(M.rows) for j in range(M.cols))

    out = []
    for k in rep.numbers():
        M = rep[k]
        if rep.parity(k) == ODD:
            br = T @ M + M @ T                      # {T, odd}
            op, cl = "{", "}"
        else:
            br = T @ M - M @ T                      # [T, even]
            op, cl = "[", "]"
        out.append(TCheck(f"{op}T, {rep.title(k)}{cl} = 0",
                          is_zero_matrix(br), br))
    return out


def anticenter_casimir_relations(rep, T, C2, chi, is_zero=is_zero_scalar):
    """The monic anticenter identity, checked on the FULL matrices:

        [C_2, chi] = 0          C_2 is even, chi is even as a matrix
        T = C_2 . chi           the identity itself
        T^2 = C_2^2             central -- Gorelik's theorem, explicitly

    This is an OPERATOR identity, not an eigenvalue relation, and the
    distinction is the whole reason to check it numerically.  On the irreps
    (N = 1) C_2 is scalar and 'T = gamma.chi with a constant gamma' would also
    fit the data; on the Matryoshka N > 1 it does not, because C_2 there has a
    nilpotent off-diagonal part -- and the identity still holds.  Nothing here
    is assumed from theory.

    T^2 = C_2^2 follows from the first two lines (chi^2 = 1 and C_2 commutes
    with chi), but it is checked independently rather than deduced: it is the
    statement that the anticenter squares into the centre, and it should fail
    loudly if the normalisation of either factor drifts.

    A rank coincidence, worth stating wherever this identity is used: it holds
    because sl(2|1) has exactly TWO positive odd roots, so the atypicality
    polynomial is quadratic -- the same degree as C_2.  The monic scale of T
    (Gorelik) generalises to sl(m|n); T = C_2 . chi does not.

    Returns a list of TCheck records -- data only, never printed."""
    def is_zero_matrix(M):
        return all(is_zero(M[i, j])
                   for i in range(M.rows) for j in range(M.cols))

    r_comm = C2 @ chi - chi @ C2
    r_ident = T - C2 @ chi
    r_square = T @ T - C2 @ C2

    return [
        TCheck("[C_2, chi] = 0", is_zero_matrix(r_comm), r_comm),
        TCheck("T = C_2 . chi", is_zero_matrix(r_ident), r_ident),
        TCheck("T^2 = C_2^2", is_zero_matrix(r_square), r_square),
    ]


# -- 3. verification ---------------------------------------------------------

def verify(rep, f=f, is_zero=is_zero_scalar):
    """Nonzero residuals of the structure-constant check (empty = all hold)."""
    return rep.check(f, is_zero=is_zero)


# -- small local display (report.py will own this properly later) ------------

def _show_matrix(title, M):
    print(f"{title} =")
    print(M)
    print()


def _combo(rep, terms):
    parts = []
    for k, c in terms.items():
        name = rep.title(k)
        if c == 1:
            parts.append(name)
        elif c == -1:
            parts.append(f"-{name}")
        else:
            parts.append(f"{c} {name}")
    return " + ".join(parts).replace("+ -", "- ")


def _bracket(rep, i, j, terms):
    op, cl = ("{", "}") if (PARITY[i] == ODD and PARITY[j] == ODD) else ("[", "]")
    return f"{op}{rep.title(i)}, {rep.title(j)}{cl} = {_combo(rep, terms)}"


def _null_pairs(rep):
    """Every canonical pair i<=j whose (super-)bracket must vanish, i.e. every
    pair the default check() sweeps that is NOT in INDEPENDENT. These are the
    relations with no right-hand side: Y central, the null mixed brackets, the
    odd squares (via the diagonal anticommutator {a,a}=2a^2), and the vanishing
    odd-odd pairs. For the superalgebra these are genuine constraints with no
    Lie-algebra analogue, so they are listed explicitly, not left implicit."""
    nums = rep.numbers()
    indep = set(INDEPENDENT)
    out = []
    for p in range(len(nums)):
        for q in range(p, len(nums)):
            i, j = nums[p], nums[q]
            if (i, j) not in indep:
                out.append((i, j))
    return out


def _null_label(rep, i, j):
    """Render a vanishing bracket, e.g. '[Y, e] = 0' or '{u, u} = 2 u^2 = 0'."""
    both_odd = PARITY[i] == ODD and PARITY[j] == ODD
    op, cl = ("{", "}") if both_odd else ("[", "]")
    ti, tj = rep.title(i), rep.title(j)
    lhs = f"{op}{ti}, {tj}{cl}"
    if i == j and both_odd:                       # {a,a} = 2 a^2
        return f"{lhs} = 2 {ti}^2 = 0"
    return f"{lhs} = 0"


def _casimir_report(rep, chi, N, a, y):
    """Killing metric + quadratic/cubic Casimir constructors and their
    commutation checks. Only run when --casimirs is given."""
    # -- the lower-index Killing metric --------------------------------------
    print()
    g = killing_metric(rep, chi, N)
    titles = [rep.title(k) for k in rep.numbers()]

    def _row(r, cols):
        return "(" + " ".join(f"{g[r, c]!s:>4}" for c in cols) + ")"

    even_c, odd_c = range(0, 4), range(4, 8)
    print(f"Lower-index Killing metric   "
          f"g_ab=(1/2N)STr(ab+ba), g_ij=(1/2N)STr(ij-ji),  N={N}")
    print(f"  even [{' '.join(titles[:4])}]:  "
          + "  ".join(_row(r, even_c) for r in range(4)))
    print(f"  odd  [{' '.join(titles[4:])}]:  "
          + "  ".join(_row(r, odd_c) for r in range(4, 8)))

    # -- the upper-index Killing metric (inverse) ----------------------------
    print()
    g_up = upper_killing_metric(g, rep)
    print(f"Upper-index Killing metric (inverse)   g^ij = (g_lower)^-1_ij")
    
    def _row_up(r, cols):
        return "(" + " ".join(f"{g_up[r, c]!s:>4}" for c in cols) + ")"
    
    print(f"  even [{' '.join(titles[:4])}]:  "
          + "  ".join(_row_up(r, even_c) for r in range(4)))
    print(f"  odd  [{' '.join(titles[4:])}]:  "
          + "  ".join(_row_up(r, odd_c) for r in range(4, 8)))

    # Test inverse property: g^ij g_jk = delta^i_k
    print()
    print("Testing metric inverse property: g^ij g_jk = delta^i_k")
    failures = verify_metric_inverse(g, g_up)
    if not failures:
        print("  All components verified. \u2713")
    else:
        print(f"  {len(failures)} component(s) failed:")
        for i, k, got, expected in failures:
            print(f"    g^{titles[i]}{titles[k]} * g_{k}? should be "
                  f"{expected}, got {got}")

    # -- the quadratic Casimir operator C_2 --------------------------------
    print()
    print("Quadratic Casimir operator (direct, no metric):")
    print("    C_2 = h^2 - Y^2 + 2(ef + fe) + 2 N_2        eigenvalue mu^2 - lambda^2")
    c2_result = casimir_quadratic_direct(rep)

    # For comparison with the metric section: this module's own form carries
    # its Dynkin index, so contracting with it gives C_2 / (4 mu).
    c2_contracted = casimir_quadratic(rep, g_up)['eigenvalue']
    if c2_contracted is not None:
        print(f"    (1/2) g^AB M_A M_B  with this module's own form = "
              f"{c2_contracted}   [= C_2 / (4(a+1)); not used]")
    
    print()
    print("C_2 even sector (generators Y, e, f, h only):")
    print(c2_result['even_sector'])
    
    print()
    print("C_2 odd sector (generators u, v, w, x only):")
    print(c2_result['odd_sector'])
    
    print()
    print("C_2 total = even + odd:")
    print(c2_result['total'])
    
    print()
    if c2_result['is_scalar_multiple']:
        print(f"✓ C_2 is a scalar multiple of identity with eigenvalue: "
              f"{c2_result['eigenvalue']}")
    else:
        print("✗ C_2 is NOT a scalar multiple of identity")
        # Show a few nonzero off-diagonal elements. For the Matryoshka these
        # live in the [i, i+d] block band (d = 4(a+1)), not the top-left
        # corner, so scan the whole matrix rather than a fixed 4x4 window.
        C2 = c2_result['total']
        found = []
        for i in range(C2.rows):
            for j in range(C2.cols):
                if i != j and simplify(C2[i, j]) != 0:
                    found.append((i, j, C2[i, j]))
                    break                          # one per row is enough
            if len(found) >= 3:
                break
        if found:
            print("  Some nonzero off-diagonal elements:")
            for i, j, val in found:
                print(f"    C_2[{i},{j}] = {val}")

    # -- C_2 commutes with every generator (the defining Casimir property) ---
    # Being scalar is sufficient but not necessary: for N > 1 the Matryoshka
    # C_2 is not proportional to the identity (it is not block-diagonal), yet a
    # genuine Casimir must still commute with all generators. This is the real
    # test in that regime.
    print()
    print("Testing that C_2 commutes with every generator: [C_2, M_k] = 0")
    casimir_checks = casimir_commutes(rep, c2_result['total'],"C2")
    for c in casimir_checks:
        mark = "\u2713" if c.ok else "FAILED"
        print(f"    {c.label:<16} {mark}")
    print()
    if all(c.ok for c in casimir_checks):
        print("C_2 commutes with all generators: it is a Casimir operator. \u2713")
    else:
        for c in casimir_checks:
            if not c.ok:
                print(f"    {c.label} residual =")
                print(c.residual)
                print()

    # -- the master equation:  d_ABC = (constant tensor) . Tr(Y) -------------
    print()
    print("=" * 80)
    print("Cubic tensor of THIS module, and its proportionality to Tr(Y)")
    print("=" * 80)
    print()

    d_lower = cubic_d(rep, chi, N)
    tr_Y = simplify(trace(rep[Y]) / N)
    tr_Y_formula = 4 * (a + 1) * (y - 1)
    mark_tr = "\u2713" if is_zero_scalar(tr_Y - tr_Y_formula) else "FAILED"
    print(f"Tr(Y) per layer = {tr_Y}      formula 4(a+1)(y-1) = "
          f"{tr_Y_formula}   {mark_tr}")
    ratio, d_failures = cubic_d_proportionality(d_lower, tr_Y)
    print(f"d_ABC(R) = (constant tensor) . Tr(Y),  constant tensor = "
          f"-d_ABC(R(0,0))/4")
    print(f"  proportionality constant  Tr(Y)/Tr(Y)_ref = {ratio}  "
          f"[= A(R) = -(a+1)(y-1)]")
    if not d_failures:
        print(f"  all {len(d_lower)} components proportional. \u2713")
    else:
        print(f"  {len(d_failures)} of {len(d_lower)} component(s) FAILED:")
        for k, got, expected in d_failures[:10]:
            name = "".join(titles[i] for i in k)
            print(f"    d_{name}: got {got}, expected {expected}")
    print()

    # -- the FIXED cubic Casimir ---------------------------------------------
    print()
    print("=" * 80)
    print("Cubic Casimir operator (FIXED):  C_3 = 8 . Σ d^ABC(ref) M_A M_B M_C")
    print("    d^ABC(ref) comes from the y = 0 anchor R(0,0), normalised by")
    print("    (a_ref+1)^2 so that any y = 0 anchor gives the same tensor.")
    print("    One fixed element of U(g), evaluated here -- the module under")
    print("    study never has its own form inverted.")
    print("    Eigenvalue: lambda (mu^2 - lambda^2) = lambda . C_2   [monic]")
    print("=" * 80)
    print()

    anchor_failures = reference_anchor_check((0, 1))
    if not anchor_failures:
        print("  anchors R(0,0) and R(1,1/2) give the identical fixed tensor "
              "(all 512 components). \u2713")
    else:
        print(f"  {len(anchor_failures)} component(s) DIFFER between anchors:")
        for anchor_a, k, got, expected in anchor_failures[:10]:
            name = "".join(titles[i] for i in k)
            print(f"    anchor a={anchor_a}  d^{name}: {got} vs {expected}")
    print()

    d_upper = reference_cubic_tensor()
    C3 = casimir_cubic(rep, d_upper, coef=8, debug=False)
    
    print("C_3 total =")
    print(C3)
    print()
    
    # Check if C_3 is diagonal
    if (0):
        C3_diag_00 = simplify(C3[0, 0])
        C3_diag_33 = simplify(C3[3, 3]) if C3.rows > 3 else 0
        diff = simplify(C3_diag_00 - C3_diag_33)
        
        print(f"C_3[0,0] = {C3_diag_00}")
        print(f"C_3[3,3] = {C3_diag_33}")
        print(f"Difference = {diff}")
        print()
    
    # Test that C_3 commutes with all generators
    print("Testing that C_3 commutes with every generator: [C_3, M_k] = 0")
    c3_checks = casimir_commutes(rep, C3,"C3")
    for c in c3_checks:
        mark = "\u2713" if c.ok else "FAILED"
        print(f"    {c.label:<16} {mark}")
    print()
    if all(c.ok for c in c3_checks):
        print("C_3 commutes with all generators: it is a Casimir operator. \u2713")
    else:
        for c in c3_checks:
            if not c.ok:
                print(f"    {c.label} residual =")
                print(c.residual)
                print()

    # Hand the computed objects back so main() can reuse them: the anticenter
    # identity needs C_2, and the export/verification block would otherwise
    # rebuild the metric and the d-tensor from scratch a second time.
    return {
        'g_lower': g, 'g_upper': g_up,
        'd_lower': d_lower, 'd_upper': d_upper,
        'c2': c2_result, 'C3': C3,
    }



def main(a, N=1, casimirs=False, zetaH=False, case=None):
    if zetaH:
        base, dy, is0 = Rzeta(a), dy_zeta, zeta_is_zero(a)
        b = AL**2
    else:
        base, dy, is0 = Rsl21(a), dy_chevalley, is_zero_scalar
        b = B
    rep = matryoshka(base, N, dy)

    y = 2 * b - a
    basis = ("zeta-Hermitian basis" if zetaH else "Chevalley basis")
    if N == 1:
        print(f"sl(2|1) Kac module R(a={a}, b={b})  [y={y}]   "
              f"(dimension {rep.dim} = 4(a+1))   b is FORMAL, {basis}")
    else:
        print(f"sl(2|1) Matryoshka indecomposable rep MR(a={a}, b={b}, N={N})  "
              f"[y={y}]   (dimension {rep.dim} = {N}*4(a+1))   "
              f"b is FORMAL, {basis}")
    if zetaH:
        print(f"case {case}:  {CASE_REGION[case]}      "
              f"alpha = sqrt(b) = sqrt(K_1),  gamma = sqrt(b-a-1) = sqrt(K_2)")
    print("=" * 46, "\n")

    # Which Taylor jets are alive.  In the Chevalley basis every generator is
    # affine in y, so the answer is [0, 1] and the Matryoshka has two bands.
    # After the zeta gauge the entries are square roots and every order is
    # alive -- which is what makes N > 2 a real test rather than a re-run.
    if N > 1:
        alive = jet_orders(base, N, dy, is0)
        missing = [p for p in range(N) if p not in alive]
        print(f"Taylor jets (1/p!) (d/dy)^p mu on sub-diagonal p, "
              f"p = 0..{N - 1}:")
        note = ("   [every order alive: a genuine jet]" if not missing else
                "   [orders " + ",".join(map(str, missing))
                + " vanish: generators affine in y]")
        print(f"    nonzero orders: {alive}{note}")
        print()

    for num in rep.numbers():
        _show_matrix(f"{rep.title(num)}  (#{num})", rep[num])

    print("Nonzero (super-)brackets to verify:")
    for (i, j), terms in INDEPENDENT.items():
        print("    " + _bracket(rep, i, j, terms))
    print()

    nulls = _null_pairs(rep)
    print("Null (super-)brackets to verify (must vanish):")
    for (i, j) in nulls:
        print("    " + _null_label(rep, i, j))
    print()

    n_pairs = len(INDEPENDENT) + len(nulls)
    print(f"Coverage: all {n_pairs} canonical pairs i<=j checked "
          f"({len(INDEPENDENT)} nonzero + {len(nulls)} null); "
          f"reversed pairs follow by graded antisymmetry.")
    print()

    residuals = verify(rep, is_zero=is0)
    if not residuals:
        print("All sl(2|1) relations verified (nonzero and null). \u2713")
    else:
        print(f"{len(residuals)} relation(s) FAILED:\n")
        for r in residuals:
            op, cl = ("{", "}") if (PARITY[r.i] == ODD and PARITY[r.j] == ODD) \
                else ("[", "]")
            print(f"    {op}{r.title_i}, {r.title_j}{cl} residual =")
            print(r.matrix)
            print()

    # -- the grading operator chi --------------------------------------------
    print()
    chi = Rchi(a, N)
    _show_matrix("chi  (grading operator: -1,+1,+1,-1 per layer)", chi)

    print("Grading-operator relations to verify "
          "([chi,even]=0, {chi,odd}=0, chi^2=1):")
    chi_checks = chi_relations(rep, chi, is0)
    for c in chi_checks:
        mark = "\u2713" if c.ok else "FAILED"
        print(f"    {c.label:<16} {mark}")
    print()

    if all(c.ok for c in chi_checks):
        print("chi commutes with the even generators and anticommutes with the "
              "odd ones, and chi^2 = 1. \u2713")
    else:
        for c in chi_checks:
            if not c.ok:
                print(f"    {c.label} residual =")
                print(c.residual)
                print()

    # -- the reality condition ------------------------------------------------
    #
    # This is the a-posteriori proof.  The eight matrices were produced above by
    # whatever means; nothing about their derivation enters here.  What is
    # checked is what they satisfy, by explicit multiplication, entry by entry.
    if zetaH:
        print()
        Z = zeta_matrix(a, case, N)
        label = ("zeta  (per layer: 1, sign K_1, sign K_2, sign K_1 sign K_2)"
                 if N == 1 else
                 "zeta_N = J (x) zeta  (J reverses the N generations)")
        _show_matrix(label, Z)

        print("Reality conditions to verify "
              "(zeta M^dag zeta = partner, zeta^2 = 1):")
        z_checks = zeta_relations(rep, Z, case, is0)
        for c in z_checks:
            print(f"    {c.label:<34} {'✓' if c.ok else 'FAILED'}")
        print()
        if all(c.ok for c in z_checks):
            print(f"All eight matrices are zeta-Hermitian in case {case}. ✓")
            if N > 1:
                print("  (with zeta_N = J (x) zeta: a block-diagonal zeta "
                      "cannot work, the Matryoshka being block triangular)")
        else:
            for c in z_checks:
                if not c.ok:
                    print(f"    {c.label} residual =")
                    print(c.residual)
                    print()

        if N > 1:
            print()
            print("Each Taylor coefficient separately zeta-Hermitian "
                  "(the hypothesis zeta_N = J (x) zeta rests on):")
            j_checks = zeta_jet_relations(base, a, case, N, is0)
            for c in j_checks:
                print(f"    {c.label:<44} {'✓' if c.ok else 'FAILED'}")
            print()
            if not all(c.ok for c in j_checks):
                for c in j_checks:
                    if not c.ok:
                        print(f"    {c.label} sample residual =")
                        print(c.residual)
                        print()

    if casimirs:
        # The Casimir section predates the zeta ground and uses bare
        # ``simplify`` throughout, so it runs on a copy with GA eliminated in
        # favour of AL, where ``is_zero_scalar`` is already correct.  C_2, C_3
        # and T are polynomials in the generators; this changes nothing about
        # what is computed.
        crep = zeta_reduce(a)(rep) if zetaH else rep
        cas = _casimir_report(crep, chi, N, a, y)

    if casimirs:
        # -- the Gorelik anticenter element T_4: N_2 and N_4 ---------------------
        print()
        print("=" * 80)
        print("Gorelik anticenter T_4:  N_2 then N_4  (compare by eye)")
        print("=" * 80)
        print()

        N2 = anticenter_N2(crep)
        print("N_2 = UV - VU + WX - XW =")
        print(N2)
        print()

        N4 = anticenter_N4(crep)
        print("N_4 = Σ_σ sign(σ) (u v w x permuted), fully antisymmetric =")
        print(N4)
        print()

        T = (Rational(2, 3) * (N4 - N2)).applyfunc(simplify)
        print("T = (2/3)(N_4 - N_2)   [monic normalisation] =")
        print(T)
        print()

        print("Anticenter relations to verify numerically "
              "([T,even]=0, {T,odd}=0):")
        t_checks = anticenter_relations(crep, T)
        for c in t_checks:
            mark = "\u2713" if c.ok else "FAILED"
            print(f"    {c.label:<16} {mark}")
        print()
        if all(c.ok for c in t_checks):
            print("T commutes with all even generators and anticommutes with all "
                  "odd ones. \u2713")
        else:
            for c in t_checks:
                if not c.ok:
                    print(f"    {c.label} residual =")
                    print(c.residual)
                    print()

        # -- the monic identity T = C_2 . chi, on the full matrices ----------
        print()
        print("Anticenter/Casimir identity to verify numerically "
              "(operator identity, not eigenvalues):")
        tc_checks = anticenter_casimir_relations(crep, T, cas['c2']['total'], chi)
        for c in tc_checks:
            mark = "\u2713" if c.ok else "FAILED"
            print(f"    {c.label:<16} {mark}")
        print()
        if all(c.ok for c in tc_checks):
            print("T = C_2 . chi and T^2 = C_2^2 hold as operator identities. "
                  "\u2713")
            if N > 1:
                print("  (checked where C_2 is NOT scalar, so this is stronger "
                      "than an eigenvalue relation)")
        else:
            for c in tc_checks:
                if not c.ok:
                    print(f"    {c.label} residual =")
                    print(c.residual)
                    print()

        # -- Export data for analytic formula discovery -----
        print()
        print("=" * 80)
        print("EXPORT: Non-zero metric elements, d_YYY, STr(N_4)")
        print("=" * 80)
        print()

        # built once in _casimir_report and reused here
        g_lower, g_upper = cas['g_lower'], cas['g_upper']
        d_lower, d_upper = cas['d_lower'], cas['d_upper']

        # Non-zero metric values
        print(f"Non-zero metric values (lower index):")
        metric_entries = []
        for i in range(8):
            for j in range(i, 8):
                val = g_lower[i, j]
                if val != 0:
                    metric_entries.append((i, j, val))
                    titles = ['Y', 'e', 'f', 'h', 'u', 'v', 'w', 'x']
                    print(f"  g_{{{titles[i]}{titles[j]}}} = {val}")
        if not metric_entries:
            print("  (all zero)")
        print()

        print(f"Non-zero metric values (upper index):")
        metric_entries_up = []
        for i in range(8):
            for j in range(i, 8):
                val = g_upper[i, j]
                if val != 0:
                    metric_entries_up.append((i, j, val))
                    titles = ['Y', 'e', 'f', 'h', 'u', 'v', 'w', 'x']
                    print(f"  g^{{{titles[i]}{titles[j]}}} = {val}")
        if not metric_entries_up:
            print("  (all zero)")
        print()

        # Analytic formulas with automatic verification
        print()
        print("=" * 80)
        print("ANALYTIC FORMULAS & VERIFICATION")
        print("=" * 80)
        print()

        y_val = 2 * b - a
        
        # d_YYY = -6(a+1)(y-1)
        d_yyy_actual = simplify(d_lower[(0, 0, 0)])
        d_yyy_formula = -6 * (a + 1) * (y_val - 1)
        d_yyy_ratio = (simplify(d_yyy_actual / d_yyy_formula) if d_yyy_formula != 0
                       else ("0" if d_yyy_actual == 0 else "MISMATCH"))
        mark_d = "✓" if d_yyy_ratio in (1, "0") else "✗"
        print(f"d_YYY = {d_yyy_actual}")
        print(f"  Formula:  -6(a+1)(y-1) = {d_yyy_formula}")
        print(f"  Ratio: {d_yyy_ratio}  {mark_d}"
              + ("   [both sides vanish]" if d_yyy_formula == 0 else ""))
        print()
        
        # C_2 = mu^2 - lambda^2 = -(y+a)(y-a-2)   [monic; direct construction]
        c2_result = cas['c2']
        # On the Matryoshka (N > 1) C_2 is not scalar, so there is no single
        # eigenvalue; the diagonal blocks are still R(a,b), so the [0,0] entry
        # is the value the formula predicts and the nilpotent part sits off the
        # diagonal.  The operator identity above is what covers N > 1 properly.
        c2_actual = (c2_result['eigenvalue'] if c2_result['is_scalar_multiple']
                     else simplify(c2_result['total'][0, 0]))
        c2_note = "" if c2_result['is_scalar_multiple'] else "   [[0,0] entry; C_2 not scalar]"
        c2_formula = -(y_val + a) * (y_val - a - 2)
        c2_ratio = (simplify(c2_actual / c2_formula) if c2_formula != 0
                    else ("0" if c2_actual == 0 else "MISMATCH"))
        mark_c2 = "✓" if c2_ratio in (1, "0") else "✗"
        print(f"C_2 = {c2_actual}{c2_note}")
        print(f"  Formula:  mu² - lambda² = -(y+a)(y-a-2) = {c2_formula}")
        print(f"  Ratio: {c2_ratio}  {mark_c2}"
              + ("   [atypical: both sides vanish]" if c2_formula == 0 else ""))
        print()
        
        # C_3 = -1/(2(a+1)) * C_2 * (y-1)²
        C3 = cas['C3']              # already rescaled by -8(a+1)^2
        c3_actual = simplify(C3[0, 0])
        # C_3 = lambda (mu^2 - lambda^2) = lambda . C_2, the universal
        # operator built from the y = 0 anchor.  Linear in lambda, as an
        # element of U(g) must be.
        lam = y_val - 1
        c3_formula = lam * (-(y_val + a) * (y_val - a - 2))
        c3_ratio = (simplify(c3_actual / c3_formula) if c3_formula != 0
                    else ("0" if c3_actual == 0 else "MISMATCH"))
        mark_c3 = "✓" if c3_ratio in (1, "0") else "✗"
        print(f"C_3 = {c3_actual}")
        print(f"  Formula:  lambda(mu² - lambda²) = lambda·C_2 = {c3_formula}")
        print(f"  Ratio: {c3_ratio}  {mark_c3}"
              + ("   [both sides vanish]" if c3_formula == 0 else ""))
        print()

        # The pair (C_2, C_3) separates the typical blocks: lambda = C_3/C_2
        # and then mu^2 = C_2 + lambda^2.
        if c2_actual != 0:
            lam_rec = simplify(c3_actual / c2_actual)
            mu2_rec = simplify(c2_actual + lam_rec**2)
            ok_rec = (is_zero_scalar(lam_rec - lam)
                      and is_zero_scalar(mu2_rec - (a + 1)**2))
            print(f"Weight recovered from the Casimir pair:  "
                  f"lambda = C_3/C_2 = {lam_rec},  "
                  f"mu² = C_2 + lambda² = {mu2_rec}   "
                  f"{'✓' if ok_rec else '✗'}")
        else:
            print("Weight recovery from (C_2, C_3): C_2 = 0 (atypical), "
                  "the pair does not separate here.")
        print()
        
        # T[0,0] = -(mu^2 - lambda^2), the sign flip being chi = -1 on the top
        # layer.  T itself is the monic (2/3)(N_4 - N_2) built above; this is
        # the eigenvalue readout, the operator identity T = C_2 . chi having
        # already been checked on the full matrices.
        t_eigenvalue = simplify(T[0, 0])
        t_formula = (y_val + a) * (y_val - a - 2)
        t_ratio = (simplify(t_eigenvalue / t_formula) if t_formula != 0
                   else ("0" if t_eigenvalue == 0 else "MISMATCH"))
        mark_t = "✓" if t_ratio in (1, "0") else "✗"
        print(f"T[0,0] = {t_eigenvalue}")
        print(f"  Formula:  -(mu² - lambda²) = (y+a)·(y-a-2) = {t_formula}")
        print(f"  Ratio: {t_ratio}  {mark_t}"
              + ("   [atypical: both sides vanish]" if t_formula == 0 else ""))
        print()
        
        # The anticenter as a polynomial in the Casimirs: T = C_2 . chi, with
        # no cubic term.  The real check is the operator identity printed
        # above; this line echoes it on the [0,0] entry, where chi = -1.
        t_from_c2 = simplify(c2_actual * chi[0, 0])
        ok_tc = is_zero_scalar(t_eigenvalue - t_from_c2)
        print(f"T[0,0] = {t_eigenvalue}")
        print(f"  Formula:  C_2[0,0] · chi[0,0] = {t_from_c2}")
        print(f"  {'✓' if ok_tc else '✗'}"
              + ("" if c2_result['is_scalar_multiple']
                 else "   [entrywise echo only; the operator identity above is the real check]"))
        print()

        # STr(N_4)
        e = simplify(supertrace(chi, N4))
        print(f"STr(N_4) = {e}")
        print()


# -- command-line argument types -------------------------------------------
#
# There is no hard cap on a or N: the construction is valid at any size (see
# ``matryoshka``).  What there is, is a guard on the product, because the cost
# grows as the cube of the dimension 4(a+1)N and it is easy to ask by accident
# for something that will not finish.  SIZE_LIMIT is CPU, not correctness --
# edit it if you mean it.  The working grid a <= 2, N <= 4 sits well inside:
# a = 0, N = 3 (three generations of quarks) gives 3, and a = 2, N = 4 gives 12.

SIZE_LIMIT = 16                 # reject (a+1)*N > SIZE_LIMIT, i.e. dim > 64


def _weight_a(s):
    """a: a non-negative integer."""
    try:
        a = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"a must be an integer, got {s!r}") from None
    if a < 0:
        raise argparse.ArgumentTypeError("a must be a non-negative integer")
    return a


def _layers_N(s):
    """N: a positive integer."""
    try:
        n = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"N must be an integer, got {s!r}") from None
    if n < 1:
        raise argparse.ArgumentTypeError("N must be a positive integer")
    return n


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="sl21_b.py",
        description=(
            "This program computes the matrices of the sl(2/1) superalgebra "
            "with the Kac-Dynkin weight b kept as a FORMAL VARIABLE: every "
            "matrix, bracket and invariant comes out as an exact function of "
            "b.  Only a (and optionally N) are given.  There is no -b and no "
            "-y: use sl21.py for a numeric weight.  Keeping b formal is what "
            "allows real derivatives of every order in b, and the "
            "zeta-Hermitian real form, which depends on sqrt(b) and "
            "sqrt(b-a-1).  Jets are taken in y = 2b - a."
        ),
        epilog=(
            "Authors: Jean Thierry-Mieg (NLM/NIH) and Claude.  "
            "License: public domain; no rights reserved."
        ),
    )
    parser.add_argument(
        "-a", metavar="A", type=_weight_a, required=True,
        help="Kac-Dynkin weight a (non-negative integer)",
    )
    parser.add_argument(
        "-N", metavar="N", type=_layers_N, default=1,
        help="number of layers of the Matryoshka indecomposable "
             "representation, built as the Taylor jet of R(a,b) in y: the p-th "
             "sub-diagonal carries (1/p!) (d/dy)^p mu "
             "(positive integer; default 1 = plain R(a,b))",
    )
    parser.add_argument(
        "--zetaH", action="store_true",
        help="build in the zeta-Hermitian basis: conjugate by a diagonal so "
             "that the odd blocks are symmetric in alpha=sqrt(b) and "
             "gamma=sqrt(b-a-1), then verify zeta M^dag zeta = partner on the "
             "finished matrices.  Requires --case",
    )
    parser.add_argument(
        "--case", metavar="C", type=int, choices=(1, 2, 3),
        help="which region of b, for --zetaH.  "
             "1: b > a+1.  2: 0 < b < a+1 (the quark window).  3: b < 0.  "
             "The eight matrices are the same in all three; only the signs of "
             "zeta and of the conjugation change.  The atypical points "
             "b(b-a-1)=0 are out of scope",
    )
    parser.add_argument(
        "--casimirs", action="store_true",
        help="also build and check the Killing metric and the quadratic and "
             "cubic Casimir operators (off by default; with a formal b this is "
             "a good deal slower than in sl21.py)",
    )
    return parser


if __name__ == "__main__":
    parser = _build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    if args.zetaH and args.case is None:
        parser.error("--zetaH needs --case 1, 2 or 3 (which region of b); "
                     "see --help")
    if args.case is not None and not args.zetaH:
        parser.error("--case only means something with --zetaH")
    if (args.a + 1) * args.N > SIZE_LIMIT:
        parser.error(
            f"(a+1)*N = {(args.a + 1) * args.N} exceeds SIZE_LIMIT = "
            f"{SIZE_LIMIT}, i.e. dimension {4 * (args.a + 1) * args.N} > "
            f"{4 * SIZE_LIMIT}.  This is a CPU guard, not a limit of the "
            f"construction, which is valid at any a and N.  If you really "
            f"want this, edit SIZE_LIMIT near the bottom of sl21_b.py."
        )

    main(args.a, args.N, casimirs=args.casimirs,
         zetaH=args.zetaH, case=args.case)
