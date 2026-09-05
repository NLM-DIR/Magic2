"""
sl21.py -- Kac modules of the superalgebra sl(2|1).

Builds R(a, b): the 4(a+1)-dimensional Kac module with even generators
Y, e, f, h (numbers 0,1,2,3) and odd generators u, v, w, x (numbers 4,5,6,7),
then verifies every sl(2|1) (anti)commutator.  Run it as

        python sl21.py -a <int> -b <rational>     (e.g. -a 1 -b 1: the adjoint)
        python sl21.py -a <int> -y <rational>     (b derived as (a+y)/2)
        python sl21.py -h                         (help)

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
"""

import sys
import argparse
from itertools import permutations
from collections import namedtuple

from sympy import Rational, sympify, simplify

from matrix import Matrix
from algebra import Algebra, EVEN, ODD, is_zero_scalar

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


def _build(a, b):
    """Assemble the eight 4(a+1)-square generator matrices for R(a, b)."""
    s12, s13, s24, s34, m21, m31, m42, m43 = _scales(a, b)
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


def Rsl21(a, b):
    """Return R(a, b): the 4(a+1)-dimensional Kac module of sl(2|1).

    ``a`` is a non-negative integer; ``b`` may be an int or a SymPy symbol."""
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


def Rderiv(a):
    """R'(a,b): the derivative of R(a,b) with respect to b, per generator.

    R(a,b) is affine in b (b enters only linearly, through y=2b-a in Y and
    through the raising scales m21,m31,m42,m43 in u,w). Its b-derivative is
    therefore the constant coefficient of b, computed exactly by finite
    difference::

        R'(a,b) = R(a,1) - R(a,0)          (independent of b)

    Returns a dict {canonical number -> Matrix}. Only Y, u, w have nonzero
    derivative blocks; e,f,h,v,x are b-independent and map to zero."""
    r1 = Rsl21(a, 1)
    r0 = Rsl21(a, 0)
    return {k: r1[k] - r0[k] for k in r1.numbers()}


def Rmatryoshka(a, b, N):
    """Return MR(a,b,N): the Matryoshka indecomposable representation, N layers.

    MR is block N x N in the layer index. Each generator matrix is::

        [ R   0   0  ... ]
        [ R'  R   0  ... ]        R  = R(a,b)  on the diagonal blocks
        [ 0   R'  R  ... ]        R' = R'(a,b) on the first sub-diagonal
        [ ...            ]        0  everywhere else (upper and lower-by-2+)

    This is R(a, b*Id_N + n) with n the strictly-lower Jordan shift (n^N = 0):
    since R is affine in b the substitution is exact, so MR satisfies the same
    b-independent structure constants and is checked by the same f. N=1 gives
    back R(a,b) unchanged."""
    if not (1 <= N <= 4):
        raise ValueError("N must be a positive integer < 5")
    base = Rsl21(a, b)
    dR = Rderiv(a)
    d = base.dim                                   # block size = 4(a+1)
    D = N * d                                       # full size
    nums = base.numbers()

    g = Algebra(name=f"MR(a={a},b={b},N={N})")
    for k in nums:
        R, Rp = base[k], dR[k]

        def entry(i, j, R=R, Rp=Rp, d=d):
            I, J = i // d, j // d                  # block row / column
            bi, bj = i % d, j % d                  # position within the block
            if I == J:
                return R[bi, bj]                   # diagonal: R(a,b)
            if I == J + 1:
                return Rp[bi, bj]                  # sub-diagonal: R'(a,b)
            return 0                               # everything else empty

        g.add(Matrix(D, D, entry), base.title(k), base.parity(k), number=k)
    return g


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


def chi_relations(rep, chi):
    """The chi relations as plain data (never printed).

    For every generator g: the super-bracket of chi (odd) with g, i.e.
    [chi, g] = chi g - g chi when g is even and {chi, g} = chi g + g chi when g
    is odd; each must be the zero matrix. Plus the involution chi^2 = 1. Zero-
    testing goes through is_zero_scalar (invariant 4), same policy as check()."""
    D = chi.rows

    def is_zero_matrix(M):
        return all(is_zero_scalar(M[i, j])
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


def casimir_commutes(rep, C2, nam):
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
        return all(is_zero_scalar(M[i, j])
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


def anticenter_relations(rep, T):
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
        return all(is_zero_scalar(M[i, j])
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


def anticenter_casimir_relations(rep, T, C2, chi):
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
        return all(is_zero_scalar(M[i, j])
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

def verify(rep, f=f):
    """Nonzero residuals of the structure-constant check (empty = all hold)."""
    return rep.check(f)


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



def main(a, b, N=1, casimirs=False):
    rep = Rmatryoshka(a, b, N)

    y = 2 * b - a
    if N == 1:
        print(f"sl(2|1) Kac module R(a={a}, b={b})  [y={y}]   "
              f"(dimension {rep.dim} = 4(a+1))")
    else:
        print(f"sl(2|1) Matryoshka indecomposable rep MR(a={a}, b={b}, N={N})  "
              f"[y={y}]   (dimension {rep.dim} = {N}*4(a+1))")
    print("=" * 46, "\n")

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

    residuals = verify(rep)
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
    chi_checks = chi_relations(rep, chi)
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

    if casimirs:
        cas = _casimir_report(rep, chi, N, a, y)

    if casimirs:
        # -- the Gorelik anticenter element T_4: N_2 and N_4 ---------------------
        print()
        print("=" * 80)
        print("Gorelik anticenter T_4:  N_2 then N_4  (compare by eye)")
        print("=" * 80)
        print()

        N2 = anticenter_N2(rep)
        print("N_2 = UV - VU + WX - XW =")
        print(N2)
        print()

        N4 = anticenter_N4(rep)
        print("N_4 = Σ_σ sign(σ) (u v w x permuted), fully antisymmetric =")
        print(N4)
        print()

        T = (Rational(2, 3) * (N4 - N2)).applyfunc(simplify)
        print("T = (2/3)(N_4 - N_2)   [monic normalisation] =")
        print(T)
        print()

        print("Anticenter relations to verify numerically "
              "([T,even]=0, {T,odd}=0):")
        t_checks = anticenter_relations(rep, T)
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
        tc_checks = anticenter_casimir_relations(rep, T, cas['c2']['total'], chi)
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

A_MAX = 10          # largest a accepted (dimension 4(a+1) = 44)
N_MAX = 4           # largest Matryoshka depth accepted
TOO_HEAVY = "valid but this would make the calculation too heavy, sorry"


def _weight_a(s):
    """a: a non-negative integer, at most A_MAX."""
    try:
        a = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"a must be an integer, got {s!r}") from None
    if a < 0:
        raise argparse.ArgumentTypeError("a must be a non-negative integer")
    if a > A_MAX:
        raise argparse.ArgumentTypeError(f"a = {a} is {TOO_HEAVY} (max {A_MAX})")
    return a


def _layers_N(s):
    """N: a positive integer, at most N_MAX."""
    try:
        n = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"N must be an integer, got {s!r}") from None
    if n < 1:
        raise argparse.ArgumentTypeError("N must be a positive integer")
    if n > N_MAX:
        raise argparse.ArgumentTypeError(f"N = {n} is {TOO_HEAVY} (max {N_MAX})")
    return n


def _fraction(s):
    """b or y: an exact rational such as 3, -2, 3/2, -7/3, (-7/3) or 1.5.

    Surrounding parentheses and blanks are ignored, so that a negative fraction
    can be written either as  -y -7/3  or  -y "(-7/3)"."""
    text = s.strip()
    while text.startswith("(") and text.endswith(")"):
        text = text[1:-1].strip()
    text = text.replace(" ", "")
    try:
        return Rational(text)
    except (TypeError, ValueError, SyntaxError):
        raise argparse.ArgumentTypeError(
            f"expected a rational number like 3/2 or -7/3, got {s!r}") from None


def _protect_negative_fractions(argv):
    """Rewrite  -b -7/3  as  -b=-7/3  so argparse does not mistake the value for
    an option; argparse already accepts plain negative integers like -7."""
    out = []
    skip = False
    for i, tok in enumerate(argv):
        if skip:
            skip = False
            continue
        if tok in ("-b", "-y") and i + 1 < len(argv) \
                and argv[i + 1].startswith("-"):
            out.append(f"{tok}={argv[i + 1]}")
            skip = True
        else:
            out.append(tok)
    return out


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="sl21.py",
        description=(
            "This program computes the matrices of the sl(2/1) superalgebra: "
            "please specify the Kac-Dynkin weights a,b (Kac convention, "
            "b=(a+y)/2, atypical at b=0 and b=a+1) or equivalently a,y "
            "where y=2b-a."
        ),
        epilog=(
            "Authors: Jean Thierry-Mieg (NLM/NIH) and Claude.  "
            "License: public domain; no rights reserved."
        ),
    )
    parser.add_argument(
        "-a", metavar="A", type=_weight_a, required=True,
        help=f"Kac-Dynkin weight a (non-negative integer, at most {A_MAX})",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-b", metavar="B", type=_fraction,
        help="Kac-Dynkin weight b = (a+y)/2 (integer or fraction, e.g. 3/2 "
             "or -7/3); atypical at b=0 and b=a+1",
    )
    group.add_argument(
        "-y", metavar="Y", type=_fraction,
        help="weight y = 2b - a (integer or fraction, e.g. 3/2 or -7/3); "
             "b is derived as (a+y)/2",
    )
    parser.add_argument(
        "-N", metavar="N", type=_layers_N, default=1,
        help="number of layers of a Matryoshka indecomposable representation "
             f"(positive integer, at most {N_MAX}; default 1 = plain R(a,b))",
    )
    parser.add_argument(
        "--casimirs", action="store_true",
        help="also build and check the Killing metric and the quadratic and "
             "cubic Casimir operators (slower; off by default)",
    )
    return parser


if __name__ == "__main__":
    parser = _build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args(_protect_negative_fractions(sys.argv[1:]))

    a = args.a
    b = args.b if args.b is not None else Rational(args.a + args.y, 2)

    main(a, b, args.N, casimirs=args.casimirs)
