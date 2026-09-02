"""
sl21.py -- Kac modules of the superalgebra sl(2|1).

Builds R(a, b): the 4(a+1)-dimensional Kac module with even generators
Y, e, f, h (numbers 0,1,2,3) and odd generators u, v, w, x (numbers 4,5,6,7),
then verifies every sl(2|1) (anti)commutator.  Run it as

        python sl21.py -a <int> -b <int>          (e.g. -a 1 -b 1: the adjoint)
        python sl21.py -a <int> -y <rational>     (y = (a+b)/2, b derived from y)
        python sl21.py -h                         (help)

``a`` and ``b`` are the Kac-Dynkin weights (Kac's own naming); y = (a+b)/2 is
the module's central-charge eigenvalue on the top layer.

The module splits into four sl(2) layers, block-diagonal for e,f,h and scalar
for Y (eigenvalues y, y-1, y-1, y-2 with y=(a+b)/2):

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
"""

import sys
import argparse
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
    carrying all the b-dependence.  m21,m43 vanish at b=-3a; m31,m42 at b=a+4."""
    D = 4 * (a + 1)
    s12, s13, s24, s34 = 1, 1, -1, 1
    m21 = Rational(1, D) * (b + 3 * a)
    m43 = Rational(1, D) * (b + 3 * a)
    m31 = Rational(1, D) * ((a + 4) - b)
    m42 = Rational(1, D) * (b - (a + 4))
    return (s12, s13, s24, s34, m21, m31, m42, m43)


def _build(a, b):
    """Assemble the eight 4(a+1)-square generator matrices for R(a, b)."""
    s12, s13, s24, s34, m21, m31, m42, m43 = _scales(a, b)
    y = (a + b) * Rational(1, 2)

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

    R(a,b) is affine in b (b enters only linearly, through y=(a+b)/2 in Y and
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


def main(a, b, N=1):
    rep = Rmatryoshka(a, b, N)

    y = Rational(a + b, 2)
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


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="sl21.py",
        description=(
            "This program computes the matrices of the sl(2/1) superalgebra: "
            "please specify the Kac-Dynkin weights a,b or a,y where "
            "y=(a+b)/2."
        ),
        epilog=(
            "Authors: Jean Thierry-Mieg (NLM/NIH) and Claude.  "
            "License: public domain; no rights reserved."
        ),
    )
    parser.add_argument(
        "-a", metavar="A", type=int, required=True,
        help="Kac-Dynkin weight a (non-negative integer)",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-b", metavar="B", type=str,
        help="Kac-Dynkin weight b (integer or fraction, e.g. 3/2)",
    )
    group.add_argument(
        "-y", metavar="Y", type=str,
        help="weight y = (a+b)/2 (integer or fraction, e.g. 3/2)",
    )
    parser.add_argument(
        "-N", metavar="N", type=int, default=1,
        help="number of layers of a Matryoshka indecomposable representation "
             "with N layers (positive integer < 5; default 1 = plain R(a,b))",
    )
    return parser


if __name__ == "__main__":
    parser = _build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    a = args.a
    if args.b is not None:
        b = Rational(args.b)
    else:
        y = Rational(args.y)
        b = 2 * y - a

    N = args.N
    if not (1 <= N <= 4):
        parser.error("N must be a positive integer smaller than 5")

    main(a, b, N)