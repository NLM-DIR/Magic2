"""
sl21.py -- Kac modules of the superalgebra sl(2|1).

Builds R(a, b): the 4(a+1)-dimensional Kac module with even generators
Y, e, f, h (numbers 0,1,2,3) and odd generators u, v, w, x (numbers 4,5,6,7),
then verifies every sl(2|1) (anti)commutator.  Run it as

        python sl21.py [a] [b]            (defaults a=1, b=1: the adjoint)

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

from sympy import Rational, sympify

from matrix import Matrix
from algebra import Algebra, EVEN, ODD

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
        put(Xmat, Lt, Ls, _shape(kind, js, 'x'), sc)

    # raising blocks carry u, w  (all the b-dependence)
    for (Lt, Ls, kind, js, sc) in [(0, 1, 'down', a + 1, m21),
                                   (0, 2, 'up',   a - 1, m31),
                                   (1, 3, 'up',   a,     m42),
                                   (2, 3, 'down', a,     m43)]:
        if d[Lt] == 0 or d[Ls] == 0:
            continue
        put(Umat, Lt, Ls, _shape(kind, js, 'u'), sc)
        put(Wmat, Lt, Ls, _shape(kind, js, 'w'), sc)

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


def main(a, b):
    rep = Rsl21(a, b)

    print(f"sl(2|1) Kac module R({a},{b})   (dimension {rep.dim} = 4(a+1))")
    print("=" * 46, "\n")

    for num in rep.numbers():
        _show_matrix(f"{rep.title(num)}  (#{num})", rep[num])

    print("Nonzero (super-)brackets to verify:")
    for (i, j), terms in INDEPENDENT.items():
        print("    " + _bracket(rep, i, j, terms))
    print()

    residuals = verify(rep)
    if not residuals:
        print("All sl(2|1) relations verified. \u2713")
    else:
        print(f"{len(residuals)} relation(s) FAILED:\n")
        for r in residuals:
            op, cl = ("{", "}") if (PARITY[r.i] == ODD and PARITY[r.j] == ODD) \
                else ("[", "]")
            print(f"    {op}{r.title_i}, {r.title_j}{cl} residual =")
            print(r.matrix)
            print()


if __name__ == "__main__":
    a = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    b = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    main(a, b)

    
