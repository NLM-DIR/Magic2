

"""
sl2.py -- first usage of the library.
 
Builds R[a], the (a+1)-dimensional irreducible representation of sl(2), whose
three even generators e, f, h carry the canonical numbers 1, 2, 3, then verifies
the Chevalley commutators
 
        [h, e] = 2 e ,   [h, f] = -2 f ,   [e, f] = h
 
with everything else zero. Run it as
 
        python sl2.py [a]          (a defaults to 2)
 
It prints the three matrices, the formal nonzero commutators, and either
"verified" or the residual matrices for whatever failed.
 
Construction (this is the standard irrep; the three matrices are not free --
once f is the sub-diagonal of ones and h is the given diagonal, the requirement
[e, f] = h pins e down completely):
 
    h[k,k]   = a - 2k                 diagonal a, a-2, ..., -a   (top = a)
    f[k+1,k] = 1                      ones just below the diagonal
    e[k,k+1] = (k+1)(a-k)             just above the diagonal; symmetric,
                                      rising then falling, ends = a
"""
 
import sys
 
from matrix import Matrix
from algebra import Algebra, EVEN
 
# Canonical numbers, fixed once for the whole family.
# Y is the central identity (gl(2) = sl(2) + scalars); it sits at 0.
Y, E, F, H = 0, 1, 2, 3
 
# The independent nonzero brackets, written exactly as stated. Everything else
# (including the reversed pairs, via antisymmetry) follows from these.
INDEPENDENT = {
    (H, E): {E: 2},     # [h, e] =  2 e
    (H, F): {F: -2},    # [h, f] = -2 f
    (E, F): {H: 1},     # [e, f] =  h
}
 
 
# -- 1. construct the representation R[a] ------------------------------------
 
def Rsl2(a):
    """Return R[a]: the (a+1)-dimensional irrep of sl(2)."""
    if a < 0:
        raise ValueError("a must be a non-negative integer")
    n = a + 1
 
    e = Matrix(n, n, lambda i, j: (i + 1) * (a - i) if j == i + 1 else 0)
    f = Matrix(n, n, lambda i, j: 1 if i == j + 1 else 0)
    h = Matrix(n, n, lambda i, j: (a - 2 * i) if i == j else 0)
 
    g = Algebra(name=f"R[{a}]")
    g.add(e, "e", EVEN)                        # -> canonical number 1
    g.add(f, "f", EVEN)                        # -> canonical number 2
    g.add(h, "h", EVEN)                        # -> canonical number 3
    g.add(Matrix.one(n), "Y", EVEN, number=0)  # central identity -> number 0
    return g
 
 
# -- 2. the structure constants f(i, j, k), initialized explicitly -----------
 
def structure_constants():
    """Build f(i, j, k) from INDEPENDENT, adding the antisymmetric partners so
    the check is correct for either ordering of a pair ([y,x] = -[x,y])."""
    table = {}
    for (i, j), terms in INDEPENDENT.items():
        table[(i, j)] = dict(terms)
        table[(j, i)] = {k: -c for k, c in terms.items()}
 
    def f(i, j, k):
        return table.get((i, j), {}).get(k, 0)
 
    return f
 
 
# -- 3. verification ---------------------------------------------------------
 
def verify(rep, f):
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
 
 
def main(a):
    rep = Rsl2(a)
    f = structure_constants()
 
    print(f"sl(2) representation R[{a}]   (dimension {rep.dim})")
    print("=" * 40, "\n")
 
    for num in rep.numbers():
        _show_matrix(f"{rep.title(num)}  (#{num})", rep[num])
 
    print("Nonzero commutators to verify:")
    for (i, j), terms in INDEPENDENT.items():
        print(f"    [{rep.title(i)}, {rep.title(j)}] = {_combo(rep, terms)}")
    if Y in rep.numbers():
        print(f"    [{rep.title(Y)}, .] = 0   ({rep.title(Y)} is central)")
    print()
 
    residuals = verify(rep, f)
    if not residuals:
        print("All Chevalley relations verified. \u2713")
    else:
        print(f"{len(residuals)} relation(s) FAILED:\n")
        for r in residuals:
            print(f"    [{r.title_i}, {r.title_j}}} residual =")
            print(r.matrix)
            print()
 
 
if __name__ == "__main__":
    a = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    main(a)
 
