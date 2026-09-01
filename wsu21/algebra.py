"""
algebra.py -- the M[] registry, the super-bracket, and the structure-constant
check, plus the R[a, b] family of representations.
 
One ``Algebra`` instance is ONE representation: a table of up to ~12 generator
matrices, each carrying a title and a parity bit. Generators are numbered
0, 1, 2, ... in the order they are added. That numbering is the abstract
algebra's labelling, so if every rep adds its generators in the same order the
structure-constant function ``f(i, j, k)`` is universal -- the same ``f`` checks
every rep (design: parity is a bit in the registry, independent of the entry
type; see MATRIX_CONTEXT sec. 4 and 6).
 
The super-bracket takes its sign from the stored parities, never by hand:
 
    [A, B}  =  A @ B  -  (-1)^(|A|*|B|) * B @ A
             =  A @ B  -  B @ A     if A or B is even   (commutator)
             =  A @ B  +  B @ A     if both are odd      (anticommutator)
 
``check(f)`` returns the nonzero residuals of  [M_i, M_j} - Σ_k f(i,j,k) M_k  as
plain data; it never prints. Formatting is report.py's job (invariant 2).
 
Zero-testing is done HERE, at a layer that knows the scalars are SymPy: a
residual entry counts as zero only after ``is_zero`` / ``simplify`` clears it,
never on a bare structural ``==`` (invariant 4). matrix.py stays blind to all
of this.
"""
 
from collections import namedtuple
 
from sympy import sympify, simplify
 
from matrix import Matrix
 
# -- parity ------------------------------------------------------------------
 
EVEN = 0
ODD = 1
 
# One reported discrepancy: which pair, the two titles, and the leftover matrix.
Residual = namedtuple("Residual", "i j title_i title_j matrix")
 
_Gen = namedtuple("_Gen", "M title parity")
 
 
# -- scalar zero-testing (SymPy-aware; the policy invariant 4 demands) -------
 
def is_zero_scalar(x):
    """True iff the SymPy scalar ``x`` is algebraically zero.
 
    Fast path on SymPy's ``.is_zero`` assumption; falls back to ``simplify``
    for the cases it leaves undecided. Never a bare ``==`` on the raw form."""
    x = sympify(x)
    decided = x.is_zero
    if decided is True:
        return True
    if decided is False:
        return False
    return simplify(x) == 0
 
 
def _matrix_is_zero(M, is_zero):
    return all(is_zero(M[i, j]) for i in range(M.rows) for j in range(M.cols))
 
 
# -- one representation ------------------------------------------------------
 
class Algebra:
    """A single representation: the table of generator matrices, keyed by their
    canonical number 1, 2, 3, ... (1-based, in the order they are added)."""
 
    def __init__(self, name=""):
        self.name = name
        self.dim = None          # common d of the d x d generator matrices
        self._gens = {}          # canonical number -> _Gen (insertion-ordered)
 
    # registration --------------------------------------------------------
 
    def add(self, M, title, parity, number=None):
        """Register a generator; return its canonical number.
 
        Numbers default to 1, 2, 3, ... in add order, but may be given
        explicitly -- e.g. a central element at 0. The number is the universal
        label ``f(i, j, k)`` uses, so assign it the same way in every rep.
        All generators in a rep must be square and share one dimension d."""
        if parity not in (EVEN, ODD):
            raise ValueError("parity must be EVEN (0) or ODD (1)")
        if M.rows != M.cols:
            raise ValueError(
                f"generator {title!r} must be square, got {M.rows}x{M.cols}"
            )
        if self.dim is None:
            self.dim = M.rows
        elif M.rows != self.dim:
            raise ValueError(
                f"generator {title!r} has dimension {M.rows}, "
                f"but this rep is {self.dim}-dimensional"
            )
        if number is None:
            number = max(self._gens) + 1 if self._gens else 1
        elif number in self._gens:
            raise ValueError(
                f"canonical number {number} is already "
                f"{self._gens[number].title!r}"
            )
        self._gens[number] = _Gen(M, title, parity)
        return number
 
    # access --------------------------------------------------------------
 
    def numbers(self):
        """The canonical numbers, ascending (e.g. [0, 1, 2, 3])."""
        return sorted(self._gens)
 
    def __getitem__(self, i):
        return self._gens[i].M
 
    def title(self, i):
        return self._gens[i].title
 
    def parity(self, i):
        return self._gens[i].parity
 
    def __len__(self):
        return len(self._gens)
 
    def __iter__(self):
        return iter(self.numbers())
 
    def __repr__(self):
        tag = f"{self.name!r}: " if self.name else ""
        return f"Algebra({tag}{len(self)} generators, dim {self.dim})"
 
    # the super-bracket ---------------------------------------------------
 
    def bracket(self, i, j):
        """The super-bracket [M_i, M_j}; sign taken from the stored parities."""
        A, B = self[i], self[j]
        if self.parity(i) == ODD and self.parity(j) == ODD:
            return A @ B + B @ A            # both odd -> anticommutator
        return A @ B - B @ A                # otherwise -> commutator
 
    # the structure-constant check ---------------------------------------
 
    def check(self, f, pairs=None, is_zero=is_zero_scalar):
        """Residuals of  [M_i, M_j} - Σ_k f(i,j,k) M_k  that are not zero.
 
        ``f(i, j, k)`` returns the scalar structure constant. By default every
        pair with ``i <= j`` is checked (the non-redundant set; the diagonal is
        kept because an odd generator's self-anticommutator is a real relation).
        Pass ``pairs`` to override, or a custom ``is_zero`` for other scalars.
 
        Returns a list of ``Residual`` records -- data only, never printed."""
        nums = self.numbers()
        d = self.dim
        if pairs is None:
            pairs = [(nums[p], nums[q])
                     for p in range(len(nums)) for q in range(p, len(nums))]
 
        out = []
        for i, j in pairs:
            combo = Matrix.zero(d, d)
            for k in nums:
                c = f(i, j, k)
                if c == 0:                 # a skip is always safe, even if only
                    continue               # structural; adding 0*M_k is a no-op
                combo = combo + c * self[k]
            residual = self.bracket(i, j) - combo
            if not _matrix_is_zero(residual, is_zero):
                out.append(
                    Residual(i, j, self.title(i), self.title(j), residual)
                )
        return out
 
 
# -- the R[a, b] family ------------------------------------------------------
 
class Family:
    """A parametrised family of representations, indexed ``R[a, b]``.
 
    Wrap a builder ``(a, b) -> Algebra``; then ``R[a, b]`` builds that rep once
    and caches it, so several reps can be held side by side cheaply. Works as a
    decorator::
 
        @Family
        def R(a, b):
            g = Algebra(name=f"R[{a},{b}]")
            g.add(..., "H", EVEN)
            g.add(..., "E", EVEN)
            ...
            return g
 
        rep = R[1, 2]        # an Algebra; R[3, 5] is another, independent one
    """
 
    def __init__(self, builder):
        self._builder = builder
        self._cache = {}
 
    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if key not in self._cache:
            self._cache[key] = self._builder(*key)
        return self._cache[key]
 
    def clear(self):
        """Drop cached reps (e.g. after changing the builder)."""
        self._cache.clear()
 
    def __repr__(self):
        return f"Family({self._builder.__name__}, {len(self._cache)} cached)"
 
