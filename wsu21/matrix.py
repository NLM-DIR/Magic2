"""
matrix.py -- exact dense matrices, blind to the scalar type T.
 
The matrix layer depends ONLY on the ring contract of its entries:
 
    a + b        __add__  / __radd__
    a - b        __sub__  / __neg__
    a * b        __mul__  / __rmul__
    a == b       __eq__
    repr(a)      __repr__
 
It never inspects, imports, or names the scalar type. Entries may be SymPy
scalars today and Grassmann elements tomorrow; this file does not change
(design invariant 1).
 
Two small ring facts are relied on, both standard:
 
  * A ring embeds the integers, so the Python literals ``0`` and ``1`` coerce
    to the additive / multiplicative identity through the entry type's
    ``__radd__`` / ``__rmul__``. This is the only reason ``zero`` and ``one``
    can exist without naming T.
  * Multiplication may be non-commutative (Grassmann). Every operation below
    therefore keeps left/right order exactly: ``c * A`` puts the scalar on the
    left of each entry, ``A * c`` on the right, and ``A @ B`` forms products in
    the order ``A[i,k] * B[k,j]``.
 
Equality here is *structural* -- it uses the entry type's ``==`` and nothing
more. For SymPy scalars that can miss algebraic zeros, so trustworthy
zero-testing is the caller's job at a layer that knows T: normalise first, e.g.
 
    residual.applyfunc(simplify) == Matrix.zero(n, n)
 
(design invariant 4). matrix.py deliberately does not simplify anything.
 
All real presentation (titles, nonzero filtering, LaTeX) lives in report.py;
``__repr__`` here is a plain developer view only (design invariant 2).
"""
 
from functools import reduce
import operator
 
 
class Matrix:
    __slots__ = ("rows", "cols", "_data")
 
    # -- construction --------------------------------------------------------
 
    def __init__(self, rows, cols, entry):
        """Build a ``rows`` x ``cols`` matrix by addressing each cell.
 
        ``entry(i, j)`` is called for every ``0 <= i < rows`` and
        ``0 <= j < cols`` and returns that cell's value (of type T). This
        callback *is* the cell-addressing initialiser, e.g.::
 
            Matrix(10, 10, lambda i, j: 1 if i == j else 0)
        """
        if rows < 0 or cols < 0:
            raise ValueError("dimensions must be non-negative")
        if not callable(entry):
            raise TypeError("entry must be callable: entry(i, j) -> cell")
        self.rows = rows
        self.cols = cols
        self._data = [[entry(i, j) for j in range(cols)] for i in range(rows)]
 
    @classmethod
    def _wrap(cls, rows, cols, grid):
        """Adopt an already-built grid (list of rows) without re-running a
        callback. Internal fast path; ``grid`` is taken as-is, not copied, so
        callers must hand over a freshly built structure."""
        m = object.__new__(cls)
        m.rows = rows
        m.cols = cols
        m._data = grid
        return m
 
    @classmethod
    def zero(cls, n, m):
        """The ``n`` x ``m`` zero matrix (entries are the ring's additive id)."""
        return cls._wrap(n, m, [[0] * m for _ in range(n)])
 
    @classmethod
    def one(cls, n):
        """The ``n`` x ``n`` identity matrix."""
        return cls._wrap(
            n, n, [[1 if i == j else 0 for j in range(n)] for i in range(n)]
        )
 
    # -- access --------------------------------------------------------------
 
    def __getitem__(self, key):
        try:
            i, j = key
        except (TypeError, ValueError):
            raise TypeError("index a Matrix as A[i, j]") from None
        return self._data[i][j]
 
    def applyfunc(self, f):
        """Return a new Matrix with ``f`` applied to every entry.
 
        The one hook the T-aware layers use to normalise before a zero-test
        (``applyfunc(simplify)``) or to expand/collect. Still blind to T: it
        only maps a caller-supplied function over the cells."""
        return Matrix._wrap(
            self.rows, self.cols, [[f(x) for x in row] for row in self._data]
        )
 
    # -- arithmetic ----------------------------------------------------------
 
    def __add__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        self._require_same_shape(other, "+")
        a, b = self._data, other._data
        grid = [[a[i][j] + b[i][j] for j in range(self.cols)]
                for i in range(self.rows)]
        return Matrix._wrap(self.rows, self.cols, grid)
 
    def __sub__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        self._require_same_shape(other, "-")
        a, b = self._data, other._data
        grid = [[a[i][j] - b[i][j] for j in range(self.cols)]
                for i in range(self.rows)]
        return Matrix._wrap(self.rows, self.cols, grid)
 
    def __neg__(self):
        return Matrix._wrap(
            self.rows, self.cols, [[-x for x in row] for row in self._data]
        )
 
    def __matmul__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        if self.cols != other.rows:
            raise ValueError(
                f"shape mismatch for @: {self._shape()} @ {other._shape()}"
            )
        n = self.cols
        a, b = self._data, other._data
 
        def cell(i, j):
            terms = [a[i][k] * b[k][j] for k in range(n)]
            return reduce(operator.add, terms) if terms else 0
 
        grid = [[cell(i, j) for j in range(other.cols)]
                for i in range(self.rows)]
        return Matrix._wrap(self.rows, other.cols, grid)
 
    def __mul__(self, scalar):
        # A * c  -- scalar on the RIGHT of each entry. A * B is rejected so the
        # non-commutative-safe matrix product always goes through @.
        if isinstance(scalar, Matrix):
            raise TypeError(
                "use @ for matrix multiplication; * is scalar multiplication"
            )
        return Matrix._wrap(
            self.rows, self.cols, [[x * scalar for x in row] for row in self._data]
        )
 
    def __rmul__(self, scalar):
        # c * A  -- scalar on the LEFT of each entry.
        return Matrix._wrap(
            self.rows, self.cols, [[scalar * x for x in row] for row in self._data]
        )
 
    # -- comparison ----------------------------------------------------------
 
    def __eq__(self, other):
        # Structural equality via the entry type's ==. See module docstring:
        # algebraic zero-testing belongs to the T-aware caller.
        if not isinstance(other, Matrix):
            return NotImplemented
        if (self.rows, self.cols) != (other.rows, other.cols):
            return False
        a, b = self._data, other._data
        return all(a[i][j] == b[i][j]
                   for i in range(self.rows) for j in range(self.cols))
 
    __hash__ = None  # matrices are values, not keys; equality is structural
 
    # -- display (developer view; real presentation lives in report.py) ------
 
    def __repr__(self):
        if self.rows == 0 or self.cols == 0:
            return f"Matrix({self.rows}, {self.cols}, [])"
        cells = [[repr(x) for x in row] for row in self._data]
        width = [max(len(cells[i][j]) for i in range(self.rows))
                 for j in range(self.cols)]
        lines = [
            "[ " + "  ".join(cells[i][j].rjust(width[j])
                             for j in range(self.cols)) + " ]"
            for i in range(self.rows)
        ]
        return "\n".join(lines)
 
    # -- internal helpers ----------------------------------------------------
 
    def _shape(self):
        return f"{self.rows}x{self.cols}"
 
    def _require_same_shape(self, other, op):
        if (self.rows, self.cols) != (other.rows, other.cols):
            raise ValueError(
                f"shape mismatch for {op}: {self._shape()} {op} {other._shape()}"
            )
 
