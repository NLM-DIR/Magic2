# Exact matrix representations of sl(2|1): mathematics and usage

This manual describes what the programs compute and how to run them.
Installation is covered in `README.md`.

The programs build the representation matrices of the Lie superalgebra
sl(2|1) in exact arithmetic (rational numbers, or symbolic expressions in the
weight b). They then verify every (super-)commutator and the classical
invariants: the grading operator, the Casimir operators and the Gorelik
anticenter. Every statement marked *verified* below is checked by the programs
themselves, on the full matrices, each time they run.

1. The library
2. sl(2): `sl2.py`
3. sl(2|1) and its Kac modules: `sl21.py`
4. The grading operator χ
5. Indecomposable representations: the Matryoshka and its jets
6. Invariant forms and the master equation for the cubic tensor
7. The Casimir operators C₂ and C₃
8. The Gorelik anticenter T
9. Formal b: `sl21_b.py`
10. The ζ-Hermitian real forms
11. Command-line reference
12. Reference values

---

## 1. The library

`matrix.py` provides an exact matrix class. Matrices are built by addressing
cells, `Matrix(n, m, lambda i, j: ...)`, and support `+`, `-`, `@` (matrix
product) and scalar multiplication. The class never inspects its entries. It
only adds, multiplies and compares them, so the entries can be integers,
rationals, square roots or symbolic expressions (all handled by SymPy).

`algebra.py` holds a list of generator matrices, each with a canonical number,
a title and a parity (even or odd). It computes the super-bracket

    [A, B} = AB − (−1)^{|A||B|} BA

which is a commutator unless both matrices are odd, in which case it is an
anticommutator. It then checks the bracket against structure constants
f(i, j, k):

    residual(i, j) = [M_i, M_j} − Σ_k f(i,j,k) M_k

Every pair i ≤ j is checked, including those whose bracket should vanish. A
residual is declared zero only after symbolic simplification, so an
expression such as `sqrt(b)**2 − b` is correctly recognised as zero.

## 2. sl(2): `sl2.py`

`sl2.py` builds the (a+1)-dimensional irreducible representation R[a] of
sl(2), with highest weight a, in the integer basis

    h[k,k]   = a − 2k          diagonal a, a−2, …, −a
    f[k+1,k] = 1               ones below the diagonal
    e[k,k+1] = (k+1)(a−k)      fixed by [e,f] = h

and verifies [h,e] = 2e, [h,f] = −2f, [e,f] = h. This integer basis, with f
made of ones, is used for every sl(2) layer below, so all sl(2|1) matrices
stay rational.

## 3. sl(2|1) and its Kac modules: `sl21.py`

### Generators and structure constants

sl(2|1) has four even generators, Y (hypercharge), e, f, h (sl(2)), and four
odd generators u, v, w, x, with canonical numbers 0…7 in that order. The
nonzero brackets are

    [h,e] = 2e     [h,f] = −2f     [e,f] = h
    [Y,u] = u      [Y,v] = −v      [Y,w] = w      [Y,x] = −x
    [h,u] = −u     [h,v] = v       [h,w] = w      [h,x] = −x
    [e,u] = w      [f,w] = u       [e,x] = v      [f,v] = x
    {u,v} = (Y+h)/2      {w,x} = (h−Y)/2      {v,w} = −e      {u,x} = f

Y is central in the even part. The odd generators form two sl(2) doublets in
the natural basis: (u, w) with Y = +1, and (x, v) with Y = −1. In each
doublet e raises the lower member to the upper one and f lowers it back.

### The Kac modules R(a, b)

R(a,b) is the Kac module induced from the (a+1)-dimensional sl(2)
representation of the top layer. It has dimension 4(a+1). We use the
Kac–Dynkin labels in Kac's convention:

    a ≥ 0 an integer,      b = (a+y)/2,      equivalently  y = 2b − a

where y is the eigenvalue of Y on the top layer. b may be any rational (in
`sl21.py`) or a formal symbol (in `sl21_b.py`, §9).

As a representation of sl(2) ⊕ Y, the module splits into four layers:

    L1 = R[a]      Y = y        dim a+1       (top)
    L2 = R[a+1]    Y = y−1      dim a+2
    L3 = R[a−1]    Y = y−1      dim a         (empty when a = 0)
    L4 = R[a]      Y = y−2      dim a+1       (bottom)

The matrices are ordered layer by layer. e, f, h are block diagonal (the
matrices of §2 on each layer) and Y is diagonal and constant on each layer.
The odd generators map between layers. On each block they are sl(2) doublet
intertwiners in the integer basis, multiplied by a scale:

* The lowering operators v and x map L1 → L2, L3 and L2, L3 → L4. They do not
  depend on b.
* The raising operators u and w map back up and carry all the b-dependence,
  through four scales proportional to b (blocks L2 → L1 and L4 → L3) or to
  a+1−b (blocks L3 → L1 and L4 → L2).

The signs are chosen so that v has non-negative entries on the chain
L1 → L2 → L4. In the fundamental case a = 0 the basis is then exactly

    L,  vL,  fvL,  vfvL

with v = E₂₁ + E₄₃ (E_ij is the matrix unit). Explicitly, for a = 0,

    Y = diag(2b, 2b−1, 2b−1, 2b−2)      e = E₂₃     f = E₃₂     h = diag(0, 1, −1, 0)
    u = b E₁₂ + (b−1) E₃₄               v = E₂₁ + E₄₃
    w = −b E₁₃ + (b−1) E₂₄              x = E₃₁ − E₄₂

For a ≥ 1, the relation v² = 0 forces a minus sign on one of the two paths
L1 → L4. It is placed on the block L3 → L4.

### Typicality

With the ρ-shifted weights

    μ = a + 1,      λ = y − 1 = 2b − a − 1

the two positive odd roots give the factors μ + λ = 2b and μ − λ = 2(a+1−b).
The module is **atypical** exactly when one of them vanishes:

    b = 0      or      b = a + 1

These are the zeros of the four raising scales. There the module is still a
representation, but it becomes reducible and indecomposable. Landmarks:

| module                 | a | b   | y  |
|------------------------|---|-----|----|
| atypical R(0,0)        | 0 | 0   | 0  |
| atypical fundamental   | 0 | 1   | 2  |
| atypical fundamental   | 1 | 0   | −1 |
| adjoint                | 1 | 1   | 1  |

## 4. The grading operator χ

χ is diagonal, constant on each layer, with values

    χ = −1 on L1,   +1 on L2 and L3,   −1 on L4,          χ² = 1

It commutes with the even generators and anticommutes with the odd ones
(*verified*). The supertrace is STr(M) = Tr(χ M).

## 5. Indecomposable representations: the Matryoshka and its jets

Replace the weight y by y·1 + n, where n is the nilpotent shift of order N
(n^N = 0). Every generator μ_k(y) then becomes an N×N block matrix, block
lower triangular, whose p-th sub-diagonal carries the Taylor coefficient

    block (I, J) = (1/p!) (d/dy)^p μ_k ,      p = I − J ≥ 0

(p = 0 is the diagonal). Because n is nilpotent the series terminates, so the
substitution is exact. The result satisfies the same y-independent structure
constants as R(a,b). It is an indecomposable representation of dimension
4(a+1)·N, called the N-layer Matryoshka. N = 1 is R(a,b) itself.

In the basis of §3 every generator is affine in b, so only the orders p = 0
and p = 1 survive and the Matryoshka has two bands. `sl21.py` builds it from
the exact derivative R′ = R(a,1) − R(a,0) = dR/db. `sl21_b.py` takes true
derivatives of every order with respect to y (so the Y block on the first
sub-diagonal is the identity). The two conventions differ by the similarity
diag(1, ½, ¼, …) acting on the generation index. In the ζ-Hermitian basis
(§10) the entries are square roots of b, and every order p is alive.

On the Matryoshka, C₂, C₃ and T are no longer multiples of the identity. The
identities of §7 and §8 are nevertheless *verified* as operator identities on
the full matrices.

## 6. Invariant forms and the master equation for the cubic tensor

### The invariant form

On any module, the supertrace defines an invariant form on the generators:

    g_AB = (1/2N) STr(M_A M_B + M_B M_A)      (A, B even)
    g_AB = (1/2N) STr(M_A M_B − M_B M_A)      (A, B odd)

sl(2|1) is simple, so this form is unique up to a scalar:

    g_AB(R) = I(R) · g_AB(ref),       I(R) = μ = a + 1   (Dynkin index)

Here ref is the anchor module R(0,0). Its form is half the Killing form. In
the basis (Y, e, f, h, u, v, w, x) it reads

    g_YY = −2    g_ef = g_fe = 1    g_hh = 2    g_uv = −g_vu = −1    g_xw = −g_wx = −1

The odd block is a symplectic form ω on the four-dimensional odd space. With
the pairing (u,v), (x,w) its two blocks are equal.

### The cubic tensor and the master equation

    d_ABC = (1/2) STr( M_A [[M_B, M_C]] )

where [[B, C]] = BC − CB if B and C are both odd, and BC + CB otherwise.
The cubic tensor of every Kac module is one fixed tensor times one number,
the ordinary trace of the hypercharge:

    d_ABC(R) = Tr(Y)_R · c_ABC,       c_ABC = −d_ABC(R(0,0)) / 4

    Tr(Y) = 4 (a+1)(y−1) = 4 μ λ

This is *verified* on all 512 components. All the module dependence of the
cubic tensor sits in Tr(Y), exactly as for a gauge anomaly. So d_ABC is
additive over composite modules, and vanishes on any module with Tr(Y) = 0.
The self-conjugate line b = (a+1)/2 is such a case, and it contains the
adjoint R(1,1). (The supertrace STr(Y) vanishes identically, so the ordinary
trace is the relevant one.) 46 of the 512 components of c_ABC are nonzero.
The module-dependent numbers are:

    anomaly coefficient     A(R)  = −μ λ = −Tr(Y)/4
    cubic constant          d_YYY = STr(Y³) = −6 μ λ

## 7. The Casimir operators C₂ and C₃

C₂ and C₃ are fixed elements of the universal enveloping algebra U(g). They
are written down once and evaluated in whatever module is under study.

### Quadratic Casimir

    C₂ = h² − Y² + 2(ef + fe) + 2 N₂,        N₂ = uv − vu + xw − wx

This is an explicit polynomial in the generators and involves no metric. Its
eigenvalue on R(a,b) is the atypicality polynomial:

    C₂ = μ² − λ² = −4 b (b − a − 1)

It vanishes exactly at the atypical weights b = 0 and b = a+1. This is Kac's
typicality condition written as a Casimir. The normalisation is monic in the
atypicality polynomial, so C₂ = 4 on the adjoint.

### Cubic Casimir

    C₃ = 8 · Σ d^ABC(ref) M_A M_B M_C

Here d^ABC(ref) is the cubic tensor of the anchor R(0,0), with indices raised
by the inverse of its form. It is computed once and never re-read from the
module under study. The master equation (§6) is what makes this legitimate.
Any anchor with y = 0 gives the same tensor after multiplication by (a+1)²
(*verified* between R(0,0) and R(1,½)). The eigenvalue is

    C₃ = λ (μ² − λ²) = λ · C₂

so the pair (C₂, C₃) recovers the weight: λ = C₃/C₂ and μ² = C₂ + λ². The two
Casimirs therefore separate the typical modules.

Both operators commute with all eight generators, and this is *verified* on
the Matryoshka as well. There they are not scalar, and commutation (not
scalarity) is the defining property.

Contracting with a module's own form instead of the anchor's drags in the
Dynkin index. For example, (1/2) g^AB(R) M_A M_B = C₂ / (4μ), and the cubic
analogue becomes quadratic in λ. The programs report these contractions as
diagnostics only.

## 8. The Gorelik anticenter T

    T = (2/3) (N₄ − N₂)

Here N₂ is as in §7, and N₄ is the fully antisymmetrised product of the four
odd generators:

    N₄ = Σ_{σ ∈ S₄} sign(σ) M_σ(1) M_σ(2) M_σ(3) M_σ(4),     reference order (u, v, x, w)

T commutes with the even generators and anticommutes with the odd ones: it
belongs to the anticenter of U(g). The factor 2/3 is Gorelik's monic
normalisation. The Harish-Chandra image of the anticenter is the product over
the positive odd roots α of ⟨Λ+ρ | α⟩, which here is μ² − λ², and N₄ − N₂ has
eigenvalue (3/2)(μ² − λ²). The following are *verified* as operator
identities on the full matrices, including the Matryoshka:

    T = C₂ · χ,       [C₂, χ] = 0,       T² = C₂²

T² being a polynomial in the Casimirs, it is central (Gorelik's theorem).
T = C₂·χ is special to sl(2|1): there are exactly two positive odd roots, so
the atypicality polynomial has the degree of C₂. For sl(m|n) the anticenter
has degree m·n and no such relation holds, while Gorelik's monic
normalisation still applies.

### Coordinate-invariant form

N₄ alone depends on the choice of odd basis: under M_i → A_ij M_j it picks up
det A. Let g^ij be the inverse of the symplectic form ω of §6, and define

    Q₂ = g^ij M_i M_j,
    Q₄ = (g^ij g^kl − g^ik g^jl + g^il g^jk) M_i M_j M_k M_l = N₄ / Pf(ω)

Both are invariant, and T = (2/3)(β Q₄ + α Q₂). A rescaling of the form
changes α and β separately but not the ratio α²/β = 1. With the anchor's
form, Pf(ω) = 1 and Q₂ = N₂, Q₄ = N₄, so the Chevalley expression above
carries no hidden coefficient. With the Killing form (twice the anchor's),
N₄ − N₂ = 4 Q₄ − 2 Q₂.

## 9. Formal b: `sl21_b.py`

`sl21_b.py` builds the same modules with b kept as a SymPy symbol, declared
positive, so that every matrix and every invariant is an exact function of b.
There is no `-b` option. Two things require the symbol:

* derivatives of every order, used for the jets of §5;
* the ζ-Hermitian real forms of §10, whose entries are √b and √(b−a−1).

The anchors of the cubic tensor stay numeric. `--casimirs` is exact here too,
but slower than in `sl21.py`, because the metric and the 512 components of the
d-tensor are symbolic.

## 10. The ζ-Hermitian real forms

### The gauge

A diagonal similarity D, constant on each sl(2) layer up to the standard sl(2)
normalisation inside the layer, redistributes the b-dependence symmetrically
between raising and lowering operators. The odd entries then become the square
roots of the two Kac conditions:

    α = √b = √⟨Λ+ρ | β₁⟩,         γ = √(b − a − 1) = √⟨Λ+ρ | β₂⟩

(β₁, β₂ are the positive odd roots). The layer constants of D are

    L1: a+1      L2: α      L3: γ √(a(a+1))      L4: αγ

For a = 0 the result is

    u = α E₁₂ + γ E₃₄          v = α E₂₁ + γ E₄₃
    w = −α E₁₃ + γ E₂₄         x = α E₃₁ − γ E₄₂

D is singular at the atypical points, which are therefore excluded.

### Reality condition

Whether α and γ are real or imaginary depends on the region of b. The
programs treat three cases:

| case | region        | sign K₁ = b | sign K₂ = b−a−1 | ζ per layer (L1…L4) |
|------|---------------|-------------|-----------------|---------------------|
| 1    | b > a+1       | +           | +               | (1, 1, 1, 1)        |
| 2    | 0 < b < a+1   | +           | −               | (1, 1, −1, −1)      |
| 3    | b < 0         | −           | −               | (1, −1, −1, 1)      |

In general ζ = diag(1, sign K₁, sign K₂, sign K₁·sign K₂) per layer. At a = 0
the layer L3 is empty, and in case 2 this reads diag(1, 1, 1, −1) on the four
states. Complex conjugation flips the sign of each imaginary radical. The
same eight matrices then satisfy, in all three cases,

    ζ Y† ζ = Y     ζ h† ζ = h     ζ e† ζ = f     ζ u† ζ = v     ζ w† ζ = −x

together with the conjugate relations and ζ² = 1 (*verified*). The sign in
ζ w† ζ = −x is forced by the structure constants. Take the adjoint of
[e,u] = w: this gives [v, f] = ζ w† ζ, so ζ w† ζ = −[f,v] = −x. The Hermitian
combinations are e, f = (λ₁ ± iλ₂)/2, u, v = (λ₆ ± iλ₇)/2, and
w = (λ₄ + iλ₅)/2, x = −(λ₄ − iλ₅)/2, where the λ's are ζ-Hermitian.

For the Matryoshka a block-diagonal ζ cannot work, since the adjoint of a
lower-triangular matrix is upper triangular. The metric is instead
ζ_N = J ⊗ ζ, with J the reversal of the N generations. It pairs generation p
with generation N+1−p. The condition then holds band by band, because each
Taylor coefficient is separately ζ-Hermitian (*verified*).

### The block rescaling M (`--rescale`)

    M = diag( 1·Id_{a+1},  α·Id_{a+2},  γ·Id_a,  1·Id_{a+1} )

M commutes with e, f, h and Y. With `--rescale`, every generator is printed as
M μ M⁻¹. Since M is not unitary, the reality condition holds with the
transported metric

    η = M⁻† ζ M⁻¹,        η⁻¹ μ† η = ± partner

For the Matryoshka, M is replaced by its own jet M_N. For a = 0:

    u = E₁₂ + αγ E₃₄               v = b E₂₁ + (γ/α) E₄₃
    w = −E₁₃ + αγ E₂₄              x = b E₃₁ − (γ/α) E₄₂

## 11. Command-line reference

Each program prints its help when called with no arguments or with `-h`. On
Windows, type `py` instead of `python3`.

### `sl2.py`

    python3 sl2.py -a A            A = 0 … 20

### `sl21.py` (numeric weight)

    python3 sl21.py -a A -b B [-N N] [--casimirs]
    python3 sl21.py -a A -y Y [-N N] [--casimirs]

| option      | meaning                                                         |
|-------------|-----------------------------------------------------------------|
| `-a A`      | integer 0 … 10                                                  |
| `-b B`      | rational, e.g. `3/2`, `-7/3`                                    |
| `-y Y`      | alternative to `-b`; b = (a+y)/2                                |
| `-N N`      | Matryoshka depth 1 … 4 (default 1)                              |
| `--casimirs`| also the invariant form, master equation, C₂, C₃, T identities  |

### `sl21_b.py` (formal b)

    python3 sl21_b.py -a A [-N N] [--zetaH --case C [--rescale]] [--casimirs]

| option       | meaning                                                          |
|--------------|------------------------------------------------------------------|
| `-a A`, `-N N` | a ≥ 0, N ≥ 1, with (a+1)·N ≤ 16 (dimension ≤ 64)             |
| `--zetaH`    | ζ-Hermitian basis (§10); requires `--case`                       |
| `--case C`   | 1: b > a+1; 2: 0 < b < a+1; 3: b < 0                             |
| `--rescale`  | with `--zetaH`, print M μ M⁻¹ and check with η                   |
| `--casimirs` | as above, symbolic in b                                          |

### What a run prints

1. The generator matrices. For N > 1, `sl21_b.py` also reports which jet
   orders are nonzero.
2. The list of nonzero and null brackets checked, then
   `All sl(2|1) relations verified (nonzero and null). ✓`, or the residual
   matrices of whatever failed.
3. χ and its relations.
4. With `--zetaH`: ζ (or η), and the reality conditions.
5. With `--casimirs`: the invariant form and its inverse, Tr(Y), the master
   equation, C₂, C₃ and T with their eigenvalues next to the closed formulas,
   the weight recovered from (C₂, C₃), and the identities T = C₂χ, T² = C₂².

Any line containing `FAILED` or `✗` indicates a relation that does not hold.
The exception is the line stating that C₂ is not a multiple of the identity
for N > 1, which is expected (§5). Output can be redirected to a file:
`python3 sl21.py -a 1 -b 1 --casimirs > R11.txt`.

### Examples

    python3 sl21.py -a 1 -b 1 --casimirs          # the adjoint
    python3 sl21.py -a 0 -b 1/2 -N 3 --casimirs   # 3-layer Matryoshka of the fundamental
    python3 sl21_b.py -a 1                        # R(1,b) as functions of b
    python3 sl21_b.py -a 0 --zetaH --case 2 -N 2  # real form, 2-layer jet
    python3 sl21_b.py -a 0 --zetaH --case 1 --rescale

## 12. Reference values

Closed forms on R(a,b):

| quantity | in (μ, λ)  | in (a, b)             | zeros                     |
|----------|------------|-----------------------|---------------------------|
| Tr(Y)    | 4μλ        | 4(a+1)(2b−a−1)        | b = (a+1)/2               |
| I(R)     | μ          | a+1                   | never                     |
| A(R)     | −μλ        | −(a+1)(2b−a−1)        | b = (a+1)/2               |
| d_YYY    | −6μλ       | −6(a+1)(2b−a−1)       | b = (a+1)/2               |
| C₂       | μ² − λ²    | −4b(b−a−1)            | b = 0, a+1                |
| C₃       | λ C₂       | (2b−a−1) C₂           | b = 0, (a+1)/2, a+1       |
| T        | C₂ χ       |                       | b = 0, a+1                |

Values printed by `sl21.py --casimirs`. T[0,0] = −C₂, because χ = −1 on the
top layer.

| a | b   | y    | μ | λ    | Tr(Y) | A(R) | d_YYY | C₂    | C₃      | T[0,0] |
|---|-----|------|---|------|-------|------|-------|-------|---------|--------|
| 0 | 0   | 0    | 1 | −1   | −4    | 1    | 6     | 0     | 0       | 0      |
| 0 | 1   | 2    | 1 | 1    | 4     | −1   | −6    | 0     | 0       | 0      |
| 1 | 0   | −1   | 2 | −2   | −16   | 4    | 24    | 0     | 0       | 0      |
| 1 | 1   | 1    | 2 | 0    | 0     | 0    | 0     | 4     | 0       | −4     |
| 1 | 1/2 | 0    | 2 | −1   | −8    | 2    | 12    | 3     | −3      | −3     |
| 0 | 9/8 | 9/4  | 1 | 5/4  | 5     | −5/4 | −15/2 | −9/16 | −45/64  | 9/16   |
| 2 | 3/2 | 1    | 3 | 0    | 0     | 0    | 0     | 9     | 0       | −9     |
| 2 | 3/4 | −1/2 | 3 | −3/2 | −18   | 9/2  | 27    | 27/4  | −81/8   | −27/4  |
| 3 | 3/4 | −3/2 | 4 | −5/2 | −40   | 10   | 60    | 39/4  | −195/8  | −39/4  |

The first three rows are atypical (C₂ = C₃ = T = 0). On the adjoint R(1,1),
Tr(Y) = 0 and the whole cubic tensor vanishes, not just C₃.

---

Authors: Jean Thierry-Mieg (NLM/NIH) and Claude. License: public domain; no
rights reserved.
