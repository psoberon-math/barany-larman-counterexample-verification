# Exact verification of the affine colorful Tverberg counterexample

Ancillary code for Section 3, Table 1, of Pablo Soberón's manuscript
**The topological Bárány–Larman conjecture for prime numbers**.

The program verifies that the thirteen points below, with color partition
$(3,3,3,3,1)$, have no four pairwise disjoint nonempty rainbow subsets whose
convex hulls intersect. It uses Sarkaria's transformation and exhaustive
**exact-arithmetic** checking.

## Run the verification

Requires **Python 3.9 or later**, with no third-party packages.

```sh
python3 verify_affine_sarkaria.py
```

The coordinates and coloring are embedded in the source. The program does
not read or write other files or access the network. It prints its summary
only after the complete scan; runtime depends on the machine.

Successful output is:

```text
Assignments checked: 55,296
Nonempty rainbow partitions: 54,648
Lifted matrices of rank 12: 55,296
Origin captures: 0
PASS: the affine colorful Tverberg counterexample is verified.
```

Exit status `0` means the complete verification succeeded. A detected
verification error prints `NOT VERIFIED` to standard error and exits with
status `1`; invalid command-line arguments give status `2`. An unexpected
lower rank stops the calculation without certifying the example. All checks
remain active with `python3 -O`: none relies on an `assert` statement.

## Coordinates, colors, and notation

Point indices are **0 through 12**, as in Table 1. Write $b_i$ for the point
of index $i$. The colors are $X_1,\ldots,X_5$, and the proposed parts are
$B_1,\ldots,B_4$. Thus colors and part labels are distinct notions.

| Index $i$ | Color | Coordinates of $b_i$ |
|---:|:---:|:---|
| 0 | $X_1$ | $(13,60,-64)$ |
| 1 | $X_1$ | $(355,256,429)$ |
| 2 | $X_1$ | $(-85,-35,-309)$ |
| 3 | $X_2$ | $(17,59,-62)$ |
| 4 | $X_2$ | $(174,-631,165)$ |
| 5 | $X_2$ | $(18,61,-63)$ |
| 6 | $X_3$ | $(12,-114,-79)$ |
| 7 | $X_3$ | $(95,132,14)$ |
| 8 | $X_3$ | $(-583,0,371)$ |
| 9 | $X_4$ | $(14,80,-66)$ |
| 10 | $X_4$ | $(17,60,-61)$ |
| 11 | $X_4$ | $(25,56,-64)$ |
| 12 | $X_5$ | $(-72,17,-212)$ |


## Why the computation proves the claim

### 1. Enumerating the rainbow assignments

An assignment $\varphi:\{0,\ldots,12\}\to\{1,2,3,4\}$ places $b_i$ in
$B_{\varphi(i)}$. It is rainbow exactly when the labels on each color class
are distinct.

Fix

$$
\varphi(0)=1,\qquad \varphi(1)=2,\qquad \varphi(2)=3.
$$

Every full rainbow assignment has exactly one such representative under
permutations of the four part labels. Each of the **remaining three** triples
has $4\cdot3\cdot2=24$ injective assignments, and the singleton has four
choices. Therefore the program checks

$$
4\cdot24^3=55\,296
$$

assignments. It checks their validity, uniqueness, and total count. It makes
no restriction on part sizes. The list includes 648 assignments with an
empty part, which cannot capture the origin and are harmless to check.
The other 54,648 assignments have four nonempty parts.

Any collection of four disjoint partial rainbow subsets extends to a full
assignment: add each unused point to a part missing its color. There are
enough available labels because each color class has at most three points.
Adding points only enlarges convex hulls. Consequently, excluding every full
assignment excludes all partial witnesses as well. Boundary intersections
are included; no positivity or general-position assumption is used here.

**Enumeration wording in the draft.** The count $55\,296$ already removes
the symmetry of permuting the four parts by fixing the first triple. Without
that normalization, the count is $4\cdot24^4=1\,327\,104$. Thus the sentence
saying that symmetry reduction was not implemented should instead say that
no **further** reduction is used after fixing the first triple's labels.
Strictly, 55,296 counts assignments, including the empty-part cases, rather
than partitions into four nonempty sets.

### 2. Sarkaria's transformation

For $r=4$ and $d=3$, choose the integer simplex vertices

$$
v_1=(1,0,0),\quad v_2=(0,1,0),\quad
v_3=(0,0,1),\quad v_4=(-1,-1,-1).
$$

Their only linear relation, up to scaling, is $v_1+v_2+v_3+v_4=0$.
This is a particular choice of the simplex allowed in Section 3.
For an assignment $\varphi$, form the thirteen lifts

$$
q_i=(b_i,1)\otimes v_{\varphi(i)}\in\mathbb Z^{12}.
$$

Sarkaria's equivalence [Sar92, BO96] is

$$
\bigcap_{j=1}^4\operatorname{conv}(B_j)\ne\varnothing
\quad\Longleftrightarrow\quad
0\in\operatorname{conv}\{q_0,\ldots,q_{12}\}.
$$

For this choice of the $v_j$, a nonnegative dependence among the lifts says
that the four vectors $\sum_{\varphi(i)=j}\lambda_i(b_i,1)$ are equal.
If $\sum_i\lambda_i=1$, their last coordinates are all $1/4$; multiplying
each part's coefficients by four gives convex representations of a common
point. Conversely, divide the coefficients of four convex representations
of a common point by four to obtain a convex representation of the origin
by the lifts. Zero coefficients are allowed in both directions.

### 3. The rank and sign test

Let $M_\varphi$ be the $12\times13$ integer matrix with columns $q_i$.
For **each** assignment the program verifies that
$\operatorname{rank}M_\varphi=12$. Thus the kernel is one-dimensional.
It computes a nonzero primitive integer generator
$z=(z_0,\ldots,z_{12})$ and substitutes it into the original matrix to verify
$M_\varphi z=0$ exactly.

The origin belongs to the convex hull of the columns precisely when the
kernel contains a nonzero coordinatewise nonnegative vector. Since the
kernel is one-dimensional, this occurs exactly when $z$ or $-z$ is
coordinatewise nonnegative. The exclusion condition is therefore simply

$$
\min_i z_i<0<\max_i z_i.
$$

This condition holds for all 55,296 assignments. The uniqueness of the
linear dependence is a verified rank condition, not an assumption about
the coordinates. Unexpected lower rank is not interpreted as an exclusion:
it would require a different feasibility test.

## Reading and checking the program

The source contains type annotations, mathematical docstrings, examples,
and comments describing the exact elimination and its checks.

| Function | Purpose |
|---|---|
| `rainbow_assignments()` | Enumerate the normalized color-respecting assignments. |
| `sarkaria_matrix(points, assignment)` | Form the twelve tensor coordinates of each lifted point. |
| `kernel_vector(matrix)` | Check full row rank and compute an exact integer dependence. |
| `origin_in_convex_hull(matrix)` | Apply the sign test, with zero coefficients allowed. |
| `verify()` | Check every assignment, the enumeration, and the final counts. |
| `main()` | Print a result only after successful completion, or report an error. |

Elimination uses the fraction-free Bareiss update with integer arithmetic.
Every division is checked for zero remainder. Back substitution uses
Python's exact `Fraction` type. After clearing denominators, the resulting
integer vector is checked in the **unmodified input matrix**. No
floating-point arithmetic enters the verification.

The small examples in the docstrings can be run separately:

```sh
python3 -m doctest -v verify_affine_sarkaria.py
```

There are 12 such checks. They include origin containment in the interior
and on the boundary, exclusion, rank-deficient input, row swaps, a skipped
pivot column, tensor coordinate order, and assignment counts. These are
sanity checks, **not a substitute for running the full verification**.

The complete scan and all 12 documented checks were rerun successfully
on 28 September 2026 using CPython 3.13.5. The complete scan also succeeded
with optimization enabled (`python3 -O verify_affine_sarkaria.py`).

## Scope

This program certifies only the finite affine counterexample with four
parts in dimension three and the specified coloring. It does not search
for the points or explain their discovery. It does not verify the
prime-case positive theorem, the nonlinear planar example, or affine
independence of every four original points. The rank checks here concern
the **lifted** matrices. The program is not a general convex-hull solver
for rank-deficient matrices.

## References

The manuscript supplies the coordinates and the computational formulation.
Its references for the transformation are:

- **[Sar92]** Karanbir S. Sarkaria, *Tverberg's theorem via number fields*,
  Israel Journal of Mathematics **79** (1992), 317–320.
- **[BO96]** Imre Bárány and Shmuel Onn, *Colourful linear programming*,
  Integer Programming and Combinatorial Optimization (Vancouver, BC, 1996),
  pp. 1–15.
