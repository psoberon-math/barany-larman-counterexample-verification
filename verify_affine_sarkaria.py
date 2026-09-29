#!/usr/bin/env python3
"""Exact verification of the affine colorful Tverberg counterexample.

Companion to Section 3, Table 1, of Pablo Soberon's manuscript
"The topological Barany-Larman conjecture for prime numbers".

Run the complete verification:
    python3 verify_affine_sarkaria.py
Run the small examples in the function documentation (not the full proof):
    python3 -m doctest -v verify_affine_sarkaria.py

Python 3.9 or later; standard library only. No input files are needed.

NOTATION
    b_i = POINTS[i], for i = 0,...,12 (the manuscript's point indices).
    X_1,...,X_5 are the four consecutive triples and the final singleton.
    B_1,...,B_4 are the proposed parts; assignment[i] is the label of B_j.
    VECTORS[j - 1] = v_j, for j = 1,...,4.
    M has columns q_i = (b_i,1) tensor v_assignment[i] in Z^12.

Sarkaria's equivalence identifies a common point of the four convex hulls
with 0 in conv{q_i}. For every assignment, we check rank(M) = 12, compute
an integer generator z of its one-dimensional kernel, and check Mz = 0
in the ORIGINAL matrix. Both positive and negative entries of z exclude
origin containment. Zero entries are allowed, so boundary intersections
are not discarded. Lower rank is an error, NEVER a successful exclusion.

All arithmetic is exact: fraction-free integer elimination (Bareiss),
rational back substitution, then integer substitution and sign checks.
All verification checks remain active under python3 -O.

The supplied four-page draft has inconsistent color labels in Table 1.
This file follows its stated "four triples and a singleton" profile, with
coordinates and indices unchanged. README.md records the discrepancies
and the part-label symmetry already used in the count 55,296.
"""

from __future__ import annotations

from argparse import ArgumentParser
from collections.abc import Iterator, Sequence
from fractions import Fraction
from itertools import permutations, product
from math import gcd, lcm
import sys


Point = tuple[int, int, int]
Assignment = tuple[int, ...]
IntegerMatrix = list[list[int]]

# Table 1 coordinates, in their original order. Colors follow the stated
# profile (3,3,3,3,1), NOT the duplicated X_2 labels in the supplied draft.
POINTS: tuple[Point, ...] = (
    (13, 60, -64),       #  0: X_1
    (355, 256, 429),     #  1: X_1
    (-85, -35, -309),    #  2: X_1
    (17, 59, -62),       #  3: X_2
    (174, -631, 165),    #  4: X_2
    (18, 61, -63),       #  5: X_2
    (12, -114, -79),     #  6: X_3
    (95, 132, 14),       #  7: X_3
    (-583, 0, 371),      #  8: X_3
    (14, 80, -66),       #  9: X_4
    (17, 60, -61),       # 10: X_4
    (25, 56, -64),       # 11: X_4
    (-72, 17, -212),     # 12: X_5 (singleton)
)
COLOR_CLASSES = ((0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11), (12,))
PART_LABELS = (1, 2, 3, 4)
# An integer simplex with v_1 + v_2 + v_3 + v_4 = 0 as its only relation.
VECTORS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, -1, -1))
EXPECTED_ASSIGNMENTS = 4 * 24**3       # 55,296, including empty parts
EXPECTED_NONEMPTY = 54_648


class VerificationError(RuntimeError):
    """A failed check or unsupported rank: no certification may be issued."""


def rainbow_assignments() -> Iterator[Assignment]:
    """Yield all full rainbow assignments after fixing the first triple.

    Part labels are 1,...,4. Fix b_0,b_1,b_2 in B_1,B_2,B_3. This removes
    the symmetry of permuting the four parts. Each of X_2,X_3,X_4 has 24
    injections into the four labels; X_5 has four choices. The resulting
    4*24^3 assignments include five-point parts and empty parts.

    Any disjoint partial rainbow partition extends to a full assignment:
    put unused points into parts missing their color. Convex hulls only
    enlarge. Thus checking this list also rules out all partial witnesses.

    >>> next(rainbow_assignments())
    (1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1)
    >>> sum(1 for _ in rainbow_assignments())
    55296
    >>> sum(len(set(phi)) == 4 for phi in rainbow_assignments())
    54648
    """
    injections = tuple(permutations(PART_LABELS, 3))
    for labels_X2, labels_X3, labels_X4 in product(injections, repeat=3):
        for label_X5 in PART_LABELS:
            yield (1, 2, 3) + labels_X2 + labels_X3 + labels_X4 + (label_X5,)


def sarkaria_matrix(
    points: Sequence[Point], assignment: Sequence[int]
) -> IntegerMatrix:
    """Return the matrix with columns (b_i,1) tensor v_assignment[i].

    This function constructs lifts for any nonempty list of integer
    triples with one part label per point. It does not impose a coloring;
    the enumeration and verify() enforce the rainbow restriction.

    The tensor coordinate order is
        (x*v_j[0], x*v_j[1], x*v_j[2],
         y*v_j[0], y*v_j[1], y*v_j[2],
         z*v_j[0], z*v_j[1], z*v_j[2],
           v_j[0],   v_j[1],   v_j[2]).
    For the counterexample, the output has 12 rows and 13 columns.

    >>> M = sarkaria_matrix([(2, 3, 5)], [4])
    >>> [row[0] for row in M]
    [-2, -2, -2, -3, -3, -3, -5, -5, -5, -1, -1, -1]
    """
    if not points or len(points) != len(assignment):
        raise ValueError("Need at least one point and one part label per point.")
    if any(len(p) != 3 or any(type(t) is not int for t in p) for p in points):
        raise ValueError("Coordinates must be integer triples.")
    if any(type(j) is not int or j not in PART_LABELS for j in assignment):
        raise ValueError("Part labels must be integers in {1,2,3,4}.")
    columns = [
        tuple(a * b for a in (*point, 1) for b in VECTORS[label - 1])
        for point, label in zip(points, assignment)
    ]
    return [list(row) for row in zip(*columns)]


def kernel_vector(matrix: Sequence[Sequence[int]]) -> tuple[int, ...]:
    """Compute a primitive integer generator of a one-dimensional kernel.

    The input must be an m-by-(m+1) integer matrix with m > 0. Full row
    rank is VERIFIED by elimination, not assumed. If rank is lower, raise
    VerificationError: one vector's signs cannot decide a larger kernel.

    The original matrix is not modified. Row swaps and skipped zero
    columns are allowed. Every Bareiss division is checked for zero
    remainder. Exact Fraction back substitution produces a rational
    dependence; clearing denominators and dividing by the gcd makes it
    primitive. Its sole free coordinate is positive. Finally check the
    nonzero integer vector in the original equations.

    >>> kernel_vector([[2, 0, 1], [0, 3, 1]])
    (-3, -2, 6)
    >>> kernel_vector([[0, 1, 2], [1, 0, 3]])  # row swap
    (-3, -2, 1)
    >>> kernel_vector([[1, 2, 0], [0, 0, 1]])  # skipped column
    (-2, 1, 0)
    """
    m = len(matrix)
    n = len(matrix[0]) if m else 0
    if not m or n != m + 1 or any(len(row) != n for row in matrix):
        raise ValueError("Expected a nonempty m-by-(m+1) matrix.")
    if any(type(t) is not int for row in matrix for t in row):
        raise ValueError("Matrix entries must be integers.")
    a = [list(row) for row in matrix]
    pivots: list[int] = []
    previous = 1

    # Fraction-free elimination. Each step is a row operation over Q:
    # new row k = (pivot * old row k - entry * pivot row) / previous.
    # Its exact integer divisions prevent roundoff or tolerance decisions.
    for col in range(n):
        row = len(pivots)
        pivot_row = next((k for k in range(row, m) if a[k][col] != 0), None)
        if pivot_row is None:
            continue
        a[row], a[pivot_row] = a[pivot_row], a[row]
        pivot = a[row][col]
        for k in range(row + 1, m):
            entry = a[k][col]
            for j in range(col + 1, n):
                numerator = pivot * a[k][j] - entry * a[row][j]
                quotient, remainder = divmod(numerator, previous)
                if remainder != 0:
                    raise VerificationError("Nonexact elimination division.")
                a[k][j] = quotient
            a[k][col] = 0
        previous = pivot
        pivots.append(col)
        if len(pivots) == m:
            break

    if len(pivots) != m:
        raise VerificationError(
            f"Rank {len(pivots)}, expected {m}; no conclusion for this matrix."
        )

    # Exactly one coordinate is free. Setting it to 1 fixes the scale.
    free = next(j for j in range(n) if j not in pivots)
    x = [Fraction(0) for _ in range(n)]
    x[free] = Fraction(1)
    for row, col in reversed(list(enumerate(pivots))):
        tail = sum((a[row][j] * x[j] for j in range(col + 1, n)), Fraction(0))
        x[col] = -tail / a[row][col]

    # Convert to a nonzero primitive integer dependence, then independently
    # substitute it into the unmodified input matrix (not the echelon form).
    denominator = lcm(*(q.denominator for q in x))
    z = [q.numerator * (denominator // q.denominator) for q in x]
    divisor = gcd(*z)
    if divisor == 0:
        raise VerificationError("The computed dependence is zero.")
    result = tuple(q // divisor for q in z)
    if any(sum(a * b for a, b in zip(row, result)) for row in matrix):
        raise VerificationError("Computed dependence failed exact substitution.")
    return result


def origin_in_convex_hull(matrix: Sequence[Sequence[int]]) -> bool:
    """Test whether 0 lies in the convex hull of the columns.

    This is an exact test ONLY for m-by-(m+1) matrices of rank m, which
    kernel_vector() checks. A nonzero kernel vector z yields a convex
    combination of 0 iff z or -z is coordinatewise nonnegative. Zeros
    are allowed. If z has both signs, every nonzero multiple does too.

    >>> origin_in_convex_hull([[1, 0, -1], [0, 1, -1]])  # interior
    True
    >>> origin_in_convex_hull([[1, 0, 1], [0, 1, 1]])    # outside
    False
    >>> origin_in_convex_hull([[0, 1, 0], [0, 0, 1]])    # boundary
    True
    >>> try:  # lower rank must stop the test, not give False
    ...     origin_in_convex_hull([[1, 2, 3], [2, 4, 6]])
    ... except VerificationError as exc:
    ...     print(exc)
    Rank 1, expected 2; no conclusion for this matrix.
    """
    z = kernel_vector(matrix)
    return not (min(z) < 0 < max(z))


def verify() -> tuple[int, int]:
    """Check every assignment and return (total, nonempty) only on success.

    In addition to the rank, residual, and sign checks, check that all
    assignments are distinct, respect the coloring, have the fixed first
    triple, and attain both expected counts. No part-size pruning is used.
    Any failed check aborts the scan without printing PASS.
    """
    if len(POINTS) != 13 or len(set(POINTS)) != 13:
        raise VerificationError("Expected thirteen distinct indexed points.")
    if tuple(i for color in COLOR_CLASSES for i in color) != tuple(range(13)):
        raise VerificationError("Color classes do not partition indices 0,...,12.")
    if tuple(map(len, COLOR_CLASSES)) != (3, 3, 3, 3, 1):
        raise VerificationError("Expected color profile (3,3,3,3,1).")

    seen: set[Assignment] = set()
    nonempty = 0
    for assignment in rainbow_assignments():
        if len(assignment) != 13 or assignment[:3] != (1, 2, 3):
            raise VerificationError(f"Invalid normalized assignment: {assignment}.")
        if assignment in seen:
            raise VerificationError(f"Repeated assignment: {assignment}.")
        if any(len({assignment[i] for i in c}) != len(c) for c in COLOR_CLASSES):
            raise VerificationError(f"Non-rainbow assignment: {assignment}.")
        seen.add(assignment)
        try:
            matrix = sarkaria_matrix(POINTS, assignment)
            if origin_in_convex_hull(matrix):
                raise VerificationError("Origin CAPTURED; not an exclusion.")
        except (VerificationError, ArithmeticError, ValueError) as exc:
            raise VerificationError(f"Assignment {assignment}: {exc}") from exc
        nonempty += len(set(assignment)) == 4

    checked = len(seen)
    if checked != EXPECTED_ASSIGNMENTS:
        raise VerificationError(
            f"Checked {checked} assignments, expected {EXPECTED_ASSIGNMENTS}."
        )
    if nonempty != EXPECTED_NONEMPTY:
        raise VerificationError(
            f"Found {nonempty} nonempty partitions, expected {EXPECTED_NONEMPTY}."
        )
    return checked, nonempty


def main() -> int:
    """Run the full check; return 0 on success and 1 on a verification error."""
    parser = ArgumentParser(
        description="Verify the 13-point affine counterexample by exact Sarkaria lifts."
    )
    parser.parse_args()
    try:
        checked, nonempty = verify()
    except (VerificationError, ArithmeticError, ValueError) as exc:
        print(f"NOT VERIFIED: {exc}", file=sys.stderr)
        return 1
    print(f"Assignments checked: {checked:,}")
    print(f"Nonempty rainbow partitions: {nonempty:,}")
    print(f"Lifted matrices of rank 12: {checked:,}")
    print("Origin captures: 0")
    print("PASS: the affine colorful Tverberg counterexample is verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
