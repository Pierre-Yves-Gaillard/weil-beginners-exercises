# https://claude.ai/chat/7c1796ac-d088-493e-b106-f1df0431e5d9 (me)
# https://claude.ai/share/6a1f0822-ce0f-4b35-ae22-831a6971b522 (public)
#!/usr/bin/env python3
"""
Generate a LaTeX table summarizing the "nearest-integer" variant of the
extended Euclidean algorithm applied to two integers a0 > a1 > 0.

  Forward pass:
    q1 = round(a0/a1),  x0=1, y0=0, x1=0, y1=1
    a_i = a_{i-2} - q_{i-1} a_{i-1}
    x_i = x_{i-2} - q_{i-1} x_{i-1}
    y_i = y_{i-2} - q_{i-1} y_{i-1}
    q_i = round(a_{i-1}/a_i)
  continuing as long as a_i != 0. We stop when a_{n+1} = 0, so
  a_n = gcd(a0, a1) and a0*x_n + a1*y_n = a_n.

  Here "round" means round to the nearest integer, with an exact tie
  (i.e. a value of the form n + 1/2 for an integer n) broken by
  ROUNDING TO THE INTEGER OF LEAST ABSOLUTE VALUE:
      n + 1/2  is rounded to  n    if n >= 0,
      n + 1/2  is rounded to  n+1  if n <  0.
  (E.g. round(3/2) = 1, round(5/2) = 2, round(-3/2) = -1,
  round(-5/2) = -2.) This is implemented EXACTLY, using Python's
  fractions.Fraction -- no floating-point arithmetic is used anywhere
  in this computation, so there is no risk of floating-point error
  changing the result. Because q_i is now defined by rounding rather
  than by floor(), the a_i, x_i, y_i, q_i can be negative from i = 2
  onward (a0, a1 are still assumed positive).

  Backward pass (Bezout coefficients u_i, v_i for 0 <= i <= n-1),
  unchanged from before:
    u_{n-1} = 0,  v_{n-1} = 1
    u_{i-1} = v_i
    v_{i-1} = u_i - q_i * v_i         (for 1 <= i <= n-1)
  satisfying a_i*u_i + a_{i+1}*v_i = a_n, so in particular
    a0*u0 + a1*v0 = a_n = gcd(a0, a1).

HOW TO USE
----------
Just edit A0 and A1 below to whatever values you want, then run the
script (e.g. press F5 in IDLE, or run "python3 ext_euclid.py" in a
terminal). It will:
  - print the LaTeX table to the screen, and
  - write it to a file called ext_euclid.tex in the same folder as
    this script.
"""

import math
import os
from fractions import Fraction

# ----------------------------------------------------------------
# EDIT THESE TWO VALUES (must satisfy A0 > A1 > 0):
A0 = 3
A1 = 2
# Optional caption/label for the LaTeX table (set to None to omit):
CAPTION = None
LABEL = None
# Name of the output file (change extension to .txt if you prefer):
OUTPUT_FILENAME = "ext_euclid.tex"
# ----------------------------------------------------------------


def round_least_abs_on_tie(numerator: int, denominator: int) -> int:
    """
    Return the integer nearest to numerator/denominator. An exact tie
    (a value of the form n + 1/2 for an integer n) is broken by
    choosing whichever of n, n+1 has the smaller absolute value:
        n + 1/2 -> n    if n >= 0
        n + 1/2 -> n+1  if n <  0
    Exact: uses Python's Fraction, so there is no floating-point
    error, and it works correctly regardless of the signs of
    numerator/denominator.
    """
    frac = Fraction(numerator, denominator)
    n = frac.numerator // frac.denominator      # floor(frac); Fraction always
                                                 # normalizes to a positive
                                                 # denominator, so // is a
                                                 # true floor division here.
    remainder = frac - n                        # in [0, 1)
    half = Fraction(1, 2)
    if remainder < half:
        return n
    elif remainder > half:
        return n + 1
    else:
        # exact tie: frac = n + 1/2
        return n if n >= 0 else n + 1


def extended_euclid_table(a0: int, a1: int):
    """
    Compute the forward pass (a_i, x_i, y_i, q_i) of the nearest-
    integer extended Euclidean algorithm for a0 > a1 > 0, where each
    q_i is round(a_{i-1}/a_i), ties broken toward least absolute
    value (see round_least_abs_on_tie above).

    Returns a list of dicts with keys "i", "a", "x", "y", "q", for
    i = 0, ..., n, where a_{n+1} = 0. The "q" value for i = 0 is
    None (there is no q0); all other rows have an integer q.
    """
    if not (a0 > a1 > 0):
        raise ValueError("Require a0 > a1 > 0")

    rows = []

    # i = 0
    a_prev2, x_prev2, y_prev2 = a0, 1, 0
    rows.append({"i": 0, "a": a_prev2, "x": x_prev2, "y": y_prev2, "q": None})

    # i = 1
    q_prev = round_least_abs_on_tie(a_prev2, a1)   # q1 = round(a0/a1)
    a_prev1, x_prev1, y_prev1 = a1, 0, 1
    rows.append({"i": 1, "a": a_prev1, "x": x_prev1, "y": y_prev1, "q": q_prev})

    i = 2
    while True:
        a_i = a_prev2 - q_prev * a_prev1
        x_i = x_prev2 - q_prev * x_prev1
        y_i = y_prev2 - q_prev * y_prev1

        if a_i == 0:
            break

        q_i = round_least_abs_on_tie(a_prev1, a_i)   # round(a_{i-1}/a_i)
        rows.append({"i": i, "a": a_i, "x": x_i, "y": y_i, "q": q_i})

        # shift window
        a_prev2, x_prev2, y_prev2 = a_prev1, x_prev1, y_prev1
        a_prev1, x_prev1, y_prev1 = a_i, x_i, y_i
        q_prev = q_i
        i += 1

    return rows


def add_backward_uv(rows):
    """
    Given the forward-pass rows (as produced by extended_euclid_table),
    compute u_i, v_i for 0 <= i <= n-1 by backward induction, and add
    them (as "u", "v" keys) to each row in place. Row n gets u = v = None
    (there is no u_n, v_n in this scheme).
    """
    n = rows[-1]["i"]          # last index (a_n = gcd, a_{n+1} = 0)

    # q_by_index[i] = q_i, needed for the backward recursion
    q_by_index = {row["i"]: row["q"] for row in rows}

    u = {n - 1: 0, n: None}
    v = {n - 1: 1, n: None}

    for i in range(n - 1, 0, -1):
        q_i = q_by_index[i]
        u[i - 1] = v[i]
        v[i - 1] = u[i] - q_i * v[i]

    for row in rows:
        i = row["i"]
        row["u"] = u.get(i)
        row["v"] = v.get(i)


def to_latex_table(rows, caption=None, label=None):
    """
    Format the rows (with forward-pass and backward-pass data already
    filled in) as a standalone LaTeX table (table + tabular).
    """
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"    \centering")
    lines.append(r"    \begin{tabular}{|r|r|r|r|r|r|r|}")
    lines.append(r"    \hline")
    lines.append(r"    i & a  &  x &  y & q &  u &  v \\ \hline")
    for row in rows:
        q_str = "" if row["q"] is None else str(row["q"])
        u_str = "" if row["u"] is None else str(row["u"])
        v_str = "" if row["v"] is None else str(row["v"])
        lines.append(
            f"    {row['i']} & {row['a']} & {row['x']} & {row['y']} & "
            f"{q_str} & {u_str} & {v_str} \\\\ \\hline"
        )
    lines.append(r"    \end{tabular}")
    if caption:
        lines.append(f"    \\caption{{{caption}}}")
    if label:
        lines.append(f"    \\label{{{label}}}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def main():
    rows = extended_euclid_table(A0, A1)
    add_backward_uv(rows)
    latex = to_latex_table(rows, caption=CAPTION, label=LABEL)

    print(latex)

    # Sanity checks.
    n = rows[-1]["i"]
    a_n = rows[-1]["a"]
    x_n, y_n = rows[-1]["x"], rows[-1]["y"]
    assert A0 * x_n + A1 * y_n == a_n
    assert math.gcd(A0, A1) == abs(a_n)

    u0 = rows[0]["u"]
    v0 = rows[0]["v"]
    assert A0 * u0 + A1 * v0 == a_n

    for row in rows[:-1]:   # i = 0, ..., n-1
        i = row["i"]
        a_i, u_i, v_i = row["a"], row["u"], row["v"]
        a_ip1 = rows[i + 1]["a"]
        assert a_i * u_i + a_ip1 * v_i == a_n

    # Write the LaTeX table to a file next to this script.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, OUTPUT_FILENAME)
    with open(output_path, "w") as f:
        f.write(latex + "\n")

    print(f"\n[Table written to {output_path}]")


if __name__ == "__main__":
    main()
