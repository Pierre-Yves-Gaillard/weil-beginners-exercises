# Claude
#!/usr/bin/env python3
"""
Generate a LaTeX table summarizing the extended Euclidean algorithm
applied to two integers a0 > a1 > 0, using the notation:

    q1 = floor(a0/a1),  x0=1, y0=0, x1=0, y1=1
    a_i = a_{i-2} - q_{i-1} a_{i-1}
    x_i = x_{i-2} - q_{i-1} x_{i-1}
    y_i = y_{i-2} - q_{i-1} y_{i-1}
    q_i = floor(a_{i-1}/a_i)

continuing as long as a_i != 0.

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

# ----------------------------------------------------------------
# EDIT THESE TWO VALUES (must satisfy A0 > A1 > 0):
A0 = 309
A1 = 186
# Optional caption/label for the LaTeX table (set to None to omit):
CAPTION = None
LABEL = None
# Name of the output file (change extension to .txt if you prefer):
OUTPUT_FILENAME = "ext_euclid.tex"
# ----------------------------------------------------------------


def extended_euclid_table(a0: int, a1: int):
    """
    Compute the rows of the extended Euclidean algorithm table for
    a0 > a1 > 0.

    Returns a list of tuples (i, a_i, x_i, y_i, q_i) for i = 0, ..., n,
    where a_{n+1} = 0. q_0 is represented as None (blank in the table).
    """
    if not (a0 > a1 > 0):
        raise ValueError("Require a0 > a1 > 0")

    rows = []

    # i = 0
    a_prev2, x_prev2, y_prev2 = a0, 1, 0
    rows.append((0, a_prev2, x_prev2, y_prev2, None))

    # i = 1
    q_prev = a_prev2 // a1          # q1 = floor(a0/a1)
    a_prev1, x_prev1, y_prev1 = a1, 0, 1
    rows.append((1, a_prev1, x_prev1, y_prev1, q_prev))

    i = 2
    while True:
        a_i = a_prev2 - q_prev * a_prev1
        x_i = x_prev2 - q_prev * x_prev1
        y_i = y_prev2 - q_prev * y_prev1

        if a_i == 0:
            break

        q_i = a_prev1 // a_i         # floor(a_{i-1}/a_i)
        rows.append((i, a_i, x_i, y_i, q_i))

        # shift window
        a_prev2, x_prev2, y_prev2 = a_prev1, x_prev1, y_prev1
        a_prev1, x_prev1, y_prev1 = a_i, x_i, y_i
        q_prev = q_i
        i += 1

    return rows


def to_latex_table(rows, a0, a1, caption=None, label=None):
    """
    Format the rows as a standalone LaTeX table (table + tabular),
    in the style requested.
    """
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"    \centering")
    lines.append(r"    \begin{tabular}{|r|r|r|r|r|}")
    lines.append(r"    \hline")
    lines.append(r"    i & a  &  x &  y & q \\ \hline")
    for (i, a, x, y, q) in rows:
        q_str = "" if q is None else str(q)
        lines.append(f"    {i} & {a} & {x} & {y} & {q_str} \\\\ \\hline")
    lines.append(r"    \end{tabular}")
    if caption:
        lines.append(f"    \\caption{{{caption}}}")
    if label:
        lines.append(f"    \\label{{{label}}}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def main():
    rows = extended_euclid_table(A0, A1)
    latex = to_latex_table(rows, A0, A1, caption=CAPTION, label=LABEL)

    print(latex)

    # Sanity check: a0*x_n + a1*y_n = a_n = gcd(a0, a1)
    n = rows[-1][0]
    a_n, x_n, y_n = rows[-1][1], rows[-1][2], rows[-1][3]
    assert A0 * x_n + A1 * y_n == a_n
    assert math.gcd(A0, A1) == a_n

    # Write the LaTeX table to a file next to this script.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, OUTPUT_FILENAME)
    with open(output_path, "w") as f:
        f.write(latex + "\n")

    print(f"\n[Table written to {output_path}]")


if __name__ == "__main__":
    main()
