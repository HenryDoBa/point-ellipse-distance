# Closest and farthest points on an ellipse

## Run the interactive application

Use Python 3.10 or newer in a desktop terminal (not an inline notebook plot).

```bash
python3 -m pip install -r requirements.txt
python3 ellipse_distance.py --seed 42
```

Omit `--seed 42` for a different random initial example each time.
Click inside the plot window to give it keyboard focus.

- **Enter or Space:** calculate and draw both extrema and their line segments.
- **N:** generate a new random ellipse and query point; clear the old solution.
- **S:** save the current figure as `example_result.png` in the working directory.
- **Q:** close the window.

On Windows, use `python` instead of `python3` if necessary.
If no window opens, use a local graphical environment with a Matplotlib GUI
backend (for example, TkAgg with Tk installed). The `--save` mode below works
without a display.

To reproduce the included figure:

```bash
python3 ellipse_distance.py --seed 42 --save example_result.png
```

The provided PNG is an export of the actual application's result figure.
If a literal desktop screenshot is required, run the interactive application,
press Enter, and capture the plot window using the operating system screenshot tool.

## Mathematical method

Reference: *Point-ellipse distance*, December 10, 2024, supplied with the assignment:
https://cv.inf.elte.hu/wp-content/uploads/2024/12/PontEllipseDist.pdf

Write the ellipse as

    x^T Q x + 2 g^T x + F = 0,
    Q = [[A, C/2], [C/2, B]],  g = [D/2, E/2]^T.

For the query x0, the Lagrangian for squared distance is

    J(x, lambda) = ||x - x0||^2 + lambda (x^T Q x + 2 g^T x + F).

Stationarity gives

    (I + lambda Q)x = x0 - lambda g.

The inverse formula in the reference produces x = [P1, P2]^T / P3, where

    P1 = (EC/4 - DB/2) lambda^2 + (B u0 - D/2 - C v0/2) lambda + u0
    P2 = (DC/4 - EA/2) lambda^2 + (A v0 - E/2 - C u0/2) lambda + v0
    P3 = (AB - C^2/4) lambda^2 + (A+B) lambda + 1.

Substituting into the ellipse equation and clearing the denominator gives

    A P1^2 + B P2^2 + C P1 P2 + D P1 P3 + E P2 P3 + F P3^2 = 0.

This is a polynomial of degree at most four. Each admissible real root gives
a stationary point. Evaluate the distances to all candidates and select the
smallest and largest. The ellipse is compact and smooth, so both global
extrema exist and satisfy the Lagrange condition.

### Equivalent coordinate system used in the code

The generated ellipse has center c, rotation matrix R, and semi-axes a,b.
Set s=max(a,b), r=(a/s,b/s), q=R^T(x0-c)/s, and x=c+s R y.
Translation, rotation, and uniform scaling preserve distance ordering.
In these coordinates:

    d1*y1^2 + d2*y2^2 = 1,  with di=1/ri^2.
    P1 = q1*(1+lambda*d2)
    P2 = q2*(1+lambda*d1)
    P3 = (1+lambda*d1)*(1+lambda*d2)
    d1*P1^2 + d2*P2^2 - P3^2 = 0.

Thus the code applies the same derivation to a simpler, centered diagonal
conic. Recover each candidate as yi=qi/(1+lambda*di), then x=c+s R y.
The multiplier belongs to this normalized coordinate system; it need not
have the same numerical value as the original-coordinate multiplier.
`Ellipse.conic()` returns the six general conic coefficients if needed.

`numpy.roots` computes polynomial roots using floating-point linear algebra.
This is an algebraic candidate method, not numerical optimization: there is
no gradient descent, numerical minimizer, or sampled search in the solver.
It is not a symbolic radical formula for the quartic.

### Singular cases and numerical scope

The inverse in the PDF is undefined when det(I+lambda*Q)=0. Simply discarding
such roots can miss a true extremum, especially for a query on an ellipse axis.
The code handles axis queries analytically: add the two endpoints on the
query axis, then consider the singular multiplier for the other axis and
recover its free coordinate from the ellipse constraint. At the center this
recovers all four endpoints. Circle queries use the radial formula; at a
circle's center every boundary point is both closest and farthest, and two
representative points are returned.

For ties, the function returns one minimizing/maximizing representative.
The candidate list retains distinct stationary points in nondegenerate cases.

Like all floating-point polynomial methods, this implementation has limits
for extremely eccentric ellipses, nearly repeated roots, or extreme input
magnitudes. The GUI generates moderate semi-axes in [1,4] and local query
offsets in [-6,6]. Small coordinate tolerances identify axis/circle cases;
constraint checks reject invalid roots. This is an educational implementation,
not an arbitrary-precision geometry library.

## Files and checks

- `ellipse_distance.py`: solver and keyboard-driven visualization.
- `example_result.png`: reproducible example result (seed 42).
- `test_solver.py`: analytical examples, degeneracies, 200 random cases, and a keyboard event check.
- `requirements.txt`: runtime dependencies.
- `HUONG_DAN_VI.md`: short Vietnamese explanation.

Run the checks:

```bash
MPLBACKEND=Agg python3 -m unittest -v test_solver.py
```

On Windows PowerShell, first run `$env:MPLBACKEND="Agg"`, then
`python -m unittest -v test_solver.py`.

Tests check the ellipse equation and Lagrange stationarity, and compare
extreme distances to an independent dense boundary sample with a sampling
error bound. Sampling is used only in the tests, never by the algorithm.
All five test methods passed in the preparation environment, including the
200 random cases. The keyboard callback was exercised with Matplotlib's
headless canvas; an actual desktop window should be checked on the user's PC.

Seed 42 result (rounded):

    Query:    ( 6.80329241,  2.88919018)
    Closest:  ( 4.25793906,  1.38373352), distance = 2.95723240
    Farthest: (-2.14096903, -1.73787343), distance = 10.07022991
