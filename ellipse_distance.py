"""Point-ellipse extrema using the Lagrange/quartic method in PontEllipseDist.pdf.

Run: python ellipse_distance.py
Keys: Enter/Space = solve, N = new random case, S = save figure, Q = quit.
No numerical optimizer or sampled search is used by the solver.
"""
from dataclasses import dataclass
import argparse
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as Poly


@dataclass
class Ellipse:
    center: np.ndarray
    axes: np.ndarray  # positive semi-axis lengths
    angle: float     # radians

    def __post_init__(self):
        self.center = np.asarray(self.center, dtype=float)
        self.axes = np.asarray(self.axes, dtype=float)
        if (self.center.shape != (2,) or self.axes.shape != (2,)
                or not np.all(np.isfinite(self.center))
                or not np.all(np.isfinite(self.axes))
                or np.any(self.axes <= 0) or not np.isfinite(self.angle)):
            raise ValueError('Finite center, angle and two positive semi-axes required.')

    @property
    def rotation(self):
        c, s = np.cos(self.angle), np.sin(self.angle)
        return np.array([[c, -s], [s, c]])

    def conic(self):
        """Return A,B,C,D,E,F for Au²+Bv²+Cuv+Du+Ev+F=0."""
        Q = self.rotation @ np.diag(1 / self.axes**2) @ self.rotation.T
        linear = -2 * Q @ self.center
        return np.array([Q[0, 0], Q[1, 1], 2*Q[0, 1], *linear,
                         self.center @ Q @ self.center - 1])

    def curve(self, count=600):
        t = np.linspace(0, 2*np.pi, count)
        return self.center + (self.rotation @
                             (self.axes[:, None]*[np.cos(t), np.sin(t)])).T


def stationary_points(ellipse, point):
    """All distinct stationary candidates (representatives for a centered circle).

    Translate, rotate and uniformly scale before applying the PDF's formula.
    These transformations preserve which point is closest/farthest.
    Axis cases are solved separately because I+lambda*Q may be singular.
    """
    point = np.asarray(point, dtype=float)
    if point.shape != (2,) or not np.all(np.isfinite(point)):
        raise ValueError('The query must be a finite 2D point.')
    scale = float(np.max(ellipse.axes))
    r = ellipse.axes / scale
    q = ellipse.rotation.T @ (point - ellipse.center) / scale
    d = 1/r**2
    tol = 1e-12
    candidates = []

    def add(y):
        y = np.asarray(y, dtype=float)
        if abs(np.sum((y/r)**2)-1) > 1e-6:
            return
        if not any(np.linalg.norm(y-z) < 1e-8 for z in candidates):
            candidates.append(y)

    if abs(r[0]-r[1]) < tol:
        # Circle: every boundary point is optimal when query equals center.
        norm = np.linalg.norm(q)
        direction = q/norm if norm > tol else np.array([1., 0.])
        add(r[0]*direction)
        add(-r[0]*direction)
    elif np.min(np.abs(q)) < tol:
        # Query lies on a principal axis. At a singular multiplier one
        # coordinate is free; enforce the ellipse constraint to recover it.
        j = int(np.argmin(np.abs(q)))
        i = 1-j
        for sign in (-1, 1):
            y = np.zeros(2)
            y[i] = sign*r[i]
            add(y)
        lam = -1/d[j]
        y = np.zeros(2)
        y[i] = q[i]/(1+lam*d[i])
        remainder = 1-(y[i]/r[i])**2
        if remainder >= -tol:
            for sign in (-1, 1):
                y[j] = sign*r[j]*np.sqrt(max(0., remainder))
                add(y.copy())
    else:
        # In principal coordinates D=E=C=0 and F=-1.
        # P1=q1(1+lambda*d2), P2=q2(1+lambda*d1),
        # P3=(1+lambda*d1)(1+lambda*d2).
        # Substitution into the constraint gives the PDF's quartic:
        # d1*P1² + d2*P2² - P3² = 0.
        h1, h2 = Poly([1., d[0]]), Poly([1., d[1]])
        p1, p2, p3 = q[0]*h2, q[1]*h1, h1*h2
        quartic = d[0]*p1*p1 + d[1]*p2*p2 - p3*p3
        coeff = quartic.coef / np.max(np.abs(quartic.coef))
        # Eigenvalue-based polynomial root calculation, NOT optimization.
        for root in np.roots(coeff[::-1]):
            if abs(root.imag) <= 1e-8*(1+abs(root.real)):
                lam = root.real
                denominator = 1+lam*d
                if np.min(np.abs(denominator)) > 1e-14:
                    add(q/denominator)
    if len(candidates) < 2:
        raise ArithmeticError('Insufficient valid roots; input is numerically ill-conditioned.')
    return ellipse.center + (ellipse.rotation @
                             (scale*np.asarray(candidates)).T).T


def closest_farthest(ellipse, point):
    """Return closest point, farthest point, their distances and all candidates."""
    candidates = stationary_points(ellipse, point)
    distances = np.linalg.norm(candidates-np.asarray(point), axis=1)
    lo, hi = np.argmin(distances), np.argmax(distances)
    return candidates[lo], candidates[hi], distances[lo], distances[hi], candidates


class Demo:
    def __init__(self, seed=None):
        import matplotlib.pyplot as plt
        self.plt = plt
        self.rng = np.random.default_rng(seed)
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.fig.subplots_adjust(bottom=.13, top=.86)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.new_case()

    def new_case(self):
        self.ellipse = Ellipse(self.rng.uniform(-2, 2, 2),
                               self.rng.uniform(1, 4, 2),
                               self.rng.uniform(0, np.pi))
        self.point = self.ellipse.center+self.rng.uniform(-6, 6, 2)
        self.solved = False
        self.draw()

    def solve(self):
        self.result = closest_farthest(self.ellipse, self.point)
        self.solved = True
        self.draw()
        near, far, dn, df, _ = self.result
        print(f'Query: {self.point}\nClosest: {near}, distance={dn:.8f}'
              f'\nFarthest: {far}, distance={df:.8f}')

    def draw(self):
        ax = self.ax
        ax.clear()
        e = self.ellipse
        curve = e.curve()
        ax.plot(*curve.T, color='#2463a5', lw=2.5, label='Ellipse')
        ax.fill(*curve.T, color='#2463a5', alpha=.05)
        ax.scatter(*e.center, marker='+', s=75, c='#64748b', label='Center')
        ax.scatter(*self.point, marker='*', s=180, c='#ee9b18',
                   edgecolors='#774400', label='Query point', zorder=5)
        title = 'Press Enter or Space to compute both extrema'
        if self.solved:
            near, far, dn, df, _ = self.result
            for p, color, name in [(near, '#13865f', 'Closest'),
                                    (far, '#cc4554', 'Farthest')]:
                ax.plot([self.point[0], p[0]], [self.point[1], p[1]],
                        '--', color=color, lw=1.8)
                ax.scatter(*p, s=80, c=color, zorder=6,
                           label=f'{name}: ({p[0]:.3f}, {p[1]:.3f})')
            title = f'Minimum distance = {dn:.6f}    |    Maximum distance = {df:.6f}'
        self.fig.suptitle('Closest and farthest points on an ellipse',
                          fontsize=18, fontweight='bold', y=.97)
        ax.set_title(title, fontsize=11, pad=14)
        ax.set_xlabel('u')
        ax.set_ylabel('v')
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(alpha=.2)
        ax.margins(.18)
        ax.legend(loc='best', fontsize=9, framealpha=.95)
        # Figure-level status is replaced, not accumulated across key presses.
        for artist in list(self.fig.texts)[1:]:
            artist.remove()
        self.fig.text(.5, .055,
                      f'Center=({e.center[0]:.2f}, {e.center[1]:.2f})   '
                      f'Semi-axes=({e.axes[0]:.2f}, {e.axes[1]:.2f})   '
                      f'Rotation={np.degrees(e.angle):.1f} degrees',
                      ha='center', fontsize=10)
        self.fig.text(.5, .025, 'Enter / Space: solve     N: new case     '
                      'S: save image     Q: quit', ha='center', fontsize=10)
        self.fig.canvas.draw_idle()

    def on_key(self, event):
        if event.key in ('enter', ' '):
            self.solve()
        elif event.key in ('n', 'N'):
            self.new_case()
        elif event.key in ('s', 'S'):
            self.fig.savefig('example_result.png', dpi=180)
            print('Saved example_result.png')
        elif event.key in ('q', 'Q'):
            self.plt.close(self.fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--save', type=Path,
                        help='Solve and export a PNG without opening a window.')
    args = parser.parse_args()
    if args.save:
        import matplotlib
        matplotlib.use('Agg')
    demo = Demo(args.seed)
    if args.save:
        demo.solve()
        demo.fig.savefig(args.save, dpi=180)
    else:
        demo.plt.show()


if __name__ == '__main__':
    main()
