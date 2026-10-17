"""Independent checks. Dense sampling is a test oracle, never the solver."""
import unittest
import numpy as np
from ellipse_distance import Ellipse, closest_farthest, Demo


class SolverTests(unittest.TestCase):
    def check_case(self, e, p):
        near, far, dn, df, candidates = closest_farthest(e, p)
        y = (e.rotation.T @ (candidates-e.center).T).T
        np.testing.assert_allclose(np.sum((y/e.axes)**2, axis=1), 1, atol=1e-6)
        # Lagrange condition: displacement is parallel to the surface normal.
        normals = (e.rotation @ (y/e.axes**2).T).T
        delta = candidates-p
        cross = delta[:, 0]*normals[:, 1]-delta[:, 1]*normals[:, 0]
        np.testing.assert_allclose(cross, 0, atol=1e-6)
        sampled = np.linalg.norm(e.curve(30001)-p, axis=1)
        self.assertLessEqual(dn, sampled.min()+1e-6)
        self.assertGreaterEqual(df, sampled.max()-1e-6)
        # Distance is 1-Lipschitz along the curve. This bounds sampling error.
        bound = 2*np.pi*max(e.axes)/30000 + 1e-6
        self.assertLess(abs(dn-sampled.min()), bound)
        self.assertLess(abs(df-sampled.max()), bound)

    def test_known_center(self):
        _, _, dn, df, _ = closest_farthest(Ellipse([0, 0], [4, 2], 0), [0, 0])
        self.assertAlmostEqual(dn, 2)
        self.assertAlmostEqual(df, 4)

    def test_circle(self):
        for p in ([0, 0], [1, 0], [3, 4]):
            e = Ellipse([0, 0], [2, 2], .7)
            self.check_case(e, np.array(p))
            _, _, dn, df, _ = closest_farthest(e, p)
            self.assertAlmostEqual(dn, abs(np.linalg.norm(p)-2))
            self.assertAlmostEqual(df, np.linalg.norm(p)+2)

    def test_axes_boundary_and_singular_cases(self):
        e = Ellipse([2, -1], [4, 2], .73)
        for local in ([0, 0], [1, 0], [3, 0], [4, 0], [0, 2],
                      [0, 1], [8, 0], [0, -8]):
            self.check_case(e, e.center+e.rotation@local)
        t = .37
        p = e.center+e.rotation@(e.axes*[np.cos(t), np.sin(t)])
        self.check_case(e, p)
        self.assertLess(closest_farthest(e, p)[2], 1e-7)

    def test_random_cases(self):
        rng = np.random.default_rng(782)
        for _ in range(200):
            e = Ellipse(rng.uniform(-10, 10, 2), rng.uniform(.5, 5, 2),
                        rng.uniform(-np.pi, np.pi))
            p = e.center+rng.uniform(-12, 12, 2)
            self.check_case(e, p)

    def test_keyboard_callback(self):
        from matplotlib.backend_bases import KeyEvent
        import matplotlib.pyplot as plt
        demo = Demo(42)
        self.assertFalse(demo.solved)
        demo.fig.canvas.callbacks.process('key_press_event',
            KeyEvent('key_press_event', demo.fig.canvas, key='enter'))
        self.assertTrue(demo.solved)
        demo.fig.canvas.callbacks.process('key_press_event',
            KeyEvent('key_press_event', demo.fig.canvas, key='n'))
        self.assertFalse(demo.solved)
        plt.close(demo.fig)


if __name__ == '__main__':
    unittest.main(verbosity=2)
