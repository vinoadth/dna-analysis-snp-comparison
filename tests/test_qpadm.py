import unittest

import numpy as np

from dna_compare.comparisons.qpadm import f4, qpadm_weights


class QpAdmTests(unittest.TestCase):
    def test_f4_identities(self):
        a = np.linspace(0.1, 0.9, 200)
        b = np.linspace(0.9, 0.1, 200)
        c = np.linspace(0.2, 0.8, 200)
        d = np.linspace(0.05, 0.7, 200)
        self.assertAlmostEqual(f4(a, a, c, d), 0.0, places=12)
        self.assertAlmostEqual(f4(a, b, c, c), 0.0, places=12)

    def test_recovers_known_mixture(self):
        n = 4000
        grid = np.linspace(0.05, 0.95, n)
        s0 = grid
        s1 = np.clip(grid[::-1] * 0.7 + 0.1, 0, 1)
        s2 = np.clip((grid * 0.4 + 0.3) ** 1.0, 0, 1)
        rights = [np.clip((grid + offset) % 1.0, 0, 1) for offset in (0.02, 0.11, 0.23, 0.37, 0.51, 0.68)]
        truth = np.array([0.20, 0.35, 0.45])
        target = truth[0] * s0 + truth[1] * s1 + truth[2] * s2
        raw, normed, rss = qpadm_weights(target, [s0, s1, s2], rights)
        self.assertTrue(np.allclose(raw, truth, atol=1e-6), raw)
        self.assertTrue(np.allclose(normed, truth, atol=1e-6), normed)
        self.assertLess(rss, 1e-12)

    def test_recovers_five_source_mixture(self):
        n = 5000
        grid = np.linspace(0.05, 0.95, n)
        sources = [
            grid,
            np.clip(grid[::-1] * 0.7 + 0.1, 0, 1),
            np.clip((grid * 0.4 + 0.3), 0, 1),
            np.clip(np.sin(grid * 3) * 0.4 + 0.4, 0, 1),
            np.clip(np.cos(grid * 5) * 0.35 + 0.45, 0, 1),
        ]
        rights = [np.clip((grid + offset) % 1.0, 0, 1) for offset in (0.02, 0.11, 0.23, 0.37, 0.51, 0.68, 0.81)]
        truth = np.array([0.40, 0.30, 0.15, 0.10, 0.05])
        target = sum(w * src for w, src in zip(truth, sources))
        raw, normed, rss = qpadm_weights(target, sources, rights)
        self.assertTrue(np.allclose(raw, truth, atol=1e-6), raw)
        self.assertTrue(np.allclose(normed, truth, atol=1e-6), normed)
        self.assertLess(rss, 1e-12)


if __name__ == "__main__":
    unittest.main()
