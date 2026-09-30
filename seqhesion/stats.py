"""Small statistics kept here so the product does not need scikit-learn."""
import numpy as np


def adjusted_rand(a, b):
    """Adjusted Rand index of two labelings (Hubert & Arabie), as sklearn computes it."""
    a, b = np.asarray(a), np.asarray(b)
    n = a.size
    if n < 2:
        return 1.0
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    pair, counts = np.unique(ai.astype(np.int64) * (bi.max() + 1) + bi, return_counts=True)
    comb = lambda x: (x * (x - 1) // 2).astype(np.float64).sum()  # noqa: E731
    s_ij = comb(counts)
    s_a = comb(np.bincount(ai))
    s_b = comb(np.bincount(bi))
    total = n * (n - 1) / 2
    expected = s_a * s_b / total
    top = (s_a + s_b) / 2
    if top == expected:
        return 1.0
    return float((s_ij - expected) / (top - expected))
