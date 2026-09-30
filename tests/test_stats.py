import numpy as np
import pytest

from seqhesion.stats import adjusted_rand


def test_identical_partitions_score_one_whatever_the_labels():
    assert adjusted_rand([0, 0, 1, 1, 2], [5, 5, 3, 3, 9]) == 1.0


def test_matches_the_textbook_value():
    # by hand: sum C(n_ij, 2) = 2; rows (3, 3) give 6; columns (2, 2, 2) give 3;
    # expected = 6 * 3 / C(6, 2) = 1.2; ARI = (2 - 1.2) / ((6 + 3) / 2 - 1.2)
    a = [0, 0, 0, 1, 1, 1]
    b = [0, 0, 1, 1, 2, 2]
    assert adjusted_rand(a, b) == pytest.approx((2 - 1.2) / ((6 + 3) / 2 - 1.2))


def test_symmetric():
    rng = np.random.default_rng(0)
    a, b = rng.integers(0, 5, 50), rng.integers(0, 7, 50)
    assert adjusted_rand(a, b) == pytest.approx(adjusted_rand(b, a))
