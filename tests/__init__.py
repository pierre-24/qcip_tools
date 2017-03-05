import math


def float_almost_equals(a, b, threshold=1e-3):
    return math.fabs(a - b) < threshold


def array_almost_equals(a, b, threshold=1e-3):
    """Check if two arrays containing float number are almost equals"""

    if len(a) != len(b):
        raise RuntimeError()

    return all(float_almost_equals(a[i], b[i], threshold) for i in range(len(a)))
